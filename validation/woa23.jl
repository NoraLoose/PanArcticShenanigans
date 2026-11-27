# The World Ocean Atlas 2023 (WOA23; Reagan et al. 2024) temperature and salinity climatologies on their own 0.25° grid,
# a validation dataset: read a climatology, build an Oceananigans latitude-longitude grid whose cells are WOA's cells,
# and the fields of one month, so that the model's code (diagnostics/) runs on them
# (woa23_ocean_diagnostics!, for analysis/validation_diagnostics.jl).
# Needs validation_config.jl (the cache, <validation_dir>/WOA23/raw) and validation/glorys.jl (level_interfaces).
#
# Each WOA23 climatology averages one period (a decade, e.g. "A5B4" for 2005–2014) and comes as an annual (00), 12
# monthly (01–12, down to 1500 m) and 4 seasonal (13–16: Jan–Mar, …, Oct–Dec, down to 5500 m) means of the
# objectively analysed field ("_an"); at 0.25°, 1985–1994 and 1995–2004 have no monthly ones. The files are large (0.4–0.9 GB, all statistics, the whole globe), so each is
# downloaded once, its analysed field north of `south` kept in a small cache file, and the download deleted.
#
# WOA's temperature is in situ; the model's and GLORYS's are potential temperatures, so it is converted with TEOS-10
# (GibbsSeaWater: absolute salinity from WOA's practical salinity, then potential temperature at 0 dbar). Salinity stays
# practical, as in GLORYS. WOA has no sea floor: a column's bottom is the bottom of the deepest level with data in the
# annual climatology, in whole cells.

using NumericalEarth, Dates, JLD2
using NumericalEarth.DataWrangling: NCDatasets
using .NCDatasets: Dataset, dimnames
using GibbsSeaWater: gsw_p_from_z, gsw_sa_from_sp, gsw_pt0_from_t
using Downloads: Downloads
using Printf

const woa_url = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA"
woa_dir = raw_dir("WOA23")

woa_filename(variable, decade, period) = "woa23_$(decade)_$(variable)$(lpad(period, 2, '0'))_04.nc"   # 04: 0.25°
woa_file_url(variable, decade, period) =
    join((woa_url, variable == "t" ? "temperature" : "salinity", "netcdf", decade, "0.25", woa_filename(variable, decade, period)), "/")

# Whether NCEI has this climatology: at 0.25°, some decades have no monthly ones (1985–1994 and 1995–2004 have only the
# annual and seasonal climatologies)
function woa_available(variable, decade, period; south)
    cache = joinpath(woa_dir, decade, replace(woa_filename(variable, decade, period), ".nc" => "_an_$(south)N.jld2"))
    isfile(cache) && return true
    response = Downloads.request(woa_file_url(variable, decade, period); method = "HEAD", throw = false)
    return response isa Downloads.Response && response.status == 200
end

"""
    read_woa(variable, decade, period; south)

WOA23's analysed field of `variable` ("t" in situ temperature, °C, or "s" practical salinity) for `decade` (e.g.
"A5B4") and `period` (0 annual, 1–12 monthly, 13–16 seasonal), north of `south` (°N): (; lon, lat, depth, data), the
data (lon, lat, depth) from the surface down, NaN where WOA has none. Downloads and caches it on first use.
"""
function read_woa(variable, decade, period; south)
    file  = woa_filename(variable, decade, period)
    cache = joinpath(mkpath(joinpath(woa_dir, decade)), replace(file, ".nc" => "_an_$(south)N.jld2"))
    isfile(cache) && return jldopen(f -> (; lon = f["lon"], lat = f["lat"], depth = f["depth"], data = f["data"]), cache)

    path = joinpath(woa_dir, decade, file)
    if !isfile(path)
        url = woa_file_url(variable, decade, period)
        @info "Downloading $url"
        for attempt in 1:5   # NCEI's server stalls now and then
            try
                Downloads.download(url, path * ".part")
                break
            catch error
                attempt == 5 && rethrow()
                @warn "Download failed ($(sprint(showerror, error)[1:min(end, 120)])); retrying in $(60attempt) s"
                sleep(60attempt)
            end
        end
        mv(path * ".part", path; force = true)
    end

    ds  = Dataset(path)
    var = ds["$(variable)_an"]
    dimnames(var)[1:3] == ("lon", "lat", "depth") || error("unexpected dimensions $(dimnames(var)) in $path")
    lon, lat, depth = Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(ds["depth"][:])
    rows = findall(≥(south), lat)
    data = Float32.(coalesce.(var[:, rows, :, 1], NaN32))
    close(ds)

    jldsave(cache; lon, lat = lat[rows], depth, data)
    rm(path)
    return (; lon, lat = lat[rows], depth, data)
end

"""
    woa_grid(decade; south)

WOA's cells north of `south` for `decade`: an `ImmersedBoundaryGrid` on the CPU, periodic in longitude, whose sea floor is
the bottom of the deepest level with data in the annual climatology (whole cells). Returns the grid and the depths of
its levels (m, from the surface down).
"""
function woa_grid(decade; south)
    annual = read_woa("t", decade, 0; south)
    lon, lat, depth = annual.lon, annual.lat, annual.depth
    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    zf = level_interfaces(depth)
    deepest = [something(findlast(isfinite, annual.data[i, j, :]), 0) for i in axes(annual.data, 1), j in axes(annual.data, 2)]
    bottom  = [k == 0 ? 0.0 : -zf[k + 1] for k in deepest]

    underlying = LatitudeLongitudeGrid(CPU(); size = (length(lon), length(lat), length(depth)),
                                       longitude = (lon[1] - Δλ / 2, lon[1] - Δλ / 2 + 360),
                                       latitude  = vcat(lat .- Δφ / 2, min(90, lat[end] + Δφ / 2)),
                                       z = reverse(-zf), topology = (Periodic, Bounded, Bounded))
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom)), depth
end

# A field of `decade` for `month` (0: annual) on every level: the monthly climatology down to 1500 m (if `monthly`;
# otherwise the seasonal one there too), the seasonal one below, and the annual one wherever those have no data, so that
# every month covers the same cells
function woa_month(variable, decade, month; south, monthly)
    annual = read_woa(variable, decade, 0; south)
    field  = copy(annual.data)
    month == 0 && return field
    season = 13 + (month - 1) ÷ 3
    for period in (monthly ? (season, month) : (season,))   # the season, then the month
        (; depth, data) = read_woa(variable, decade, period; south)
        levels = 1:length(depth)
        depth == annual.depth[levels] || error("WOA23 $decade period $period has other levels than the annual climatology")
        fill = isfinite.(data)
        view(field, :, :, levels)[fill] .= data[fill]
    end
    return field
end

"""
    woa_fields(grid, decade, month; south, monthly = true, within = trues(size(grid)[1:2]))

WOA's potential temperature (°C, at 0 dbar) and practical salinity of `decade` for `month` (1–12, or 0 for the annual
climatology) on the levels of `woa_grid`, in Oceananigans' order (bottom first); NaN where WOA has no data and outside
the columns `within` (a Bool matrix; the conversion of the temperature is slow, so only those columns are converted).
With `monthly = false`, the month's season stands for it (decades without monthly climatologies).
"""
function woa_fields(grid, decade, month; south, monthly = true, within = trues(size(grid)[1:2]))
    t = woa_month("t", decade, month; south, monthly)
    S = woa_month("s", decade, month; south, monthly)
    depth = read_woa("t", decade, 0; south).depth
    Nx, Ny, Nz = size(t)
    λ = [λnode(i, 1, 1, grid, Center(), Center(), Center()) for i in 1:Nx]
    φ = [φnode(1, j, 1, grid, Center(), Center(), Center()) for j in 1:Ny]
    θ = fill(NaN, Nx, Ny, Nz)
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        (within[i, j] && isfinite(t[i, j, k]) && isfinite(S[i, j, k])) || continue
        p  = gsw_p_from_z(-depth[k], φ[j])
        SA = gsw_sa_from_sp(Float64(S[i, j, k]), p, λ[i], φ[j])
        θ[i, j, k] = gsw_pt0_from_t(SA, Float64(t[i, j, k]), p)
    end
    return reverse(θ, dims = 3), reverse(Float64.(S), dims = 3)
end

"""
    woa23_ocean_diagnostics!(decades; settings, set_names)

WOA23's heat, salt and freshwater content (diagnostics/ocean_content.jl) and T, S profiles (diagnostics/ts_profiles.jl)
for each decade of `decades` (code => [first year, last year]) and each region set of `set_names`, over the columns
inside the model domain (domain.jl): the 12 monthly climatologies (monthly down to 1500 m, seasonal below; the seasonal
ones for decades without monthly ones at 0.25°) and the annual one ("annual/…"), kept in
<validation_dir>/WOA23/<family>/<set>/<decade>.jld2 with the decade's regions, levels, "years" and "monthly". Decades
already there are skipped. Needs domain.jl, diagnostics/ocean_content.jl, diagnostics/ts_profiles.jl and
validation_config.jl.
"""
function woa23_ocean_diagnostics!(decades; settings, set_names)
    for (decade, years) in sort(collect(decades), by = last)
        all(isfile(period_file("WOA23", family, set, decade)) for set in set_names, family in ocean_families) && continue
        grid, depths = woa_grid(decade; south = φ_south)
        within = within_model_domain(grid)
        sets   = [region_masks(grid; set, within) for set in set_names]
        @info @sprintf("WOA23 %s (%d–%d): %d × %d columns, %d levels; %d columns in the model domain", decade, years...,
                       size(grid)[1:2]..., length(depths), count(within))
        for regions in sets, family in ocean_families
            save_dataset_regions("WOA23", family, regions, nothing; settings = family_file_settings(family, settings))
        end

        monthly = all(woa_available(v, decade, m; south = φ_south) for v in ("t", "s"), m in 1:12)
        monthly || @info "WOA23 $decade has no monthly climatologies at 0.25°: each month gets its season's"

        T, S = CenterField(grid), CenterField(grid)
        diagnostics(month) = begin   # month 0: the annual climatology
            θ, s = woa_fields(grid, decade, month; south = φ_south, monthly, within)
            set!(T, θ); set!(S, s)
            [Dict(family => ocean_family(family, grid, T, S, regions; settings) for family in ocean_families)
             for regions in sets]
        end
        annual = diagnostics(0)
        months = map(1:12) do month
            result = diagnostics(month)
            @info "WOA23 $decade, month $month"
            result
        end

        for (n, regions) in enumerate(sets), family in ocean_families
            save_dataset_period("WOA23", family, regions.set, decade, [m[n][family] for m in months]; dates = 1:12,
                                regions, grid, level_depths = reverse(depths), years, monthly,
                                (Symbol("annual/$name") => value for (name, value) in pairs(annual[n][family]))...)
        end
    end
end
