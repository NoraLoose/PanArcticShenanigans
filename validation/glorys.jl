# GLORYS12 daily means on their own 1/12° grid, a validation dataset: read a variable, and build an Oceananigans
# latitude-longitude grid whose cells are GLORYS's cells, with GLORYS's levels and sea floor (deptho, as partial bottom
# cells), so the model's code runs on it: regional_diagnostics.jl (glorys_regional_diagnostics!) and transports.jl
# (glorys_transports!), both for analysis/validation_diagnostics.jl.
# Needs `data_dir`, where the GLORYS downloads are (paths.jl: paths.toml or your paths_local.toml).

using NumericalEarth, Dates
using CopernicusMarine   # enables the GLORYS download extension
using NumericalEarth.DataWrangling: BoundingBox, Metadatum, NCDatasets, Downloads, default_download_directory
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily, GLORYSStatic, GLORYS_dataset_variable_names
using .NCDatasets: Dataset, dimnames
using Printf
using Oceananigans.BoundaryConditions: fill_halo_regions!

glorys_dir = isempty(data_dir) ? default_download_directory(GLORYSDaily()) : data_dir

# A GLORYS variable over `region` (one day, or static), with NaN on land: (; lon, lat, depth, data)
function read_glorys(name, dataset, region; date = nothing)
    metadata = isnothing(date) ? Metadatum(name; dataset, region, dir = glorys_dir) :
                                 Metadatum(name; dataset, date, region, dir = glorys_dir)
    path = Downloads.download(metadata)
    ds   = Dataset(path)
    var  = ds[GLORYS_dataset_variable_names[name]]
    lon, lat = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:])
    dims = dimnames(var)
    if "depth" in dims
        dims[1:3] == ("longitude", "latitude", "depth") || error("unexpected dimensions $dims in $path")
        data, depth = Float64.(coalesce.(var[:, :, :, 1], NaN)), Float64.(ds["depth"][:])
    else
        dims[1:2] == ("longitude", "latitude") || error("unexpected dimensions $dims in $path")
        data, depth = Float64.(coalesce.(ndims(var) == 3 ? var[:, :, 1] : var[:, :], NaN)), nothing
    end
    close(ds)
    return (; lon, lat, depth, data)
end

# Interfaces of GLORYS's levels (m below the surface), halfway between the level depths, the top one at the surface
function level_interfaces(depth)
    zf = vcat(0.0, (depth[1:end-1] .+ depth[2:end]) ./ 2)
    return vcat(zf, depth[end] + (depth[end] - zf[end]))
end

"""
    glorys_grid(region, date; sample = :temperature)

GLORYS's cells over `region` (a `BoundingBox`): an `ImmersedBoundaryGrid` on the CPU with GLORYS's sea floor (deptho) as
partial bottom cells, periodic in longitude if the region goes all the way round. Returns the grid, the levels kept
(every level with ocean somewhere in the region, numbered from the surface down) and their depths (m). The levels are
read from the `sample` variable on `date`.
"""
function glorys_grid(region, date; sample = :temperature)
    static  = read_glorys(:depth, GLORYSStatic(), region)
    sampled = read_glorys(sample, GLORYSDaily(), region; date)
    (static.lon ≈ sampled.lon && static.lat ≈ sampled.lat) || error("GLORYS's deptho and $sample are on different grids")

    lon, lat = static.lon, static.lat
    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    deptho = [isfinite(d) && d > 0 ? d : 0.0 for d in static.data]
    zf     = level_interfaces(sampled.depth)
    levels = 1:findlast(<(maximum(deptho)), zf[1:end-1])

    periodic  = isapprox(length(lon) * Δλ, 360; atol = Δλ / 2)
    west      = lon[1] - Δλ / 2
    east      = periodic ? west + 360 : lon[end] + Δλ / 2
    latitude  = vcat(lat .- Δφ / 2, min(90, lat[end] + Δφ / 2))   # faces; the last cell is cut at the pole
    underlying = LatitudeLongitudeGrid(CPU(); size = (length(lon), length(lat), length(levels)),
                                       longitude = (west, east), latitude,
                                       z = reverse(-zf[1:last(levels) + 1]),
                                       topology = (periodic ? Periodic : Bounded, Bounded, Bounded))
    grid = ImmersedBoundaryGrid(underlying, PartialCellBottom(-deptho))
    return grid, levels, sampled.depth[levels]
end

# A GLORYS variable of one day on the levels of glorys_grid, in Oceananigans' order (bottom first); NaN on land
glorys_level_data(name, region, levels, date) =
    reverse(read_glorys(name, GLORYSDaily(), region; date).data[:, :, levels], dims = 3)

"""
    glorys_regional_diagnostics!(days; settings, set_names)

GLORYS's regional diagnostics (regional_diagnostics.jl) of each day of `days` and each region set of `set_names`, over
the columns inside the model domain (domain.jl), kept in <validation_dir>/GLORYS/diagnostics/<set>/<day>.jld2 (with
the regions once in regions.jld2). Days already there are skipped, so runs over the same days share them. Needs
domain.jl, regional_diagnostics.jl and validation_config.jl.
"""
function glorys_regional_diagnostics!(days; settings, set_names)
    region = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))
    grid, levels, depths = glorys_grid(region, first(days))
    within = within_model_domain(grid)
    sets   = [region_masks(grid; set, within) for set in set_names]
    @info @sprintf("GLORYS: %d × %d columns, %d levels; %d columns in the model domain", size(grid)[1:2]...,
                   length(levels), count(within))
    for regions in sets
        save_dataset_regions("GLORYS", regions, grid; level_depths = reverse(depths), settings)
    end

    T, S = CenterField(grid), CenterField(grid)
    for date in days
        missing_sets = filter(regions -> !isfile(period_file("GLORYS", regions.set, date)), sets)
        isempty(missing_sets) && continue
        set!(T, glorys_level_data(:temperature, region, levels, date))
        set!(S, glorys_level_data(:salinity,    region, levels, date))
        for regions in missing_sets
            result = regional_diagnostics(grid, T, S, regions; settings...)
            save_dataset_period("GLORYS", regions.set, date, [result]; dates = [date])
        end
    end
end

# ## Transports through the sections of a run's config

# GLORYS's cells over a box around one section: the grid, the box and the levels kept (from the surface down)
function glorys_section_grid(lons, lats, date)
    region = BoundingBox(longitude = (minimum(lons) - 1, maximum(lons) + 1), latitude = (minimum(lats) - 1, maximum(lats) + 1))
    grid, levels, _ = glorys_grid(region, date; sample = :u_velocity)
    return grid, region, levels
end

# GLORYS's u, v, T, S of one day on the section's grid; the levels flipped to Oceananigans' order (bottom first). GLORYS
# has every variable at cell centres, so the velocity on a face is the mean of the two cells on either side (a face
# with land on either side is closed and carries nothing).
function glorys_section_fields(grid, region, levels, date)
    level_data(name) = replace(glorys_level_data(name, region, levels, date), NaN => 0.0)
    uᶜ, vᶜ = level_data(:u_velocity), level_data(:v_velocity)
    Nλ, Nφ, Nz = size(uᶜ)

    u, v = XFaceField(grid), YFaceField(grid)
    T, S = CenterField(grid), CenterField(grid)
    uᶠ = zeros(Nλ + 1, Nφ, Nz); uᶠ[2:Nλ, :, :] .= (uᶜ[1:end-1, :, :] .+ uᶜ[2:end, :, :]) ./ 2
    vᶠ = zeros(Nλ, Nφ + 1, Nz); vᶠ[:, 2:Nφ, :] .= (vᶜ[:, 1:end-1, :] .+ vᶜ[:, 2:end, :]) ./ 2
    set!(u, uᶠ); set!(v, vᶠ)
    set!(T, level_data(:temperature)); set!(S, level_data(:salinity))
    foreach(fill_halo_regions!, (u, v, T, S))
    return u, v, T, S
end

"""
    glorys_transports!(days, config)

GLORYS's daily-mean transports through the sections of `config` ([transports.sections]), on GLORYS's own 1/12° grid
with the same code as the model's (transports.jl): for each section, an Oceananigans grid of GLORYS's cells over a box
around it, with GLORYS's levels and sea floor; then section_faces and section_transports give the faces, signs and
formulas the model uses (the same staircase rule, positive into the Arctic, the config's reference temperature and
salinity, NumericalEarth's reference density and heat capacity). GLORYS's daily means do not contain the covariance of
velocity and tracers within the day.

Kept in glorys_transport_dir(config) (validation_config.jl), a folder per set of sections and reference values, one
file per day: "date" and "transports/<section>_<quantity>" (in the units of transport_units). Days already there are
skipped. Downloads, once, small boxes around each section: daily uo, vo, thetao, so, and the static deptho.
"""
function glorys_transports!(days, config)
    dir = glorys_transport_dir(config)
    missing_days = filter(date -> !isfile(joinpath(dir, "$date.jld2")), days)
    isempty(missing_days) && return @info "GLORYS transports: all $(length(days)) days already in $dir"

    values = Dict(date => Dict{String, Float64}() for date in missing_days)
    for (name, lons, lats) in transport_sections(config)
        grid, region, levels = glorys_section_grid(lons, lats, first(missing_days))
        sᵘ, sᵛ = section_faces(grid, lons, lats)
        u_faces, v_faces, deepest = open_faces(grid, sᵘ, sᵛ)
        @info @sprintf("GLORYS, section %s: %d open faces at the surface, deepest %.0f m, %d levels",
                       name, length(u_faces) + length(v_faces), deepest, length(levels))
        for date in missing_days
            u, v, T, S = glorys_section_fields(grid, region, levels, date)
            for (key, field) in section_transports(grid, u, v, T, S, [(name, sᵘ, sᵛ)]; reference_values(config)...)
                compute!(field)
                values[date][key] = field[1, 1, 1]
            end
        end
    end

    for date in missing_days
        jldopen(joinpath(dir, "$date.jld2"), "w") do file
            file["date"] = string(date)
            for (key, value) in values[date]
                file["transports/$key"] = value
            end
        end
    end
    @info "GLORYS transports of $(length(missing_days)) days written to $dir"
end
