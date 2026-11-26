# UDASH, the Unified Database for Arctic and Subarctic Hydrography (Behrendt et al., 2018, ESSD 10, 1119; data:
# Behrendt, Sumata, Rabe and Schauer, 2017, PANGAEA, doi:10.1594/PANGAEA.872931, CC-BY 3.0): temperature and salinity
# profiles 1980–2015 from ships, drifting stations, submarines, floats and moorings, a validation dataset of monthly
# mean profiles by region, as in hycom-validation's UDASH scripts.
#
# The yearly files (ArcticOcean_phys_oce_<year>.txt, space-separated, one row per measurement) are looked for in
# <validation_dir>/UDASH/raw and downloaded from PANGAEA if missing: one zip per year (3–90 MB; 778 MB for all years),
# kept on tape there, so the first request starts loading it and the download waits until it is ready (minutes).
# A copy you already have can be linked into that folder instead.
#
# For each profile (one whitespace-separated file per year, one row per measurement):
#   1. Quality control as hycom-validation: the profile's position flag good (0); temperature and salinity flags 0, 1
#      or 2 (good, averaged duplicate, suspect gradient), dropping inversions, outliers, missing values and profiles that
#      failed the statistical screening (+10).
#   2. In-situ temperature to potential temperature with TEOS-10 (GibbsSeaWater), with S = 34.9 where salinity is missing
#      (XBT, MBT; the error is < 0.01 °C).
#   3. Interpolated onto GLORYS's levels down to 2000 m, a level only where the observations around it are at most
#      max(20 m, 0.25 z) apart; levels above the shallowest observation get its value if it is at most 12 m deep.
#   4. Assigned to a region through the GLORYS cell it lies in (regions.jl on GLORYS's grid and sea floor, over the
#      model domain), so the regions are those of the gridded datasets.
# Then, for each month, the mean of the profiles in each region on each level, with the number of profiles. The mean is
# over profiles, not measurements, so a finely sampled cast does not outweigh a coarse one; it is not area-weighted,
# and only as representative as the profiles are, so read it together with the counts. No heat or freshwater content.
#
# udash_regional_diagnostics! keeps them in <validation_dir>/UDASH/diagnostics/<set>/<year>.jld2. Needs
# validation_config.jl, domain.jl, regions.jl, regional_diagnostics.jl and validation/glorys.jl.

using GibbsSeaWater: gsw_sa_from_sp, gsw_pt0_from_t
using Dates, Printf
using Downloads: Downloads

const udash_good_qc            = (0, 1, 2)
const udash_salinity_fallback  = 34.9
const udash_gap_minimum        = 20.0    # m
const udash_gap_fraction       = 0.25
const udash_surface_maximum    = 12.0    # m
const udash_maximum_depth      = 2000.0  # m

const udash_url = "https://hs.pangaea.de/Projects/UDASH"

# The file of `year` in <validation_dir>/UDASH/raw, downloaded from PANGAEA and unzipped if it is not there
function udash_file(year)
    name = "ArcticOcean_phys_oce_$year"
    file = joinpath(raw_dir("UDASH"), "$name.txt")
    isfile(file) && return file
    zip = joinpath(raw_dir("UDASH"), "$name.zip")
    download_from_tape("$udash_url/$name.zip", zip)
    run(`unzip -o -q $zip -d $(raw_dir("UDASH"))`)
    isfile(file) || error("$zip has no $name.txt")
    rm(zip)
    return file
end

# Download `url` to `path` from a server that answers with an error page (HTTP 503) while it loads the file from tape:
# ask again every `wait` seconds until the file comes
function download_from_tape(url, path; wait = 30, attempts = 120)
    for attempt in 1:attempts
        response = Downloads.request(url; output = path * ".part", throw = false)
        if response isa Downloads.Response && response.status == 200
            mv(path * ".part", path; force = true)
            return path
        end
        attempt == 1 && @info "$(basename(url)) is being loaded from PANGAEA's tape; asking again every $wait s"
        sleep(wait)
    end
    rm(path * ".part"; force = true)
    error("$url did not arrive within $(attempts * wait ÷ 60) minutes")
end

"""
    read_udash_year(path)

The good profiles of one UDASH file: a vector of `(; lon, lat, date, depth, θ, S)`, θ the potential temperature and S
the practical salinity (NaN where not good), on the profile's measured depths.
"""
function read_udash_year(path)
    profiles = NamedTuple{(:lon, :lat, :date, :depth, :θ, :S), Tuple{Float64, Float64, Date, Vector{Float64}, Vector{Float64}, Vector{Float64}}}[]
    id, lon, lat, date, good_position = "", NaN, NaN, Date(0), false
    depth, θ, S = Float64[], Float64[], Float64[]

    function finish!()
        good_position && !isempty(depth) && push!(profiles, (; lon, lat, date, depth = copy(depth), θ = copy(θ), S = copy(S)))
        empty!(depth); empty!(θ); empty!(S)
    end

    for (n, line) in enumerate(eachline(path))
        n == 1 && continue   # header
        f = split(line)
        length(f) ≥ 15 || continue
        if f[1] != id
            finish!()
            id   = f[1]
            lon  = parse(Float64, f[7])
            lat  = parse(Float64, f[8])
            date = Date(f[6][1:10])
            good_position = parse(Int, f[11]) == 0
        end
        p, z = parse(Float64, f[9]), parse(Float64, f[10])
        (good_position && z ≥ 0) || continue
        T  = parse(Float64, f[12]); good_T = T > -900 && parse(Int, f[13]) in udash_good_qc
        Sₚ = parse(Float64, f[14]); good_S = Sₚ > -900 && parse(Int, f[15]) in udash_good_qc
        (good_T || good_S) || continue
        SA = gsw_sa_from_sp(good_S ? Sₚ : udash_salinity_fallback, p, lon, lat)
        push!(depth, z)
        push!(θ, good_T ? gsw_pt0_from_t(SA, T, p) : NaN)
        push!(S, good_S ? Sₚ : NaN)
    end
    finish!()
    return profiles
end

"""
    udash_levels(depth, values, levels)

One variable of one profile on `levels` (m, from the surface down; NaN where not covered), as hycom-validation's
`to_levels`: linear between observations at most max(20 m, 0.25 z) apart, and the shallowest value up to the surface if
it is at most 12 m deep.
"""
function udash_levels(depth, values, levels)
    ok = isfinite.(depth) .& isfinite.(values)
    z, x = depth[ok], values[ok]
    out = fill(NaN, length(levels))
    isempty(z) && return out
    order = sortperm(z)
    z, x = z[order], x[order]
    keep = [true; diff(z) .> 0]   # one value per depth (the first)
    z, x = z[keep], x[keep]
    gap(zᵢ) = max(udash_gap_minimum, udash_gap_fraction * zᵢ)

    if length(z) == 1
        k = argmin(abs.(levels .- z[1]))
        abs(levels[k] - z[1]) ≤ gap(z[1]) / 2 && (out[k] = x[1])
    else
        for (k, zᵢ) in enumerate(levels)
            z[1] ≤ zᵢ ≤ z[end] || continue
            b = clamp(searchsortedfirst(z, zᵢ), 2, length(z))
            if zᵢ == z[b] || zᵢ == z[b - 1]
                out[k] = zᵢ == z[b] ? x[b] : x[b - 1]
            elseif z[b] - z[b - 1] ≤ gap(zᵢ)
                out[k] = x[b - 1] + (zᵢ - z[b - 1]) / (z[b] - z[b - 1]) * (x[b] - x[b - 1])
            end
        end
    end
    z[1] ≤ udash_surface_maximum && (out[levels .< z[1]] .= x[1])
    return out
end

# The regions of each set on GLORYS's cells over the model domain, and a lookup from (λ, φ) to the cell
function udash_region_lookup(set_names)
    region = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))
    static = read_glorys(:depth, GLORYSStatic(), region)
    lon, lat = static.lon, static.lat
    λ = [lon[i] for i in eachindex(lon), j in eachindex(lat)]
    φ = [lat[j] for i in eachindex(lon), j in eachindex(lat)]
    depth = [isfinite(d) && d > 0 ? d : 0.0 for d in static.data]
    depth[.!in_model_domain.(λ, φ)] .= 0
    sets = [region_masks(λ, φ, depth; set, periodic_longitude = true) for set in set_names]

    # Region number of each cell in each set (0: none)
    numbers = map(sets) do regions
        number = zeros(Int, size(depth))
        for (r, code) in enumerate(regions.codes)
            number[regions.masks[code]] .= r
        end
        number
    end
    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    function cell(λᵢ, φᵢ)
        j = round(Int, (φᵢ - lat[1]) / Δφ) + 1
        1 ≤ j ≤ length(lat) || return nothing
        return CartesianIndex(mod(round(Int, (λᵢ - lon[1]) / Δλ), length(lon)) + 1, j)
    end
    return sets, numbers, cell
end

# GLORYS's level depths down to udash_maximum_depth (m, from the surface down), from a small box of one day
function udash_target_levels()
    sample = read_glorys(:temperature, GLORYSDaily(), BoundingBox(longitude = (0, 1), latitude = (75, 76)); date = Date(1993, 1, 1))
    return filter(≤(udash_maximum_depth), sample.depth)
end

"""
    udash_regional_diagnostics!(options, years; settings, set_names)

UDASH's monthly mean profiles by region for each year of `years` and each region set of `set_names`:
<validation_dir>/UDASH/diagnostics/<set>/<year>.jld2 with "temperature", "salinity" (level, region, month), the levels
from the bottom up, "n_temperature", "n_salinity" (the number of profiles), and "dates" the first of each month. Years
already there are skipped.
"""
function udash_regional_diagnostics!(options, years; settings, set_names)
    years = filter(year -> !all(isfile(period_file("UDASH", set, year)) for set in set_names), collect(years))
    isempty(years) && return @info "UDASH: all years already computed"

    levels = udash_target_levels()
    sets, numbers, cell = udash_region_lookup(set_names)
    foreach(regions -> save_dataset_regions("UDASH", regions, nothing; level_depths = reverse(levels), settings), sets)
    L = length(levels)

    for year in years
        profiles = read_udash_year(udash_file(year))
        sums = [(θ = zeros(L, length(r.codes), 12), S = zeros(L, length(r.codes), 12),
                 nθ = zeros(Int, L, length(r.codes), 12), nS = zeros(Int, L, length(r.codes), 12)) for r in sets]
        used = 0
        for profile in profiles
            I = cell(profile.lon, profile.lat)
            isnothing(I) && continue
            θ = udash_levels(profile.depth, profile.θ, levels)
            S = udash_levels(profile.depth, profile.S, levels)
            m = month(profile.date)
            counted = false
            for (s, number) in enumerate(numbers)
                r = number[I]
                r == 0 && continue
                counted = true
                for k in 1:L
                    isfinite(θ[k]) && (sums[s].θ[k, r, m] += θ[k]; sums[s].nθ[k, r, m] += 1)
                    isfinite(S[k]) && (sums[s].S[k, r, m] += S[k]; sums[s].nS[k, r, m] += 1)
                end
            end
            used += counted
        end
        @info @sprintf("UDASH %d: %d profiles, %d in a region", year, length(profiles), used)

        for (s, regions) in enumerate(sets)
            (; θ, S, nθ, nS) = sums[s]
            mean(total, count) = ifelse.(count .> 0, total ./ max.(count, 1), NaN)
            results = [(; temperature = mean(θ[end:-1:1, :, m], nθ[end:-1:1, :, m]),   # levels from the bottom up
                          salinity    = mean(S[end:-1:1, :, m], nS[end:-1:1, :, m]),
                          n_temperature = nθ[end:-1:1, :, m], n_salinity = nS[end:-1:1, :, m]) for m in 1:12]
            save_dataset_period("UDASH", regions.set, year, results; dates = [Date(year, m, 1) for m in 1:12])
        end
    end
end
