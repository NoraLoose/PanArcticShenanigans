# GLORYS's daily-mean transports through the sections of [transports], computed on GLORYS's own 1/12° grid with the
# same code as the model's (transports.jl): the reference for the model's arctic_transports.jld2.
#
#   julia +1.13 --project=. analysis/glorys_transports.jl [config.toml]
#
# For each section, an Oceananigans latitude-longitude grid whose cells are GLORYS's cells is built over a box around
# the section, with GLORYS's levels and sea floor (deptho, as partial bottom cells). GLORYS has every variable at cell
# centres, so the velocity on a face is the mean of the two cells on either side (a face with land on either side is
# closed and carries nothing). Then section_faces and section_transports of transports.jl give the faces, signs and
# formulas the model uses: the same staircase rule, positive into the Arctic, the same reference temperature and
# salinity ([transports]), and NumericalEarth's reference density and heat capacity (as the model's online transports).
# GLORYS's daily means do not contain the covariance of velocity and tracers within the day.
#
# Downloads, once, small boxes around each section: daily uo, vo, thetao, so for every day from start_date to
# end_date, and the static deptho. Needs no GPU and no model run.
#
# Writes glorys_transports.jld2 in output_dir: "dates" (the days, as strings) and "transports/<section>_<quantity>"
# (one value per day, in the units of transport_units), read by plotting/gateway_transports.jl.

include("../config.jl")
include("../transports.jl")

using NumericalEarth, Dates, JLD2
using CopernicusMarine   # enables the GLORYS download extension
using Oceananigans.BoundaryConditions: fill_halo_regions!
using NumericalEarth.DataWrangling: BoundingBox, Metadatum, NCDatasets, Downloads, default_download_directory
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily, GLORYSStatic, GLORYS_dataset_variable_names
using .NCDatasets: Dataset, dimnames

glorys_dir = isempty(data_dir) ? default_download_directory(GLORYSDaily()) : data_dir
start_date = config["time"]["start_date"]
end_date   = config["time"]["end_date"]
days       = collect(start_date:Day(1):end_date - Day(1))   # GLORYS's daily means over the run

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

# GLORYS's cells over a box around one section: the grid, and the levels kept (from the surface down)
function glorys_section_grid(lons, lats)
    region = BoundingBox(longitude = (minimum(lons) - 1, maximum(lons) + 1), latitude = (minimum(lats) - 1, maximum(lats) + 1))
    static = read_glorys(:depth, GLORYSStatic(), region)
    sample = read_glorys(:u_velocity, GLORYSDaily(), region; date = first(days))
    (static.lon ≈ sample.lon && static.lat ≈ sample.lat) || error("GLORYS's deptho and uo are on different grids")

    lon, lat = static.lon, static.lat
    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    deptho = [isfinite(d) && d > 0 ? d : 0.0 for d in static.data]
    zf     = level_interfaces(sample.depth)
    levels = 1:findlast(<(maximum(deptho)), zf[1:end-1])   # every level with ocean somewhere in the box

    underlying = LatitudeLongitudeGrid(CPU(); size = (length(lon), length(lat), length(levels)),
                                       longitude = (lon[1] - Δλ / 2, lon[end] + Δλ / 2),
                                       latitude  = (lat[1] - Δφ / 2, lat[end] + Δφ / 2),
                                       z = reverse(-zf[1:last(levels) + 1]),
                                       topology = (Bounded, Bounded, Bounded))
    grid = ImmersedBoundaryGrid(underlying, PartialCellBottom(-deptho))
    return grid, region, levels
end

# GLORYS's u, v, T, S of one day on the section's grid; the levels flipped to Oceananigans' order (bottom first)
function glorys_fields(grid, region, levels, date)
    level_data(name) = reverse(replace(read_glorys(name, GLORYSDaily(), region; date).data[:, :, levels], NaN => 0.0), dims = 3)
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

results = Dict{String, Vector{Float64}}()

for (name, lons, lats) in transport_sections(config)
    grid, region, levels = glorys_section_grid(lons, lats)
    sᵘ, sᵛ = section_faces(grid, lons, lats)
    u_faces, v_faces, deepest = open_faces(grid, sᵘ, sᵛ)
    @info @sprintf("GLORYS, section %s: %d open faces at the surface, deepest %.0f m, %d levels",
                   name, length(u_faces) + length(v_faces), deepest, length(levels))

    for date in days
        u, v, T, S = glorys_fields(grid, region, levels, date)
        transports = section_transports(grid, u, v, T, S, [(name, sᵘ, sᵛ)]; reference_values(config)...)
        for (key, field) in transports
            compute!(field)
            push!(get!(results, key, Float64[]), field[1, 1, 1])
        end
    end

    @info @sprintf("    volume (Sv): %s", join([@sprintf("%.2f", Q) for Q in results["$(name)_volume"]], ", "))
end

file = joinpath(output_dir, "glorys_transports.jld2")
jldopen(file, "w") do f
    f["dates"] = string.(days)
    f["units"] = transport_units
    for (key, values) in results
        f["transports/$key"] = values
    end
end
@info "Wrote $file"
