# GLORYS's net transport through the model's two open edges, computed on GLORYS's own grid (its velocities, levels and
# land mask, before any interpolation), next to what the boundary data carry after interpolation onto the model grid
# (the `*_glorys_Sv` columns of a transport monitor). A difference points to the interpolation, the inpainting or the
# bathymetry rather than to GLORYS itself.
#
#   julia +1.13 --project=. plotting/glorys_section_transport.jl [config.toml] [monitor csv]
#
# Reads the GLORYS daily files already downloaded for the run (u and v for each date). Transports in Sv; Atlantic
# positive into the domain, Pacific positive out of it (as in the monitor).

include("../setup.jl")

using Printf
using NumericalEarth.DataWrangling: NCDatasets
using .NCDatasets: Dataset, dimnames
using Oceananigans.Grids: λnode, φnode
using NumericalEarth.DataWrangling: metadata_path

monitor_file = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(output_dir, "transport_monitor.csv")

cpu_grid = Oceananigans.Architectures.on_architecture(CPU(), underlying_grid)
R = cpu_grid.radius
edge_corners(j) = [(λnode(i, j, 1, cpu_grid, Face(), Face(), Center()), φnode(i, j, 1, cpu_grid, Face(), Face(), Center()))
                   for i in 1:Nx+1]
atlantic_corners = edge_corners(1)
pacific_corners  = edge_corners(Ny + 1)

all_φ = [last.(atlantic_corners); last.(pacific_corners)]
φ_band = (minimum(all_φ) - 1, maximum(all_φ) + 1)
@info @sprintf("Edges between %.1f°N and %.1f°N", extrema(all_φ)...)

wrap(Δλ) = mod(Δλ + 180, 360) - 180

# Depth-integrated velocity (m² s⁻¹) on GLORYS's grid, over the latitude band of the edges; land and below-bottom
# cells (missing in the file) count as zero.
function native_depth_integral(name, date)
    path = metadata_path(Metadatum(name; date, dataset, region, dir_kw...))
    isfile(path) || error("$path not found; run panarctic.jl for these dates first")
    ds  = Dataset(path)
    var = ds[name === :u_velocity ? "uo" : "vo"]
    dimnames(var)[1:3] == ("longitude", "latitude", "depth") || error("unexpected dimensions $(dimnames(var))")
    lon = Float64.(ds["longitude"][:])
    lat = Float64.(ds["latitude"][:])
    z   = Float64.(ds["depth"][:])
    band = findall(φ -> φ_band[1] ≤ φ ≤ φ_band[2], lat)
    u = coalesce.(var[:, band, :, 1], 0f0)
    close(ds)
    interfaces = vcat(0, (z[1:end-1] .+ z[2:end]) ./ 2, z[end] + (z[end] - z[end-1]) / 2)
    Δz = diff(interfaces)
    U = dropdims(sum(Float64.(u) .* reshape(Δz, 1, 1, :), dims = 3), dims = 3)
    return lon, lat[band], U
end

# Bilinear sample of a (lon, lat) field on a regular grid, periodic in longitude
function sample(U, lon, lat, λ, φ)
    Nλ = length(lon)
    fi = mod(λ - lon[1], 360) / (lon[2] - lon[1])
    fj = (φ - lat[1]) / (lat[2] - lat[1])
    i₀ = floor(Int, fi); j₀ = clamp(floor(Int, fj), 0, length(lat) - 2)
    a = fi - i₀; b = clamp(fj - j₀, 0, 1)
    i₁ = mod(i₀, Nλ) + 1; i₂ = mod(i₀ + 1, Nλ) + 1
    j₁ = j₀ + 1; j₂ = j₀ + 2
    return (1 - a) * (1 - b) * U[i₁, j₁] + a * (1 - b) * U[i₂, j₁] + (1 - a) * b * U[i₁, j₂] + a * b * U[i₂, j₂]
end

# Transport across a polyline of corners, positive towards increasing j (left of the direction of increasing i)
function section_transport(corners, lon, lat, Uᵉ, Uⁿ)
    T = 0.0
    for k in 1:length(corners) - 1
        (λa, φa), (λb, φb) = corners[k], corners[k + 1]
        Δλ = wrap(λb - λa)
        λm, φm = λa + Δλ / 2, (φa + φb) / 2
        dx = R * cosd(φm) * deg2rad(Δλ)
        dy = R * deg2rad(φb - φa)
        T += sample(Uᵉ, lon, lat, λm, φm) * (-dy) + sample(Uⁿ, lon, lat, λm, φm) * dx
    end
    return T
end

function read_monitor(path)
    isfile(path) || return nothing
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(l)]
    return Dict(String(h) => [r[k] for r in rows] for (k, h) in enumerate(header))
end
monitor = read_monitor(monitor_file)

# GLORYS daily means are centred at noon; the run starts at noon on dates[start_index], so snapshot n is at monitor
# hour 24 (n - start_index).
function monitor_value(column, n)
    isnothing(monitor) && return NaN
    k = findfirst(h -> h ≈ 24 * (n - start_index), monitor["hour"])
    return isnothing(k) ? NaN : monitor[column][k]
end

@printf("\n%-12s %26s %26s %26s\n", "", "Atlantic (into domain)", "Pacific (out of domain)", "net inflow")
@printf("%-12s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", "date", "native", "model", "diff",
        "native", "model", "diff", "native", "model", "diff")
for (n, date) in enumerate(dates)
    lon, lat, Uᵉ = native_depth_integral(:u_velocity, date)
    _,   _,   Uⁿ = native_depth_integral(:v_velocity, date)
    A = section_transport(atlantic_corners, lon, lat, Uᵉ, Uⁿ) / 1e6
    P = section_transport(pacific_corners,  lon, lat, Uᵉ, Uⁿ) / 1e6
    Am, Pm = monitor_value("atlantic_glorys_Sv", n), monitor_value("pacific_glorys_Sv", n)
    @printf("%-12s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", Date(date), A, Am, Am - A,
            P, Pm, Pm - P, A - P, Am - Pm, (Am - Pm) - (A - P))
end
println("\nnative: GLORYS's own grid; model: the boundary data on the model grid (from $monitor_file, NaN if absent)")
