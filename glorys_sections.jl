# GLORYS's net transport through the model's open edges on GLORYS's own grid: its velocities, levels and land mask,
# before any interpolation onto the model grid. Needs setup.jl (grid, dates, dataset, region, dir_kw).
#
# Each edge is the polyline through the corners of its boundary faces; the transport is positive towards increasing j
# (Atlantic: into the domain; Pacific: out of it).

using Oceananigans.Grids: λnode, φnode
using NumericalEarth.DataWrangling: NCDatasets, metadata_path
using .NCDatasets: Dataset, dimnames

section_grid = Oceananigans.Architectures.on_architecture(CPU(), underlying_grid)
edge_corners(j) = [(λnode(i, j, 1, section_grid, Face(), Face(), Center()), φnode(i, j, 1, section_grid, Face(), Face(), Center()))
                   for i in 1:Nx+1]
atlantic_corners = edge_corners(1)
pacific_corners  = edge_corners(Ny + 1)

edge_latitudes = [last.(atlantic_corners); last.(pacific_corners)]
section_band = (minimum(edge_latitudes) - 1, maximum(edge_latitudes) + 1)

wrap_longitude(Δλ) = mod(Δλ + 180, 360) - 180

# Depth-integrated velocity (m² s⁻¹) on GLORYS's grid over the latitude band of the edges; missing cells count as zero
function native_depth_integral(name, date)
    path = metadata_path(Metadatum(name; date, dataset, region, dir_kw...))
    isfile(path) || error("$path not found; the GLORYS files for these dates must be downloaded first")
    ds  = Dataset(path)
    var = ds[name === :u_velocity ? "uo" : "vo"]
    dimnames(var)[1:3] == ("longitude", "latitude", "depth") || error("unexpected dimensions $(dimnames(var))")
    lon = Float64.(ds["longitude"][:])
    lat = Float64.(ds["latitude"][:])
    z   = Float64.(ds["depth"][:])
    band = findall(φ -> section_band[1] ≤ φ ≤ section_band[2], lat)
    u = coalesce.(var[:, band, :, 1], 0f0)
    close(ds)
    interfaces = vcat(0, (z[1:end-1] .+ z[2:end]) ./ 2, z[end] + (z[end] - z[end-1]) / 2)
    U = dropdims(sum(Float64.(u) .* reshape(diff(interfaces), 1, 1, :), dims = 3), dims = 3)
    return lon, lat[band], U
end

# Bilinear sample of a (lon, lat) field on a regular grid, periodic in longitude
function sample_lonlat(U, lon, lat, λ, φ)
    Nλ = length(lon)
    fi = mod(λ - lon[1], 360) / (lon[2] - lon[1])
    fj = (φ - lat[1]) / (lat[2] - lat[1])
    i₀ = floor(Int, fi); j₀ = clamp(floor(Int, fj), 0, length(lat) - 2)
    a = fi - i₀; b = clamp(fj - j₀, 0, 1)
    i₁ = mod(i₀, Nλ) + 1; i₂ = mod(i₀ + 1, Nλ) + 1
    j₁ = j₀ + 1; j₂ = j₀ + 2
    return (1 - a) * (1 - b) * U[i₁, j₁] + a * (1 - b) * U[i₂, j₁] + (1 - a) * b * U[i₁, j₂] + a * b * U[i₂, j₂]
end

# Transport (m³ s⁻¹) across a polyline of corners, positive to the left of the direction of increasing i
function section_transport(corners, lon, lat, Uᵉ, Uⁿ)
    R = section_grid.radius
    T = 0.0
    for k in 1:length(corners) - 1
        (λa, φa), (λb, φb) = corners[k], corners[k + 1]
        Δλ = wrap_longitude(λb - λa)
        λm, φm = λa + Δλ / 2, (φa + φb) / 2
        dx = R * cosd(φm) * deg2rad(Δλ)
        dy = R * deg2rad(φb - φa)
        T += sample_lonlat(Uᵉ, lon, lat, λm, φm) * (-dy) + sample_lonlat(Uⁿ, lon, lat, λm, φm) * dx
    end
    return T
end

# (Atlantic into the domain, Pacific out of it), m³ s⁻¹, for each date
function native_edge_transports(dates)
    A = zeros(length(dates)); P = zeros(length(dates))
    for (n, date) in enumerate(dates)
        lon, lat, Uᵉ = native_depth_integral(:u_velocity, date)
        _,   _,   Uⁿ = native_depth_integral(:v_velocity, date)
        A[n] = section_transport(atlantic_corners, lon, lat, Uᵉ, Uⁿ)
        P[n] = section_transport(pacific_corners,  lon, lat, Uᵉ, Uⁿ)
    end
    return A, P
end
