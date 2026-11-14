# Volume, heat and freshwater transport through the Arctic gateways of Shu et al. (2023, GMD 16, 2539): Fram Strait,
# the Barents Sea Opening, Bering Strait and Davis Strait. Positive into the Arctic. Heat relative to 0 °C (−1.9 °C at
# Bering Strait), freshwater relative to 34.8, as in their Table 3. Included by panarctic_nested.jl.
#
# Each section is a polyline of (longitude, latitude) points, mapped to a staircase of model faces through the cell
# corners nearest the line, so it stays in place when the resolution changes. The setup log lists each section's faces
# and warns when an end is in water. Hourly snapshots go to <output_dir>/gateways.csv; plot_gateways.jl takes daily
# means and compares them with the mooring estimates. glorys_gateways (below) puts GLORYS through the same sections.

using Oceananigans.Grids: inactive_node
using NumericalEarth.NestedModels: child_frame_sources, parent_ocean_sources
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶠᶜ, Δzᶠᶜᶜ

const ρ₀cₚ = 1026 * 3991.9   # J m⁻³ K⁻¹
const Sʳᵉᶠ = 34.8

# (name, points, reference temperature in °C); each line runs so that the Arctic is on its left
const gateway_lines = [
    ("fram",   [(-22.0, 78.83), (13.0, 78.83)],              0.0),   # 78°50′N mooring line, Greenland to Spitsbergen
    ("bso",    [(16.5, 76.6), (19.0, 74.45), (20.6, 70.1)],  0.0),   # Sørkapp, Bjørnøya, Norway
    ("bering", [(-170.3, 66.2), (-167.6, 65.55)],            -1.9),  # Cape Dezhnev to Cape Prince of Wales
    ("davis",  [(-63.0, 66.6), (-52.5, 66.6)],               0.0),   # Baffin Island to Greenland
]

great_circle(λ₁, φ₁, λ₂, φ₂) = 2asind(sqrt(sind((φ₂ - φ₁) / 2)^2 + cosd(φ₁) * cosd(φ₂) * sind((λ₂ - λ₁) / 2)^2))

# Corners visited by the line through `points`, joined by unit steps of i or j
function corner_path(points, corner_λ, corner_φ)
    Nx, Ny = size(corner_λ) .- 1
    function nearest_corner(λ, φ, near)
        candidates = isnothing(near) ? CartesianIndices(corner_λ) :
                     CartesianIndices((max(1, near[1] - 4):min(Nx + 1, near[1] + 4), max(1, near[2] - 4):min(Ny + 1, near[2] + 4)))
        best = argmin(c -> great_circle(λ, φ, corner_λ[c], corner_φ[c]), candidates)
        return (best[1], best[2])
    end
    spacing = minimum(great_circle(corner_λ[i, j], corner_φ[i, j], corner_λ[i+1, j], corner_φ[i+1, j]) for i in 1:Nx, j in (1, Ny ÷ 2, Ny + 1))
    path = [nearest_corner(points[1]..., nothing)]
    for p in 1:length(points) - 1
        (λa, φa), (λb, φb) = points[p], points[p + 1]
        n = ceil(Int, 8 * great_circle(λa, φa, λb, φb) / spacing)
        for s in 1:n
            α = s / n
            c = nearest_corner((1 - α) * λa + α * λb, (1 - α) * φa + α * φb, path[end])
            while c != path[end]
                i, j = path[end]
                push!(path, i != c[1] ? (i + sign(c[1] - i), j) : (i, j + sign(c[2] - j)))
            end
        end
    end
    return path
end

# Faces crossed by the path, each once, as (:u or :v, i, j, sign) with sign × velocity the transport to the left.
#   corner (i, j) → (i + 1, j): v-face (i, j), +1      corner (i, j) → (i, j + 1): u-face (i, j), −1
function path_faces(path)
    signs = Dict{Tuple{Symbol, Int, Int}, Int}()
    for (a, b) in zip(path[1:end-1], path[2:end])
        face = b[2] == a[2] ? (:v, min(a[1], b[1]), a[2]) : (:u, a[1], min(a[2], b[2]))
        signs[face] = get(signs, face, 0) + (b[2] == a[2] ? b[1] - a[1] : a[2] - b[2])
    end
    return [(f..., s) for (f, s) in signs if s != 0]
end

wet_face(kind, i, j, k, g) = kind === :v ? !inactive_node(i, j, k, g, Center(), Face(), Center()) :
                                           !inactive_node(i, j, k, g, Face(), Center(), Center())

# Volume (m³ s⁻¹), heat (W) and freshwater (m³ s⁻¹) through a section, from CPU copies of u, v, T, S. A face where
# the velocity is not finite (land in a parent dataset) carries nothing; where T or S is not, it carries volume only.
function section_transport(faces, Tʳᵉᶠ, u, v, T, S)
    g  = v.grid
    Nz = size(g, 3)
    volume = heat = freshwater = 0.0
    for (kind, i, j, sgn) in faces, k in 1:Nz
        wet_face(kind, i, j, k, g) || continue
        q = kind === :v ? sgn * v[i, j, k] * Δxᶜᶠᶜ(i, j, k, g) * Δzᶜᶠᶜ(i, j, k, g) :
                          sgn * u[i, j, k] * Δyᶠᶜᶜ(i, j, k, g) * Δzᶠᶜᶜ(i, j, k, g)
        isfinite(q) || continue
        i₀, j₀ = kind === :v ? (i, j - 1) : (i - 1, j)            # the cell on the other side of the face
        wet₀ = !inactive_node(i₀, j₀, k, g, Center(), Center(), Center())
        wet₁ = !inactive_node(i, j, k, g, Center(), Center(), Center())
        Tᶠ = (wet₀ * T[i₀, j₀, k] + wet₁ * T[i, j, k]) / (wet₀ + wet₁)
        Sᶠ = (wet₀ * S[i₀, j₀, k] + wet₁ * S[i, j, k]) / (wet₀ + wet₁)
        volume += q
        if isfinite(Tᶠ) && isfinite(Sᶠ)
            heat       += ρ₀cₚ * q * (Tᶠ - Tʳᵉᶠ)
            freshwater += q * (1 - Sᶠ / Sʳᵉᶠ)
        end
    end
    return volume, heat, freshwater
end

# The four sections as staircases of faces on `grid` (a CPU grid)
function gateway_sections(grid)
    Nx, Ny, Nz = size(grid)
    corner_λ = [λnode(i, j, 1, grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]
    corner_φ = [φnode(i, j, 1, grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]
    return map(gateway_lines) do (name, points, Tʳᵉᶠ)
        path  = corner_path(points, corner_λ, corner_φ)
        faces = path_faces(path)
        wet   = count(f -> any(k -> wet_face(f[1], f[2], f[3], k, grid), 1:Nz), faces)
        @info @sprintf("Gateway %-6s %3d faces (%d wet), %.2f°E %.2f°N to %.2f°E %.2f°N", name, length(faces), wet,
                       corner_λ[path[1]...], corner_φ[path[1]...], corner_λ[path[end]...], corner_φ[path[end]...])
        ends = [path_faces(path[1:2]); path_faces(path[end-1:end])]
        any(f -> any(k -> wet_face(f[1], f[2], f[3], k, grid), 1:Nz), ends) &&
            @warn "Gateway $name: an end of the line is in water at this resolution, so water can pass round it"
        (; name, faces, Tʳᵉᶠ)
    end
end

gateway_header(sections) = join(["hour"; [s.name * suffix for s in sections for suffix in ("_volume_Sv", "_heat_TW", "_freshwater_mSv")]], ',')

# One CSV row: the hour, then volume (Sv), heat (TW) and freshwater (mSv) for each section, from fields on any architecture
function gateway_row(sections, hour, u, v, T, S)
    u, v, T, S = map(f -> on_architecture(CPU(), f), (u, v, T, S))
    values = Float64[hour]
    for s in sections
        Q, H, F = section_transport(s.faces, s.Tʳᵉᶠ, u, v, T, S)
        append!(values, (Q / 1e6, H / 1e12, F / 1e3))
    end
    return values
end

write_gateway_row(file, values) = open(io -> println(io, join((@sprintf("%.4f", x) for x in values), ',')), file, "a")

function add_gateway_monitor!(simulation, child, output_dir; schedule = TimeInterval(1hours))
    sections = gateway_sections(on_architecture(CPU(), child.grid))
    file = joinpath(output_dir, "gateways.csv")
    open(io -> println(io, gateway_header(sections)), file, "w")

    function gateway_monitor(sim)
        values = gateway_row(sections, time(sim) / hours, child.velocities.u, child.velocities.v, child.tracers.T, child.tracers.S)
        write_gateway_row(file, values)
        @info "         gateways (Sv): " * join((@sprintf("%s %+.2f", s.name, values[3n - 1]) for (n, s) in enumerate(sections)), ", ")
        return nothing
    end

    add_callback!(simulation, gateway_monitor, schedule; name = :gateway_monitor)
    return file
end

# GLORYS through the same sections: its daily means, interpolated onto `grid` and rotated into the grid's frame as
# the initial state is, sampled at the center of each day of the run (noon) and written to
# <output_dir>/glorys_gateways.csv. Run by glorys_gateways_run.jl.
function glorys_gateways(grid, dataset, dates, output_dir; days, dir = default_download_directory(dataset))
    sections = gateway_sections(on_architecture(CPU(), grid))
    region   = BoundingBox(grid.underlying_grid; padding = default_horizontal_padding(dataset))
    glorys   = PrescribedOcean(region, dates, dataset; architecture = Oceananigans.Architectures.architecture(grid), dir)
    sources  = child_frame_sources(grid, parent_ocean_sources(glorys))

    u, v, T, S = XFaceField(grid), YFaceField(grid), CenterField(grid), CenterField(grid)
    file = joinpath(output_dir, "glorys_gateways.csv")
    open(io -> println(io, gateway_header(sections)), file, "w")
    for n in 0:ceil(Int, days) - 1
        t = (n + 1/2) * Oceananigans.Units.day
        glorys.clock.time = t
        Oceananigans.TimeSteppers.update_state!(glorys)   # loads the days around t into memory
        for (field, name) in ((u, :u), (v, :v), (T, :T), (S, :S))
            sample_parent!(field, boundary_location(name), getproperty(sources, name), t)
        end
        values = gateway_row(sections, t / hours, u, v, T, S)
        write_gateway_row(file, values)
        @info @sprintf("GLORYS day %d: ", n) * join((@sprintf("%s %+.2f", s.name, values[3m - 1]) for (m, s) in enumerate(sections)), ", ") * " Sv"
    end
    return file
end
