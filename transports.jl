# Volume, heat, salt and freshwater transports through sections (by default Fram Strait, the Barents Sea Opening,
# Bering Strait and Davis Strait; set in [transports] in the config) and through the two open boundaries. Included by
# panarctic.jl (online) and analysis/transports.jl (offline), so both use the same faces and the same formulas.
#
# A section is a list of (longitude, latitude) waypoints joined by great circles, so its definition does not depend on
# the grid. `section_faces` finds the u and v faces of the model grid that the section crosses: each face gets a sign,
# +1 where positive u (or v) flows to the right of the section, walking from the first waypoint to the last, and -1
# where it flows to the left. A transport is
#
#     Σ over faces and levels of  sign × velocity × face area × (α + β c),
#
# with c the tracer averaged from the two cells on either side of the face. The face areas are the model's own (partial
# bottom cells, z⋆ stretching), and only faces between two wet cells count, so the volume transport is the one the model
# uses in its continuity equation. The open boundaries are sections too: the rows of v faces at the domain's southern
# (Atlantic) and northern (Pacific) edges, positive into the domain.
#
# Online, the transports are computed every time step and averaged like the other output, so they include the
# covariance of velocity and tracers within the averaging window. Offline, from the saved averages of u, v, T and S,
# they miss that part (the average of u T is not the average of u times the average of T).

using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: λnode, φnode, inactive_node, immersed_peripheral_node, immersed_inactive_node
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Printf

# Density and heat capacity that convert volume fluxes to salt and heat transports. The defaults are NumericalEarth's for
# `ocean_simulation` (1020 kg m⁻³, TEOS-10 heat capacity); panarctic.jl passes the model's own values.
const default_reference_density = 1020.0                 # kg m⁻³
const default_heat_capacity     = 3991.86795711963       # J kg⁻¹ K⁻¹

transport_settings(config) = get(config, "transports", Dict{String, Any}())

# Sections from the config, sorted by name: [(name, longitudes, latitudes), ...]
function transport_sections(config)
    sections = get(transport_settings(config), "sections", Dict{String, Any}())
    return [(name, Float64.(s["lon"]), Float64.(s["lat"])) for (name, s) in sort(collect(sections), by = first)]
end

#####
##### Faces crossed by a section
#####

unit_vector(λ, φ) = (cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ))
cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
dot3(a, b)   = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

"""
    section_faces(grid, lons, lats)

Signs of the u and v faces crossed by the section through the waypoints `lons`, `lats`: integer arrays `sᵘ`
(`(Nx+1) × Ny`, like u) and `sᵛ` (`Nx × (Ny+1)`, like v), with `sᵘ[i, j]` the sign of the u face between cells
`(i-1, j)` and `(i, j)` (0 if the section does not cross it), and likewise `sᵛ[i, j]` for the v face between `(i, j-1)`
and `(i, j)`.

For each segment, the cells to the right of its great circle and between its two ends are marked. The faces between a
marked cell and one to the left of the great circle are the faces the segment crosses, with sign +1 where the positive
velocity points into the marked cell. Cells on the edge of the domain are never marked, so the domain edges are not
section faces.
"""
function section_faces(grid, lons, lats)
    grid = on_architecture(CPU(), grid)
    Nx, Ny = size(grid, 1), size(grid, 2)
    ccc = (Center(), Center(), Center())
    x = [unit_vector(λnode(i, j, 1, grid, ccc...), φnode(i, j, 1, grid, ccc...)) for i in 1:Nx, j in 1:Ny]

    sᵘ = zeros(Int, Nx + 1, Ny)
    sᵛ = zeros(Int, Nx, Ny + 1)

    for n in 1:length(lons)-1
        a = unit_vector(lons[n],   lats[n])
        b = unit_vector(lons[n+1], lats[n+1])
        normal = cross3(a, b)
        side = dot3.(x, Ref(normal))   # < 0: right of the great circle through a and b

        # Between the two ends: in the wedge spanned by a and b, on their side of the globe (not the antipodal wedge)
        between(p) = dot3(cross3(a, p), cross3(p, b)) ≥ 0 && dot3(p, a .+ b) > 0
        right = (side .< 0) .& between.(x)
        right[[1, Nx], :] .= false
        right[:, [1, Ny]] .= false

        # A face is on the section if its unmarked cell is left of the great circle, not just beyond the segment's ends.
        for j in 1:Ny, i in 2:Nx
            f = right[i, j] - right[i-1, j]
            f != 0 && side[f > 0 ? i-1 : i, j] ≥ 0 && (sᵘ[i, j] += f)
        end

        for j in 2:Ny, i in 1:Nx
            f = right[i, j] - right[i, j-1]
            f != 0 && side[i, f > 0 ? j-1 : j] ≥ 0 && (sᵛ[i, j] += f)
        end
    end

    return sᵘ, sᵛ
end

# A face is open if it is not on or below the bathymetry. Only the immersed boundary counts, so the faces on the open
# boundaries are open too (`peripheral_node` would also close every face on the domain's edges).
@inline open_face(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

const fcc = (Face(), Center(), Center())
const cfc = (Center(), Face(), Center())

# Faces of a section that are open at the surface, and the depth of the deepest one
function open_faces(grid, sᵘ, sᵛ)
    grid = on_architecture(CPU(), grid)
    Nz = size(grid, 3)
    u_faces = [(i, j) for (i, j) in Tuple.(findall(!=(0), sᵘ)) if open_face(i, j, Nz, grid, fcc...)]
    v_faces = [(i, j) for (i, j) in Tuple.(findall(!=(0), sᵛ)) if open_face(i, j, Nz, grid, cfc...)]
    depth(i, j, Δz, ℓ) = sum(Δz(i, j, k, grid) for k in 1:Nz if open_face(i, j, k, grid, ℓ...))
    depths = vcat([depth(i, j, Δzᶠᶜᶜ, fcc) for (i, j) in u_faces], [depth(i, j, Δzᶜᶠᶜ, cfc) for (i, j) in v_faces])
    return u_faces, v_faces, maximum(depths, init = 0.0)
end

"""
    gateway_sections(grid, config)

Faces of every section in the config on `grid`, as `[(name, sᵘ, sᵛ), ...]`. Logs each section's number of open faces
and deepest face, and warns when a section end is in the ocean (the section may then miss part of the strait).
"""
function gateway_sections(grid, config)
    cpu_grid = on_architecture(CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    ccc = (Center(), Center(), Center())
    x = [unit_vector(λnode(i, j, 1, cpu_grid, ccc...), φnode(i, j, 1, cpu_grid, ccc...)) for i in 1:Nx, j in 1:Ny]

    sections = []
    for (name, lons, lats) in transport_sections(config)
        sᵘ, sᵛ = section_faces(cpu_grid, lons, lats)
        u_faces, v_faces, deepest = open_faces(cpu_grid, sᵘ, sᵛ)
        isempty(u_faces) && isempty(v_faces) && error("Section $name crosses no open faces of the grid: check its waypoints.")

        for (λ, φ) in ((first(lons), first(lats)), (last(lons), last(lats)))
            i, j = Tuple(argmax(dot3.(x, Ref(unit_vector(λ, φ)))))   # nearest cell centre
            inactive_node(i, j, Nz, cpu_grid, ccc...) ||
                @warn "Section $name ends in the ocean at ($λ°E, $φ°N), cell ($i, $j): move the waypoint onto land."
        end

        @info @sprintf("Section %s: %d open faces at the surface (%d u, %d v), deepest %.0f m",
                       name, length(u_faces) + length(v_faces), length(u_faces), length(v_faces), deepest)
        push!(sections, (name, sᵘ, sᵛ))
    end

    return sections
end

"""
    boundary_sections(grid)

The open boundaries of panarctic.jl as sections, `[(name, sᵘ, sᵛ), ...]`: the v faces at `j = 1` (Atlantic) and
`j = Ny + 1` (Pacific), signed so that positive transport is into the domain. Keep in step with the open boundary
conditions in panarctic.jl.
"""
function boundary_sections(grid)
    Nx, Ny = size(grid, 1), size(grid, 2)
    no_u_faces = zeros(Int, Nx + 1, Ny)
    south = zeros(Int, Nx, Ny + 1)
    north = zeros(Int, Nx, Ny + 1)
    south[:, 1]      .= +1   # v > 0 flows into the domain
    north[:, Ny + 1] .= -1   # v < 0 flows into the domain
    return [("atlantic_boundary", no_u_faces, south), ("pacific_boundary", copy(no_u_faces), north)]
end

#####
##### Transports
#####

@inline tracerᶠᶜᶜ(i, j, k, grid, ::Nothing) = zero(grid)
@inline tracerᶜᶠᶜ(i, j, k, grid, ::Nothing) = zero(grid)
@inline tracerᶠᶜᶜ(i, j, k, grid, c) = ℑxᶠᵃᵃ(i, j, k, grid, c)
@inline tracerᶜᶠᶜ(i, j, k, grid, c) = ℑyᵃᶠᵃ(i, j, k, grid, c)

# Flux through u face (i, j, k) and v face (i, j, k): zero unless the face is on the section and open
@inline function u_flux(i, j, k, grid, sᵘ, u, c, α, β)
    s = @inbounds sᵘ[i, j]
    F = @inbounds s * u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid) * (α + β * tracerᶠᶜᶜ(i, j, k, grid, c))
    return ifelse((s != 0) & open_face(i, j, k, grid, fcc...), F, zero(F))
end

@inline function v_flux(i, j, k, grid, sᵛ, v, c, α, β)
    s = @inbounds sᵛ[i, j]
    F = @inbounds s * v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) * (α + β * tracerᶜᶠᶜ(i, j, k, grid, c))
    return ifelse((s != 0) & open_face(i, j, k, grid, cfc...), F, zero(F))
end

# Flux through the faces of cell (i, j, k) on the section: its west and south faces, and in the last column and row
# also its east and north faces, so that every face, including those on the domain's edges, is counted once.
@inline function section_flux(i, j, k, grid, sᵘ, sᵛ, u, v, c, α, β)
    F = u_flux(i, j, k, grid, sᵘ, u, c, α, β) + v_flux(i, j, k, grid, sᵛ, v, c, α, β)
    F += ifelse(i == grid.Nx, u_flux(i + 1, j, k, grid, sᵘ, u, c, α, β), zero(F))
    F += ifelse(j == grid.Ny, v_flux(i, j + 1, k, grid, sᵛ, v, c, α, β), zero(F))
    return F
end

"""
    section_transports(grid, u, v, T, S, sections; reference_temperature = 0, reference_salinity = 34.8,
                       reference_density = $default_reference_density, heat_capacity = $default_heat_capacity)

Volume, heat, salt and freshwater transports through each section of `sections` (from `gateway_sections` or
`boundary_sections`), as a
`Dict` of scalar fields named `"<section>_volume"` etc., computed from `u`, `v`, `T`, `S`.
Units: volume and freshwater in Sv, heat in TW, salt in kg s⁻¹. Heat is relative to `reference_temperature`,
freshwater relative to `reference_salinity`: (S_ref - S) / S_ref times the volume flux.
"""
function section_transports(grid, u, v, T, S, sections;
                            reference_temperature = 0, reference_salinity = 34.8,
                            reference_density = default_reference_density, heat_capacity = default_heat_capacity)
    FT = eltype(grid)
    ρc = reference_density * heat_capacity
    Tᵣ, Sᵣ = reference_temperature, reference_salinity

    # (tracer, α, β): flux density = velocity × (α + β × tracer)
    quantities = (volume     = (nothing, 1e-6,             0.0),                       # Sv
                  heat       = (T,       -ρc * Tᵣ * 1e-12, ρc * 1e-12),                # TW
                  salt       = (S,       0.0,              reference_density / 1000),  # kg s⁻¹ (S in g kg⁻¹)
                  freshwater = (S,       1e-6,             -1e-6 / Sᵣ))                # Sv

    arch = architecture(grid)
    Nx, Ny = size(grid, 1), size(grid, 2)
    transports = Dict{String, Any}()

    for (name, sᵘ, sᵛ) in sections
        # Sum over the box of cells around the section's faces (a face on the east or north edge belongs to the last cell)
        faces = vcat(findall(!=(0), sᵘ), findall(!=(0), sᵛ))
        i₁, i₂ = extrema(min(I[1], Nx) for I in faces)
        j₁, j₂ = extrema(min(I[2], Ny) for I in faces)
        signs = map(s -> on_architecture(arch, FT.(s)), (sᵘ, sᵛ))

        for (quantity, (c, α, β)) in pairs(quantities)
            flux = KernelFunctionOperation{Center, Center, Center}(section_flux, grid, signs..., u, v, c, FT(α), FT(β))
            flux = Field(flux; indices = (i₁:i₂, j₁:j₂, :))
            transports["$(name)_$(quantity)"] = Field(Reduction(sum!, flux, dims = (1, 2, 3)))
        end
    end

    return transports
end

transport_quantities = ("volume", "heat", "salt", "freshwater")
transport_units      = Dict("volume" => "Sv", "heat" => "TW", "salt" => "kg s⁻¹", "freshwater" => "Sv")

# Reference temperature and salinity from the config, as keyword arguments for `section_transports`
function reference_values(config)
    settings = transport_settings(config)
    return (reference_temperature = Float64(get(settings, "reference_temperature", 0)),
            reference_salinity    = Float64(get(settings, "reference_salinity", 34.8)))
end

# Record how the transports were computed in a JLD2 file: units, reference values, and per section the waypoints and
# the signed faces, so the transports can be traced back to the grid they were computed on.
function save_transport_metadata!(file, config, sections; reference_density, heat_capacity)
    file["metadata/units"] = transport_units
    file["metadata/reference_density"] = reference_density
    file["metadata/heat_capacity"]     = heat_capacity
    for (key, value) in pairs(reference_values(config))
        file["metadata/$key"] = value
    end
    waypoints = Dict(name => (lons, lats) for (name, lons, lats) in transport_sections(config))
    # Faces as rows (i, j, sign); the full sign arrays would be grid-sized
    face_table(s) = Int32[f(I) for I in findall(!=(0), s), f in (I -> I[1], I -> I[2], I -> s[I])]
    for (name, sᵘ, sᵛ) in sections
        haskey(waypoints, name) && ((file["sections/$name/lon"], file["sections/$name/lat"]) = waypoints[name])
        file["sections/$name/u_faces"] = face_table(sᵘ)
        file["sections/$name/v_faces"] = face_table(sᵛ)
    end
    return nothing
end
