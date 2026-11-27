# Heat, salt and freshwater content of the regions of regions.jl, for the model's daily averages
# (analysis/ocean_diagnostics.jl) and for the gridded datasets: it needs only the grid (an ImmersedBoundaryGrid, whose
# cell thicknesses follow the free surface on the model's z-star grid), T and S, so it runs on any grid.
#
# For each region:
#   heat                  ρ₀ cₚ ∫ (T − T_ref) dV                   J
#   salt                  ρ₀ ∫ S dV / 1000                        kg  (S in g/kg)
#   freshwater            ∫ (S_ref − S) / S_ref dV                m³
# each from the surface down to every depth of `depth_bounds` (a cell across a depth counts with its part above it) and
# over the whole water column, and
#   freshwater_isohaline  ∫ (S_ref − S) / S_ref dV from the surface down to the first cell with S ≥ S_ref (Aagaard and
#                         Carmack, 1989): fresher water below saltier is left out, as in hycom-validation
# T_ref, S_ref and the depths are those of validation.toml [ocean_content] (validation_config.jl), the same for every run;
# ρ₀ and cₚ are NumericalEarth's (transports.jl), as for the transports.

@isdefined(region_index) || include("diagnostic_files.jl")
@isdefined(default_reference_density) || include("../transports.jl")   # default_reference_density, default_heat_capacity

"""
    ocean_content(grid, T, S, regions; reference_temperature, reference_salinity, depth_bounds,
                  reference_density = default_reference_density, heat_capacity = default_heat_capacity)

The contents above for one time, on the CPU: `T` and `S` are fields on `grid`, `regions` a
set of `region_masks(grid; set)`. Cells where T or S is not finite count as land. Returns a `NamedTuple` of arrays:
(depth bound, region) for heat, salt and freshwater, (region) for freshwater_isohaline, the regions in the order of
`regions.codes`.
"""
function ocean_content(grid, T, S, regions; reference_temperature, reference_salinity, depth_bounds,
                       reference_density = default_reference_density, heat_capacity = default_heat_capacity)
    R, B = length(regions.codes), length(depth_bounds)
    contents = (heat = zeros(B, R), salt = zeros(B, R), freshwater = zeros(B, R), freshwater_isohaline = zeros(R))
    add_content!(contents, grid, Array(interior(T)), Array(interior(S)), region_index(grid, regions),
                 Float64.(depth_bounds), Float64(reference_temperature), Float64(reference_salinity),
                 Float64(reference_density), Float64(heat_capacity))
    return contents
end

# The sums over every ocean column of every region (a function of its own, so that it compiles for the concrete types)
function add_content!(contents, grid, θ, s, region, depth_bounds, T_ref, S_ref, ρ₀, cₚ)
    (; heat, salt, freshwater, freshwater_isohaline) = contents
    Nx, Ny, Nz = size(grid)
    ccc = (Center(), Center(), Center())

    for j in 1:Ny, i in 1:Nx
        r = region[i, j]
        r == 0 && continue
        A = Azᶜᶜᶜ(i, j, 1, grid)
        top = 0.0              # depth of the top of cell k below the surface
        above_isohaline = true
        for k in Nz:-1:1
            (inactive_node(i, j, k, grid, ccc...) || !isfinite(θ[i, j, k]) || !isfinite(s[i, j, k])) && break
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᵏ, Sᵏ = θ[i, j, k], s[i, j, k]

            for (b, depth) in enumerate(depth_bounds)
                h = clamp(depth - top, zero(Δz), Δz)   # the part of the cell above `depth`
                heat[b, r]       += ρ₀ * cₚ * (Tᵏ - T_ref) * A * h
                salt[b, r]       += ρ₀ * Sᵏ / 1000 * A * h
                freshwater[b, r] += (S_ref - Sᵏ) / S_ref * A * h
            end

            above_isohaline &= Sᵏ < S_ref
            above_isohaline && (freshwater_isohaline[r] += (S_ref - Sᵏ) / S_ref * A * Δz)
            top += Δz
        end
    end
    return nothing
end

# The settings written with the contents
content_file_settings(settings; reference_density = default_reference_density, heat_capacity = default_heat_capacity) =
    (; depth_bounds = settings.depth_bounds,   # m; Inf is the whole water column
       settings.reference_temperature, settings.reference_salinity, reference_density, heat_capacity)
