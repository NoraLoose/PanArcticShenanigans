# Heat, salt and freshwater content and mean temperature and salinity profiles of the regions of regions.jl, for the
# model's daily averages (analysis/regional_diagnostics.jl): it needs only the grid (an ImmersedBoundaryGrid, whose cell
# thicknesses follow the free surface on the model's z-star grid), T and S, so it runs on any grid.
#
# For each region:
#   heat                  ρ₀ cₚ ∫ (T − T_ref) dV                   J
#   salt                  ρ₀ ∫ S dV / 1000                        kg  (S in g/kg)
#   freshwater            ∫ (S_ref − S) / S_ref dV                m³
# each from the surface down to every depth of `depth_bounds` (a cell across a depth counts with its part above it) and
# over the whole water column, and
#   freshwater_isohaline  ∫ (S_ref − S) / S_ref dV from the surface down to the first cell with S ≥ S_ref (Aagaard and
#                         Carmack, 1989): fresher water below saltier is left out, as in hycom-validation
#   temperature, salinity volume-weighted mean on each level (°C, g/kg)
#   volume                volume of each level (m³)
# T_ref, S_ref, the depths and the region sets are those of validation.toml [regional_diagnostics] (validation_config.jl),
# the same for every run; ρ₀ and cₚ are NumericalEarth's (transports.jl), as for the transports.

include("regions.jl")
include("transports.jl")   # default_reference_density, default_heat_capacity

using Oceananigans.Grids: inactive_node
using Oceananigans.Operators: Δzᶜᶜᶜ, Azᶜᶜᶜ
using JLD2

"""
    regional_diagnostics(grid, T, S, regions; reference_temperature, reference_salinity, depth_bounds,
                         reference_density = default_reference_density, heat_capacity = default_heat_capacity)

The diagnostics above for one time, on the CPU: `T` and `S` are fields (or arrays of the interior) on `grid`, `regions`
a set of `region_masks(grid; set)`. Cells where T or S is not finite count as land. Returns a `NamedTuple` of arrays:
(depth bound, region) for heat, salt and freshwater, (region) for freshwater_isohaline, (level, region) for
temperature, salinity and volume, the regions in the order of `regions.codes` and the levels from the bottom up.
"""
function regional_diagnostics(grid, T, S, regions; reference_temperature, reference_salinity, depth_bounds,
                              reference_density = default_reference_density, heat_capacity = default_heat_capacity)
    Nx, Ny, Nz = size(grid)
    region = zeros(Int, Nx, Ny)
    for (r, code) in enumerate(regions.codes)
        region[regions.masks[code]] .= r
    end

    R, B = length(regions.codes), length(depth_bounds)
    sums = (heat = zeros(B, R), salt = zeros(B, R), freshwater = zeros(B, R), freshwater_isohaline = zeros(R),
            ∑TV = zeros(Nz, R), ∑SV = zeros(Nz, R), volume = zeros(Nz, R))
    add_columns!(sums, grid, Array(interior(T)), Array(interior(S)), region, Float64.(depth_bounds),
                 Float64(reference_temperature), Float64(reference_salinity), Float64(reference_density),
                 Float64(heat_capacity))

    (; heat, salt, freshwater, freshwater_isohaline, ∑TV, ∑SV, volume) = sums
    temperature = ifelse.(volume .> 0, ∑TV ./ volume, NaN)
    salinity    = ifelse.(volume .> 0, ∑SV ./ volume, NaN)
    return (; heat, salt, freshwater, freshwater_isohaline, temperature, salinity, volume)
end

# The sums over every ocean column of every region (a function of its own, so that it compiles for the concrete types)
function add_columns!(sums, grid, θ, s, region, depth_bounds, T_ref, S_ref, ρ₀, cₚ)
    (; heat, salt, freshwater, freshwater_isohaline, ∑TV, ∑SV, volume) = sums
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
            V = A * Δz

            ∑TV[k, r]    += Tᵏ * V
            ∑SV[k, r]    += Sᵏ * V
            volume[k, r] += V

            for (b, depth) in enumerate(depth_bounds)
                h = clamp(depth - top, zero(Δz), Δz)   # the part of the cell above `depth`
                heat[b, r]       += ρ₀ * cₚ * (Tᵏ - T_ref) * A * h
                salt[b, r]       += ρ₀ * Sᵏ / 1000 * A * h
                freshwater[b, r] += (S_ref - Sᵏ) / S_ref * A * h
            end

            above_isohaline &= Sᵏ < S_ref
            above_isohaline && (freshwater_isohaline[r] += (S_ref - Sᵏ) / S_ref * V)
            top += Δz
        end
    end
    return nothing
end

# Mean depth (m, positive) of each level of `grid` below the surface at rest, from the bottom up
level_depths(grid) = -collect(znodes(grid.underlying_grid, Center()))

# ## Files
#
# The regions (set, codes, names, masks, areas), the levels and the settings; and the diagnostics of a list of times
# (`results`, one `regional_diagnostics` per time), each with time as the last dimension, with `dates` (strings) and any
# `extra` entries: a run's file holds both (save_regional_diagnostics).

function write_regions!(file, regions, grid; level_depths = nothing, settings,
                        reference_density = default_reference_density, heat_capacity = default_heat_capacity)
    if !isnothing(grid)
        Nx, Ny = size(grid)[1:2]
        area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
        file["regions/area"] = [sum(area[regions.masks[code]]) for code in regions.codes]   # m²
        for code in regions.codes
            file["regions/masks/$code"] = Array(regions.masks[code])
        end
    end
    file["regions/set"]   = regions.set
    file["regions/codes"] = regions.codes
    file["regions/names"] = regions.names
    isnothing(level_depths) || (file["level_depths"] = level_depths)   # m, from the bottom up, as the profiles' levels
    write_settings!(file, settings; reference_density, heat_capacity)
end

function write_settings!(file, settings; reference_density = default_reference_density,
                         heat_capacity = default_heat_capacity)
    file["depth_bounds"] = settings.depth_bounds   # m; Inf is the whole water column
    file["reference_temperature"] = settings.reference_temperature
    file["reference_salinity"]    = settings.reference_salinity
    file["reference_density"]     = reference_density
    file["heat_capacity"]         = heat_capacity
end

function write_results!(file, results; dates, extra...)
    file["dates"] = string.(dates)
    for name in keys(first(results))
        file[string(name)] = stack(getfield.(results, name))   # time last
    end
    for (name, value) in extra
        file[string(name)] = value
    end
end

"""
    save_regional_diagnostics(path, results, regions, grid; dates, level_depths, settings, extra...)

A run's diagnostics of one region set: regions, settings and results in the JLD2 file `path`.
"""
function save_regional_diagnostics(path, results, regions, grid; dates, level_depths, settings, extra...)
    jldopen(path, "w") do file
        write_regions!(file, regions, grid; level_depths, settings)
        write_results!(file, results; dates, extra...)
    end
    @info "Wrote $path"
end
