# Mean temperature and salinity profiles of the regions of regions.jl, for the model's daily averages
# (analysis/ocean_diagnostics.jl) and for the gridded datasets: it needs only the grid (an ImmersedBoundaryGrid, whose
# cell thicknesses follow the free surface on the model's z-star grid), T and S, so it runs on any grid. The observed
# profiles (validation/) are written in the same form.
#
# For each region and level:
#   temperature, salinity  volume-weighted mean (°C, g/kg)
#   volume                 volume of the region's cells on that level (m³)

@isdefined(region_index) || include("diagnostic_files.jl")

"""
    ts_profiles(grid, T, S, regions)

The profiles above for one time, on the CPU: `T` and `S` are fields on `grid`, `regions` a
set of `region_masks(grid; set)`. Cells where T or S is not finite count as land. Returns a `NamedTuple` of (level,
region) arrays temperature, salinity and volume, the regions in the order of `regions.codes` and the levels from the
bottom up; a level without ocean in a region has NaN.
"""
function ts_profiles(grid, T, S, regions)
    Nz, R = size(grid, 3), length(regions.codes)
    sums = (∑TV = zeros(Nz, R), ∑SV = zeros(Nz, R), volume = zeros(Nz, R))
    add_profiles!(sums, grid, Array(interior(T)), Array(interior(S)), region_index(grid, regions))
    (; ∑TV, ∑SV, volume) = sums
    temperature = ifelse.(volume .> 0, ∑TV ./ volume, NaN)
    salinity    = ifelse.(volume .> 0, ∑SV ./ volume, NaN)
    return (; temperature, salinity, volume)
end

# The sums over every ocean column of every region (a function of its own, so that it compiles for the concrete types)
function add_profiles!(sums, grid, θ, s, region)
    (; ∑TV, ∑SV, volume) = sums
    Nx, Ny, Nz = size(grid)
    ccc = (Center(), Center(), Center())

    for j in 1:Ny, i in 1:Nx
        r = region[i, j]
        r == 0 && continue
        A = Azᶜᶜᶜ(i, j, 1, grid)
        for k in Nz:-1:1
            (inactive_node(i, j, k, grid, ccc...) || !isfinite(θ[i, j, k]) || !isfinite(s[i, j, k])) && break
            V = A * Δzᶜᶜᶜ(i, j, k, grid)
            ∑TV[k, r]    += θ[i, j, k] * V
            ∑SV[k, r]    += s[i, j, k] * V
            volume[k, r] += V
        end
    end
    return nothing
end
