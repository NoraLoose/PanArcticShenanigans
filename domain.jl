# The model domain: its map projection and extent, which do not depend on the grid spacing or the config. Included by
# setup.jl (the model grid); within_model_domain tells which columns of another grid (e.g. a dataset's) lie inside the
# domain, without a run config or the bathymetry.

using Oceananigans
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid, geographic_to_conformal_conic

# Polar stereographic (ConformalConicGrid, λ₀ = −40°). Open boundaries on the Atlantic (south) and Pacific (north) edges.
const domain_projection = (central_longitude  = -40,   # λ₀: mid-Atlantic
                           standard_parallel  = 90,
                           latitude_of_origin = 90)
const domain_x = (-3280e3, 3616e3)   # west–east projected extent (m)
const domain_y = (-4100e3, 2700e3)   # south–north projected extent (m)
const φ_south  = 43                  # southern limit of the data read for the domain (the grid reaches 43.6°N)

# The domain's map projection, from a coarse grid over the domain
const domain_mapping = ConformalConicGrid(CPU(); size = (4, 4, 1), x = domain_x, y = domain_y, z = (-1, 0),
                                          domain_projection...).conformal_mapping

# Whether the point (λ, φ) (degrees) lies inside the model domain
function in_model_domain(λ, φ)
    x, y = geographic_to_conformal_conic(domain_mapping, λ, φ)
    return domain_x[1] ≤ x ≤ domain_x[2] && domain_y[1] ≤ y ≤ domain_y[2]
end

"""
    within_model_domain(grid)

For each column of `grid` (e.g. GLORYS's or WOA's cells), whether its centre lies inside the model domain: a Bool
matrix (Nx, Ny).
"""
function within_model_domain(grid)
    Nx, Ny = size(grid)[1:2]
    ccc = (Center(), Center(), Center())
    return [in_model_domain(λnode(i, j, 1, grid, ccc...), φnode(i, j, 1, grid, ccc...)) for i in 1:Nx, j in 1:Ny]
end
