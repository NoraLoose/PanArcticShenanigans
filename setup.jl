# Shared setup for panarctic.jl and the plotting scripts (plotting/): settings, grid and dates.
# Included with `include("setup.jl")`; the config file is the first command-line argument (default config.toml).

include("config.jl")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: GEBCO2026   # not exported by NumericalEarth
using Dates, CUDA
using CopernicusMarine   # enables the GLORYS download extension

arch = config["architecture"] == "GPU" ? GPU() : CPU()

# ## Grid
#
# Polar stereographic (ConformalConicGrid, λ₀ = −40°). Open boundaries on the Atlantic
# (j = 1) and Pacific (j = Ny) edges; east/west edges are land. `spacing_km` is the largest
# spacing (at the pole; cells shrink ~15% toward the Atlantic edge).

Nz      = config["grid"]["Nz"]
depth   = config["grid"]["depth"] * meters
spacing = config["grid"]["spacing_km"] * kilometers
z       = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

x₁, x₂  = -3280e3, 3616e3   # west–east projected extent (m)
y₁, y₂  = -4100e3, 2700e3   # south–north projected extent (m)
underlying_grid = ConformalConicGrid(arch;
                                     size               = (round(Int, (x₂ - x₁) / spacing), round(Int, (y₂ - y₁) / spacing), Nz),
                                     x                  = (x₁, x₂),
                                     y                  = (y₁, y₂),
                                     central_longitude  = -40,   # λ₀: mid-Atlantic
                                     standard_parallel  = 90,
                                     latitude_of_origin = 90,
                                     z,
                                     halo = (7, 7, 7))
φ_south = 43   # southern limit of GLORYS data (grid reaches 43.6°N)

Nx, Ny = size(underlying_grid)[1:2]
@info "Grid: $Nx × $Ny × $Nz, largest spacing $(spacing / 1e3) km"

# GEBCO 2026 (15″; includes IBCAO for the Arctic). First run downloads ~7.5 GB, and regridding a new grid
# needs ~55 GB of memory (the global file is read in full); the result is cached, so later runs are quick.
# Regrid on a CPU copy of the grid: the whole file would otherwise go to the GPU (~28 GB, too much for 40 GB cards).
# The cache does not depend on the architecture, so CPU and GPU runs share it.
bathymetry = Metadatum(:bottom_height; dataset = GEBCO2026(), dir_kw...)
bottom_height = regrid_bathymetry(on_architecture(CPU(), underlying_grid), bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(Array(interior(bottom_height, :, :, 1)));
                            active_cells_map = true)

# ## Dates
#
# The run goes from 00:00 on start_date to 00:00 on end_date (or to end_date's time of day, if it has one).
# Model time counts from 00:00 on start_date, as in nested_ocean_model, which starts the ocean and sea ice
# from GLORYS on start_date.

dataset    = GLORYSDaily()
start_date = DateTime(config["time"]["start_date"])
end_date   = DateTime(config["time"]["end_date"])
end_date > start_date || error("end_date ($end_date) must be after start_date ($start_date).")

seconds(Δ) = Dates.value(Millisecond(Δ)) / 1000
stop_time  = seconds(end_date - start_date)

# GLORYS daily means over the run and the region they are read from, for plotting/glorys_free_surface.jl.
dates   = floor(start_date, Day) : Day(1) : floor(end_date, Day)
region  = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))

# GLORYS free surface on the model grid, for plotting/glorys_free_surface.jl.
glorys_free_surface() = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid,
                                        inpainting=100, time_indices_in_memory=length(dates))
