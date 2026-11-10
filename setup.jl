# Shared setup for panarctic.jl and the plotting scripts (plotting/): settings, grid and GLORYS dates/region.
# Included with `include("setup.jl")`; the config file is the first command-line argument (default config.toml).

include("config.jl")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: GEBCO2026, all_dates   # not exported by NumericalEarth
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

# ## GLORYS dates and region
#
# The run goes from noon on start_date to noon on end_date (or to end_date's time of day, if it has one).
# GLORYS daily means (centred at noon) are loaded with one extra day before and after the run as margin,
# where GLORYS has them (it starts on 1993-01-01). Model time counts from 00:00 on dates[1].

dataset    = GLORYSDaily()
available  = all_dates(dataset, :temperature)
start_date = DateTime(config["time"]["start_date"])
end_date   = config["time"]["end_date"]
run_start  = start_date + Hour(12)
run_end    = end_date isa Date ? DateTime(end_date) + Hour(12) : DateTime(end_date)
run_end > run_start || error("end_date ($end_date) must be after noon on start_date ($(Date(start_date))).")

last_mean  = ceil(run_end - Hour(12), Day)   # date of the first daily mean at or after the end of the run
dates      = max(start_date - Day(1), first(available)) : Day(1) : min(last_mean + Day(1), last(available))
start_index = findfirst(==(start_date), dates)
isnothing(start_index) && error("start_date $start_date is outside the GLORYS record ($(first(available)) to $(last(available))).")
last_mean ≤ last(available) || error("The run ends after the last GLORYS daily mean ($(last(available))).")

seconds(Δ) = Dates.value(Millisecond(Δ)) / 1000
start_time = seconds(run_start - dates[1])
stop_time  = seconds(run_end   - dates[1])

region  = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))

# GLORYS free surface on the model grid. Used by panarctic.jl and plotting/glorys_free_surface.jl,
# so both read the same inpainted data from the cache.
glorys_free_surface() = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid,
                                        inpainting=100, time_indices_in_memory=length(dates))
