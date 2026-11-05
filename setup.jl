# Shared setup for panarctic.jl and the plotting scripts (plotting/): settings, grid and GLORYS dates/region.
# Included with `include("setup.jl")`; the config file is the first command-line argument (default config.toml).

include("config.jl")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using NumericalEarth.DataWrangling: all_dates
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

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...)
bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

# ## GLORYS dates and region
#
# Daily means (centred at noon), with one extra day before the start date and after the end of the run as
# margin, where GLORYS has them (it starts on 1993-01-01). Time counts from 00:00 on dates[1], and the run
# starts at noon on the start date, dates[start_index].

dataset    = GLORYSDaily()
available  = all_dates(dataset, :temperature)
start_date = DateTime(config["time"]["start_date"])
run_days   = ceil(Int, config["time"]["stop_hours"] / 24)
dates      = max(start_date - Day(1), first(available)) : Day(1) : min(start_date + Day(max(run_days, 1) + 1), last(available))
start_index = findfirst(==(start_date), dates)
isnothing(start_index) && error("start_date $start_date is outside the GLORYS record ($(first(available)) to $(last(available))).")
start_time  = (start_index - 1) * days + 12hours
last_noon   = (length(dates) - 1) * days + 12hours
start_time + config["time"]["stop_hours"] * hours ≤ last_noon ||
    error("The run ends after the last GLORYS daily mean ($(last(available))); use an earlier start_date or a shorter run.")

region  = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))

# GLORYS free surface on the model grid. Used by panarctic.jl and plotting/glorys_free_surface.jl,
# so both read the same inpainted data from the cache.
glorys_free_surface() = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid,
                                        inpainting=100, time_indices_in_memory=length(dates))
