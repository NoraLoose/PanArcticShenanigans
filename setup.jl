# Shared setup for panarctic.jl and the plotting scripts (plotting/): settings, grid and GLORYS dates/region.
# Included with `include("setup.jl")`; the config file is the first command-line argument (default config.toml).

include("config.jl")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using NumericalEarth.DataWrangling: GEBCO2026, all_dates   # not exported by NumericalEarth
using Dates, CUDA, Printf
using CopernicusMarine   # enables the GLORYS download extension

arch = config["architecture"] == "GPU" ? GPU() : CPU()

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

φ_south = 43   # southern limit of GLORYS data (grid reaches 43.6°N)
region  = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))

# ## Grid
#
# Polar stereographic (ConformalConicGrid, λ₀ = −40°). Open boundaries on the Atlantic
# (j = 1) and Pacific (j = Ny) edges; east/west edges are land. `spacing_km` is the largest
# spacing (at the pole; cells shrink ~15% toward the Atlantic edge).

# [grid] bathymetry = "gebco" (default) or "glorys"; vertical = "exponential" (default, Nz levels to `depth`) or
# "glorys" (GLORYS's own 50 levels; Nz and depth are then ignored). See glorys_bathymetry.jl.
bathymetry_source = get(config["grid"], "bathymetry", "gebco")
vertical_levels   = get(config["grid"], "vertical", "exponential")
bathymetry_source in ("gebco", "glorys") || error("[grid] bathymetry must be \"gebco\" or \"glorys\"")
vertical_levels in ("exponential", "glorys") || error("[grid] vertical must be \"exponential\" or \"glorys\"")
(bathymetry_source == "glorys" || vertical_levels == "glorys") && include("glorys_bathymetry.jl")

spacing = config["grid"]["spacing_km"] * kilometers
if vertical_levels == "glorys"
    z_faces = glorys_vertical_faces()
    Nz      = length(z_faces) - 1
    depth   = -first(z_faces)
    z       = Oceananigans.Grids.MutableVerticalDiscretization(z_faces)
    @info @sprintf("Vertical: GLORYS's %d levels, top cell %.2f m, bottom at %.0f m", Nz, -z_faces[end-1], depth)
else
    Nz      = config["grid"]["Nz"]
    depth   = config["grid"]["depth"] * meters
    z       = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)
end

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

Nx, Ny = size(underlying_grid)[1:2]
@info "Grid: $Nx × $Ny × $Nz, largest spacing $(spacing / 1e3) km"

# GEBCO 2026 (15″; includes IBCAO for the Arctic). First run downloads ~7.5 GB, and regridding a new grid
# needs ~55 GB of memory (the global file is read in full); the result is cached, so later runs are quick.
# Or GLORYS's own sea floor ([grid] bathymetry = "glorys"), regridded on the CPU at every start (glorys_bathymetry.jl).
if bathymetry_source == "glorys"
    bottom_height = Field{Center, Center, Nothing}(underlying_grid)
    set!(bottom_height, Array(interior(glorys_bottom_height(underlying_grid; minimum_depth = 15,
                                                                            interpolation_passes = 2,
                                                                            major_basins = 1))))
else
    bathymetry = Metadatum(:bottom_height; dataset = GEBCO2026(), dir_kw...)
    bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                      minimum_depth = 15,
                                      interpolation_passes = 2,
                                      major_basins = 1)
end

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

# GLORYS free surface on the model grid. Used by panarctic.jl and plotting/glorys_free_surface.jl,
# so both read the same inpainted data from the cache.
glorys_free_surface() = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid,
                                        inpainting=100, time_indices_in_memory=length(dates))
