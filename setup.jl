# Shared setup for panarctic.jl and the plotting scripts (plotting/): settings, grid and dates.
# Included with `include("setup.jl")`; the config file is the first command-line argument (default config.toml).

include("config.jl")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: GEBCO2026, ETOPO2022   # not exported by NumericalEarth
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

# Bathymetry ([grid] bathymetry), regridded on a CPU copy of the grid (the global file would otherwise go to the GPU:
# ~28 GB for GEBCO, too much for 40 GB cards):
# - "GEBCO2026" (default): 15″, includes IBCAO for the Arctic. First run downloads ~7.5 GB, and regridding a new grid
#   needs ~55 GB of memory (the global file is read in full); the result is cached, so later runs are quick.
# - "ETOPO": ETOPO 2022 (60″, NumericalEarth's default); cached.
# - "GLORYS": GLORYS's own sea floor (glorys_bathymetry.jl), so the gateways and sills have the depths GLORYS's
#   transports were computed over; recomputed each run.
# The cache does not depend on the architecture, so CPU and GPU runs share it.
include("glorys_bathymetry.jl")   # glorys_bottom_height; included here, not in the `if` below (world age)
bathymetry_kw = (; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
cpu_grid      = on_architecture(CPU(), underlying_grid)
bathymetry    = get(config["grid"], "bathymetry", "GEBCO2026")
bottom_height = if bathymetry == "GEBCO2026"
    regrid_bathymetry(cpu_grid, Metadatum(:bottom_height; dataset = GEBCO2026(), dir_kw...); bathymetry_kw...)
elseif bathymetry == "ETOPO"
    regrid_bathymetry(cpu_grid, Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...); bathymetry_kw...)
elseif bathymetry == "GLORYS"
    glorys_bottom_height(cpu_grid; bathymetry_kw...)
else
    error("[grid] bathymetry must be \"GEBCO2026\", \"ETOPO\" or \"GLORYS\", not \"$bathymetry\"")
end

# Bottom cells ([grid] bottom): "partial" (default) cuts the bottom cell at the sea floor; "full" uses whole cells
# (the sea floor moves to the nearest cell face). With partial cells, [grid] partial_cell_fix (panarctic.jl) removes
# the bottom jets they otherwise allow (Oceananigans issue #6532).
bottom_kind = get(config["grid"], "bottom", "partial")
bottom_kind in ("partial", "full") || error("[grid] bottom must be \"partial\" or \"full\", not \"$bottom_kind\"")
Bottom = bottom_kind == "full" ? GridFittedBottom : PartialCellBottom
grid = ImmersedBoundaryGrid(underlying_grid, Bottom(Array(interior(bottom_height, :, :, 1))); active_cells_map = true)
@info "Bathymetry: $bathymetry, $bottom_kind bottom cells"

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
