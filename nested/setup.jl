# Settings, packages, paths, the grid and the run's dates for the nested pan-Arctic model: included by
# panarctic_nested.jl and by glorys_gateways_run.jl. Reads the config file named in ARGS (default config.toml).

using TOML
config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using NumericalEarth.DataWrangling: GEBCO2026
using CopernicusMarine   # enables the GLORYS download extension
using CUDA, Dates, Printf


arch       = config["architecture"] == "GPU" ? GPU() : CPU()
data_dir   = config["paths"]["data_dir"]
dir_kw     = isempty(data_dir) ? (;) : (; dir = data_dir)
output_dir = mkpath(config["paths"]["output_dir"])

# ## Grid: as panarctic.jl (λ₀ = −40°, pole inside; Atlantic edge j = 1, Pacific edge j = Ny)

Nz      = config["grid"]["Nz"]
depth   = config["grid"]["depth"] * meters
spacing = config["grid"]["spacing_km"] * kilometers
z       = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

x₁, x₂ = -3280e3, 3616e3
y₁, y₂ = -4100e3, 2700e3
underlying_grid = ConformalConicGrid(arch;
                                     size               = (round(Int, (x₂ - x₁) / spacing), round(Int, (y₂ - y₁) / spacing), Nz),
                                     x                  = (x₁, x₂),
                                     y                  = (y₁, y₂),
                                     central_longitude  = -40,
                                     standard_parallel  = 90,
                                     latitude_of_origin = 90,
                                     z,
                                     halo = (7, 7, 7))

Nx, Ny = size(underlying_grid)[1:2]
@info "Grid: $Nx × $Ny × $Nz, largest spacing $(spacing / 1e3) km"

# Regridded on a CPU copy of the grid, as in setup.jl; the result is cached
bathymetry_kw = (; minimum_depth = 15, interpolation_passes = 2, major_basins = 1)
cpu_grid = on_architecture(CPU(), underlying_grid)
bottom_height = if config["grid"]["bathymetry"] == "GEBCO2026"
    regrid_bathymetry(cpu_grid, Metadatum(:bottom_height; dataset = GEBCO2026(), dir_kw...); bathymetry_kw...)
else
    regrid_bathymetry(cpu_grid; bathymetry_kw...)
end

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(Array(interior(bottom_height, :, :, 1)));
                            active_cells_map = true)

# GLORYS daily means from start_date to one day past the end of the run
start_date = DateTime(config["time"]["start_date"])
run_days   = config["time"]["days"]
dates      = (start_date, start_date + Day(ceil(Int, run_days) + 1))
