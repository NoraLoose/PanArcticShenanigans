# The pan-Arctic ocean nested in GLORYS with NumericalEarth's nested_ocean_model (PR #676), on the same polar
# stereographic grid as panarctic.jl. Ocean only by default; with [coupling] coupled = true, also sea ice nested in
# GLORYS and the JRA55 atmosphere, radiation and rivers.
#
# The nested model packages the open boundaries of panarctic.jl (oblique radiation on u, v; tracer reservoirs;
# Flather) and adds Chapman on η, a sponge on T, S and the baroclinic velocity, and the seafloor blended to GLORYS's
# near the open edges. The initial state is GLORYS interpolated onto the grid.
#
# Run from this folder (its own environment: Oceananigans 0.117, NumericalEarth pinned in Project.toml):
#
#   julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate()'   # once
#   julia +1.13 --project=. panarctic_nested.jl [config file]   # default config.toml
#
# Hourly diagnostics go to the log and to <output_dir>/nested_monitor.csv; surface fields every 3 hours and η every
# hour to <output_dir>/nested_surface.jld2 and nested_free_surface.jld2.

using TOML
config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: ConformalConicGrid
using NumericalEarth.DataWrangling: GEBCO2026
using CopernicusMarine   # enables the GLORYS download extension
using CUDA, Dates, Printf

include("monitor.jl")
include("net_inflow_control.jl")
include("glorys_volume_target.jl")

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

# ## Ocean (and sea ice) nested in GLORYS
#
# GLORYS daily means from start_date to one day past the end of the run. Only the Atlantic and Pacific edges are open
# (the default would open every side that touches water).

start_date = DateTime(config["time"]["start_date"])
run_days   = config["time"]["days"]
dates      = (start_date, start_date + Day(ceil(Int, run_days) + 1))
Δt         = config["time"]["dt_minutes"] * minutes
stop_time  = run_days * Oceananigans.Units.day   # qualified, so a `days` defined in the REPL cannot change it

parent_dataset  = GLORYSDaily()
zero_net_inflow    = config["boundaries"]["zero_net_inflow"]
net_inflow_control = get(config["boundaries"], "net_inflow_control", false)
net_inflow_control && !zero_net_inflow && error("[boundaries] net_inflow_control = true needs zero_net_inflow = true")
edge_inflow        = get(config["boundaries"], "edge_inflow_Sv", nothing)
edge_targets       = isnothing(edge_inflow) ? nothing : (south = 1e6 * edge_inflow["atlantic"], north = 1e6 * edge_inflow["pacific"])
total_inflow       = get(config["boundaries"], "net_inflow_Sv", nothing)
follow_glorys      = get(config["boundaries"], "follow_glorys_volume", false)
net_inflow         = follow_glorys ? Ref(0.0) : isnothing(total_inflow) ? nothing : 1e6 * total_inflow
@info "Open boundaries: zero_net_inflow = $zero_net_inflow, net_inflow = $net_inflow, edge_transports = $edge_targets (m³/s, into the domain)"
coupled         = config["coupling"]["coupled"]

ocean = nested_ocean_simulation(grid, parent_dataset; dates, open_sides = (:south, :north), zero_net_inflow,
                                net_inflow, edge_transports = edge_targets,
                                Δt, stop_time, dir_kw...)
child = ocean.model.child

if coupled
    sea_ice = nested_sea_ice_simulation(ocean, parent_dataset; dates, dir_kw...)

    jra55 = (; dataset = MultiYearJRA55(), start_date = dates[1], end_date = dates[2], dir_kw...)
    atmosphere = JRA55PrescribedAtmosphere(arch; jra55...)
    if !config["coupling"]["pressure"]
        # The whole record in memory, so the uniform pressure is not overwritten by a reload
        atmosphere = JRA55PrescribedAtmosphere(arch; jra55..., time_indices_in_memory = length(atmosphere.velocities.u.times))
        fill!(parent(atmosphere.pressure.data), 101325)
        @warn "JRA55 sea-level pressure set to 101325 Pa everywhere ([coupling] pressure = false)"
    end
    radiation = JRA55PrescribedRadiation(arch; jra55...)
    land      = JRA55PrescribedLand(child.grid; jra55...)

    simulation = Simulation(EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation); Δt, stop_time)
else
    simulation = ocean
end

# ## Diagnostics and output

parent_η = parent_free_surface(grid, parent_dataset, dates; dir_kw...)
add_nested_monitor!(simulation, child, output_dir; exchanger = ocean.model.exchanger, parent_η)
net_inflow_control && add_net_inflow_control!(simulation, child, ocean.model.exchanger, output_dir)
follow_glorys && add_glorys_volume_target!(simulation, net_inflow, child, parent_η,
                                           config["boundaries"]["volume_relaxation_days"] * Oceananigans.Units.day)

surface = (; u = child.velocities.u, v = child.velocities.v, T = child.tracers.T, S = child.tracers.S)
simulation.output_writers[:surface] = JLD2Writer(child, surface; indices = (:, :, Nz), schedule = TimeInterval(3hours),
                                                 filename = joinpath(output_dir, "nested_surface"),
                                                 overwrite_files = true)
simulation.output_writers[:free_surface] = JLD2Writer(child, (; η = child.free_surface.displacement),
                                                      schedule = TimeInterval(1hours),
                                                      filename = joinpath(output_dir, "nested_free_surface"),
                                                      overwrite_files = true)

if coupled
    function ice_progress(sim)
        h, ℵ = sea_ice.model.ice_thickness, sea_ice.model.ice_concentration
        @info @sprintf("         ice: max h %.2f m, mean ℵ %.3f", maximum(h), sum(interior(ℵ)) / length(interior(ℵ)))
        return nothing
    end
    add_callback!(simulation, ice_progress, TimeInterval(6hours))
end

run!(simulation)
