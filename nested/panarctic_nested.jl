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

include("setup.jl")
include("monitor.jl")
include("gateways.jl")
include("net_inflow_control.jl")
include("glorys_volume_target.jl")

# ## Ocean (and sea ice) nested in GLORYS
#
# Only the Atlantic and Pacific edges are open (the default would open every side that touches water).

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
net_inflow_depth_exponent = get(config["boundaries"], "net_inflow_depth_exponent", 0)
@info "Open boundaries: zero_net_inflow = $zero_net_inflow, net_inflow = $net_inflow (pin ∝ H^$net_inflow_depth_exponent), edge_transports = $edge_targets (m³/s, into the domain)"
coupled         = config["coupling"]["coupled"]

ocean = nested_ocean_simulation(grid, parent_dataset; dates, open_sides = (:south, :north), zero_net_inflow,
                                net_inflow, net_inflow_depth_exponent, edge_transports = edge_targets,
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
add_nested_monitor!(simulation, child, output_dir; exchanger = ocean.model.exchanger, parent_η,
                    depth_exponent = net_inflow_depth_exponent)
add_gateway_monitor!(simulation, child, output_dir)
net_inflow_control && add_net_inflow_control!(simulation, child, ocean.model.exchanger, output_dir)
follow_glorys && add_glorys_volume_target!(simulation, net_inflow, child, parent_η,
                                           config["boundaries"]["volume_relaxation_days"] * Oceananigans.Units.day;
                                           tendency = get(config["boundaries"], "follow_glorys_tendency", false))

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

# PANARCTIC_DOWNLOAD_ONLY=1: build everything (which downloads all the GLORYS and JRA55 data for the run's dates)
# and stop before the first time step. Use with PANARCTIC_ARCHITECTURE=CPU to leave the GPUs free.
if get(ENV, "PANARCTIC_DOWNLOAD_ONLY", "0") == "1"
    @info "PANARCTIC_DOWNLOAD_ONLY=1: the data for $(dates[1]) to $(dates[2]) are downloaded; not running"
else
    run!(simulation)
end
