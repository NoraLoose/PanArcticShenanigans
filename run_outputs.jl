# Monitors, callbacks and output writers for one run. Included by panarctic.jl after the simulation is built, and again
# by rerun.jl for each rerun (new output_dir, same model). Every callback and writer added here is named and listed in
# run_callback_names / run_writer_names, so rerun.jl can remove them before adding them again.

run_callback_names = Symbol[]
run_writer_names   = Symbol[]

function add_run_callback!(func, schedule, name)
    add_callback!(simulation, func, schedule; name)
    push!(run_callback_names, name)
    return nothing
end

function add_run_writer!(name, writer)
    ocean.output_writers[name] = writer
    push!(run_writer_names, name)
    return nothing
end

# ## Progress and timing

wall_time = Ref(time_ns())

function progress(sim)
    model = ocean_model_of(sim)
    T = model.tracers.T
    S = model.tracers.S
    u, v, w = model.velocities
    η = model.free_surface.displacement
    ice = ocean_only ? "" : @sprintf(" max(h): %.2f m,", maximum(sim.model.sea_ice.model.ice_thickness))
    msg = @sprintf("time: %s, iter: %d, extrema(T, S): (%.1f, %.1f) °C (%.1f, %.1f) psu, extrema(u) (%.2e, %.2e, %.2e) max|η|: %.3f m,%s wall: %s",
                   prettytime(sim), iteration(sim),
                   extrema(T)..., extrema(S)..., maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(abs, η), ice,
                   prettytime(1e-9 * (time_ns() - wall_time[])))
    @info msg
    wall_time[] = time_ns()
    return nothing
end

add_run_callback!(progress, TimeInterval(1hours), :panarctic_progress)
add_run_callback!(step_timer, IterationInterval(1), :panarctic_step_timer)

function ice_progress(sim)
    uᵢ, vᵢ = sim.model.sea_ice.model.velocities
    @info @sprintf("    ice: max |u| = %.2f m/s, max |v| = %.2f m/s", maximum(abs, uᵢ), maximum(abs, vᵢ))
    return nothing
end
ocean_only || add_run_callback!(ice_progress, TimeInterval(1hours), :panarctic_ice_progress)

# ## Open-boundary transports (open_boundary_transport.jl)

reset_transport_state!()
monitor_file = start_transport_monitor()
add_run_callback!(transport_monitor, TimeInterval(1hours), :panarctic_transport_monitor)
add_run_callback!(update_surface_freshwater!, IterationInterval(1), :panarctic_surface_freshwater)

# Gateway transports (Greenland–Scotland passages, Svinøy, Lofoten slope, Barents Sea Opening, Fram, Davis and Bering
# Straits), daily means near the Atlantic edge and of the barotropic flow
get(get(config, "diagnostics", Dict()), "sections", false) && include("section_monitor.jl")

# Depth-integrated momentum budget by depth class: which term moves water onto or off the shelves
get(get(config, "diagnostics", Dict()), "tendencies", false) && include("tendency_budget.jl")

# ## Output: hourly fields near the surface

u, v, w = ocean.model.velocities
𝒱 = @at((Center, Center, Center), sqrt(u^2 + v^2))
ocean_outputs = merge(ocean.model.tracers, (; 𝒱))

add_run_writer!(:surface, JLD2Writer(ocean.model, ocean_outputs;
                                     filename = ocean_surface_file,
                                     indices = (:, :, grid.Nz-2),
                                     schedule = TimeInterval(1hours),
                                     overwrite_files = true))

if !ocean_only
    h = sea_ice.model.ice_thickness
    ℵ = sea_ice.model.ice_concentration
    he = h * ℵ
    sea_ice_outputs = (; he, uᵢ = sea_ice.model.velocities.u, vᵢ = sea_ice.model.velocities.v)
    sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                  filename = sea_ice_surface_file,
                                                  schedule = TimeInterval(1hours),
                                                  overwrite_files = true)
end
