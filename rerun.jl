# Run the ocean-only simulation again in the same Julia session, from the same initial state, without rebuilding the
# model: no data loading, no balancing and no compilation. Use it in the REPL after panarctic.jl has run (to the end,
# or stopped with Ctrl-C) with [initialization] ocean_only = true (or a rest test):
#
#     empty!(ARGS); push!(ARGS, "my_next_run.toml"); include("rerun.jl")
#
# The new config may change [paths] (give each run its own output_dir), stop_hours in [time] (up to the last GLORYS day
# the first run loaded) and [diagnostics]. Everything else ([grid], [ocean], [boundaries], [initialization],
# [atmosphere], start_date, dt_minutes) is built into the model or its initial state, so rerun.jl stops if it differs:
# use a fresh session for those.

(@isdefined(initial_state) && @isdefined(simulation)) ||
    error("rerun.jl needs a session in which panarctic.jl has built an ocean-only model")
ocean_only || error("rerun.jl supports ocean-only runs ([initialization] ocean_only = true, or a rest test)")

phase_timer[] = time_ns()
empty!(phase_times)

previous = (; config, output_dir, ocean_surface_file, sea_ice_surface_file)
include("config.jl")   # the new config, output_dir and file names

let changed = String[]
    for key in union(keys(previous.config), keys(config))
        key in ("paths", "diagnostics", "time") && continue
        get(previous.config, key, nothing) == get(config, key, nothing) || push!(changed, "[$key]")
    end
    for key in ("start_date", "dt_minutes")
        previous.config["time"][key] == config["time"][key] || push!(changed, "$key in [time]")
    end
    new_stop_hours = config["time"]["stop_hours"]
    too_long = start_time + new_stop_hours * hours > last_noon
    if !isempty(changed) || too_long
        global config               = previous.config
        global output_dir           = previous.output_dir
        global ocean_surface_file   = previous.ocean_surface_file
        global sea_ice_surface_file = previous.sea_ice_surface_file
        isempty(changed) || error("rerun.jl: $(join(changed, ", ")) differ from the model in this session; " *
                                  "start a fresh session for this config")
        error("rerun.jl: stop_hours = $new_stop_hours goes past the GLORYS days loaded in this session " *
              "(at most $(Int((last_noon - start_time) / hours)) hours); start a fresh session for a longer run")
    end
end

for name in run_callback_names
    delete!(simulation.callbacks, name)
end
for name in run_writer_names
    delete!(ocean.output_writers, name)
end

# Oceananigans 0.116.1 restores the open-boundary radiation buffers with copyto! from host to device offset arrays,
# which falls back to scalar indexing on the GPU. The buffers are small (one row per edge), so allow it here.
CUDA.@allowscalar Oceananigans.restore_prognostic_state!(ocean.model, initial_state)
simulation.stop_time = start_time + config["time"]["stop_hours"] * hours

include("run_outputs.jl")
@info "Rerun from the initial state: model time $(prettytime(ocean.model.clock.time)), iteration " *
      "$(iteration(simulation)), to $(prettytime(simulation.stop_time)); output in $output_dir"
finish_phase!("rerun set-up: restore the initial state, monitors and writers")

start_run_timing!()
run!(simulation)
report_run_timing(simulation)
