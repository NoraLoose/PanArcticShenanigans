# Wall-clock timing: how long each set-up phase takes, the first time step (mostly compilation), and the stepping
# speed. Printed as each phase finishes and summarised when the run ends. Included once by panarctic.jl, before
# setup.jl, so the set-up of the grid and bathymetry is timed too; rerun.jl uses the same functions.

phase_timer = Ref(time_ns())
phase_times = Pair{String, Float64}[]

seconds_since(t) = 1e-9 * (time_ns() - t)
format_seconds(s) = s < 120 ? string(round(s, digits = 1), " s") : string(round(s / 60, digits = 1), " min")

function finish_phase!(label)
    seconds = seconds_since(phase_timer[])
    push!(phase_times, label => seconds)
    @info "Timing: $label took $(format_seconds(seconds))"
    phase_timer[] = time_ns()
    return nothing
end

# Times inside run!: initialization (iteration 0: update_state!, first outputs and callbacks), the first time step
# (compiles the time step when the model type is new in this session), and every later step
run_timing = Dict{Symbol, Any}()

function start_run_timing!()
    empty!(run_timing)
    run_timing[:start] = time_ns()
    return nothing
end

function step_timer(sim)
    n = iteration(sim)
    if n == 0 && !haskey(run_timing, :initialized)
        run_timing[:initialized] = time_ns()
        @info "Timing: initialization inside run! took $(format_seconds(1e-9 * (run_timing[:initialized] - run_timing[:start])))"
    elseif n == 1 && !haskey(run_timing, :first_step)
        run_timing[:first_step] = time_ns()
        run_timing[:first_step_model_time] = time(sim)
        since = get(run_timing, :initialized, run_timing[:start])
        @info "Timing: first time step took $(format_seconds(1e-9 * (run_timing[:first_step] - since))) " *
              "(compilation when this model type is new in the session)"
    end
    return nothing
end

function report_run_timing(sim)
    total = seconds_since(run_timing[:start])
    lines = ["Timing summary:"]
    for (label, seconds) in phase_times
        push!(lines, "   $(rpad(label, 60)) $(format_seconds(seconds))")
    end
    if haskey(run_timing, :first_step)
        initialization = 1e-9 * (get(run_timing, :initialized, run_timing[:start]) - run_timing[:start])
        first_step = 1e-9 * (run_timing[:first_step] - get(run_timing, :initialized, run_timing[:start]))
        stepping = seconds_since(run_timing[:first_step])
        model_hours = (time(sim) - run_timing[:first_step_model_time]) / 3600
        push!(lines, "   $(rpad("run!: initialization", 60)) $(format_seconds(initialization))")
        push!(lines, "   $(rpad("run!: first time step", 60)) $(format_seconds(first_step))")
        push!(lines, "   $(rpad("run!: later steps ($(round(model_hours, digits = 1)) model hours)", 60)) " *
                     "$(format_seconds(stepping))" *
                     (model_hours > 0 ? ", $(round(stepping / model_hours, digits = 2)) s per model hour" : ""))
    end
    push!(lines, "   $(rpad("run! in total", 60)) $(format_seconds(total))")
    @info join(lines, "\n")
    return nothing
end
