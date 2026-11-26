# Everything needed to compare a run with the validation datasets, in one command (see validation/README.md):
#
#   julia +1.13 --project=. analysis/compare_with_validation.jl [config.toml]
#   julia +1.13 --project=. analysis/compare_with_validation.jl config_a.toml config_b.toml ...   (several runs)
#
# Runs only the steps whose results are missing or out of date, each with its own script (which can also be run alone):
#   1. the validation datasets' regional diagnostics, if not there yet (analysis/validation_diagnostics.jl WOA23, TSAR,
#      and UDASH for the run's years)
#   2. GLORYS on the run's days: regional diagnostics and transports (analysis/validation_diagnostics.jl GLORYS --run)
#   3. the run's regional diagnostics, if missing or older than its output (analysis/regional_diagnostics.jl)
#   4. the figures (in figure_dir): plotting/regional_diagnostics.jl, plotting/gateway_transports.jl, and
#      analysis/tsubouchi_comparison.jl if the run overlaps Tsubouchi et al.'s Sep 2005 - Aug 2006
# With several runs, steps 1-3 are done for each, and the figures show them all together (in a folder
# compare_<experiment>_<experiment>... in the figure folder). Needs no GPU. Each step is a separate Julia process, so a
# step that fails leaves the results of the others.

include("../config.jl")
include("../transports.jl")
include("../validation_config.jl")

using Dates

repo = joinpath(@__DIR__, "..")
config_files = filter(endswith(".toml"), ARGS)
figures_wanted = !("--no-figures" in ARGS)

# Several runs: steps 1-3 for each (this script, without figures), then the figures of all of them
if length(config_files) > 1
    for file in config_files
        run(`$(Base.julia_cmd()) --project=$repo $(@__FILE__) $file --no-figures`)
    end
    run(`$(Base.julia_cmd()) --project=$repo $(joinpath(repo, "plotting", "regional_diagnostics.jl")) $config_files`)
    @info "Figures of $(length(config_files)) runs in the figure folder's compare_... folder"
    exit()
end
julia_script(script, args...) = `$(Base.julia_cmd()) --project=$repo $(joinpath(repo, script)) $args`

function step(description, script, args...)
    @info "$description: $script $(join(args, " "))"
    run(julia_script(script, args...))
end

set_names = region_set_names(validation)
run_days  = Date(config["time"]["start_date"]):Day(1):Date(DateTime(config["time"]["end_date"]) - Millisecond(1))
has(dataset, set, period) = isfile(period_file(dataset, set, period))

# 1. Datasets that do not depend on the run
haskey(validation["datasets"], "WOA23") &&
    !all(has("WOA23", set, decade) for set in set_names, decade in keys(dataset_options("WOA23")["decades"])) &&
    step("WOA23 climatologies", "analysis/validation_diagnostics.jl", "WOA23")
haskey(validation["datasets"], "TSAR") && "TSAR" in set_names && isempty(period_files("TSAR", "TSAR")) &&
    step("TSAR observed profiles", "analysis/validation_diagnostics.jl", "TSAR")

# UDASH: the run's years
run_years = unique(year.(run_days))
haskey(validation["datasets"], "UDASH") && !all(has("UDASH", set, y) for set in set_names, y in run_years) &&
    step("UDASH profiles of the run's years", "analysis/validation_diagnostics.jl", "UDASH", "--run", config_file)

# 2. GLORYS on the run's days
glorys_done = all(has("GLORYS", set, day) for set in set_names, day in run_days) &&
              all(isfile(joinpath(glorys_transport_dir(config), "$day.jld2")) for day in run_days)
haskey(validation["datasets"], "GLORYS") && !glorys_done &&
    step("GLORYS on the run's days", "analysis/validation_diagnostics.jl", "GLORYS", "--run", config_file)

# 3. The run's regional diagnostics
outdated(file) = !isfile(file) || mtime(file) < mtime(ocean_output_file)
isfile(ocean_output_file) || error("No model output $ocean_output_file: run panarctic.jl $config_file first")
any(outdated(regional_output_file(set)) for set in set_names) &&
    step("The run's regional diagnostics", "analysis/regional_diagnostics.jl", config_file)

# 4. Figures
figures_wanted || exit()
step("Regional figures", "plotting/regional_diagnostics.jl", config_file)
isfile(transport_output_file) ? step("Gateway transports", "plotting/gateway_transports.jl", config_file) :
                                @info "No online transports ($transport_output_file): no gateway figure"
first(run_days) < Date(2006, 9, 1) && last(run_days) ≥ Date(2005, 9, 5) &&
    step("Against Tsubouchi et al. (2018)", "analysis/tsubouchi_comparison.jl", config_file)

@info "Figures in $figure_dir"
