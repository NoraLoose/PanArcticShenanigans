# Diagnostics of a validation dataset (validation.toml [datasets]), computed with the same code and settings as a run's,
# each on the dataset's own grid: regional diagnostics (regional_diagnostics.jl, validation.toml [regional_diagnostics])
# over the columns inside the model domain (domain.jl), kept in <validation_dir>/<dataset>/diagnostics/<set>/; for
# GLORYS with --run also the transports through the run's sections (transports.jl, the config's [transports]), kept in
# <validation_dir>/GLORYS/transports/. Both are shared by all runs; what is already there is skipped (the regional
# diagnostics are computed again if the settings changed).
#
#   julia +1.13 --project=. analysis/validation_diagnostics.jl WOA23
#   julia +1.13 --project=. analysis/validation_diagnostics.jl TSAR
#   julia +1.13 --project=. analysis/validation_diagnostics.jl UDASH [--run config.toml]   (all years, or the run's)
#   julia +1.13 --project=. analysis/validation_diagnostics.jl GLORYS --run config.toml
#   julia +1.13 --project=. analysis/validation_diagnostics.jl GLORYS --from 1993-01-01 --to 1993-03-01
#
# GLORYS: the days of a run (--run: its config's start_date to end_date), with the transports through its sections; or
# the regional diagnostics only, from --from up to, not including, --to.
# Needs no GPU and no model run; downloads what is missing (validation/<dataset>.jl).

include("../validation_config.jl")   # validation, validation_dir, the paths and settings

dataset = isempty(ARGS) || startswith(first(ARGS), "--") ?
          error("Which dataset? One of $(join(dataset_names, ", ")): analysis/validation_diagnostics.jl <dataset>") :
          first(ARGS)
options = dataset_options(dataset)

using Dates

include("../domain.jl")
include("../regional_diagnostics.jl")
include("../validation/glorys.jl")
include("../validation/woa23.jl")
include("../validation/tsar.jl")
include("../validation/udash.jl")

settings  = regional_settings(validation)
set_names = region_set_names(validation)

# The run's config (--run), if any
run_file   = command_line_option("run")
run_config = isnothing(run_file) ? nothing : TOML.parsefile(run_file)

# The days of GLORYS's daily means: those of the run, or --from up to --to
function glorys_days()
    if !isnothing(run_config)
        from, to = DateTime(run_config["time"]["start_date"]), DateTime(run_config["time"]["end_date"])
    else
        from, to = command_line_option("from"), command_line_option("to")
        (isnothing(from) || isnothing(to)) && error("GLORYS needs --run config.toml, or --from and --to")
        from, to = DateTime(from), DateTime(to)
    end
    return collect(Date(from):Day(1):Date(to - Millisecond(1)))
end

if dataset == "GLORYS"
    days = glorys_days()
    glorys_regional_diagnostics!(days; settings, set_names)
    isnothing(run_config) ? @info("No --run config: no transports (they are through a run's sections)") :
                            glorys_transports!(days, run_config)
elseif dataset == "WOA23"
    woa23_regional_diagnostics!(options["decades"]; settings, set_names)
elseif dataset == "UDASH"
    years = isnothing(run_config) ? (first(options["years"]):last(options["years"])) : unique(year.(glorys_days()))
    udash_regional_diagnostics!(options, years; settings, set_names)
elseif dataset == "TSAR"
    tsar_regional_diagnostics!(options; settings, set_names)
else
    error("No diagnostics for dataset $dataset yet: add validation/$(lowercase(dataset)).jl and a branch here")
end
