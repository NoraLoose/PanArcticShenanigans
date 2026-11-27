# The runs of a plotting script of the regional diagnostics (plotting/ocean_content.jl, plotting/ts_profiles.jl): every
# config file on the command line, their colours (--colors red,royalblue: Makie's colour names or hex like "#1f77b4", in
# the order of the config files; otherwise a palette apart from WOA23's colours) and the figure folder (the run's own, or for several runs a folder
# compare_<experiment>_<experiment>... in the figure folder). The runs are told apart by their experiment names.

include("../config.jl")
include("../validation_config.jl")

using CairoMakie, Dates, JLD2, TOML

config_files = isempty(ARGS) ? [config_file] : filter(endswith(".toml"), ARGS)   # config.jl reads the first
configs      = [TOML.parsefile(file) for file in config_files]
experiments  = [get(c, "experiment", splitext(basename(f))[1]) for (c, f) in zip(configs, config_files)]
allunique(experiments) || error("The runs need different experiment names: $experiments")
colors_option = command_line_option("colors")
run_colors   = isnothing(colors_option) ? [:red, :darkgoldenrod, :saddlebrown, :darkcyan, :deeppink, :olive] :
                                          [Makie.to_color(String(strip(c))) for c in split(colors_option, ",")]
length(configs) ≤ length(run_colors) ||
    error("$(length(configs)) runs, but $(length(run_colors)) colours: give one per run with --colors")

figures = length(configs) == 1 ? figure_dir : mkpath(joinpath(paths["figure_dir"], "compare_" * join(experiments, "_")))
title   = join(experiments, ", ")

# A run's diagnostics file `name`_<set>.jld2 (e.g. "ocean_content"), with its dates as DateTime, its start, and its
# averaging window
function read_run(c, experiment, name, set)
    base = isempty(get(paths, "diagnostics_dir", "")) ? paths["output_dir"] : paths["diagnostics_dir"]   # as config.jl
    file = joinpath(base, get(c, "experiment", ""), "$(name)_$(lowercase(set)).jld2")
    isfile(file) || error("No $file: run analysis/ocean_diagnostics.jl for $experiment first")
    d = read_entries(file)
    return merge(d, (; dates = DateTime.(d.dates), start = DateTime(c["time"]["start_date"]),
                       window = Millisecond(round(Int, d.average_hours * 3600_000)), experiment))
end

read_runs(name, set) = [read_run(c, e, name, set) for (c, e) in zip(configs, experiments)]

# One panel per region, up to 5 in a row
columns(n) = min(n, 5)
panel(fig, r, n) = fig[fldmod1(r, columns(n))...]
figure_size(n) = (300 * columns(n), 280 * cld(n, columns(n)) + 60)

# All entries of a JLD2 file (groups flattened, "regions/codes" → regions_codes) as a NamedTuple
function read_entries(path)
    entries = Pair{Symbol, Any}[]
    jldopen(path) do file
        function walk(group, prefix)
            for key in keys(group)
                value = group[key]
                name = isempty(prefix) ? key : "$(prefix)_$key"
                value isa JLD2.Group ? walk(value, name) : push!(entries, Symbol(name) => value)
            end
        end
        walk(file, "")
    end
    return (; entries...)
end

# ## The validation datasets' diagnostics (analysis/validation_diagnostics.jl)

# A dataset's diagnostics `family` of a set: its regions and settings (regions.jld2), and each period's diagnostics, with
# the period's own regions and levels where it has them; `nothing` without any. `same_settings(stored)` checks the
# settings against the runs'.
function read_dataset(dataset, family, set; same_settings = d -> true)
    haskey(validation["datasets"], dataset) || return nothing
    isfile(regions_file(dataset, family, set)) || return nothing
    base = read_entries(regions_file(dataset, family, set))
    same_settings(base) ||
        @warn "$dataset's $family ($set) was computed with other settings than the runs': compute it again"
    periods = [merge(base, read_entries(period_file(dataset, family, set, splitext(f)[1])))
               for f in sort(period_files(dataset, family, set))]
    return isempty(periods) ? nothing : periods
end

# Daily periods (GLORYS) within [first_day, last_day), joined along time: dates, levels and the fields `names`, time last
function join_days(periods, first_day, last_day, names)
    days = filter(p -> first_day ≤ DateTime(only(p.dates)) < last_day, periods)
    isempty(days) && return nothing
    along_time(name) = cat((getfield(p, name) for p in days)...; dims = ndims(getfield(first(days), name)))
    level_depths = hasproperty(first(days), :level_depths) ? first(days).level_depths : nothing
    return (; dates = [DateTime(only(p.dates)) for p in days], level_depths, (name => along_time(name) for name in names)...)
end

# Each run's span (its first average's start to its last one's end), and that of all runs
run_spans(runs) = [(minimum(r.dates), maximum(r.dates) + r.window) for r in runs]

woa_colors = Makie.wong_colors()[[1, 3, 4, 6, 2, 5]]
woa_label(d) = "WOA23 $(first(d.years))–$(last(d.years))" * (d.monthly ? "" : " (seasonal)")
