# One or more runs' regional diagnostics (arctic_regional_diagnostics_<set>.jld2, from analysis/regional_diagnostics.jl),
# one panel per region, for each region set of validation.toml. Writes to `figure_dir` (for several runs, to a folder
# compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   regional_heat_content_<set>.png, regional_freshwater_content_<set>.png, regional_salt_content_<set>.png
#       heat content, freshwater content down to the isohaline, and salt content over the run; heat and salt over the
#       whole water column, or down to `content_depth`. Each run's averages are drawn at the middle of their windows.
#   regional_temperature_profiles_<set>.png, regional_salinity_profiles_<set>.png
#       the mean profiles on the run's first and last day (with several runs, each run's last day)
#
#   julia +1.13 --project=. analysis/regional_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/regional_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/regional_diagnostics.jl config_a.toml config_b.toml ...   (several runs)
#   julia +1.13 --project=. plotting/regional_diagnostics.jl config_a.toml config_b.toml --colors red,royalblue
#
# The runs are told apart by their experiment names; their time axis is days since the earliest start. --colors gives
# their colours, in the order of the config files (Makie's colour names, or hex like "#1f77b4").

include("../config.jl")
include("../validation_config.jl")

using CairoMakie, Dates, JLD2, TOML

content_depth = Inf   # heat and salt content from the surface down to this depth (m): one of validation.toml's
                      # depth_bounds, or Inf for the whole water column

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

settings_of(d) = (d.reference_temperature, d.reference_salinity, d.depth_bounds)

# The runs: every config file on the command line (config.jl reads the first)
config_files = isempty(ARGS) ? [config_file] : filter(endswith(".toml"), ARGS)
configs      = [TOML.parsefile(file) for file in config_files]
experiments  = [get(c, "experiment", splitext(basename(f))[1]) for (c, f) in zip(configs, config_files)]
allunique(experiments) || error("The runs need different experiment names: $experiments")
colors_at    = findfirst(==("--colors"), ARGS)
run_colors   = isnothing(colors_at) ? [:red, :darkgoldenrod, :saddlebrown, :darkcyan, :deeppink, :olive] :
                                      [Makie.to_color(String(strip(c))) for c in split(ARGS[colors_at + 1], ",")]
length(configs) ≤ length(run_colors) ||
    error("$(length(configs)) runs, but $(length(run_colors)) colours: give one per run with --colors")

# Figures: the run's figure folder, or for several runs a folder named after them
figures = length(configs) == 1 ? figure_dir : mkpath(joinpath(paths["figure_dir"], "compare_" * join(experiments, "_")))
title   = join(experiments, ", ")

# A run's diagnostics of a set, with its dates as DateTime, its start, and its averaging window
function read_run(c, experiment, set)
    base = isempty(get(paths, "diagnostics_dir", "")) ? paths["output_dir"] : paths["diagnostics_dir"]   # as config.jl
    file = joinpath(base, get(c, "experiment", ""), "arctic_regional_diagnostics_$(lowercase(set)).jld2")
    isfile(file) || error("No $file: run analysis/regional_diagnostics.jl for $experiment first")
    d = read_entries(file)
    return merge(d, (; dates = DateTime.(d.dates), start = DateTime(c["time"]["start_date"]),
                       window = Millisecond(round(Int, d.average_hours * 3600_000)), experiment))
end

columns(n) = min(n, 5)
panel(fig, r, n) = fig[fldmod1(r, columns(n))...]
figure_size(n) = (300 * columns(n), 280 * cld(n, columns(n)) + 60)

for set in region_set_names(validation)
    runs = [read_run(c, e, set) for (c, e) in zip(configs, experiments)]
    run  = first(runs)
    all(settings_of(r) == settings_of(run) && r.regions_codes == run.regions_codes for r in runs) ||
        error("The runs' regional diagnostics ($set) were computed with different settings or regions")
    codes, names = run.regions_codes, run.regions_names
    R = length(codes)

    start = minimum(r.start for r in runs)   # the time axis: days since the earliest start
    days(dates) = Dates.value.(dates .- start) ./ 86400_000

    # The days whose profiles are drawn: one run's first and last day; with several runs, each run's last day
    profile_days(r) = length(runs) == 1 ? ((Date(first(r.dates)), :dash), (Date(last(r.dates)), :solid)) :
                                          ((Date(last(r.dates)), :solid),)
    run_label(r, day) = length(runs) == 1 ? "model $day" : "$(r.experiment) $day"

    bound(d) = findfirst(==(content_depth), d.depth_bounds)
    isnothing(bound(run)) && error("content_depth = $content_depth is not one of the depth bounds $(run.depth_bounds)")
    depth_label = isinf(content_depth) ? "whole water column" : "0–$(round(Int, content_depth)) m"
    contents = [("heat",       "heat content, $depth_label (ZJ)", d -> d.heat[bound(run), :, :] ./ 1e21),
                ("freshwater", "freshwater content down to $(run.reference_salinity) g/kg (km³)",
                                                                   d -> d.freshwater_isohaline ./ 1e9),
                ("salt",       "salt content, $depth_label (Gt)",  d -> d.salt[bound(run), :, :] ./ 1e12)]

    for (name, label, value) in contents
        fig = Figure(size = figure_size(R))
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = "days since $(Date(start))")
            for (k, run) in enumerate(runs)
                scatterlines!(ax, days(run.dates .+ run.window / 2), value(run)[r, :]; color = run_colors[k],
                              markersize = 5, label = length(runs) == 1 ? "model" : run.experiment)
            end
            r == 1 && length(runs) > 1 && axislegend(ax; position = :lb, framevisible = false, labelsize = 10)
        end
        Label(fig[0, :], "$title, $set regions: $label"; fontsize = 18)
        file = joinpath(figures, "regional_$(name)_content_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end

    for (name, label, field) in (("temperature", "temperature (°C)", :temperature), ("salinity", "salinity (g/kg)", :salinity))
        fig = Figure(size = figure_size(R) .+ (0, 100))
        legend = Pair{String, Any}[]   # the figure's legend: every line drawn in any panel, once
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = label, ylabel = "depth (m)", yreversed = true)
            for (k, run) in enumerate(runs), (day, style) in profile_days(run)
                lines!(ax, getfield(run, field)[:, r, findfirst(==(day), Date.(run.dates))], run.level_depths;
                       color = run_colors[k], linestyle = style)
                entry = run_label(run, day)
                any(==(entry) ∘ first, legend) || push!(legend, entry => LineElement(color = run_colors[k], linestyle = style))
            end
        end
        Legend(fig[cld(R, columns(R)) + 1, :], last.(legend), first.(legend); orientation = :horizontal,
               framevisible = false, labelsize = 11)
        Label(fig[0, :], "$title, $set regions: mean $label"; fontsize = 18)
        file = joinpath(figures, "regional_$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
