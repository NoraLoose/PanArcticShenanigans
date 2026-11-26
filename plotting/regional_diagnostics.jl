# One or more runs' regional diagnostics (arctic_regional_diagnostics_<set>.jld2, from analysis/regional_diagnostics.jl)
# against those of the validation datasets in <validation_dir> (analysis/validation_diagnostics.jl), one panel per
# region, for each region set of validation.toml. Writes to `figure_dir` (for several runs, to a folder
# compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   regional_heat_content_<set>.png, regional_freshwater_content_<set>.png, regional_salt_content_<set>.png
#       heat content, freshwater content down to the isohaline, and salt content; heat and salt over the whole water
#       column, or down to `content_depth`. Each run's averages are drawn at the middle of their windows, GLORYS's daily
#       means at noon of their day, and each WOA23 decade's monthly climatology as a step over every month of the years
#       it stands for (seasonal for the decades without monthly ones). The regions of the different grids differ
#       slightly in volume (coasts, sea floor), so their contents differ by an offset as well.
#   regional_temperature_profiles_<set>.png, regional_salinity_profiles_<set>.png
#       the mean profiles of the run and GLORYS on the run's first and last day (with several runs, each run's last day),
#       the annual climatology of each WOA23 decade, the observed annual-mean profiles (TSAR) of the runs' years, and
#       UDASH's mean of the profiles in the region over the runs' months (with the largest number of profiles on a level;
#       levels without any left out)
# A dataset without diagnostics for the set, or for the run's days, is left out (with a warning).
#
#   julia +1.13 --project=. analysis/regional_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/regional_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/regional_diagnostics.jl config_a.toml config_b.toml ...   (several runs)
#
# The runs are told apart by their experiment names; their time axis is days since the earliest start.

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
run_colors   = [:red, :darkgoldenrod, :saddlebrown, :darkcyan, :deeppink, :olive]   # apart from WOA23's colours
length(configs) ≤ length(run_colors) || error("At most $(length(run_colors)) runs in one figure")

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

# A dataset's diagnostics of a set: its regions and settings (regions.jld2), and each period's diagnostics, with the
# period's own regions and levels where it has them; `nothing` without any
function read_dataset(dataset, set, run)
    haskey(validation["datasets"], dataset) || return nothing
    isfile(regions_file(dataset, set)) || return nothing
    base = read_entries(regions_file(dataset, set))
    settings_of(base) == settings_of(run) ||
        @warn "$dataset's diagnostics ($set) were computed with other settings than the run's: compute them again"
    periods = [merge(base, read_entries(period_file(dataset, set, splitext(f)[1]))) for f in sort(period_files(dataset, set))]
    return isempty(periods) ? nothing : periods
end

# Daily periods (GLORYS) within [first_day, last_day), joined along time: dates and the diagnostics with time last
function join_days(periods, first_day, last_day)
    days = filter(p -> first_day ≤ DateTime(only(p.dates)) < last_day, periods)
    isempty(days) && return nothing
    along_time(name) = cat((getfield(p, name) for p in days)...; dims = ndims(getfield(first(days), name)))
    return (; dates = [DateTime(only(p.dates)) for p in days], level_depths = first(days).level_depths,
              heat = along_time(:heat), salt = along_time(:salt), freshwater_isohaline = along_time(:freshwater_isohaline),
              temperature = along_time(:temperature), salinity = along_time(:salinity))
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
    spans = [(minimum(r.dates), maximum(r.dates) + r.window) for r in runs]   # each run's averages
    first_day, last_day = minimum(first.(spans)), maximum(last.(spans))

    glorys_periods = read_dataset("GLORYS", set, run)
    glorys = isnothing(glorys_periods) ? nothing : join_days(glorys_periods, first_day, last_day)
    isnothing(glorys) && @warn "No GLORYS diagnostics ($set) for the runs' days: analysis/validation_diagnostics.jl GLORYS --run <config>"
    woa  = something(read_dataset("WOA23", set, run), [])
    isempty(woa) && @warn "No WOA23 diagnostics ($set): analysis/validation_diagnostics.jl WOA23"
    tsar = read_dataset("TSAR", set, run)
    udash_periods = read_dataset("UDASH", set, run)

    # UDASH's monthly mean profiles over the run's months, weighted by their numbers of profiles: (; temperature,
    # salinity, n_temperature, n_salinity, level_depths, months) with the profiles' count per level, or nothing
    function udash_over_run()
        isnothing(udash_periods) && return nothing
        months = [(p, m) for p in udash_periods for (m, d) in enumerate(p.dates)
                  if any(from < DateTime(d) + Month(1) && DateTime(d) < to for (from, to) in spans)]
        isempty(months) && return nothing
        function combine(field, count)   # each month's mean weighted by its number of profiles
            n = sum(getfield(p, count)[:, :, m] for (p, m) in months)
            total = sum(replace(getfield(p, field)[:, :, m], NaN => 0.0) .* getfield(p, count)[:, :, m] for (p, m) in months)
            return ifelse.(n .> 0, total ./ max.(n, 1), NaN), n
        end
        θ, nθ = combine(:temperature, :n_temperature)
        S, nS = combine(:salinity, :n_salinity)
        label = join(unique(Dates.format(DateTime(p.dates[m]), "u yyyy") for (p, m) in months), ", ")
        return (; temperature = θ, salinity = S, n_temperature = nθ, n_salinity = nS,
                  level_depths = first(udash_periods).level_depths, label)
    end
    udash = udash_over_run()
    isnothing(udash) && @info "No UDASH profiles ($set) for the run's months: analysis/validation_diagnostics.jl UDASH --run $config_file"
    woa_colors = Makie.wong_colors()[[1, 3, 4, 6, 2, 5]]

    glorys_index(day) = isnothing(glorys) ? nothing : findfirst(==(day), Date.(glorys.dates))
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

    # A WOA decade's months through the years it stands for, within each run: (segment ends in days, month)
    function woa_steps(d)
        steps = Tuple{Float64, Float64, Int}[]
        for (from, to) in spans, year in first(d.years):last(d.years), month in 1:12
            a, b = max(DateTime(year, month), from), min(DateTime(year, month) + Month(1), to)
            a < b && push!(steps, (days(a), days(b), month))
        end
        return steps
    end
    woa_label(d) = "WOA23 $(first(d.years))–$(last(d.years))" * (d.monthly ? "" : " (seasonal)")

    for (name, label, value) in contents
        fig = Figure(size = figure_size(R))
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = "days since $(Date(start))")
            for (k, d) in enumerate(woa)
                steps = woa_steps(d)
                isempty(steps) && continue
                vw = value(d)
                linesegments!(ax, [Point2f(t, vw[r, m]) for (a, b, m) in steps for t in (a, b)];
                              color = woa_colors[k], linewidth = 3, label = woa_label(d))
            end
            for (k, run) in enumerate(runs)
                lines!(ax, days(run.dates .+ run.window / 2), value(run)[r, :]; color = run_colors[k],
                       label = length(runs) == 1 ? "model" : run.experiment)
            end
            isnothing(glorys) ||
                scatter!(ax, days(glorys.dates .+ Hour(12)), value(glorys)[r, :]; color = :black, marker = :diamond,
                         label = "GLORYS")
            r == 1 && axislegend(ax; position = :lb, framevisible = false, merge = true, labelsize = 10)
        end
        Label(fig[0, :], "$title, $set regions: $label"; fontsize = 18)
        file = joinpath(figures, "regional_$(name)_content_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end

    run_years = unique(y for (from, to) in spans for y in year(from):year(to - Millisecond(1)))   # the runs' years only
    for (name, label, field) in (("temperature", "temperature (°C)", :temperature), ("salinity", "salinity (g/kg)", :salinity))
        fig = Figure(size = figure_size(R) .+ (0, 150))
        legend = Pair{String, Any}[]   # the figure's legend: every line drawn in any panel, once
        function entry!(label, element)
            any(==(label) ∘ first, legend) || push!(legend, label => element)
        end
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = label, ylabel = "depth (m)", yreversed = true)
            for (k, d) in enumerate(woa)
                lines!(ax, getfield(d, Symbol("annual_$field"))[:, r], d.level_depths; color = woa_colors[k], linewidth = 1.5)
                entry!(woa_label(d), LineElement(color = woa_colors[k], linewidth = 1.5))
            end
            if !isnothing(tsar)
                for p in tsar, (n, y) in enumerate(p.dates)
                    parse(Int, y) in run_years || continue
                    lines!(ax, getfield(p, field)[:, r, n], p.level_depths; color = :darkgreen, linewidth = 2)
                    entry!("observed $y (TSAR)", LineElement(color = :darkgreen, linewidth = 2))
                end
            end
            if !isnothing(udash)
                n = getfield(udash, Symbol("n_$field"))[:, r]
                if any(n .> 0)
                    observed = n .> 0
                    scatterlines!(ax, getfield(udash, field)[observed, r], udash.level_depths[observed]; color = :purple,
                                  markersize = 5, linewidth = 1.5)
                    text!(ax, 0.98, 0.02; text = "UDASH: $(maximum(n)) profiles", space = :relative,
                          align = (:right, :bottom), fontsize = 10, color = :purple)
                    entry!("UDASH $(udash.label), mean of the profiles", [LineElement(color = :purple, linewidth = 1.5),
                                                                          MarkerElement(color = :purple, marker = :circle, markersize = 5)])
                end
            end
            for (k, run) in enumerate(runs), (day, style) in profile_days(run)
                lines!(ax, getfield(run, field)[:, r, findfirst(==(day), Date.(run.dates))], run.level_depths;
                       color = run_colors[k], linestyle = style)
                entry!(run_label(run, day), LineElement(color = run_colors[k], linestyle = style))
                g = glorys_index(day)
                isnothing(g) && continue
                lines!(ax, getfield(glorys, field)[:, r, g], glorys.level_depths; color = :black, linestyle = style)
                entry!("GLORYS $day", LineElement(color = :black, linestyle = style))
            end
        end
        Legend(fig[cld(R, columns(R)) + 1, :], last.(legend), first.(legend); orientation = :horizontal, nbanks = 2,
               framevisible = false, labelsize = 11)
        Label(fig[0, :], "$title, $set regions: mean $label"; fontsize = 18)
        file = joinpath(figures, "regional_$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
