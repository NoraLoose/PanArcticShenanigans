# One or more runs' mean temperature and salinity profiles (ts_profiles_<set>.jld2, from analysis/ocean_diagnostics.jl)
# against those of the validation datasets in <validation_dir> (analysis/validation_diagnostics.jl), one panel per
# region, for each region set of validation.toml. Writes to `figure_dir` (for several runs, to a folder
# compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   temperature_profiles_<set>.png, salinity_profiles_<set>.png
#       the profiles of the run and GLORYS on the run's first and last day (with several runs, each run's last day),
#       the annual climatology of each WOA23 decade, the observed annual-mean profiles (TSAR) of the runs' years, and
#       UDASH's mean of the profiles in the region over the runs' months (with the largest number of profiles on a level;
#       levels without any left out)
# A dataset without profiles for the set, or for the runs' days, is left out (with a warning).
#
#   julia +1.13 --project=. analysis/ocean_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/ts_profiles.jl [config.toml]
#   julia +1.13 --project=. plotting/ts_profiles.jl config_a.toml config_b.toml ...   (several runs)
#   julia +1.13 --project=. plotting/ts_profiles.jl config_a.toml config_b.toml --colors red,royalblue
#
# See plotting/runs.jl for the runs and colours.

include("runs.jl")

for set in region_set_names(validation)
    runs = read_runs("ts_profiles", set)
    run  = first(runs)
    all(r.regions_codes == run.regions_codes for r in runs) || error("The runs' profiles ($set) are of different regions")
    codes, names = run.regions_codes, run.regions_names
    R = length(codes)

    spans = run_spans(runs)
    first_day, last_day = minimum(first.(spans)), maximum(last.(spans))
    glorys_periods = read_dataset("GLORYS", "ts_profiles", set)
    glorys = isnothing(glorys_periods) ? nothing : join_days(glorys_periods, first_day, last_day, (:temperature, :salinity))
    isnothing(glorys) && @warn "No GLORYS profiles ($set) for the runs' days: analysis/validation_diagnostics.jl GLORYS --run <config>"
    woa = something(read_dataset("WOA23", "ts_profiles", set), [])
    isempty(woa) && @warn "No WOA23 profiles ($set): analysis/validation_diagnostics.jl WOA23"
    glorys_index(day) = isnothing(glorys) ? nothing : findfirst(==(day), Date.(glorys.dates))
    tsar = read_dataset("TSAR", "ts_profiles", set)
    udash_periods = read_dataset("UDASH", "ts_profiles", set)
    run_years = unique(y for (from, to) in spans for y in year(from):year(to - Millisecond(1)))   # the runs' years only

    # UDASH's monthly mean profiles over the runs' months, weighted by their numbers of profiles: (; temperature,
    # salinity, n_temperature, n_salinity, level_depths, label) with the profiles' count per level, or nothing
    function udash_over_runs()
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
    udash = udash_over_runs()
    isnothing(udash) && @info "No UDASH profiles ($set) for the runs' months: analysis/validation_diagnostics.jl UDASH --run <config>"

    # The days whose profiles are drawn: one run's first and last day; with several runs, each run's last day
    profile_days(r) = length(runs) == 1 ? ((Date(first(r.dates)), :dash), (Date(last(r.dates)), :solid)) :
                                          ((Date(last(r.dates)), :solid),)
    run_label(r, day) = length(runs) == 1 ? "model $day" : "$(r.experiment) $day"

    for (name, label, field) in (("temperature", "temperature (°C)", :temperature), ("salinity", "salinity (g/kg)", :salinity))
        fig = Figure(size = figure_size(R) .+ (0, 150))
        legend = Pair{String, Any}[]   # the figure's legend: every line drawn in any panel, once
        entry!(label, element) = any(==(label) ∘ first, legend) || push!(legend, label => element)
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = label, ylabel = "depth (m)", yreversed = true)
            for (k, d) in enumerate(woa)
                lines!(ax, getfield(d, Symbol("annual_$field"))[:, r], d.level_depths; color = woa_colors[k], linewidth = 1.5)
                entry!(woa_label(d), LineElement(color = woa_colors[k], linewidth = 1.5))
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
        end
        Legend(fig[cld(R, columns(R)) + 1, :], last.(legend), first.(legend); orientation = :horizontal, nbanks = 2,
               framevisible = false, labelsize = 11)
        Label(fig[0, :], "$title, $set regions: mean $label"; fontsize = 18)
        file = joinpath(figures, "$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
