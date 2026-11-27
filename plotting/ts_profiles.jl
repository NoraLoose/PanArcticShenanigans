# One or more runs' mean temperature and salinity profiles (ts_profiles_<set>.jld2, from analysis/ocean_diagnostics.jl)
# against those of the validation datasets in <validation_dir> (analysis/validation_diagnostics.jl), one panel per
# region, for each region set of validation.toml, all averaged over one period. Writes to `figure_dir` (for several
# runs, to a folder compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   temperature_profiles_<set>.png, salinity_profiles_<set>.png
#       the mean over the period of each run's averages and of GLORYS's daily means, and of the WOA23 decades whose
#       years overlap the period: their monthly climatologies of the period's months, weighted by the period's days in
#       each month; the observed annual-mean profiles (TSAR) of the period's years, and UDASH's mean of the profiles in
#       the region over the period's months (with the largest number of profiles on a level; levels without any left out)
# The period is --from up to, not including, --to (dates); by default all of the runs' days. A dataset without profiles
# for the set or the period is left out (with a warning).
#
#   julia +1.13 --project=. analysis/ocean_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/ts_profiles.jl [config.toml]
#   julia +1.13 --project=. plotting/ts_profiles.jl config.toml --from 1993-01-01 --to 1993-02-01
#   julia +1.13 --project=. plotting/ts_profiles.jl config_a.toml config_b.toml ...   (several runs)
#   julia +1.13 --project=. plotting/ts_profiles.jl config_a.toml config_b.toml --colors red,royalblue
#
# See plotting/runs.jl for the runs and colours.

include("runs.jl")

# The mean over the last dimension, weighted by `w`, of the finite values only (NaN where there are none)
function finite_mean(a, w = ones(size(a, ndims(a))))
    w = reshape(w, ntuple(_ -> 1, ndims(a) - 1)..., :)
    weight = sum(isfinite.(a) .* w; dims = ndims(a))
    total  = sum(ifelse.(isfinite.(a), a, zero(eltype(a))) .* w; dims = ndims(a))
    return dropdims(ifelse.(weight .> 0, total ./ max.(weight, eps()), NaN); dims = ndims(a))
end

# The days of the period in each calendar month (1–12)
function days_per_month(from, to)
    days = zeros(12)
    for day in Date(from):Day(1):Date(to - Millisecond(1))
        days[month(day)] += 1
    end
    return days
end

for set in region_set_names(validation)
    runs = read_runs("ts_profiles", set)
    run  = first(runs)
    all(r.regions_codes == run.regions_codes for r in runs) || error("The runs' profiles ($set) are of different regions")
    codes, names = run.regions_codes, run.regions_names
    R = length(codes)

    # The period: --from, --to, or all of the runs' days
    spans = run_spans(runs)
    from = DateTime(something(command_line_option("from"), string(minimum(first.(spans)))))
    to   = DateTime(something(command_line_option("to"),   string(maximum(last.(spans)))))
    from < to || error("--from must be before --to")
    period = "$(Date(from)) – $(Date(to - Millisecond(1)))"

    # Each run's averages that lie within the period
    within(r) = findall(t -> from ≤ t && t + r.window ≤ to, r.dates)
    for r in runs
        isempty(within(r)) && @warn "$(r.experiment) has no averages in $period"
    end

    glorys_periods = read_dataset("GLORYS", "ts_profiles", set)
    glorys = isnothing(glorys_periods) ? nothing : join_days(glorys_periods, from, to, (:temperature, :salinity))
    isnothing(glorys) && @warn "No GLORYS profiles ($set) in $period: analysis/validation_diagnostics.jl GLORYS --run <config>"

    # The WOA23 decades whose years overlap the period, with the period's months
    years  = year(from):year(to - Millisecond(1))
    weights = days_per_month(from, to)
    months = findall(>(0), weights)
    woa = filter(d -> any(in(first(d.years):last(d.years)), years), something(read_dataset("WOA23", "ts_profiles", set), []))
    isempty(woa) && @warn "No WOA23 profiles ($set) for the years $(first(years))–$(last(years)): analysis/validation_diagnostics.jl WOA23"
    woa_mean(d, field) = finite_mean(stack(getfield(d, field)[:, :, m] for m in months), weights[months])
    month_label = join(Dates.monthabbr.(months), ", ")

    tsar = read_dataset("TSAR", "ts_profiles", set)
    udash_periods = read_dataset("UDASH", "ts_profiles", set)

    # UDASH's monthly mean profiles over the period's months, weighted by their numbers of profiles: (; temperature,
    # salinity, n_temperature, n_salinity, level_depths, label) with the profiles' count per level, or nothing
    function udash_over_period()
        isnothing(udash_periods) && return nothing
        chosen = [(p, m) for p in udash_periods for (m, d) in enumerate(p.dates)
                  if from < DateTime(d) + Month(1) && DateTime(d) < to]
        isempty(chosen) && return nothing
        function combine(field, count)   # each month's mean weighted by its number of profiles
            n = sum(getfield(p, count)[:, :, m] for (p, m) in chosen)
            total = sum(replace(getfield(p, field)[:, :, m], NaN => 0.0) .* getfield(p, count)[:, :, m] for (p, m) in chosen)
            return ifelse.(n .> 0, total ./ max.(n, 1), NaN), n
        end
        θ, nθ = combine(:temperature, :n_temperature)
        S, nS = combine(:salinity, :n_salinity)
        label = join(unique(Dates.format(DateTime(p.dates[m]), "u yyyy") for (p, m) in chosen), ", ")
        return (; temperature = θ, salinity = S, n_temperature = nθ, n_salinity = nS,
                  level_depths = first(udash_periods).level_depths, label)
    end
    udash = udash_over_period()
    isnothing(udash) && @info "No UDASH profiles ($set) in $period: analysis/validation_diagnostics.jl UDASH --run <config>"

    for (name, label, field) in (("temperature", "temperature (°C)", :temperature), ("salinity", "salinity (g/kg)", :salinity))
        fig = Figure(size = figure_size(R) .+ (0, 150))
        legend = Pair{String, Any}[]   # the figure's legend: every line drawn in any panel, once
        entry!(label, element) = any(==(label) ∘ first, legend) || push!(legend, label => element)
        for (r, code) in enumerate(codes)
            ax = Axis(panel(fig, r, R); title = "$code  $(names[r])", xlabel = label, ylabel = "depth (m)", yreversed = true)
            for (k, d) in enumerate(woa)
                lines!(ax, woa_mean(d, field)[:, r], d.level_depths; color = woa_colors[k], linewidth = 1.5)
                entry!("$(woa_label(d)), $month_label", LineElement(color = woa_colors[k], linewidth = 1.5))
            end
            if !isnothing(glorys)
                lines!(ax, finite_mean(getfield(glorys, field))[:, r], glorys.level_depths; color = :black)
                entry!("GLORYS", LineElement(color = :black))
            end
            if !isnothing(tsar)
                for p in tsar, (n, y) in enumerate(p.dates)
                    parse(Int, y) in years || continue
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
            for (k, run) in enumerate(runs)
                n = within(run)
                isempty(n) && continue
                lines!(ax, finite_mean(getfield(run, field)[:, r:r, n])[:, 1], run.level_depths; color = run_colors[k])
                entry!(length(runs) == 1 ? "model" : run.experiment, LineElement(color = run_colors[k]))
            end
        end
        Legend(fig[cld(R, columns(R)) + 1, :], last.(legend), first.(legend); orientation = :horizontal, nbanks = 2,
               framevisible = false, labelsize = 11)
        Label(fig[0, :], "$title, $set regions: mean $label, $period"; fontsize = 18)
        file = joinpath(figures, "$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
