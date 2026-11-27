# One or more runs' mean temperature and salinity profiles (ts_profiles_<set>.jld2, from analysis/ocean_diagnostics.jl),
# one panel per region, for each region set of validation.toml. Writes to `figure_dir` (for several runs, to a folder
# compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   temperature_profiles_<set>.png, salinity_profiles_<set>.png
#       the profiles on the run's first and last day (with several runs, each run's last day)
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

    # The days whose profiles are drawn: one run's first and last day; with several runs, each run's last day
    profile_days(r) = length(runs) == 1 ? ((Date(first(r.dates)), :dash), (Date(last(r.dates)), :solid)) :
                                          ((Date(last(r.dates)), :solid),)
    run_label(r, day) = length(runs) == 1 ? "model $day" : "$(r.experiment) $day"

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
        file = joinpath(figures, "$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
