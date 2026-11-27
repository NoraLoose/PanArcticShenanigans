# One or more runs' heat, freshwater and salt content (ocean_content_<set>.jld2, from analysis/ocean_diagnostics.jl)
# against those of the validation datasets in <validation_dir> (analysis/validation_diagnostics.jl), one panel per
# region, for each region set of validation.toml. Writes to `figure_dir` (for several runs, to a folder
# compare_<experiment>_<experiment>... in the figure folder), for each set (e.g. _pan, _tsar):
#   heat_content_<set>.png, freshwater_content_<set>.png, salt_content_<set>.png
#       heat and salt content over the whole water column, or down to `content_depth`, and freshwater content down to
#       the isohaline. Each run's averages are drawn at the middle of their windows, GLORYS's daily means at noon of
#       their day, and each WOA23 decade's monthly climatology as a step over every month of the years it stands for
#       (seasonal for the decades without monthly ones). The regions of the different grids differ slightly in volume
#       (coasts, sea floor), so their contents differ by an offset as well.
# A dataset without contents for the set, or for the runs' days, is left out (with a warning).
#
#   julia +1.13 --project=. analysis/ocean_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/ocean_content.jl [config.toml]
#   julia +1.13 --project=. plotting/ocean_content.jl config_a.toml config_b.toml ...   (several runs)
#   julia +1.13 --project=. plotting/ocean_content.jl config_a.toml config_b.toml --colors red,royalblue
#
# The time axis is days since the earliest start (see plotting/runs.jl for the runs and colours).

include("runs.jl")

content_depth = Inf   # heat and salt content from the surface down to this depth (m): one of validation.toml's
                      # depth_bounds, or Inf for the whole water column

settings_of(d) = (d.reference_temperature, d.reference_salinity, d.depth_bounds)

for set in region_set_names(validation)
    runs = read_runs("ocean_content", set)
    run  = first(runs)
    all(settings_of(r) == settings_of(run) && r.regions_codes == run.regions_codes for r in runs) ||
        error("The runs' contents ($set) were computed with different settings or regions")
    codes, names = run.regions_codes, run.regions_names
    R = length(codes)

    start = minimum(r.start for r in runs)   # the time axis: days since the earliest start
    days(dates) = Dates.value.(dates .- start) ./ 86400_000
    spans = run_spans(runs)
    first_day, last_day = minimum(first.(spans)), maximum(last.(spans))

    same_settings(d) = settings_of(d) == settings_of(run)
    glorys_periods = read_dataset("GLORYS", "ocean_content", set; same_settings)
    glorys = isnothing(glorys_periods) ? nothing :
             join_days(glorys_periods, first_day, last_day, (:heat, :salt, :freshwater_isohaline))
    isnothing(glorys) && @warn "No GLORYS contents ($set) for the runs' days: analysis/validation_diagnostics.jl GLORYS --run <config>"
    woa = something(read_dataset("WOA23", "ocean_content", set; same_settings), [])
    isempty(woa) && @warn "No WOA23 contents ($set): analysis/validation_diagnostics.jl WOA23"

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
        file = joinpath(figures, "$(name)_content_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
