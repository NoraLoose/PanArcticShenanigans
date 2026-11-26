# A run's regional diagnostics (arctic_regional_diagnostics_<set>.jld2, from analysis/regional_diagnostics.jl) against
# those of the validation datasets in <validation_dir> (analysis/validation_diagnostics.jl), one panel per region, for
# each region set of validation.toml. Writes to `figure_dir`, for each set (e.g. _pan, _tsar):
#   regional_heat_content_<set>.png, regional_freshwater_content_<set>.png, regional_salt_content_<set>.png
#       heat content, freshwater content down to the isohaline, and salt content; heat and salt over the whole water
#       column, or down to `content_depth`. The run's averages are drawn at the middle of their windows, GLORYS's daily
#       means at noon of their day, and each WOA23 decade's monthly climatology as a step over every month of the years
#       it stands for (seasonal for the decades without monthly ones). The regions of the different grids differ
#       slightly in volume (coasts, sea floor), so their contents differ by an offset as well.
#   regional_temperature_profiles_<set>.png, regional_salinity_profiles_<set>.png
#       the mean profiles of the run and GLORYS on the run's first and last day, and the annual climatology of each WOA23
#       decade
# A dataset without diagnostics for the set, or for the run's days, is left out (with a warning).
#
#   julia +1.13 --project=. analysis/regional_diagnostics.jl [config.toml]
#   julia +1.13 --project=. plotting/regional_diagnostics.jl [config.toml]

include("../config.jl")
include("../validation_config.jl")

using CairoMakie, Dates, JLD2

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

# The run's diagnostics of a set, with its dates as DateTime
read_run(set) = (d = read_entries(regional_output_file(set)); merge(d, (; dates = DateTime.(d.dates))))

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

# Daily periods (GLORYS) within the run, joined along time: dates and the diagnostics with time last
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
    run = read_run(set)
    codes, names = run.regions_codes, run.regions_names
    R = length(codes)

    start  = DateTime(config["time"]["start_date"])
    window = Millisecond(round(Int, run.average_hours * 3600_000))
    days(dates) = Dates.value.(dates .- start) ./ 86400_000
    first_day, last_day = minimum(run.dates), maximum(run.dates) + window   # the span of the run's averages

    glorys_periods = read_dataset("GLORYS", set, run)
    glorys = isnothing(glorys_periods) ? nothing : join_days(glorys_periods, first_day, last_day)
    isnothing(glorys) && @warn "No GLORYS diagnostics ($set) for the run's days: analysis/validation_diagnostics.jl GLORYS --run $config_file"
    woa  = something(read_dataset("WOA23", set, run), [])
    isempty(woa) && @warn "No WOA23 diagnostics ($set): analysis/validation_diagnostics.jl WOA23"
    woa_colors = Makie.wong_colors()[[1, 3, 4, 6, 2, 5]]

    # The run's first and last day, and the GLORYS day of each
    run_days = Date.(run.dates)
    glorys_index(day) = isnothing(glorys) ? nothing : findfirst(==(day), Date.(glorys.dates))

    bound(d) = findfirst(==(content_depth), d.depth_bounds)
    isnothing(bound(run)) && error("content_depth = $content_depth is not one of the depth bounds $(run.depth_bounds)")
    depth_label = isinf(content_depth) ? "whole water column" : "0–$(round(Int, content_depth)) m"
    contents = [("heat",       "heat content, $depth_label (ZJ)", d -> d.heat[bound(run), :, :] ./ 1e21),
                ("freshwater", "freshwater content down to $(run.reference_salinity) g/kg (km³)",
                                                                   d -> d.freshwater_isohaline ./ 1e9),
                ("salt",       "salt content, $depth_label (Gt)",  d -> d.salt[bound(run), :, :] ./ 1e12)]

    # A WOA decade's months through the years it stands for, within the run: (segment ends in days, month)
    function woa_steps(d)
        steps = Tuple{Float64, Float64, Int}[]
        for year in first(d.years):last(d.years), month in 1:12
            a, b = max(DateTime(year, month), first_day), min(DateTime(year, month) + Month(1), last_day)
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
            lines!(ax, days(run.dates .+ window / 2), value(run)[r, :]; color = :red, label = "model")
            isnothing(glorys) ||
                scatter!(ax, days(glorys.dates .+ Hour(12)), value(glorys)[r, :]; color = :black, marker = :diamond,
                         label = "GLORYS")
            r == 1 && axislegend(ax; position = :lb, framevisible = false, merge = true, labelsize = 10)
        end
        Label(fig[0, :], "$experiment, $set regions: $label"; fontsize = 18)
        file = joinpath(figure_dir, "regional_$(name)_content_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end

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
            for (day, style) in ((first(run_days), :dash), (last(run_days), :solid))
                lines!(ax, getfield(run, field)[:, r, findfirst(==(day), run_days)], run.level_depths;
                       color = :red, linestyle = style)
                entry!("model $day", LineElement(color = :red, linestyle = style))
                g = glorys_index(day)
                isnothing(g) && continue
                lines!(ax, getfield(glorys, field)[:, r, g], glorys.level_depths; color = :black, linestyle = style)
                entry!("GLORYS $day", LineElement(color = :black, linestyle = style))
            end
        end
        Legend(fig[cld(R, columns(R)) + 1, :], last.(legend), first.(legend); orientation = :horizontal, nbanks = 2,
               framevisible = false, labelsize = 11)
        Label(fig[0, :], "$experiment, $set regions: mean $label"; fontsize = 18)
        file = joinpath(figure_dir, "regional_$(name)_profiles_$(lowercase(set)).png")
        save(file, fig)
        @info "Wrote $file"
    end
end
