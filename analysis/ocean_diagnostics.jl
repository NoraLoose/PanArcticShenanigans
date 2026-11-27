# Heat, salt and freshwater content (ocean_content.jl) and mean T, S profiles (ts_profiles.jl) of the regions of
# regions.jl, computed offline from the time-averaged output of panarctic.jl. Each average of T and S is read once for
# both.
#
#   julia +1.13 --project=. analysis/ocean_diagnostics.jl [config.toml]
#   julia +1.13 --project=. analysis/ocean_diagnostics.jl config.toml --only content    (or --only profiles)
#
# Reads arctic_ocean.jld2 from `output_dir` and writes, for each region set of validation.toml (e.g. _pan, _tsar),
# ocean_content_<set>.jld2 and ts_profiles_<set>.jld2 to `diagnostics_dir` (paths.toml; by default output_dir), with
# "dates", the start of each averaging window. Runs on the CPU and reads one average at a time (the first output, the
# initial state, is left out). Uses the grid saved with the output, so the regions are those of the run's grid; the cell
# thicknesses follow the averaged free surface. The settings and region sets are those of validation.toml, the same for
# every run.
#
# Heat and salt content and the profiles are linear in T and S, so the averages give the averages of the diagnostics
# (up to the covariance of T and S with the cell thicknesses within each average, which is tiny). The freshwater
# content down to the isohaline is not quite linear: where the isohaline lies is found from the averaged S.

include("../config.jl")
include("../validation_config.jl")
include("../diagnostics/ocean_content.jl")
include("../diagnostics/ts_profiles.jl")

using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_zstar_scaling!
using Dates, Printf

only_at = findfirst(==("--only"), ARGS)
families = isnothing(only_at) ? ["content", "profiles"] : [ARGS[only_at + 1]]
all(in(("content", "profiles")), families) || error("--only content or --only profiles, not $(join(families, ", "))")

settings = content_settings(validation)
saved = Dict(name => FieldTimeSeries(ocean_output_file, name; backend = OnDisk()) for name in ("T", "S", "η"))
grid  = saved["T"].grid
times = saved["T"].times
averages = findall(>(0), times)   # the first output, at time 0, is the initial state, not an average
depths = level_depths(grid)   # before the free surface moves the levels
sets   = [region_masks(grid; set) for set in region_set_names(validation)]

# For each average, the contents and profiles of every region set
results = map(averages) do n
    update_zstar_scaling!(grid, saved["η"][n])
    T, S = saved["T"][n], saved["S"][n]
    result = [(content  = "content"  in families ? ocean_content(grid, T, S, regions; settings...) : nothing,
               profiles = "profiles" in families ? ts_profiles(grid, T, S, regions) : nothing) for regions in sets]
    @info @sprintf("%d / %d: %s", n, last(averages), prettytime(times[n]))
    result
end

# Each average covers the `average_hours` before its time; model time 0 is 00:00 on start_date
average = Millisecond(round(Int, config["output"]["average_hours"] * 3600_000))
dates   = [DateTime(config["time"]["start_date"]) + Millisecond(round(Int, 1000t)) - average for t in times[averages]]
extra   = (; times = times[averages], average_hours = config["output"]["average_hours"])

for (s, regions) in enumerate(sets)
    "content" in families &&
        save_diagnostics(ocean_content_file(regions.set), [r[s].content for r in results], regions, grid; dates,
                         settings = content_file_settings(settings), extra...)
    "profiles" in families &&
        save_diagnostics(ts_profiles_file(regions.set), [r[s].profiles for r in results], regions, grid; dates,
                         level_depths = depths, extra...)
end
