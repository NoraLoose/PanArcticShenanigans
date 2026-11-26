# Heat, salt and freshwater content and mean T, S profiles of the regions of regions.jl (see regional_diagnostics.jl),
# computed offline from the time-averaged output of panarctic.jl.
#
#   julia +1.13 --project=. analysis/regional_diagnostics.jl [config.toml]
#
# Reads arctic_ocean.jld2 from `output_dir` and writes arctic_regional_diagnostics_<set>.jld2 for each region set (e.g.
# _pan, _tsar) to `diagnostics_dir` (paths.toml; by default output_dir), with "dates", the start of each averaging
# window. Runs on the CPU and reads one average at a time (the first output, the initial state, is left out). Uses the
# grid saved with the output, so the regions are those of the run's grid; the cell thicknesses follow the averaged free
# surface. The settings and region sets are those of validation.toml [regional_diagnostics], the same for every run.
#
# Heat and salt content and the profiles are linear in T and S, so the averages give the averages of the diagnostics
# (up to the covariance of T and S with the cell thicknesses within each average, which is tiny). The freshwater
# content down to the isohaline is not quite linear: where the isohaline lies is found from the averaged S.

include("../config.jl")
include("../validation_config.jl")
include("../regional_diagnostics.jl")

using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_zstar_scaling!
using Dates, Printf

settings = regional_settings(validation)
saved = Dict(name => FieldTimeSeries(ocean_output_file, name; backend = OnDisk()) for name in ("T", "S", "η"))
grid  = saved["T"].grid
times = saved["T"].times
averages = findall(>(0), times)   # the first output, at time 0, is the initial state, not an average
depths = level_depths(grid)   # before the free surface moves the levels
sets   = [region_masks(grid; set) for set in region_set_names(validation)]

# For each average, the diagnostics of every region set
results = map(averages) do n
    update_zstar_scaling!(grid, saved["η"][n])
    T, S = saved["T"][n], saved["S"][n]
    result = [regional_diagnostics(grid, T, S, regions; settings...) for regions in sets]
    @info @sprintf("%d / %d: %s", n, last(averages), prettytime(times[n]))
    result
end

# Each average covers the `average_hours` before its time; model time 0 is 00:00 on start_date
average = Millisecond(round(Int, config["output"]["average_hours"] * 3600_000))
dates   = [DateTime(config["time"]["start_date"]) + Millisecond(round(Int, 1000t)) - average for t in times[averages]]

for (s, regions) in enumerate(sets)
    save_regional_diagnostics(regional_output_file(regions.set), getindex.(results, s), regions, grid; dates,
                              level_depths = depths, settings, times = times[averages],
                              average_hours = config["output"]["average_hours"])
end
