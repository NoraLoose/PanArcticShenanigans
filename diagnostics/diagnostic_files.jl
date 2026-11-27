# What the regional diagnostics (ocean_content.jl, ts_profiles.jl) have in common: the regions as an index per column,
# the model levels, and their files. A file holds the regions (set, codes, names, masks, areas), the levels and the
# settings, and the diagnostics of a list of times (`results`, one NamedTuple of arrays per time), each array with time as
# its last dimension, with `dates` (strings) and any `extra` entries.

@isdefined(region_masks) || include("regions.jl")

using Oceananigans.Grids: inactive_node
using Oceananigans.Operators: Δzᶜᶜᶜ, Azᶜᶜᶜ
using JLD2

# The region of each column of `grid`: its index in `regions.codes`, 0 outside every region
function region_index(grid, regions)
    Nx, Ny = size(grid)[1:2]
    region = zeros(Int, Nx, Ny)
    for (r, code) in enumerate(regions.codes)
        region[regions.masks[code]] .= r
    end
    return region
end

# Mean depth (m, positive) of each level of `grid` below the surface at rest, from the bottom up
level_depths(grid) = -collect(znodes(grid.underlying_grid, Center()))

function write_regions!(file, regions, grid; level_depths = nothing)
    if !isnothing(grid)
        Nx, Ny = size(grid)[1:2]
        area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
        file["regions/area"] = [sum(area[regions.masks[code]]) for code in regions.codes]   # m²
        for code in regions.codes
            file["regions/masks/$code"] = Array(regions.masks[code])
        end
    end
    file["regions/set"]   = regions.set
    file["regions/codes"] = regions.codes
    file["regions/names"] = regions.names
    isnothing(level_depths) || (file["level_depths"] = level_depths)   # m, from the bottom up, as the profiles' levels
end

function write_results!(file, results; dates, extra...)
    file["dates"] = string.(dates)
    for name in keys(first(results))
        file[string(name)] = stack(getfield.(results, name))   # time last
    end
    for (name, value) in extra
        file[string(name)] = value
    end
end

"""
    save_diagnostics(path, results, regions, grid; dates, level_depths = nothing, settings = (;), extra...)

A run's diagnostics of one region set: regions, levels, settings (written as entries of their own) and results in the
JLD2 file `path`.
"""
function save_diagnostics(path, results, regions, grid; dates, level_depths = nothing, settings = (;), extra...)
    jldopen(path, "w") do file
        write_regions!(file, regions, grid; level_depths)
        foreach(((name, value),) -> file[string(name)] = value, pairs(settings))
        write_results!(file, results; dates, extra...)
    end
    @info "Wrote $path"
end

# ## Validation datasets (analysis/validation_diagnostics.jl; paths in validation_config.jl)
#
# A dataset keeps, per family of diagnostics and region set, its regions once (save_dataset_regions) and its
# diagnostics in a file per period (save_dataset_period).

# The families computed from gridded T and S, and what each needs: its function, and the settings written with it
const ocean_families = ("ocean_content", "ts_profiles")

ocean_family(family, grid, T, S, regions; settings) =
    family == "ocean_content" ? ocean_content(grid, T, S, regions; settings...) : ts_profiles(grid, T, S, regions)

family_file_settings(family, settings) = family == "ocean_content" ? content_file_settings(settings) : (;)

"""
    save_dataset_regions(dataset, family, regions, grid; level_depths = nothing, settings = (;))

A validation dataset's regions of one set for one family (regions.jld2), with the family's `settings`, unless there
are already some made with the same settings; made with other settings, the dataset's diagnostics of that family and
set are out of date, and are removed. `grid = nothing` for a dataset whose cells differ between periods (its masks,
areas and levels are in each period's file) or that has no cells (profiles that are already regional means).
"""
function save_dataset_regions(dataset, family, regions, grid; level_depths = nothing, settings = (;))
    path = regions_file(dataset, family, regions.set)
    if isfile(path)
        stored = jldopen(file -> (; (name => (haskey(file, string(name)) ? file[string(name)] : nothing)
                                         for name in keys(settings))...), path)
        stored == settings && return
        @warn "The settings of $dataset's $family ($(regions.set)) have changed: computing them again"
        dir = diagnostics_dir(dataset, family, regions.set)
        foreach(f -> rm(joinpath(dir, f)), readdir(dir))
    end
    jldopen(path, "w") do file
        write_regions!(file, regions, grid; level_depths)
        foreach(((name, value),) -> file[string(name)] = value, pairs(settings))
    end
end

"""
    save_dataset_period(dataset, family, set, period, results; dates, regions = nothing, grid = nothing,
                        level_depths = nothing, extra...)

A validation dataset's diagnostics of one family and region set over one period (a day, a decade, ...):
`<period>.jld2`. A dataset whose cells differ between periods (e.g. WOA's sea floor, between decades) gives the
period's `regions`, `grid` and `level_depths`, which are then kept in the period's file instead of regions.jld2.
"""
function save_dataset_period(dataset, family, set, period, results; dates, regions = nothing, grid = nothing,
                             level_depths = nothing, extra...)
    path = period_file(dataset, family, set, period)
    jldopen(path, "w") do file
        if !isnothing(regions)
            Nx, Ny = size(grid)[1:2]
            area = [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny]
            file["regions/area"] = [sum(area[regions.masks[code]]) for code in regions.codes]   # m²
            for code in regions.codes
                file["regions/masks/$code"] = Array(regions.masks[code])
            end
            file["level_depths"] = level_depths
        end
        write_results!(file, results; dates, extra...)
    end
    @info "Wrote $path"
end
