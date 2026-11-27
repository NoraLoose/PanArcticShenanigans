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
