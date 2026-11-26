# The ocean regions of regions.jl on the model grid, one map per region set: regions_pan.png and regions_tsar.png in
# `figure_dir`.
#
#   julia +1.13 --project=. plotting/regions.jl [config.toml]
#
# Builds the grid from the config, like panarctic.jl. Each region is labelled with its code; black lines are the borders
# between regions, dashed grey lines the 60°, 70° and 80°N parallels. Thin grey lines are the 200, 1000, 2000 and 3000 m
# isobaths; the thicker dark one is the 500 m isobath, deep_basin_depth, which bounds the deep basins EUBA and AMBA.
# Ocean outside every region of a set (TSAR) is light grey, land darker.

include("../setup.jl")
include("../regions.jl")

using CairoMakie
using Printf
using Oceananigans.Operators: Az

cpu_grid = on_architecture(CPU(), grid)
ccc  = (Center(), Center(), Center())
φ    = [φnode(i, j, 1, cpu_grid, ccc...) for i in 1:Nx, j in 1:Ny]
area = [Az(i, j, 1, cpu_grid, ccc...) for i in 1:Nx, j in 1:Ny]   # horizontal cell area
bottom_depth = -Array(interior(bottom_height, :, :, 1))
bottom_depth[bottom_depth .≤ 0] .= NaN

for set in sort(collect(keys(region_sets)))
    (; codes, names, masks) = region_masks(cpu_grid; set)
    R = length(codes)

    # Region number of each column: 0 in the ocean outside every region, NaN on land
    region_index = fill(NaN, Nx, Ny)
    region_index[isfinite.(bottom_depth)] .= 0
    for (n, code) in enumerate(codes)
        region_index[masks[code]] .= n
    end

    # Cell faces between columns of different regions (or between a region and land)
    borders = Point2f[]
    for j in 1:Ny, i in 1:Nx
        i < Nx && region_index[i, j] !== region_index[i+1, j] && append!(borders, [Point2f(i + 0.5, j - 0.5), Point2f(i + 0.5, j + 0.5)])
        j < Ny && region_index[i, j] !== region_index[i, j+1] && append!(borders, [Point2f(i - 0.5, j + 0.5), Point2f(i + 0.5, j + 0.5)])
    end

    colors = Makie.to_colormap(:tab20)[1:R]

    fig = Figure(size = (1000, 640))
    ax  = Axis(fig[1, 1], title = "Ocean regions, $set (regions.jl)", xlabel = "grid i", ylabel = "grid j",
               aspect = DataAspect())
    heatmap!(ax, region_index, colormap = cgrad(vcat(Makie.to_color(:gray92), colors), categorical = true),
             colorrange = (-0.5, R + 0.5), nan_color = :gray75)
    contour!(ax, φ, levels = [60, 70, 80], color = :gray40, linewidth = 0.8, linestyle = :dash)
    contour!(ax, bottom_depth, levels = [200, 1000, 2000, 3000], color = (:black, 0.35), linewidth = 0.6)
    contour!(ax, bottom_depth, levels = [deep_basin_depth], color = :gray15, linewidth = 1.4)
    linesegments!(ax, borders, color = :black, linewidth = 0.8)

    # Each label goes on the region's column nearest to its area-weighted centre, so it lies inside the region.
    for code in codes
        I = findall(masks[code])
        w = area[I]
        cᵢ = sum(w .* first.(Tuple.(I))) / sum(w)
        cⱼ = sum(w .* last.(Tuple.(I))) / sum(w)
        i, j = Tuple(I[argmin([(i - cᵢ)^2 + (j - cⱼ)^2 for (i, j) in Tuple.(I)])])
        text!(ax, i, j, text = code, align = (:center, :center), fontsize = 14, font = :bold, color = :black)
    end

    Legend(fig[1, 2], [PolyElement(color = c) for c in colors],
           [@sprintf("%-5s %s (%.2f M km²)", code, name, sum(area[masks[code]]) / 1e12) for (code, name) in zip(codes, names)],
           "Regions", framevisible = false, labelsize = 12)

    file = joinpath(figure_dir, "regions_$(lowercase(set)).png")
    save(file, fig)
    @info "Wrote $file"
end
