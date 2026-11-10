# Horizontal grid spacing (Δx, Δy) over the ocean. Writes grid_spacing.png to `figure_dir`.
#
#   julia +1.13 --project=. plotting/grid_spacing.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../setup.jl")

using CairoMakie

land = Array(interior(bottom_height, :, :, 1)) .≥ 0

Δx = Array(interior(compute!(Field(xspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
Δy = Array(interior(compute!(Field(yspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
Δx[land] .= NaN
Δy[land] .= NaN

fig = Figure(size = (1000, 450))
ax1 = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j", title = "Δx [km]", aspect = DataAspect())
hm1 = heatmap!(ax1, Δx, colormap = :viridis)
Colorbar(fig[1, 2], hm1)
ax2 = Axis(fig[1, 3], xlabel = "grid i", title = "Δy [km]", aspect = DataAspect())
hm2 = heatmap!(ax2, Δy, colormap = :viridis)
Colorbar(fig[1, 4], hm2)

save(joinpath(figure_dir, "grid_spacing.png"), fig)
