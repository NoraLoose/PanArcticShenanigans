# Bathymetry of the model grid. Writes bathymetry.png to `output_dir`.
#
#   julia +1.13 --project=. plotting/bathymetry.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../setup.jl")

using CairoMakie

h_bottom = Array(interior(bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN

fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
          title = "Arctic bathymetry", aspect = DataAspect())
hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-depth, 0), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "bottom height [m]")

save(joinpath(output_dir, "bathymetry.png"), fig)
