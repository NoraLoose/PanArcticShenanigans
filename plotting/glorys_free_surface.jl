# Movie of the GLORYS free surface (the boundary data) on the model grid.
# Writes arctic_free_surface.mp4 to `figure_dir`.
#
#   julia +1.13 --project=. plotting/glorys_free_surface.jl [config.toml]
#
# Uses the same config as panarctic.jl. The GLORYS data are read the same way as in panarctic.jl (from
# the download/inpainting cache), so this is quick after a model run (otherwise downloads and inpaints them).

include("../setup.jl")

using CairoMakie

land = Array(interior(bottom_height, :, :, 1)) .≥ 0

ηᵉˣᵗ = glorys_free_surface()

iter = Observable(1)
ηv = @lift begin
    η = Array(interior(ηᵉˣᵗ[$iter], :, :, 1))
    η[land] .= NaN
    η
end
fig = Figure(size = (700, 600))
ax  = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
           title = "GLORYS free surface", aspect = DataAspect())
hm  = heatmap!(ax, ηv, colormap = :balance, colorrange = (-1.5, 0.2), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "Free surface [m]")
CairoMakie.record(fig, joinpath(figure_dir, "arctic_free_surface.mp4"), 1:length(dates), framerate = 4) do i
    iter[] = i
end
