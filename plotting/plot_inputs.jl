# Figures of the model inputs: bathymetry, grid spacing, and a movie of the GLORYS free surface
# (the boundary data) on the model grid. Written to `output_dir`.
#
#   julia +1.13 --project=. plotting/plot_inputs.jl [config.toml]
#
# Uses the same config as panarctic.jl. The GLORYS data are read from the download/inpainting cache,
# so this is quick after a model run (or downloads and inpaints the data otherwise).

include("../setup.jl")

using CairoMakie

h_bottom = Array(interior(bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN
dry = isnan.(h_bottom)

# ## Bathymetry

fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
          title = "Arctic bathymetry", aspect = DataAspect())
hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-depth, 0), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "bottom height [m]")

save(joinpath(output_dir, "bathymetry.png"), fig)

# ## Grid spacing

Δx = Array(interior(compute!(Field(xspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
Δy = Array(interior(compute!(Field(yspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
Δx[dry] .= NaN
Δy[dry] .= NaN

fig = Figure(size = (1000, 450))
ax1 = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j", title = "Δx [km]", aspect = DataAspect())
hm1 = heatmap!(ax1, Δx, colormap = :viridis)
Colorbar(fig[1, 2], hm1)
ax2 = Axis(fig[1, 3], xlabel = "grid i", title = "Δy [km]", aspect = DataAspect())
hm2 = heatmap!(ax2, Δy, colormap = :viridis)
Colorbar(fig[1, 4], hm2)

save(joinpath(output_dir, "grid_spacing.png"), fig)

# ## Movie of the GLORYS free surface
#
# Same settings as in panarctic.jl, so the inpainted GLORYS data are reused from the cache.

ηᵉˣᵗ = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))

iter = Observable(1)
ηv = @lift begin
    η = Array(interior(ηᵉˣᵗ[$iter], :, :, 1))
    η[dry] .= NaN
    η
end
fig = Figure(size = (700, 600))
ax  = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
           title = "GLORYS free surface", aspect = DataAspect())
hm  = heatmap!(ax, ηv, colormap = :balance, colorrange = (-1.5, 0.2), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "Free surface [m]")
CairoMakie.record(fig, joinpath(output_dir, "arctic_free_surface.mp4"), 1:length(dates), framerate = 4) do i
    iter[] = i
end
