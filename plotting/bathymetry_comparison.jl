# GEBCO 2026 against GLORYS's own sea floor, both regridded onto the model grid with the same steps (2 interpolation
# passes, 15 m minimum depth, one basin). Whole domain and the Nordic Seas window of barotropic_map.jl, with f/H
# contours (log10 spacing): the barotropic flow follows f/H, so where these differ GLORYS's barotropic flow is not
# steady on our grid.
#
#   julia --project=. plotting/bathymetry_comparison.jl [config.toml]
#
# Writes bathymetry_comparison.png in the output directory and prints depth and wet-mask statistics.

include("../setup.jl")
@isdefined(glorys_bottom_height) || include("../glorys_bathymetry.jl")

using CairoMakie, Statistics
using Oceananigans.Architectures: on_architecture
using Oceananigans.Grids: φnode

gebco  = Array(interior(regrid_bathymetry(underlying_grid, Metadatum(:bottom_height; dataset = GEBCO2026(), dir_kw...);
                                          minimum_depth = 15, interpolation_passes = 2, major_basins = 1)))[:, :, 1]
glorys = Array(interior(glorys_bottom_height(underlying_grid; minimum_depth = 15, interpolation_passes = 2,
                                             major_basins = 1)))[:, :, 1]

cpu_grid_u = on_architecture(CPU(), underlying_grid)
φ = [φnode(i, j, 1, cpu_grid_u, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
f = 2 * 7.292e-5 .* sind.(φ)

depth_of(z) = ifelse.(z .< 0, -z, NaN)
Hb, Hg = depth_of(gebco), depth_of(glorys)
both = isfinite.(Hb) .& isfinite.(Hg)
window = (200:380, 60:240)

function report(name, mask)
    d = (Hg .- Hb)[mask .& both]
    r = (Hg ./ Hb)[mask .& both]
    @printf("%-14s wet in both %6d, GEBCO only %5d, GLORYS only %5d; GLORYS − GEBCO: mean %+6.1f m, rms %6.1f m; |relative difference| > 10%%: %4.1f%% of columns\n", name, count(mask .& both),
            count(mask .& isfinite.(Hb) .& .!isfinite.(Hg)), count(mask .& .!isfinite.(Hb) .& isfinite.(Hg)),
            mean(d), sqrt(mean(abs2, d)), 100 * count(abs.(r .- 1) .> 0.1) / length(r))
end
everywhere = trues(Nx, Ny)
in_window = falses(Nx, Ny); in_window[window...] .= true
println()
report("whole domain", everywhere)
report("Nordic Seas", in_window)

logfH(H) = log10.(f ./ H)
levels = range(-7.7, -5.0, length = 28)   # f/H from ~2e-8 (5 km deep) to ~9e-6 s⁻¹ m⁻¹ (15 m)

fig = Figure(size = (1800, 1150), fontsize = 13)
depth_kw = (colormap = :deep, colorrange = (0, 4000), nan_color = :gray80)
for (r, (is, js)) in enumerate(((1:Nx, 1:Ny), window))
    sub(a) = a[is, js]
    area_name = r == 1 ? "whole domain" : "Nordic Seas"
    ax1 = Axis(fig[r, 1], title = "GEBCO depth (m), $area_name", aspect = DataAspect())
    hm = heatmap!(ax1, is, js, sub(Hb); depth_kw...)
    ax2 = Axis(fig[r, 2], title = "GLORYS depth (m), $area_name", aspect = DataAspect())
    heatmap!(ax2, is, js, sub(Hg); depth_kw...)
    Colorbar(fig[r, 3], hm)
    if r == 2
        contour!(ax1, is, js, sub(logfH(Hb)); levels, color = :white, linewidth = 0.5)
        contour!(ax2, is, js, sub(logfH(Hg)); levels, color = :white, linewidth = 0.5)
        ax1.title = "GEBCO depth, f/H contours"; ax2.title = "GLORYS depth, f/H contours"
    end
    ax3 = Axis(fig[r, 4], title = "GLORYS / GEBCO depth − 1, $area_name", aspect = DataAspect())
    hd = heatmap!(ax3, is, js, sub(ifelse.(both, Hg ./ Hb .- 1, NaN)); colormap = :balance, colorrange = (-0.5, 0.5),
                  nan_color = :gray80)
    Colorbar(fig[r, 5], hd)
end

figure_file = joinpath(output_dir, "bathymetry_comparison.png")
save(figure_file, fig)
@info "Wrote $figure_file"
