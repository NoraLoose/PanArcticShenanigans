# Movie of SST, SSS, surface speed and ice volume from the time-averaged output of panarctic.jl,
# one frame per saved average (daily by default).
#
#   julia +1.13 --project=. plotting/surface_movie.jl [config.toml]
#
# Reads arctic_ocean.jld2 and arctic_sea_ice.jld2 from `output_dir` and writes arctic_sea.mp4
# to `figure_dir`. Runs on the CPU and does not need the grid or input data. The 3D fields are read one frame at a time.

include("../config.jl")

using Oceananigans, Oceananigans.Units
using Oceananigans.OutputReaders: OnDisk
using CairoMakie

T  = FieldTimeSeries(ocean_output_file,   "T";  backend = OnDisk())
S  = FieldTimeSeries(ocean_output_file,   "S";  backend = OnDisk())
u  = FieldTimeSeries(ocean_output_file,   "u";  backend = OnDisk())
v  = FieldTimeSeries(ocean_output_file,   "v";  backend = OnDisk())
h  = FieldTimeSeries(sea_ice_output_file, "h";  backend = OnDisk())
ℵ  = FieldTimeSeries(sea_ice_output_file, "ℵ";  backend = OnDisk())

times = T.times
Nz    = size(T.grid, 3)
n     = Observable(1)

surface(φ) = interior(φ, :, :, Nz)
speed(n)   = interior(compute!(Field(@at((Center, Center, Center), sqrt(u[n]^2 + v[n]^2)), indices = (:, :, Nz))), :, :, 1)

title = @lift "Arctic — " * prettytime(times[$n] - times[1]) * " after the start"

Tₙ = @lift surface(T[$n])
Sₙ = @lift surface(S[$n])
Uₙ = @lift speed($n)
hₙ = @lift interior(h[$n], :, :, 1) .* interior(ℵ[$n], :, :, 1)   # ice volume per area

fig = Figure(size = (1100, 1000))
fig[0, 1:4] = Label(fig, title, fontsize = 20, tellwidth = false)

ax = Axis(fig[1, 1], ylabel = "grid j", aspect = DataAspect())
hm_T = heatmap!(ax, Tₙ, colormap = :thermal, colorrange = (-2, 8), nan_color = :gray80)
ax = Axis(fig[1, 3], aspect = DataAspect())
hm_S = heatmap!(ax, Sₙ, colormap = :haline, colorrange = (32.5, 35.5), nan_color = :gray80)
ax = Axis(fig[2, 1], xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
hm_U = heatmap!(ax, Uₙ, colormap = Reverse(:solar), colorrange = (0, 0.5), nan_color = :gray80)
ax = Axis(fig[2, 3], xlabel = "grid i", aspect = DataAspect())
hm_h = heatmap!(ax, hₙ, colormap = Reverse(:blues), colorrange = (0.01, 1.5), lowclip = :gray80)
Colorbar(fig[1, 2], hm_T, label = "SST [°C]")
Colorbar(fig[1, 4], hm_S, label = "SSS [psu]")
Colorbar(fig[2, 2], hm_U, label = "Surface speed [ms⁻¹]")
Colorbar(fig[2, 4], hm_h, label = "ice volume per area [m]")

CairoMakie.record(fig, joinpath(figure_dir, "arctic_sea.mp4"), 1:length(times), framerate = 4) do i
    n[] = i
end
