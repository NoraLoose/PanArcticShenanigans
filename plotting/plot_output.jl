# Movie of SST, SSS, surface speed and ice volume from the output of panarctic.jl.
#
#   julia +1.13 --project=. plotting/plot_output.jl [config.toml]
#
# Reads arctic_ocean_surface.jld2 and arctic_sea_ice_surface.jld2 from `output_dir` and writes
# arctic_sea.mp4 there. Runs on the CPU and does not need the grid or input data.

using Pkg
Pkg.activate(dirname(@__DIR__))

using Oceananigans, Oceananigans.Units
using CairoMakie
using TOML

config_file = isempty(ARGS) ? joinpath(dirname(@__DIR__), "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
output_dir = config["paths"]["output_dir"]

To = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "T")
So = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "S")
Uo = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "𝒱")
hi = FieldTimeSeries(joinpath(output_dir, "arctic_sea_ice_surface.jld2"), "he")

times = To.times
n = Observable(length(times))

title = @lift "Arctic — hour " * string(round(Int, times[$n] / hours))

Tₙ = @lift(To[$n])
Sₙ = @lift(So[$n])
Uₙ = @lift(Uo[$n])
hₙ = @lift(hi[$n])

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

CairoMakie.record(fig, joinpath(output_dir, "arctic_sea.mp4"), 1:length(times), framerate = 4) do i
    n[] = i
end
