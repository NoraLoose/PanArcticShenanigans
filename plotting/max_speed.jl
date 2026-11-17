# Where the fastest velocity of a run is, from its monitor (monitor.jl, which records the grid point of the larger of
# max |u| and max |v| every hour): the largest speed and its model level against time (top), and the hourly locations
# of the fastest point on the seafloor depth, in grid indices, with the overall maximum marked (bottom).
# Writes max_speed.png to `figure_dir`.
#
#   julia +1.13 --project=. plotting/max_speed.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../config.jl")
include("../monitor_file.jl")

using CairoMakie, Printf

header, column = read_monitor(output_dir)
hours = column("hour")
speed = max.(column("max_u"), column("max_v"))
level = column("fastest_k")
i, j  = column("fastest_i"), column("fastest_j")
n     = argmax(speed)
@info @sprintf("Fastest: %s = %.2f m/s at hour %.0f, (i, j, k) = (%d, %d, %d), %.1f°E %.1f°N, z = %.0f m",
               column("fastest_component")[n] == 1 ? "u" : "v", speed[n], hours[n], i[n], j[n], level[n],
               column("fastest_longitude")[n], column("fastest_latitude")[n], column("fastest_z")[n])

depth = reduce(hcat, [parse.(Float64, split(line, ',')) for line in readlines(depth_file(output_dir))])   # (i, j)
depth[depth .== 0] .= NaN

fig = Figure(size = (800, 1000))
ax₁ = Axis(fig[1, 1]; xlabel = "days", ylabel = "largest |u|, |v| (m/s)", title = experiment)
lines!(ax₁, hours ./ 24, speed; color = :black)
ax₂ = Axis(fig[1, 1]; ylabel = "model level k (red)", yaxisposition = :right, ygridvisible = false)
scatter!(ax₂, hours ./ 24, level; color = :red, markersize = 3)
hidexdecorations!(ax₂)

ax₃ = Axis(fig[2, 1]; xlabel = "i", ylabel = "j", aspect = DataAspect(),
           title = "seafloor depth (m); dots: hourly location of the fastest point, colored by speed; ×: the maximum")
hm = heatmap!(ax₃, depth; colormap = :deep)
sc = scatter!(ax₃, i, j; color = speed, colormap = :hot, markersize = 6)
scatter!(ax₃, [i[n]], [j[n]]; color = :cyan, marker = :xcross, markersize = 18)
Colorbar(fig[2, 2], hm; label = "depth (m)")
Colorbar(fig[2, 3], sc; label = "speed (m/s)")

file = joinpath(figure_dir, "max_speed.png")
save(file, fig)
@info "Wrote $file"
