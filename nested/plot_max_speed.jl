# Where the fastest velocity of one nested run is, from its hourly monitor (which records the grid point of the larger
# of max |u| and max |v|): the largest speed and its model level against time (top), and the hourly locations of the
# fastest point on the seafloor depth, in grid indices, with the overall maximum marked (bottom).
#
#   julia --project=.. plot_max_speed.jl output-coupled-glorys-relax-30d      # writes <folder>/max_speed.png

using CairoMakie, Printf

folder = isempty(ARGS) ? "output-nested" : ARGS[1]
text   = readlines(joinpath(folder, "nested_monitor.csv"))
header = split(text[1], ',')
rows   = [parse.(Float64, split(line, ',')) for line in text[2:end] if !isempty(line)]
column(name) = [row[findfirst(==(name), header)] for row in rows]

hours = column("hour")
speed = max.(column("max_u"), column("max_v"))
level = column("fastest_k")
i, j  = column("fastest_i"), column("fastest_j")
n     = argmax(speed)
@info @sprintf("Fastest: %s = %.2f m/s at hour %.0f, (i, j, k) = (%d, %d, %d), %.1f°E %.1f°N, z = %.0f m",
               column("fastest_component")[n] == 1 ? "u" : "v", speed[n], hours[n], i[n], j[n], level[n],
               column("fastest_longitude")[n], column("fastest_latitude")[n], column("fastest_z")[n])

depth = reduce(hcat, [parse.(Float64, split(line, ',')) for line in readlines(joinpath(folder, "depth.csv"))])   # (i, j)
depth[depth .== 0] .= NaN

fig = Figure(size = (800, 1000))
ax₁ = Axis(fig[1, 1]; xlabel = "days", ylabel = "largest |u|, |v| (m/s)", title = basename(rstrip(folder, '/')))
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

file = joinpath(folder, "max_speed.png")
save(file, fig)
@info "Wrote $file"
