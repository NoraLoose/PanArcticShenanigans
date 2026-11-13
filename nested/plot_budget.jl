# Volume budget of one nested run, from its hourly monitor: the change in domain-mean η, the parts that came through
# the open edges and across the surface (freshwater), the residual, and GLORYS's change in mean η.
#
#   julia --project=.. plot_budget.jl output-coupled-control-30d      # writes <folder>/volume_budget.png

using CairoMakie

folder = isempty(ARGS) ? "output-nested" : ARGS[1]
text   = readlines(joinpath(folder, "nested_monitor.csv"))
header = split(text[1], ',')
rows   = [parse.(Float64, split(line, ',')) for line in text[2:end] if !isempty(line)]
column(name) = [row[findfirst(==(name), header)] for row in rows]

days   = column("hour") ./ 24
η      = column("eta_mean_cm")
parent = column("parent_eta_mean_cm")

fig = Figure(size = (800, 450))
ax  = Axis(fig[1, 1]; xlabel = "days", ylabel = "change in domain-mean η (cm)", title = basename(rstrip(folder, '/')))
hlines!(ax, 0; color = :black, linewidth = 0.5)
lines!(ax, days, parent .- parent[1]; color = :black, linewidth = 3, label = "GLORYS")
lines!(ax, days, η .- η[1]; color = :black, label = "model")
lines!(ax, days, column("edge_volume_cm"); label = "through the edges")
lines!(ax, days, column("surface_volume_cm"); label = "across the surface (freshwater)")
lines!(ax, days, column("budget_residual_cm"); linestyle = :dash, label = "residual")
axislegend(ax; position = :lb, framevisible = false)

file = joinpath(folder, "volume_budget.png")
save(file, fig)
@info "Wrote $file"
