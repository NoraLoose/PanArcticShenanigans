# Volume budget of one nested run, from its hourly monitor: the change in domain-mean η, the parts that came through
# the open edges and across the surface (freshwater), the residual, and GLORYS's change in mean η (top); the net inflow
# through the edges carried by the velocities at the edge faces (the ones that carry T and S) against the inflow that
# actually changed the volume (bottom). The bottom panel starts after the first step: the first monitor sample holds
# the interpolated initial velocities, whose net inflow (about 22 Sv) the boundary conditions remove in the first step.
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

fig = Figure(size = (800, 800))
ax  = Axis(fig[1, 1]; xlabel = "days", ylabel = "change in domain-mean η (cm)", title = basename(rstrip(folder, '/')))
hlines!(ax, 0; color = :black, linewidth = 0.5)
lines!(ax, days, parent .- parent[1]; color = :black, linewidth = 3, label = "GLORYS")
lines!(ax, days, η .- η[1]; color = :black, label = "model")
lines!(ax, days, column("edge_volume_cm"); label = "through the edges")
lines!(ax, days, column("surface_volume_cm"); label = "across the surface (freshwater)")
lines!(ax, days, column("budget_residual_cm"); linestyle = :dash, label = "residual")
axislegend(ax; position = :lb, framevisible = false)

if "edge_inflow_from_volume_Sv" in header
    ax₂ = Axis(fig[2, 1]; xlabel = "days", ylabel = "net inflow through the edges (Sv)")
    hlines!(ax₂, 0; color = :black, linewidth = 0.5)
    stepped = 2:length(days)
    lines!(ax₂, days[stepped], column("net_in_Sv")[stepped]; label = "velocities at the edge faces (carry T, S)")
    lines!(ax₂, days[stepped], column("edge_inflow_from_volume_Sv")[stepped]; label = "from the volume (moves η)")
    axislegend(ax₂; position = :lb, framevisible = false)
end

file = joinpath(folder, "volume_budget.png")
save(file, fig)
@info "Wrote $file"
