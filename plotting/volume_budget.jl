# Volume budget of a run, from its monitor (monitor.jl): the change in domain-mean η, the parts that came through the
# open sides and across the surface (freshwater), the residual (should be zero), and GLORYS's change in mean η (top);
# the net inflow through the open sides carried by the velocities at the boundary faces (the ones that carry T and S)
# against the inflow that actually changed the volume (bottom). The bottom panel starts after the first sample, which
# holds the interpolated initial velocities. Writes volume_budget.png to `figure_dir`.
#
#   julia +1.13 --project=. plotting/volume_budget.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../config.jl")
include("../monitor_file.jl")

using CairoMakie

header, column = read_monitor(output_dir)
days = column("hour") ./ 24
η, parent = column("eta_mean_cm"), column("parent_eta_mean_cm")

fig = Figure(size = (800, 800))
ax  = Axis(fig[1, 1]; xlabel = "days", ylabel = "change in domain-mean η (cm)", title = experiment)
hlines!(ax, 0; color = :black, linewidth = 0.5)
lines!(ax, days, parent .- parent[1]; color = :black, linewidth = 3, label = "GLORYS")
lines!(ax, days, η .- η[1]; color = :red, linestyle = :dash, label = "model")
lines!(ax, days, column("edge_volume_cm"); label = "through the open sides")
lines!(ax, days, column("surface_volume_cm"); label = "across the surface (freshwater)")
lines!(ax, days, column("budget_residual_cm"); linestyle = :dash, label = "residual")
axislegend(ax; position = :lb, framevisible = false)

ax₂ = Axis(fig[2, 1]; xlabel = "days", ylabel = "net inflow through the open sides (Sv)")
hlines!(ax₂, 0; color = :black, linewidth = 0.5)
stepped = 2:length(days)
lines!(ax₂, days[stepped], column("net_in_Sv")[stepped]; label = "velocities at the boundary faces (carry T, S)")
lines!(ax₂, days[stepped], column("edge_inflow_from_volume_Sv")[stepped]; label = "from the volume (moves η)")
axislegend(ax₂; position = :lb, framevisible = false)

file = joinpath(figure_dir, "volume_budget.png")
save(file, fig)
@info "Wrote $file"
