# Compare nested runs from their hourly monitors: the change in domain-mean η against GLORYS's, the net inflow
# through the open edges, and the change in η over the shelf and the deep ocean.
#
# Uses CairoMakie from the main environment (one folder up); reads only the CSV files, so no GPU is needed:
#
#   julia --project=.. plot_monitor.jl output-nested output-coupled output-coupled-zero
#
# Each argument is an output folder holding nested_monitor.csv; the folder name labels the run. Writes
# monitor_comparison.png in the current folder.

using CairoMakie

runs = isempty(ARGS) ? ["output-nested"] : ARGS

function read_monitor(folder)
    text   = readlines(joinpath(folder, "nested_monitor.csv"))
    header = split(text[1], ',')
    rows   = [parse.(Float64, split(line, ',')) for line in text[2:end] if !isempty(line)]
    column(name) = [row[findfirst(==(name), header)] for row in rows]
    return (; days    = column("hour") ./ 24,
              η       = column("eta_mean_cm"),
              parent  = column("parent_eta_mean_cm"),
              inflow  = column("net_in_Sv"),
              shelf   = column("shelf_change_cm"),
              deep    = column("deep_change_cm"))
end

monitors = [read_monitor(run) for run in runs]
labels   = basename.(rstrip.(runs, '/'))

fig = Figure(size = (800, 900))
ax_η      = Axis(fig[1, 1]; ylabel = "change in domain-mean η (cm)")
ax_inflow = Axis(fig[2, 1]; ylabel = "net inflow through the edges (Sv)")
ax_shelf  = Axis(fig[3, 1]; ylabel = "change in η (cm)", xlabel = "days",
                 title = "shelf (H < 200 m, solid) and deep ocean (H > 1000 m, dashed)")
linkxaxes!(ax_η, ax_inflow, ax_shelf)

# GLORYS from the longest run (the same in every run over their common days)
longest = monitors[argmax(length.(getproperty.(monitors, :days)))]
lines!(ax_η, longest.days, longest.parent .- longest.parent[1]; color = :black, linewidth = 3, label = "GLORYS")
hlines!(ax_inflow, 0; color = :black, linewidth = 0.5)

for (n, (m, label)) in enumerate(zip(monitors, labels))
    color = Makie.wong_colors()[n]
    lines!(ax_η, m.days, m.η .- m.η[1]; color, label)
    lines!(ax_inflow, m.days, m.inflow; color)
    lines!(ax_shelf, m.days, m.shelf; color)
    lines!(ax_shelf, m.days, m.deep; color, linestyle = :dash)
end

ylims!(ax_inflow, -20, 20)
axislegend(ax_η; position = :lb, framevisible = false)
save("monitor_comparison.png", fig)
@info "Wrote monitor_comparison.png" runs
