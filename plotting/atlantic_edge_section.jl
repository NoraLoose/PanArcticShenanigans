# Temperature across the Atlantic edge: daily means from atlantic_edge_daily.jld2 against GLORYS on the same rows
# (atlantic_edge_glorys.jld2), for one run. Both are written when [diagnostics] sections = true.
#
#   julia +1.13 --project=. plotting/atlantic_edge_section.jl <output directory> [row] [day]
#
# row: which saved row (1 = the first row of cells inside the edge; default 4), day: which daily mean (default the
# last). Writes atlantic_edge_section.png in the output directory: model T, GLORYS T and model − GLORYS against i and
# depth, and the rms difference on every saved row for every day.

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.OutputWriters: jldopen

dir = ARGS[1]
row = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 4

T = FieldTimeSeries(joinpath(dir, "atlantic_edge_daily.jld2"), "T")
g = jldopen(joinpath(dir, "atlantic_edge_glorys.jld2"))
Tᵍ, tᵍ, t₀ = g["T"], g["times"], g["start_time"]
close(g)

# Daily means end at t₀ + 1, 2, ... days (the record at t₀, if any, is not an average)
averages = [n for n in eachindex(T.times) if T.times[n] > t₀ + 12 * 3600]
day_of(n) = round(Int, (T.times[n] - t₀) / 86400)
day = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : day_of(averages[end])
nd  = averages[findfirst(n -> day_of(n) == day, averages)]
z = znodes(T.grid, Center())
wet = interior(T[nd])[:, row, :] .!= 0

# GLORYS at the centre of the model's averaging window (the window ends at T.times[n])
function glorys_at(t)
    s = clamp((t - tᵍ[1]) / (tᵍ[2] - tᵍ[1]), 0, length(tᵍ) - 1)
    n = min(floor(Int, s) + 1, length(tᵍ) - 1)
    α = s - (n - 1)
    return (1 - α) .* Tᵍ[:, :, :, n] .+ α .* Tᵍ[:, :, :, n + 1]
end
window_centre(n) = T.times[n] - 12 * 3600
masked(a) = ifelse.(wet, a, NaN)

Tm = masked(interior(T[nd])[:, row, :])
Tg = masked(glorys_at(window_centre(nd))[:, row, :])
i = 1:size(Tm, 1)

fig = Figure(size = (1500, 1000), fontsize = 14)
kw = (colormap = :thermal, colorrange = (-2, 12))
ax1 = Axis(fig[1, 1], title = @sprintf("(a) model, day %d mean, row %d", day, row), ylabel = "z (m)")
hm = heatmap!(ax1, i, z, Tm; kw...)
ax2 = Axis(fig[1, 2], title = "(b) GLORYS, same window")
heatmap!(ax2, i, z, Tg; kw...)
Colorbar(fig[1, 3], hm, label = "T (°C)")
ax3 = Axis(fig[2, 1], title = "(c) model − GLORYS", xlabel = "i", ylabel = "z (m)")
hd = heatmap!(ax3, i, z, Tm .- Tg; colormap = :balance, colorrange = (-3, 3))
Colorbar(fig[2, 3], hd, label = "ΔT (°C)")

# rms difference on each saved row, day by day
ax4 = Axis(fig[2, 2], title = "(d) rms model − GLORYS on each saved row", xlabel = "day", ylabel = "°C")
for r in axes(Tᵍ, 2)
    rms_day = map(averages) do n
        m = interior(T[n])[:, r, :]
        d = (m .- glorys_at(window_centre(n))[:, r, :])[m .!= 0]
        sqrt(mean(abs2, d))
    end
    lines!(ax4, day_of.(averages), rms_day, label = "row $r", color = r, colormap = :viridis,
           colorrange = (1, size(Tᵍ, 2)))
end
axislegend(ax4, position = :lt, nbanks = 2)

save(joinpath(dir, "atlantic_edge_section.png"), fig)
println("saved ", joinpath(dir, "atlantic_edge_section.png"))
