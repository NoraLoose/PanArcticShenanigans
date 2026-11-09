# How the model's free surface changes hour by hour: η(t) − η(0), from free_surface_hourly.jld2 (written when
# [diagnostics] sections = true). The domain-mean change (the volume change through the edges and the surface) is
# removed, so what is left is the redistribution of water. Prints, for each hour, the rms change in the window and
# where the largest change is; maps a selection of hours.
#
#   julia --project=. plotting/eta_hourly.jl <config.toml> [hour ...] [window=south|nordic|all|i₁:i₂,j₁:j₂]
#
# Default hours: 1 2 3 6 9 12 18 24 30 36 42 48; default window: south (Atlantic edge to Fram Strait). Writes
# eta_hourly_<window>.png in the output directory.

include("../setup.jl")

using Printf, Statistics, CairoMakie
using Oceananigans.Grids: λnode, φnode, static_column_depthᶜᶜᵃ

options = ARGS[2:end]
window_option = something(findfirst(a -> startswith(a, "window="), options), 0)
window_name = window_option > 0 ? options[window_option][8:end] : "south"
hour_arguments = [a for a in options if !startswith(a, "window=")]
hours_to_plot = isempty(hour_arguments) ? [1, 2, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48] : parse.(Int, hour_arguments)
irange, jrange = if window_name == "nordic"
    200:380, 60:240
elseif window_name == "south"
    120:420, 1:240
elseif window_name == "all"
    1:Nx, 1:Ny
else
    Tuple(Meta.eval(Meta.parse(x)) for x in split(window_name, ','))
end

file = joinpath(output_dir, "free_surface_hourly.jld2")
η = FieldTimeSeries(file, "η")
mgrid = η.grid
wet = try
    [static_column_depthᶜᶜᵃ(i, j, mgrid) > 0 for i in 1:Nx, j in 1:Ny]
catch
    interior(η[end])[:, :, 1] .!= 0
end
record_at(h) = argmin(abs.(η.times .- (start_time + h * hours)))
n₀ = record_at(0)
abs(η.times[n₀] - start_time) < 60 || @warn "no record at the start time; using the first one, at $(η.times[n₀] / hours) h"
η₀ = Array(interior(η[n₀]))[:, :, 1]

function change(n)
    d = ifelse.(wet, Array(interior(η[n]))[:, :, 1] .- η₀, NaN)
    return d .- mean(filter(isfinite, d))   # remove the domain-mean change
end

println("\nη(t) − η(0), domain-mean change removed; window i = $irange, j = $jrange")
@printf("%6s %12s %12s   %s\n", "hour", "rms (cm)", "max |.| (cm)", "where the largest is (i, j, °E, °N)")
for h in sort(unique(vcat(0:min(48, round(Int, (η.times[end] - start_time) / hours)), hours_to_plot)))
    n = record_at(h)
    w = change(n)[irange, jrange]
    finite = filter(isfinite, w)
    isempty(finite) && continue
    k = argmax(ifelse.(isfinite.(w), abs.(w), -Inf))
    i, j = irange[k[1]], jrange[k[2]]
    @printf("%6.1f %12.2f %12.2f   (%d, %d, %.1f°E, %.1f°N)\n", (η.times[n] - start_time) / hours,
            100 * sqrt(mean(abs2, finite)), 100 * abs(w[k]), i, j,
            λnode(i, j, 1, mgrid, Center(), Center(), Center()), φnode(i, j, 1, mgrid, Center(), Center(), Center()))
end

ncols = 4
nrows = cld(length(hours_to_plot), ncols)
fig = Figure(size = (450 * ncols + 100, 400 * nrows), fontsize = 13)
heatmaps = []
for (m, h) in enumerate(hours_to_plot)
    n = record_at(h)
    ax = Axis(fig[(m - 1) ÷ ncols + 1, (m - 1) % ncols + 1], aspect = DataAspect(),
              title = @sprintf("hour %.0f: η − η(0), mean removed (cm)", (η.times[n] - start_time) / hours))
    push!(heatmaps, heatmap!(ax, irange, jrange, 100 .* change(n)[irange, jrange]; colormap = :balance, colorrange = (-15, 15),
                              nan_color = :gray80))
end
Colorbar(fig[:, ncols + 1], heatmaps[end])
figure_file = joinpath(output_dir, "eta_hourly_$(replace(window_name, r"[^A-Za-z0-9]" => "_")).png")
save(figure_file, fig)
@info "Wrote $figure_file"
