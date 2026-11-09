# How the model's free surface changes hour by hour: η(t) − η(0), from free_surface_hourly.jld2 (written when
# [diagnostics] sections = true). The domain-mean change (the volume change through the edges and the surface) is
# removed, so what is left is the redistribution of water. Prints, for each hour, the rms change in the window and
# where the largest change is; maps a selection of hours; and the hourly means over three regions (subpolar North
# Atlantic, Nordic Seas, Arctic) with the period of the see-saw between the first two (eta_hourly_regions.png); the
# change by depth class (shelf, slope, deep) without removing the domain mean; and the same for GLORYS, day by day.
#
#   julia --project=. plotting/eta_hourly.jl <config.toml> [hour ...] [window=south|nordic|all|i₁:i₂,j₁:j₂]
#
# Default hours: 1 2 3 6 9 12 18 24 30 36 42 48; default window: south (Atlantic edge to Fram Strait). Writes
# eta_hourly_<window>.png in the output directory.

include("../setup.jl")

using Printf, Statistics, CairoMakie
using Oceananigans.Grids: λnode, φnode, static_column_depthᶜᶜᵃ
using NumericalEarth.DataWrangling: NCDatasets, metadata_path
using .NCDatasets: Dataset

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
    δ = change(n)[irange, jrange]
    finite = filter(isfinite, δ)
    isempty(finite) && continue
    k = argmax(ifelse.(isfinite.(δ), abs.(δ), -Inf))
    i, j = irange[k[1]], jrange[k[2]]
    @printf("%6.1f %12.2f %12.2f   (%d, %d, %.1f°E, %.1f°N)\n", (η.times[n] - start_time) / hours,
            100 * sqrt(mean(abs2, finite)), 100 * abs(δ[k]), i, j,
            λnode(i, j, 1, mgrid, Center(), Center(), Center()), φnode(i, j, 1, mgrid, Center(), Center(), Center()))
end

# Region means of η(t) − η(0) (domain mean removed), every hour: a basin-scale oscillation shows as a see-saw between
# the subpolar North Atlantic and the Nordic Seas. The period and decay of their difference are estimated from its
# zero crossings and successive extrema.
regions = (subpolar = (150:330, 15:90),     # Irminger Sea, Iceland Basin, Rockall (south of Greenland–Scotland)
           nordic   = (250:360, 120:215),   # Norwegian, Lofoten and Greenland basins, Iceland Sea
           arctic   = (200:330, 250:330))   # Eurasian Basin
function region_mean(a, (is, js))
    w = a[is, js]
    finite = filter(isfinite, w)
    return isempty(finite) ? NaN : mean(finite)
end
hours_all = [(η.times[n] - start_time) / hours for n in eachindex(η.times) if η.times[n] ≥ start_time - 60]
records = [n for n in eachindex(η.times) if η.times[n] ≥ start_time - 60]
series = Dict(name => [100 * region_mean(change(n), r) for n in records] for (name, r) in pairs(regions))
seesaw = series[:nordic] .- series[:subpolar]
println("\nRegion means of η(t) − η(0), domain mean removed (cm)")
@printf("%6s %10s %10s %10s %16s\n", "hour", "subpolar", "nordic", "arctic", "nordic − subpolar")
for (m, h) in enumerate(hours_all)
    @printf("%6.1f %10.2f %10.2f %10.2f %16.2f\n", h, series[:subpolar][m], series[:nordic][m], series[:arctic][m], seesaw[m])
end
crossings = [hours_all[m] + (hours_all[m+1] - hours_all[m]) * seesaw[m] / (seesaw[m] - seesaw[m+1])
             for m in 2:length(seesaw)-1 if sign(seesaw[m]) != sign(seesaw[m+1]) && seesaw[m] != 0]
extrema_at = [m for m in 2:length(seesaw)-1 if (seesaw[m] - seesaw[m-1]) * (seesaw[m+1] - seesaw[m]) < 0]
println("Zero crossings of nordic − subpolar (h): ", join((@sprintf("%.1f", c) for c in crossings), ", "))
length(crossings) ≥ 2 && @printf("Period from the crossings: %.1f h (twice the mean spacing)\n", 2 * mean(diff(crossings)))
println("Extrema of nordic − subpolar: ", join((@sprintf("%.1f h: %+.1f cm", hours_all[m], seesaw[m]) for m in extrema_at), "; "))

# The same change WITHOUT removing the domain mean, split by depth over the whole domain: the domain mean itself (the
# net volume change through the edges and the surface), and the means over the shelves (H < 200 m), the slopes
# (200–1000 m) and the deep ocean (H > 1000 m). The region means above are relative to the domain mean, so a fall on
# the large shelves shows there as a rise everywhere deep; this table separates the two.
column_depth = try
    [static_column_depthᶜᶜᵃ(i, j, mgrid) for i in 1:Nx, j in 1:Ny]
catch
    nothing
end
if isnothing(column_depth)
    @warn "the saved grid has no bathymetry: skipping the means by depth"
else
    classes = (shelf = wet .& (column_depth .< 200), slope = wet .& (200 .≤ column_depth .< 1000), deep = wet .& (column_depth .≥ 1000))
    full_change(n) = ifelse.(wet, Array(interior(η[n]))[:, :, 1] .- η₀, NaN)
    class_mean(a, mask) = any(mask) ? mean(a[mask]) : NaN
    println("\nη(t) − η(0) WITHOUT removing the domain mean (cm); share of the wet cells: ",
            join((@sprintf("%s %.0f%%", name, 100 * count(m) / count(wet)) for (name, m) in pairs(classes)), ", "))
    @printf("%6s %12s %10s %10s %10s\n", "hour", "domain mean", "shelf", "slope", "deep")
    for (m, h) in enumerate(hours_all)
        a = full_change(records[m])
        @printf("%6.1f %12.2f %10.2f %10.2f %10.2f\n", h, 100 * mean(a[wet]),
                100 * class_mean(a, classes.shelf), 100 * class_mean(a, classes.slope), 100 * class_mean(a, classes.deep))
    end

    # The same means for GLORYS, the reference: its daily means are at noon and the run starts at noon on
    # dates[start_index], so GLORYS day k is at hour 24k. GLORYS zos is sampled bilinearly at the model's cell centres;
    # both are averaged over the cells where GLORYS is defined, so the two columns cover the same water.
    λc = [λnode(i, j, 1, mgrid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    φc = [φnode(i, j, 1, mgrid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    function bilinear(F, lon, lat, x, y)    # regular (lon, lat) grid, periodic in longitude; NaN if any corner is NaN
        fi = mod(x - lon[1], 360) / (lon[2] - lon[1])
        fj = (y - lat[1]) / (lat[2] - lat[1])
        i₀ = floor(Int, fi); j₀ = clamp(floor(Int, fj), 0, length(lat) - 2)
        a, b = fi - i₀, clamp(fj - j₀, 0, 1)
        i₁, i₂ = mod(i₀, length(lon)) + 1, mod(i₀ + 1, length(lon)) + 1
        return (1 - a) * (1 - b) * F[i₁, j₀+1] + a * (1 - b) * F[i₂, j₀+1] + (1 - a) * b * F[i₁, j₀+2] + a * b * F[i₂, j₀+2]
    end
    function glorys_zos(date)
        ds  = Dataset(metadata_path(Metadatum(:free_surface; date, dataset, region, dir_kw...)))
        lon = Float64.(ds["longitude"][:])
        lat = Float64.(ds["latitude"][:])
        zos = Float64.(coalesce.(ds["zos"][:, :, 1], NaN32))
        close(ds)
        return [w ? bilinear(zos, lon, lat, x, y) : NaN for (w, x, y) in zip(wet, λc, φc)]
    end
    glorys_days = [k for k in 0:(length(dates) - start_index) if 24k ≤ hours_all[end] + 0.5]
    zos₀ = glorys_zos(dates[start_index])
    common = wet .& isfinite.(zos₀)
    common_mean(a, mask) = (sel = mask .& common .& isfinite.(a); any(sel) ? 100 * mean(a[sel]) : NaN)
    println("\nGLORYS against the model at the same hours: change since the start (cm), domain mean NOT removed")
    @printf("%6s | %8s %8s %8s %8s | %8s %8s %8s %8s\n", "hour", "GLORYS", "shelf", "slope", "deep", "model", "shelf", "slope", "deep")
    for k in glorys_days
        glorys_δ = glorys_zos(dates[start_index + k]) .- zos₀
        model_δ  = full_change(record_at(24k))
        @printf("%6d | %8.2f %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f\n", 24k,
                common_mean(glorys_δ, wet), common_mean(glorys_δ, classes.shelf),
                common_mean(glorys_δ, classes.slope), common_mean(glorys_δ, classes.deep),
                common_mean(model_δ, wet), common_mean(model_δ, classes.shelf),
                common_mean(model_δ, classes.slope), common_mean(model_δ, classes.deep))
    end
end

fig_ts = Figure(size = (1000, 450), fontsize = 14)
ax_ts = Axis(fig_ts[1, 1], xlabel = "hour", ylabel = "cm", title = "region means of η(t) − η(0), domain mean removed")
for (name, color) in ((:subpolar, :steelblue), (:nordic, :firebrick), (:arctic, :gray50))
    lines!(ax_ts, hours_all, series[name], color = color, label = String(name))
end
lines!(ax_ts, hours_all, seesaw, color = :black, linestyle = :dash, label = "nordic − subpolar")
axislegend(ax_ts, position = :lb, framevisible = false)
save(joinpath(output_dir, "eta_hourly_regions.png"), fig_ts)

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
