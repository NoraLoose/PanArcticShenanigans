# GLORYS's own volume budget over the model domain, for every date the config loads: the net transport through the two
# open edges on GLORYS's own grid, and the change of GLORYS's domain-mean free surface on the model grid (the same mean
# as the `glorys_eta_mean_m` column of the transport monitor). Their difference,
#
#     D = A d⟨ηᴳ⟩/dt − (Q_Atlantic − Q_Pacific),
#
# is the part of GLORYS's sea-level change that its edge transports do not explain: its surface freshwater flux, its
# assimilation increments, and the error of the native-section estimate. With a transport monitor from a model run, the
# model's surface freshwater flux is shown beside it.
#
#   julia --project=. plotting/glorys_volume_budget.jl [config.toml] [monitor csv]
#
# The config sets the dates (use one with the run length you want, e.g. stop_hours = 720 for 32 days of GLORYS).
# Missing daily u, v files are downloaded; η is read through the same inpainted cache as the model. Writes
# glorys_volume_budget.csv and glorys_volume_budget.png in output_dir. Transports in Sv; Atlantic positive into the
# domain, Pacific positive out of it (as in the monitor).

include("../setup.jl")

using Printf, Statistics, CairoMakie
using Downloads
using Oceananigans.Architectures: on_architecture
using Oceananigans.Grids: static_column_depthᶜᶜᵃ
using Oceananigans.Operators: Azᶜᶜᶜ

include("../glorys_sections.jl")

# ## GLORYS's transports through the edges, on its own grid

for date in dates, name in (:u_velocity, :v_velocity)
    Downloads.download(Metadatum(name; date, dataset, region, dir_kw...))   # skips files already there
end
# Net and gross (each direction separately) transport across a polyline of corners, as section_transport in
# glorys_sections.jl. A net that is a small difference of large gross transports is sensitive to any section error.
function section_net_and_gross(corners, lon, lat, Uᵉ, Uⁿ)
    R = section_grid.radius
    positive = negative = 0.0
    for k in 1:length(corners) - 1
        (λa, φa), (λb, φb) = corners[k], corners[k + 1]
        Δλ = wrap_longitude(λb - λa)
        λm, φm = λa + Δλ / 2, (φa + φb) / 2
        dx = R * cosd(φm) * deg2rad(Δλ)
        dy = R * deg2rad(φb - φa)
        T = sample_lonlat(Uᵉ, lon, lat, λm, φm) * (-dy) + sample_lonlat(Uⁿ, lon, lat, λm, φm) * dx
        T > 0 ? (positive += T) : (negative -= T)
    end
    return positive - negative, positive, negative
end

N = length(dates)
atlantic, pacific = zeros(N), zeros(N)                     # net, Sv
atlantic_gross, pacific_gross = zeros(N, 2), zeros(N, 2)   # (positive, negative) parts, Sv
for (n, date) in enumerate(dates)
    lon, lat, Uᵉ = native_depth_integral(:u_velocity, date)
    _,   _,   Uⁿ = native_depth_integral(:v_velocity, date)
    a = section_net_and_gross(atlantic_corners, lon, lat, Uᵉ, Uⁿ) ./ 1e6
    p = section_net_and_gross(pacific_corners,  lon, lat, Uᵉ, Uⁿ) ./ 1e6
    atlantic[n], atlantic_gross[n, 1], atlantic_gross[n, 2] = a
    pacific[n],  pacific_gross[n, 1],  pacific_gross[n, 2]  = p
end
net = atlantic .- pacific   # net inflow
@printf("\nGross transports (mean over the record): Atlantic %.1f Sv in, %.1f Sv out (net %.2f); Pacific %.1f Sv out, %.1f Sv in (net %.2f)\n",
        mean(atlantic_gross[:, 1]), mean(atlantic_gross[:, 2]), mean(atlantic),
        mean(pacific_gross[:, 1]), mean(pacific_gross[:, 2]), mean(pacific))

# ## GLORYS's domain-mean free surface, on the model grid (as in the transport monitor)

cpu_grid = on_architecture(CPU(), grid)
wet_area = [static_column_depthᶜᶜᵃ(i, j, cpu_grid) > 0 ? Azᶜᶜᶜ(i, j, 1, cpu_grid) : 0.0 for i in 1:Nx, j in 1:Ny]
domain_area = sum(wet_area)

ηᴳ = glorys_free_surface()
length(ηᴳ.times) == length(dates) || error("expected one GLORYS free surface per date")
mean_η = [sum(Array(interior(ηᴳ[n]))[:, :, 1] .* wet_area) / domain_area for n in eachindex(dates)]   # m

# ## The budget between consecutive daily means

Sv_per_cm_per_day = domain_area * 0.01 / 86400 / 1e6
volume_rate = [domain_area * (mean_η[n+1] - mean_η[n]) / 86400 / 1e6 for n in 1:N-1]   # Sv
mean_net    = [(net[n] + net[n+1]) / 2 for n in 1:N-1]                                  # Sv, trapezoid
D           = volume_rate .- mean_net

# The model's surface freshwater flux over the same intervals, from a transport monitor if there is one. Snapshot n
# (noon on dates[n]) is at monitor hour 24 (n - start_index).
monitor_file = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(output_dir, "transport_monitor.csv")
function interval_means(path, column)
    isfile(path) || return fill(NaN, N - 1)
    lines  = readlines(path)
    header = split(lines[1], ',')
    rows   = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(l)]
    hour   = [r[1] for r in rows]
    value  = [r[findfirst(==(column), header)] for r in rows]
    means  = fill(NaN, N - 1)
    for n in 1:N-1
        h₀, h₁ = 24 * (n - start_index), 24 * (n + 1 - start_index)
        inside = findall(h -> h₀ < h ≤ h₁, hour)
        length(inside) ≥ 20 && (means[n] = mean(value[inside]))   # only intervals the run covered
    end
    return means
end
freshwater = interval_means(monitor_file, "surface_freshwater_Sv")

# ## Table

@printf("\nDomain area %.3e m²: 1 Sv = %.3f cm/day of domain-mean sea level\n\n", domain_area, 1 / Sv_per_cm_per_day)
@printf("%-23s %9s %9s %9s %9s %10s %9s %9s\n", "interval", "Atl in", "Pac out", "net Q", "A dη/dt", "D", "D", "model F")
@printf("%-23s %9s %9s %9s %9s %10s %9s %9s\n", "", "[Sv]", "[Sv]", "[Sv]", "[Sv]", "[Sv]", "[cm/day]", "[Sv]")
for n in 1:N-1
    @printf("%-23s %9.3f %9.3f %9.3f %9.3f %10.3f %9.2f %9.3f\n", "$(Date(dates[n])) → $(Date(dates[n+1]))",
            (atlantic[n] + atlantic[n+1]) / 2, (pacific[n] + pacific[n+1]) / 2, mean_net[n], volume_rate[n],
            D[n], D[n] / Sv_per_cm_per_day, freshwater[n])
end

@printf("\nD over %d days: mean %.2f Sv, std %.2f Sv, range %.2f to %.2f Sv\n", N - 1, mean(D), std(D), extrema(D)...)
covered = findall(isfinite, freshwater)
isempty(covered) ||
    @printf("Over the %d days the run covers: mean D %.2f Sv, mean model freshwater flux %.2f Sv, difference %.2f Sv\n",
            length(covered), mean(D[covered]), mean(freshwater[covered]), mean(D[covered] .- freshwater[covered]))
println("Change of ⟨ηᴳ⟩ over the record: ", @sprintf("%.2f cm", 100 * (mean_η[end] - mean_η[1])),
        "; implied by the edge transports alone: ", @sprintf("%.2f cm", sum(mean_net) / Sv_per_cm_per_day))

# Autocorrelation of D: a peak near 7 days would point to GLORYS's weekly assimilation cycle
if N - 1 ≥ 16
    d = D .- mean(D)
    print("Autocorrelation of D at lags 1–10 days:")
    for lag in 1:10
        @printf(" %.2f", sum(d[1:end-lag] .* d[1+lag:end]) / sum(abs2, d))
    end
    println()
end

# ## CSV

csv_file = joinpath(output_dir, "glorys_volume_budget.csv")
open(csv_file, "w") do io
    println(io, "date,atlantic_native_Sv,pacific_native_Sv,net_native_Sv,glorys_eta_mean_m,",
                "interval_volume_rate_Sv,interval_net_Sv,interval_D_Sv,interval_model_freshwater_Sv,",
                "atlantic_in_Sv,atlantic_out_Sv,pacific_out_Sv,pacific_in_Sv")
    for n in 1:N
        k = min(n, N - 1)   # interval n → n+1 (the last row repeats the last interval)
        @printf(io, "%s,%.4f,%.4f,%.4f,%.5f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.3f,%.3f\n", Date(dates[n]), atlantic[n],
                pacific[n], net[n], mean_η[n], volume_rate[k], mean_net[k], D[k], freshwater[k],
                atlantic_gross[n, 1], atlantic_gross[n, 2], pacific_gross[n, 1], pacific_gross[n, 2])
    end
end

# ## Figure

day      = [(n - start_index) for n in 1:N]                  # days from the start date
midpoint = [(n - start_index) + 0.5 for n in 1:N-1]

fig = Figure(size = (1100, 1000), fontsize = 15)

ax1 = Axis(fig[1, 1]; ylabel = "Sv", title = "(a) GLORYS's transport through the edges, on its own grid")
lines!(ax1, day, atlantic; color = :steelblue, label = "Atlantic, into the domain")
lines!(ax1, day, pacific;  color = :darkorange, label = "Pacific, out of the domain")
lines!(ax1, day, net;      color = :black, label = "net inflow")
axislegend(ax1; position = :lb, framevisible = false)

ax2 = Axis(fig[2, 1]; ylabel = "cm", title = "(b) change of the domain-mean free surface")
lines!(ax2, day, 100 .* (mean_η .- mean_η[1]); color = :black, label = "GLORYS ⟨ηᴳ⟩")
lines!(ax2, day, vcat(0, cumsum(mean_net)) ./ Sv_per_cm_per_day; color = :gray50, linestyle = :dash,
       label = "implied by GLORYS's edge transports")
if !isempty(covered)
    implied = vcat(0, cumsum(mean_net .+ ifelse.(isfinite.(freshwater), freshwater, 0.0))) ./ Sv_per_cm_per_day
    lines!(ax2, day, implied; color = :firebrick, linestyle = :dot, label = "edge transports + model freshwater")
end
axislegend(ax2; position = :lb, framevisible = false)

ax3 = Axis(fig[3, 1]; xlabel = "days from $(Date(start_date))", ylabel = "Sv",
           title = "(c) D = A d⟨ηᴳ⟩/dt − net transport, daily")
barplot!(ax3, midpoint, D; color = :gray60, label = "D")
hlines!(ax3, [mean(D)]; color = :black, linestyle = :dash, label = @sprintf("mean %.2f Sv", mean(D)))
isempty(covered) || scatterlines!(ax3, midpoint[covered], freshwater[covered]; color = :firebrick,
                                  label = "model surface freshwater flux")
axislegend(ax3; position = :lb, framevisible = false)
linkxaxes!(ax1, ax2, ax3)

figure_file = joinpath(output_dir, "glorys_volume_budget.png")
save(figure_file, fig)
@info "Wrote $csv_file and $figure_file"
