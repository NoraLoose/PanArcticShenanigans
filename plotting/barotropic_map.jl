# Daily-mean depth-integrated transport |U| (m² s⁻¹) and free surface around Iceland and the Nordic Seas, model
# (barotropic_daily.jld2, written when [diagnostics] sections = true) against GLORYS, both on the model grid. η is
# drawn as contours on top of |U|: a balanced barotropic flow follows them. Model and GLORYS η are each relative to
# their own mean over the window, so a datum difference drops out. Also printed: how far the model's and GLORYS's η′
# have each moved from the common initial state (GLORYS at the start), and a map of model − GLORYS η′.
#
#   julia --project=. plotting/barotropic_map.jl <config.toml> [day ...]
#
# The config is the run's own (it sets the dates and the output directory). Days are the model's daily means, counted
# from the start (1 = the first 24 hours; default 1 2 3 4). Each model mean (noon to noon) is compared with GLORYS at
# the centre of its window, i.e. the average of the two GLORYS daily means either side. Writes barotropic_map.png in
# the output directory.

include("../setup.jl")

using Printf, Statistics, CairoMakie
using Oceananigans.Grids: λnode, φnode, static_column_depthᶜᶜᵃ
using NumericalEarth.DataWrangling: NCDatasets, metadata_path
using .NCDatasets: Dataset, dimnames

days_to_plot = length(ARGS) ≥ 2 ? parse.(Int, ARGS[2:end]) : [1, 2, 3, 4]
irange, jrange = 200:380, 60:240     # Denmark Strait to Fram Strait and the Barents Sea Opening (15 km grid)

# ## Model daily means

file = joinpath(output_dir, "barotropic_daily.jld2")
U = FieldTimeSeries(file, "U")
V = FieldTimeSeries(file, "V")
η = FieldTimeSeries(file, "η")
mgrid = U.grid
averages = [n for n in eachindex(U.times) if U.times[n] > start_time + 12hours]   # records that are 24-hour means
day_of(n) = round(Int, (U.times[n] - start_time) / 1days)

λ = [λnode(i, j, 1, mgrid, Center(), Center(), Center()) for i in irange, j in jrange]
φ = [φnode(i, j, 1, mgrid, Center(), Center(), Center()) for i in irange, j in jrange]
wet = try
    [static_column_depthᶜᶜᵃ(i, j, mgrid) > 0 for i in irange, j in jrange]
catch                                     # if the saved grid has no bathymetry, land is where η is exactly zero
    interior(η[end])[irange, jrange, 1] .!= 0
end
masked(a) = ifelse.(wet, a, NaN)

# Speed of the depth-integrated flow at cell centres (an orthogonal grid, so no rotation is needed for the magnitude)
function model_speed(n)
    u = interior(U[n])[:, :, 1]
    v = interior(V[n])[:, :, 1]
    return masked([hypot((u[i, j] + u[i+1, j]) / 2, (v[i, j] + v[i, j+1]) / 2) for i in irange, j in jrange])
end
model_η(n) = masked(interior(η[n])[irange, jrange, 1])

# ## GLORYS on its own grid, sampled at the model's cell centres

band = (minimum(φ) - 1, maximum(φ) + 1)

function glorys_fields(date)
    U₂ = map((:u_velocity, :v_velocity)) do name
        ds  = Dataset(metadata_path(Metadatum(name; date, dataset, region, dir_kw...)))
        var = ds[name === :u_velocity ? "uo" : "vo"]
        dimnames(var)[1:3] == ("longitude", "latitude", "depth") || error("unexpected dimensions $(dimnames(var))")
        lon, lat, z = Float64.(ds["longitude"][:]), Float64.(ds["latitude"][:]), Float64.(ds["depth"][:])
        rows = findall(y -> band[1] ≤ y ≤ band[2], lat)
        u = coalesce.(var[:, rows, :, 1], 0f0)
        close(ds)
        interfaces = vcat(0, (z[1:end-1] .+ z[2:end]) ./ 2, z[end] + (z[end] - z[end-1]) / 2)
        (lon, lat[rows], dropdims(sum(Float64.(u) .* reshape(diff(interfaces), 1, 1, :), dims = 3), dims = 3))
    end
    ds  = Dataset(metadata_path(Metadatum(:free_surface; date, dataset, region, dir_kw...)))
    lat = Float64.(ds["latitude"][:])
    rows = findall(y -> band[1] ≤ y ≤ band[2], lat)
    zos = Float64.(coalesce.(ds["zos"][:, rows, 1], NaN32))
    lon = Float64.(ds["longitude"][:])
    close(ds)
    return U₂[1], U₂[2], (lon, lat[rows], zos)
end

# Bilinear sample on a regular (lon, lat) grid, periodic in longitude; NaN where any of the four points is NaN
function sample(F, lon, lat, x, y)
    fi = mod(x - lon[1], 360) / (lon[2] - lon[1])
    fj = (y - lat[1]) / (lat[2] - lat[1])
    i₀ = floor(Int, fi); j₀ = clamp(floor(Int, fj), 0, length(lat) - 2)
    a, b = fi - i₀, clamp(fj - j₀, 0, 1)
    i₁, i₂ = mod(i₀, length(lon)) + 1, mod(i₀ + 1, length(lon)) + 1
    return (1 - a) * (1 - b) * F[i₁, j₀+1] + a * (1 - b) * F[i₂, j₀+1] + (1 - a) * b * F[i₁, j₀+2] + a * b * F[i₂, j₀+2]
end

glorys_cache = Dict{Int, Any}()
function glorys_on_model(k)   # k-th GLORYS date: |U| and η at the model's cell centres
    get!(glorys_cache, k) do
        (lu, au, Uᵉ), (_, _, Uⁿ), (lz, az, zos) = glorys_fields(dates[k])
        speed = masked([hypot(sample(Uᵉ, lu, au, x, y), sample(Uⁿ, lu, au, x, y)) for (x, y) in zip(λ, φ)])
        surface = masked([sample(zos, lz, az, x, y) for (x, y) in zip(λ, φ)])
        (speed, surface)
    end
end

# GLORYS at the centre of a model window: GLORYS daily means are at noon on dates[k], i.e. at (k − 1) days + 12 h
function glorys_at(t)
    s = clamp((t - 12hours) / 1days, 0, length(dates) - 1)
    k = min(floor(Int, s) + 1, length(dates) - 1)
    α = s - (k - 1)
    (s₁, η₁), (s₂, η₂) = glorys_on_model(k), glorys_on_model(k + 1)
    return (1 - α) .* s₁ .+ α .* s₂, (1 - α) .* η₁ .+ α .* η₂
end

anomaly(a) = a .- mean(filter(isfinite, a))

# ## Figure: one row per day; model |U|, GLORYS |U|, model − GLORYS |U|; η contours every 5 cm

# The model starts from GLORYS at the start time, so GLORYS there is the common initial η′
η₀ = anomaly(last(glorys_at(start_time)))

fig = Figure(size = (2000, 430 * length(days_to_plot)), fontsize = 13)
speed_kw = (colormap = :speed, colorrange = (0, 600), nan_color = :gray80)
levels = -0.5:0.05:0.5
for (r, d) in enumerate(days_to_plot)
    nidx = findfirst(n -> day_of(n) == d, averages)
    isnothing(nidx) && error("no daily mean for day $d in $file (days available: $(day_of.(averages)))")
    n  = averages[nidx]
    tc = U.times[n] - 12hours
    Sm, ηm = model_speed(n), anomaly(model_η(n))
    Sg, ηg = glorys_at(tc)
    ηg = anomaly(ηg)
    rmscm(a) = 100 * sqrt(mean(filter(isfinite, a .^ 2)))
    @printf("day %d: max |U| model %.0f m²/s, GLORYS %.0f m²/s; mean |U| model %.1f, GLORYS %.1f; rms η′ model − GLORYS %.1f cm; rms change of η′ since the start: model %.1f cm, GLORYS %.1f cm\n",
            d, maximum(filter(isfinite, Sm)), maximum(filter(isfinite, Sg)), mean(filter(isfinite, Sm)),
            mean(filter(isfinite, Sg)), rmscm(ηm .- ηg), rmscm(ηm .- η₀), rmscm(ηg .- η₀))

    ax1 = Axis(fig[r, 1], title = "day $d, model |U| (m² s⁻¹), η′ contours", ylabel = "j", aspect = DataAspect())
    hm = heatmap!(ax1, irange, jrange, Sm; speed_kw...)
    contour!(ax1, irange, jrange, ηm; levels, color = :black, linewidth = 0.6)
    ax2 = Axis(fig[r, 2], title = "day $d, GLORYS |U|, η′ contours", aspect = DataAspect())
    heatmap!(ax2, irange, jrange, Sg; speed_kw...)
    contour!(ax2, irange, jrange, ηg; levels, color = :black, linewidth = 0.6)
    Colorbar(fig[r, 3], hm)
    ax3 = Axis(fig[r, 4], title = "day $d, model − GLORYS |U|", aspect = DataAspect())
    hd = heatmap!(ax3, irange, jrange, Sm .- Sg; colormap = :balance, colorrange = (-300, 300), nan_color = :gray80)
    Colorbar(fig[r, 5], hd)
    ax4 = Axis(fig[r, 6], title = "day $d, η′ model − GLORYS (m)", aspect = DataAspect())
    he = heatmap!(ax4, irange, jrange, ηm .- ηg; colormap = :balance, colorrange = (-0.2, 0.2), nan_color = :gray80)
    Colorbar(fig[r, 7], he)
    r == length(days_to_plot) && (ax1.xlabel = ax2.xlabel = ax3.xlabel = ax4.xlabel = "i")
end

figure_file = joinpath(output_dir, "barotropic_map.png")
save(figure_file, fig)
@info "Wrote $figure_file (window i = $irange, j = $jrange; η contours every 5 cm)"
