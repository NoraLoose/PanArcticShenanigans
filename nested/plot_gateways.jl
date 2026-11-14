# Volume, heat and freshwater transport through Fram Strait, the Barents Sea Opening, Bering Strait and Davis Strait,
# positive into the Arctic:
#   - the model (<folder>/gateways.csv): hourly (thin red line) and daily means (red dots);
#   - GLORYS on its own grid (glorys_gateways_native.jl): daily means through the true lines (black diamonds) and
#     through the model's staircase (open circles);
#   - the mooring-based annual-mean ranges compiled by Shu et al. (2023, GMD 16, 2539, Table 3) (gray band);
#   - with a run file whose dates fall in September 2005 – August 2006: the monthly volume transports of Tsubouchi et
#     al. (2018, JPO 48, 2029; Fig. 5a, references/tsubouchi2018_monthly_volume.csv) as orange bars over their 30-day
#     blocks, with ±0.2 Sv for reading the figure.
# Bottom row: the sum of the four (the net inflow to the Arctic; Tsubouchi et al. 2018: −0.2 Sv) and, with
# glorys_arctic_budget.csv, GLORYS's Arctic volume budget (inflow, A d⟨η⟩/dt and their difference).
# Works without gateways.csv (GLORYS and observations only). Daily means and run means use whole days only.
#
#   julia --project=.. plot_gateways.jl output-glorys-sqrtH-4d [runs/glorys-sqrtH-4d.toml]   # writes <folder>/gateways.png

using CairoMakie, Printf, Dates, TOML

folder   = isempty(ARGS) ? "output-nested" : ARGS[1]
run_file = length(ARGS) ≥ 2 ? ARGS[2] : nothing

function read_columns(file)
    isfile(file) || return nothing
    text   = filter(line -> !startswith(line, '#'), readlines(file))
    header = [strip(name, ['(', ')', '"', ' ']) for name in split(text[1], ',')]   # also reads the first version's header
    rows   = [split(line, ',') for line in text[2:end] if !isempty(line)]
    return name -> begin
        k = findfirst(==(name), header)
        [parse(Float64, row[k]) for row in rows]
    end
end

model     = read_columns(joinpath(folder, "gateways.csv"))
glorys    = read_columns(joinpath(folder, "glorys_gateways.csv"))
modelline = read_columns(joinpath(folder, "glorys_gateways_modelline.csv"))   # GLORYS through the model's staircase
budget    = read_columns(joinpath(folder, "glorys_arctic_budget.csv"))
isnothing(model) && isnothing(glorys) && error("neither gateways.csv nor glorys_gateways.csv in $folder")

mean(x) = sum(x) / length(x)

# Model: daily means over whole days
if !isnothing(model)
    hour  = model("hour")
    day   = floor.(Int, hour ./ 24)
    days  = [d for d in unique(day) if count(==(d), day) ≥ 24]
    whole = in(days).(day)
    daily(x) = [sum(x[day .== d]) / count(==(d), day) for d in days]
    mean_over_run(x) = sum(x[whole]) / count(whole)
end

# Tsubouchi et al. (2018): 30-day blocks from 5 September 2005, the tenth 35 days; days since the run's start
tsubouchi = nothing
if !isnothing(run_file)
    start = Date(TOML.parsefile(run_file)["time"]["start_date"])
    reference = read_columns(joinpath(@__DIR__, "references", "tsubouchi2018_monthly_volume.csv"))
    lengths = [30, 30, 30, 30, 30, 30, 30, 30, 30, 35, 30, 30]
    first_day = Date(2005, 9, 5) .+ Day.(cumsum([0; lengths[1:end-1]]))
    blocks = [(Dates.value(first_day[k] - start), Dates.value(first_day[k] - start) + lengths[k]) for k in 1:12]
    run_length = isnothing(model) ? maximum(glorys("hour")) / 24 + 0.5 : maximum(hour) / 24
    overlapping = [k for k in 1:12 if blocks[k][2] > 0 && blocks[k][1] < run_length]
    if isempty(overlapping)
        @info "The run ($start) does not overlap Tsubouchi et al.'s year (September 2005 – August 2006)"
    else
        tsubouchi = (; blocks, overlapping, value = g -> reference(g))
    end
end

# (low, high) of the observed annual means; freshwater converted from 10³ km³ yr⁻¹ (1 mSv = 31.56 km³ yr⁻¹)
observed = Dict(
    "fram"   => (volume = (-4.7, 0.7),  heat = (26, 50),    freshwater = (-102.4, -68.8)),  # Schauer et al. 2008; Serreze et al. 2006 (heat: Atlantic Water > 1 °C only)
    "bso"    => (volume = (1.5, 2.5),   heat = (48, 74),    freshwater = (-5.0, -1.0)),     # Smedsrud et al. 2010 (about 2 Sv; ±0.5 assumed); Skagseth et al. 2008; Serreze et al. 2006 (−2.9 mSv)
    "bering" => (volume = (0.6, 1.05),  heat = (9.5, 19.0), freshwater = (66.5, 110.9)),    # Woodgate 2018 (heat relative to −1.9 °C)
    "davis"  => (volume = (-2.1, -1.1), heat = (1, 35),     freshwater = (-98.9, -86.8)),   # Curry et al. 2014; Cuny et al. 2005
)

quantities = (("volume", "_volume_Sv", "Sv"), ("heat", "_heat_TW", "TW"), ("freshwater", "_freshwater_mSv", "mSv"))
gateways   = ("fram", "bso", "bering", "davis")

function draw!(ax, series)   # series(g) → the values to draw for gateway g (or a sum), for each source
    if !isnothing(model)
        x = series(model)
        lines!(ax, hour ./ 24, x; color = (:red, 0.4))
        scatter!(ax, days .+ 0.5, daily(x); color = :red)
    end
    isnothing(glorys) || scatter!(ax, glorys("hour") ./ 24, series(glorys); color = :black, marker = :diamond)
    isnothing(modelline) || scatter!(ax, modelline("hour") ./ 24, series(modelline);
                                     color = :transparent, strokecolor = :black, strokewidth = 1, marker = :circle)
end

function draw_tsubouchi!(ax, values)
    isnothing(tsubouchi) && return nothing
    for k in tsubouchi.overlapping
        a, b = tsubouchi.blocks[k]
        poly!(ax, Rect(a, values[k] - 0.2, b - a, 0.4); color = (:orange, 0.35))
        lines!(ax, [a, b], [values[k], values[k]]; color = :darkorange, linewidth = 2)
    end
end

title = basename(rstrip(folder, '/')) * ": model hourly (line) and daily (red), GLORYS daily (black; open: model staircase), " *
        "observed annual range (gray)" * (isnothing(tsubouchi) ? "" : ", Tsubouchi et al. 2018 monthly (orange)")
fig = Figure(size = (1200, 1200))
Label(fig[0, 1:3], title; fontsize = 16)
for (r, g) in enumerate(gateways), (c, (quantity, suffix, unit)) in enumerate(quantities)
    ax = Axis(fig[r, c]; xlabel = r == 4 && c > 1 ? "day" : "", ylabel = c == 1 ? "$g ($unit)" : unit, title = r == 1 ? quantity : "")
    low, high = getproperty(observed[g], Symbol(quantity))
    hspan!(ax, low, high; color = (:gray, 0.3))
    hlines!(ax, 0; color = :black, linewidth = 0.5)
    quantity == "volume" && !isnothing(tsubouchi) && draw_tsubouchi!(ax, tsubouchi.value(g))
    draw!(ax, source -> source(g * suffix))
end

ax = Axis(fig[5, 1]; xlabel = "day", ylabel = "sum of the four (Sv)")
hlines!(ax, 0; color = :black, linewidth = 0.5)
hlines!(ax, -0.2; color = :gray, linestyle = :dash)
isnothing(tsubouchi) || draw_tsubouchi!(ax, sum(tsubouchi.value(g) for g in gateways))
draw!(ax, source -> sum(source(g * "_volume_Sv") for g in gateways))

if !isnothing(budget)
    ax = Axis(fig[5, 2:3]; xlabel = "day", ylabel = "GLORYS Arctic budget (Sv)")
    t = budget("hour") ./ 24
    hlines!(ax, 0; color = :black, linewidth = 0.5)
    hlines!(ax, -0.17; color = :gray, linestyle = :dash)   # −(runoff + P − E)
    scatterlines!(ax, t, budget("gateways_Sv") .+ budget("fury_hecla_Sv"); color = :black, label = "inflow (gateways + Fury and Hecla)")
    scatterlines!(ax, t, budget("dVdt_Sv"); color = :blue, label = "A d⟨η⟩/dt")
    scatterlines!(ax, t, budget("residual_Sv"); color = :purple, label = "residual (expected ≈ −0.17)")
    axislegend(ax; position = :lt, labelsize = 11)
end

# ## Printed summary

if !isnothing(model)
    for g in gateways
        @printf("%-6s mean over %d days: volume %+6.2f Sv, heat %+7.1f TW, freshwater %+7.1f mSv\n", g, length(days),
                mean_over_run(model(g * "_volume_Sv")), mean_over_run(model(g * "_heat_TW")), mean_over_run(model(g * "_freshwater_mSv")))
    end
    @printf("sum of the four volume transports over %d days: %+6.2f Sv (Tsubouchi et al. 2018: -0.2 Sv)\n", length(days),
            mean_over_run(sum(model(g * "_volume_Sv") for g in gateways)))
end
if !isnothing(glorys)
    for g in gateways
        @printf("GLORYS %-6s mean: volume %+6.2f Sv, heat %+7.1f TW, freshwater %+7.1f mSv\n", g,
                mean(glorys(g * "_volume_Sv")), mean(glorys(g * "_heat_TW")), mean(glorys(g * "_freshwater_mSv")))
    end
    @printf("GLORYS sum of the four: %+6.2f Sv\n", mean(sum(glorys(g * "_volume_Sv") for g in gateways)))
end
if !isnothing(budget)
    @printf("GLORYS Arctic budget: mean residual %+6.3f Sv (expected about -0.17), range %+6.3f to %+6.3f\n",
            mean(budget("residual_Sv")), extrema(budget("residual_Sv"))...)
end

# Means over each of Tsubouchi et al.'s blocks, from the days of the run inside it (volume, Sv)
if !isnothing(tsubouchi)
    println("\nVolume transport (Sv) over Tsubouchi et al.'s 30-day blocks: observed / GLORYS / model (days of the run inside the block)")
    for k in tsubouchi.overlapping
        a, b = tsubouchi.blocks[k]
        inside(t) = [a ≤ x < b for x in t]
        line = @sprintf("block %2d (days %4d to %4d)", k, a, b)
        for g in gateways
            G = isnothing(glorys) ? NaN : (m = inside(glorys("hour") ./ 24); any(m) ? mean(glorys(g * "_volume_Sv")[m]) : NaN)
            M = isnothing(model)  ? NaN : (m = inside(days .+ 0.5); any(m) ? mean(daily(model(g * "_volume_Sv"))[m]) : NaN)
            line *= @sprintf("  %s %+5.1f/%+5.2f/%+5.2f", g, tsubouchi.value(g)[k], G, M)
        end
        nG = isnothing(glorys) ? 0 : count(inside(glorys("hour") ./ 24))
        println(line, "  ($nG days)")
    end
end

file = joinpath(folder, "gateways.png")
save(file, fig)
@info "Wrote $file"
