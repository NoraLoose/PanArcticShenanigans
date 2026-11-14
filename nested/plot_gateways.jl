# Hourly (thin line) and daily-mean (dots) volume, heat and freshwater transport through Fram Strait, the Barents Sea Opening, Bering Strait and
# Davis Strait (gateways.csv), positive into the Arctic, against the mooring-based ranges compiled by Shu et al. (2023,
# GMD 16, 2539, Table 3). The observed ranges are multi-year annual means; the run covers days to weeks of one winter.
# The bottom row is the sum of the four volume transports: the net inflow to the Arctic, which only a change in its
# volume can balance (Tsubouchi et al. 2018: −0.2 Sv). Daily means and run means use whole days only.
#
#   julia --project=.. plot_gateways.jl output-glorys-sqrtH-4d      # writes <folder>/gateways.png

using CairoMakie, Printf

folder = isempty(ARGS) ? "output-nested" : ARGS[1]
text   = readlines(joinpath(folder, "gateways.csv"))
header = [strip(name, ['(', ')', '"', ' ']) for name in split(text[1], ',')]   # also reads the first version's header
rows   = [parse.(Float64, split(line, ',')) for line in text[2:end] if !isempty(line)]
column(name) = [row[findfirst(==(name), header)] for row in rows]

hour = column("hour")
day  = floor.(Int, hour ./ 24)
days = [d for d in unique(day) if count(==(d), day) ≥ 24]
whole = in(days).(day)
daily(x) = [sum(x[day .== d]) / count(==(d), day) for d in days]
mean_over_run(x) = sum(x[whole]) / count(whole)

# (low, high) of the observed annual means; freshwater converted from 10³ km³ yr⁻¹ (1 mSv = 31.56 km³ yr⁻¹)
observed = Dict(
    "fram"   => (volume = (-4.7, 0.7),  heat = (26, 50),    freshwater = (-102.4, -68.8)),  # Schauer et al. 2008; Serreze et al. 2006 (heat: Atlantic Water > 1 °C only)
    "bso"    => (volume = (1.5, 2.5),   heat = (48, 74),    freshwater = (-5.0, -1.0)),     # Smedsrud et al. 2010 (about 2 Sv; ±0.5 assumed); Skagseth et al. 2008; Serreze et al. 2006 (−2.9 mSv)
    "bering" => (volume = (0.6, 1.05),  heat = (9.5, 19.0), freshwater = (66.5, 110.9)),    # Woodgate 2018 (heat relative to −1.9 °C)
    "davis"  => (volume = (-2.1, -1.1), heat = (1, 35),     freshwater = (-98.9, -86.8)),   # Curry et al. 2014; Cuny et al. 2005
)

quantities = (("volume", "_volume_Sv", "Sv"), ("heat", "_heat_TW", "TW"), ("freshwater", "_freshwater_mSv", "mSv"))
gateways   = ("fram", "bso", "bering", "davis")

fig = Figure(size = (1200, 1200))
Label(fig[0, 1:3], basename(rstrip(folder, '/')) * ": hourly (line) and daily means (dots) against observed annual means (band)"; fontsize = 18)
for (r, g) in enumerate(gateways), (c, (quantity, suffix, unit)) in enumerate(quantities)
    ax = Axis(fig[r, c]; xlabel = r == 4 && c > 1 ? "day" : "", ylabel = c == 1 ? "$g ($unit)" : unit, title = r == 1 ? quantity : "")
    low, high = getproperty(observed[g], Symbol(quantity))
    hspan!(ax, low, high; color = (:gray, 0.3))
    hlines!(ax, 0; color = :black, linewidth = 0.5)
    lines!(ax, hour ./ 24, column(g * suffix); color = (:red, 0.4))
    scatter!(ax, days .+ 0.5, daily(column(g * suffix)); color = :red)
end

net = sum(column(g * "_volume_Sv") for g in gateways)
ax = Axis(fig[5, 1]; xlabel = "day", ylabel = "sum of the four (Sv)")
hlines!(ax, 0; color = :black, linewidth = 0.5)
hlines!(ax, -0.2; color = :gray, linestyle = :dash)
lines!(ax, hour ./ 24, net; color = (:red, 0.4))
scatter!(ax, days .+ 0.5, daily(net); color = :red)

for g in gateways
    @printf("%-6s mean over %d days: volume %+6.2f Sv, heat %+7.1f TW, freshwater %+7.1f mSv\n", g, length(days),
            mean_over_run(column(g * "_volume_Sv")), mean_over_run(column(g * "_heat_TW")), mean_over_run(column(g * "_freshwater_mSv")))
end
@printf("sum of the four volume transports over %d days: %+6.2f Sv (Tsubouchi et al. 2018: -0.2 Sv)\n", length(days), mean_over_run(net))

file = joinpath(folder, "gateways.png")
save(file, fig)
@info "Wrote $file"
