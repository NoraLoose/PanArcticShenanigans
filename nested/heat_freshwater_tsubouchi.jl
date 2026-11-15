# Volume, heat and freshwater through Davis, Fram, the BSO and Bering for the model and GLORYS, on Tsubouchi et al.'s
# (2018) reference values, against Tsubouchi et al.: first over their monthly blocks that the run covers (same days),
# then the run mean against their annual mean (Table 3).
#
#   julia --project=. heat_freshwater_tsubouchi.jl <output folder> [start date, default 2005-09-01]
#
# Reads <folder>/gateways.csv (model, hourly) and, if present, <folder>/glorys_gateways.csv (GLORYS on its own grid,
# daily; from glorys_gateways_native.jl). No new run needed. Everything is positive into the Arctic.
#
# Heat and freshwater per gateway depend on the reference values; the net over a closed boundary does not. Ours:
# heat relative to 0 degC (Bering -1.9 degC), freshwater relative to 34.8. Tsubouchi et al.: 1.01 degC and 34.67.
# From the volume Q, our heat H and our freshwater F:
#     heat on their reference        = H - rho0 cp (1.01 - T_ours) Q
#     freshwater on their reference  = Q - (34.8 / 34.67) (Q - F)        (F = Q - sum(v S)/34.8)
# Tsubouchi et al.'s monthly values are read from their Figs. 5a, 7a, 9a (references/tsubouchi2018_monthly_*.csv):
# about +-0.2 Sv, +-5 TW, +-5 mSv. Their "months" are 30-day blocks from 5 September 2005.

using Printf, Dates

folder = length(ARGS) ≥ 1 ? ARGS[1] : error("usage: julia --project=. heat_freshwater_tsubouchi.jl <output folder> [start date]")
start  = length(ARGS) ≥ 2 ? Date(ARGS[2]) : Date(2005, 9, 1)

const ρ₀cₚ  = 1026 * 3991.9                       # J m⁻³ K⁻¹, as in gateways.jl
const S_ours, S_tsu = 34.8, 34.67
const T_tsu  = 1.01
const T_ours = Dict("fram" => 0.0, "bso" => 0.0, "davis" => 0.0, "bering" => -1.9)   # as in gateways.jl
const gateways = ["davis", "fram", "bso", "bering"]
const references = joinpath(@__DIR__, "references")

function read_columns(file)
    lines = filter(l -> !isempty(strip(l)) && !startswith(l, "#"), readlines(file))
    header = [strip(name, ['(', ')', '"', ' ']) for name in split(lines[1], ',')]
    rows = [split(l, ',') for l in lines[2:end]]
    return name -> begin
        k = findfirst(==(name), header)
        isnothing(k) && error("$file has no column $name")
        [something(tryparse(Float64, r[k]), NaN) for r in rows]
    end
end

finite_mean(x) = (y = filter(isfinite, x); isempty(y) ? NaN : sum(y) / length(y))

# Means over days [a, b) since the run's start, converted to Tsubouchi et al.'s references: (Q Sv, H TW, F mSv)
function converted(column, g, a, b)
    day = column("hour") ./ 24
    inside = (a .≤ day .< b)
    Q = finite_mean(column(g * "_volume_Sv")[inside])
    H = finite_mean(column(g * "_heat_TW")[inside])
    F = finite_mean(column(g * "_freshwater_mSv")[inside])
    H_tsu = H - ρ₀cₚ * (T_tsu - T_ours[g]) * Q * 1e6 / 1e12
    F_tsu = 1000Q - (S_ours / S_tsu) * (1000Q - F)
    return (; Q, H = H_tsu, F = F_tsu)
end

model  = read_columns(joinpath(folder, "gateways.csv"))
glorys = isfile(joinpath(folder, "glorys_gateways.csv")) ? read_columns(joinpath(folder, "glorys_gateways.csv")) : nothing
none   = (Q = NaN, H = NaN, F = NaN)
run_days = maximum(filter(isfinite, model("hour"))) / 24

monthly = (Q = read_columns(joinpath(references, "tsubouchi2018_monthly_volume.csv")),
           H = read_columns(joinpath(references, "tsubouchi2018_monthly_heat.csv")),
           F = read_columns(joinpath(references, "tsubouchi2018_monthly_freshwater.csv")))
annual = read_columns(joinpath(references, "tsubouchi2018_table3.csv"))
row(name) = findfirst(==(name), ["davis", "fram", "bso", "bering", "net_ocean"])

const quantities = ((:Q, "volume (Sv)", "%6.2f"), (:H, "heat (TW)", "%6.0f"), (:F, "freshwater (mSv)", "%6.0f"))
fmt(f, x) = Printf.format(Printf.Format(f), x)

@printf("\n%s, %s: %.0f days. Positive into the Arctic. Heat relative to %.2f degC, freshwater to %.2f.\n",
        folder, start, run_days, T_tsu, S_tsu)

# 1. Tsubouchi et al.'s monthly blocks covered by the run (at least 10 days)
first_block = Date(2005, 9, 5)
for b in 1:12
    a = Dates.value(first_block + Day(30(b - 1)) - start)
    e = a + (b == 10 ? 35 : 30)
    covered = min(e, run_days) - max(a, 0)
    covered < 10 && continue
    println("\nBlock $b ($(first_block + Day(30(b - 1))), $(round(Int, covered)) of $(e - a) days in the run): model / GLORYS / Tsubouchi et al.")
    println("gateway   ", join([rpad(q[2], 26) for q in quantities]))
    sums = Dict(k => zeros(3) for k in (:Q, :H, :F))
    for g in gateways
        m = converted(model, g, a, e)
        G = isnothing(glorys) ? none : converted(glorys, g, a, e)
        print(rpad(g, 8))
        for (k, _, f) in quantities
            t = monthly[k](g)[b]
            sums[k] .+= [getfield(m, k), getfield(G, k), t]
            print(fmt(f, getfield(m, k)), fmt(f, getfield(G, k)), fmt(f, t), "        ")
        end
        println()
    end
    print(rpad("sum", 8))
    for (k, _, f) in quantities
        print(fmt(f, sums[k][1]), fmt(f, sums[k][2]), fmt(f, sums[k][3]), "        ")
    end
    println()
end

# 2. Run mean against the annual mean (Table 3)
println("\nRun mean against Tsubouchi et al.'s annual mean ± std over Sep 2005 - Aug 2006 (seasons differ: a check of")
println("sign and size only): model / GLORYS / Tsubouchi et al.")
println("gateway   ", join([rpad(q[2], 34) for q in quantities]))
sums = Dict(k => zeros(2) for k in (:Q, :H, :F))
for name in [gateways; "sum"]
    print(rpad(name, 8))
    for (k, _, f) in quantities
        col = Dict(:Q => "volume", :H => "heat", :F => "freshwater")[k]
        unit = Dict(:Q => "_Sv", :H => "_TW", :F => "_mSv")[k]
        r = row(name == "sum" ? "net_ocean" : name)
        mean_, std_ = annual(col * unit)[r], annual(col * "_std")[r]
        if name == "sum"
            m, G = sums[k]
        else
            m = getfield(converted(model, name, 0, run_days + 1), k)
            G = isnothing(glorys) ? NaN : getfield(converted(glorys, name, 0, run_days + 1), k)
            sums[k] .+= [m, G]
        end
        print(fmt(f, m), fmt(f, G), fmt(f, mean_), " ± ", rpad(fmt(k == :Q ? "%4.2f" : "%3.0f", std_), 7), "    ")
    end
    println()
end
println("""
Tsubouchi et al.'s "sum" is their ocean-only net through the whole boundary (also the small Fury and Hecla Strait).
Their ± is the standard deviation of the twelve monthly values (the seasonal cycle), not a measurement error.""")
