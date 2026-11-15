# Heat and freshwater through Fram, the BSO, Bering and Davis on Tsubouchi et al.'s (2018) reference values, for the
# model and for GLORYS, against Tsubouchi et al.'s Table 3 (annual means, Sep 2005 - Aug 2006).
#
#   julia --project=. heat_freshwater_tsubouchi.jl <output folder>
#
# Reads <folder>/gateways.csv (model, hourly) and, if present, <folder>/glorys_gateways.csv (GLORYS on its own grid,
# daily; from glorys_gateways_native.jl). No new run needed.
#
# Heat and freshwater per gateway depend on the reference values; the net over a closed boundary does not. Ours:
# heat relative to 0 degC (Bering -1.9 degC), freshwater relative to 34.8. Tsubouchi et al.: 1.01 degC and 34.67.
# From the volume Q, our heat H and our freshwater F (all positive into the Arctic):
#     heat on their reference        = H - rho0 cp (1.01 - T_ours) Q
#     freshwater on their reference  = Q - (34.8 / 34.67) (Q - F)
# (since F = Q - (1/34.8) sum(v S), so sum(v S) = 34.8 (Q - F)).
#
# Caveat: the run covers part of one year, Tsubouchi et al.'s values are annual means. Heat transport in particular has
# a large seasonal cycle, so this is a check of sign and size; the real test is the full year.

using Printf

folder = length(ARGS) ≥ 1 ? ARGS[1] : error("usage: julia --project=. heat_freshwater_tsubouchi.jl <output folder>")

const ρ₀cₚ  = 1026 * 3991.9                       # J m⁻³ K⁻¹, as in gateways.jl
const S_ours, S_tsu = 34.8, 34.67
const T_tsu  = 1.01
const T_ours = Dict("fram" => 0.0, "bso" => 0.0, "davis" => 0.0, "bering" => -1.9)   # as in gateways.jl
const gateways = ["davis", "fram", "bso", "bering"]

function read_columns(file)
    lines = filter(l -> !isempty(strip(l)) && !startswith(l, "#"), readlines(file))
    header = [strip(name, ['(', ')', '"', ' ']) for name in split(lines[1], ',')]
    rows = [split(l, ',') for l in lines[2:end]]
    return name -> begin
        k = findfirst(==(name), header)
        isnothing(k) && error("$file has no column $name")
        [tryparse(Float64, r[k]) for r in rows]
    end, header
end

mean_of(x) = (y = filter(v -> !isnothing(v) && isfinite(v), x); isempty(y) ? NaN : sum(y) / length(y))

# Run means (Sv, TW, mSv) on Tsubouchi et al.'s references
function converted(column, g)
    Q = mean_of(column(g * "_volume_Sv"))
    H = mean_of(column(g * "_heat_TW"))
    F = mean_of(column(g * "_freshwater_mSv"))
    H_tsu = H - ρ₀cₚ * (T_tsu - T_ours[g]) * Q * 1e6 / 1e12
    F_tsu = 1000Q - (S_ours / S_tsu) * (1000Q - F)
    return (; Q, H = H_tsu, F = F_tsu)
end

model, _ = read_columns(joinpath(folder, "gateways.csv"))
glorys_file = joinpath(folder, "glorys_gateways.csv")
glorys = isfile(glorys_file) ? first(read_columns(glorys_file)) : nothing
reference, _ = read_columns(joinpath(@__DIR__, "references", "tsubouchi2018_table3.csv"))
days = (maximum(filter(!isnothing, model("hour"))) - minimum(filter(!isnothing, model("hour")))) / 24

ref(name) = (; zip(Symbol.(["volume_Sv", "volume_std", "heat_TW", "heat_std", "freshwater_mSv", "freshwater_std"]),
                  [reference(c)[findfirst(==(name), ["davis", "fram", "bso", "bering", "net_ocean"])]
                   for c in ["volume_Sv", "volume_std", "heat_TW", "heat_std", "freshwater_mSv", "freshwater_std"]])...)

@printf("\n%s: mean over %.0f days. Positive into the Arctic. Heat relative to %.2f degC, freshwater to %.2f.\n",
        folder, days, T_tsu, S_tsu)
println("Tsubouchi et al. (2018): annual mean ± std, Sep 2005 - Aug 2006.\n")
println("gateway   volume (Sv)                   heat (TW)                     freshwater (mSv)")
println("          model  GLORYS  Tsubouchi      model  GLORYS  Tsubouchi      model  GLORYS  Tsubouchi")
sums = Dict(:m => [0.0, 0.0, 0.0], :g => [0.0, 0.0, 0.0])
for g in gateways
    m = converted(model, g)
    G = isnothing(glorys) ? (Q = NaN, H = NaN, F = NaN) : converted(glorys, g)
    r = ref(g)
    sums[:m] .+= [m.Q, m.H, m.F]; sums[:g] .+= [G.Q, G.H, G.F]
    @printf("%-7s %6.2f  %6.2f  %+5.1f ± %3.1f   %6.0f  %6.0f  %+4.0f ± %3.0f    %6.0f  %6.0f  %+4.0f ± %3.0f\n",
            g, m.Q, G.Q, r.volume_Sv, r.volume_std, m.H, G.H, r.heat_TW, r.heat_std, m.F, G.F, r.freshwater_mSv, r.freshwater_std)
end
r = ref("net_ocean")
m, G = sums[:m], sums[:g]
@printf("%-7s %6.2f  %6.2f  %+5.2f ± %4.2f %6.0f  %6.0f  %+4.0f ± %3.0f    %6.0f  %6.0f  %+4.0f ± %3.0f\n",
        "sum", m[1], G[1], r.volume_Sv, r.volume_std, m[2], G[2], r.heat_TW, r.heat_std, m[3], G[3], r.freshwater_mSv, r.freshwater_std)
println("""
The "sum" row's Tsubouchi values are their ocean-only net through the whole boundary (which also includes the small
Fury and Hecla Strait). Each gateway's heat and freshwater depend on the reference values; the sum does not, if the
boundary is closed.""")
