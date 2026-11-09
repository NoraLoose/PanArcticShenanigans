# Compare runs with different barotropic open boundary conditions from their transport_monitor.csv files.
#
#   julia +1.13 --project=. plotting/obc_comparison.jl flather hybrid total_only clamped
#
# Each argument is an output directory (its name is used as the label). Writes obc_comparison.png in the current
# directory and prints a table. Transports in Sv: Atlantic positive into the domain, Pacific positive out of it.

using CairoMakie, Printf, Statistics

function read_monitor(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(l)]
    return Dict(String(h) => [r[k] for r in rows] for (k, h) in enumerate(header))
end

runs = isempty(ARGS) ? ["flather", "hybrid", "total_only", "clamped"] : ARGS
data = [(name, read_monitor(joinpath(name, "transport_monitor.csv"))) for name in runs]
colours = Dict("flather" => :steelblue, "hybrid" => :black, "total_only" => :darkorange, "clamped" => :firebrick)
colour(name) = get(colours, name, :gray50)

rms(a) = sqrt(mean(abs2, a))
Δfrom_start(a) = 100 .* (a .- a[1])      # cm

@printf("\n%-12s %14s %14s %14s %16s %16s\n", "run", "Atlantic rms", "Pacific rms", "net std", "Δ⟨η⟩ − GLORYS", "budget residual")
@printf("%-12s %14s %14s %14s %16s %16s\n", "", "error [Sv]", "error [Sv]", "[Sv]", "day 5 [cm]", "rms [Sv]")
scores = []
for (name, d) in data
    eA = rms(d["atlantic_model_Sv"] .- d["atlantic_glorys_Sv"])
    eP = rms(d["pacific_model_Sv"] .- d["pacific_glorys_Sv"])
    sN = std(d["net_model_Sv"])
    eη = Δfrom_start(d["eta_mean_m"])[end] - Δfrom_start(d["glorys_eta_mean_m"])[end]
    r  = d["budget_residual_Sv"][2:end]
    @printf("%-12s %14.3f %14.3f %14.3f %16.2f %16.3f\n", name, eA, eP, sN, eη, rms(filter(isfinite, r)))
    push!(scores, (; name, eA, eP, sN, eη))
end

fig = Figure(size = (1400, 900), fontsize = 15)
axA = Axis(fig[1, 1], title = "(a) Atlantic edge: net transport into the domain", xlabel = "hour", ylabel = "Sv")
axP = Axis(fig[1, 2], title = "(b) Pacific edge: net transport out of the domain", xlabel = "hour", ylabel = "Sv")
axη = Axis(fig[2, 1], title = "(c) change in the domain-mean free surface", xlabel = "hour", ylabel = "cm")
axs = Axis(fig[2, 2], title = "(d) errors against GLORYS", ylabel = "Sv",
           xticks = (1:3, ["Atlantic rms error", "Pacific rms error", "net-transport std"]))
for (name, d) in data
    h = d["hour"]
    lines!(axA, h, d["atlantic_model_Sv"], color = colour(name), label = name)
    lines!(axP, h, d["pacific_model_Sv"],  color = colour(name))
    lines!(axη, h, Δfrom_start(d["eta_mean_m"]), color = colour(name))
end
h, d = data[1][2]["hour"], data[1][2]
lines!(axA, h, d["atlantic_glorys_Sv"], color = :gray, linestyle = :dash, label = "GLORYS")
lines!(axP, h, d["pacific_glorys_Sv"],  color = :gray, linestyle = :dash)
lines!(axη, h, Δfrom_start(d["glorys_eta_mean_m"]), color = :gray, linestyle = :dash)
for (k, s) in enumerate(scores)
    barplot!(axs, 1:3, [s.eA, s.eP, s.sN]; dodge = fill(k, 3), n_dodge = length(scores), color = colour(s.name))
end
Legend(fig[3, 1:2], axA, orientation = :horizontal, framevisible = false)
save("obc_comparison.png", fig)
println("\nsaved obc_comparison.png")
