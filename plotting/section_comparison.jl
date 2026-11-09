# Compare runs at the gateways from their section_monitor.csv files (written when [diagnostics] sections = true).
#
#   julia +1.13 --project=. plotting/section_comparison.jl flather hybrid total_only clamped
#
# Each argument is an output directory (its name is used as the label). Writes section_comparison.png in the current
# directory and prints a table: for each gateway, the rms error of the net transport against GLORYS (interpolated onto
# the model's faces, so the geometry is the same) and of the transport-weighted temperature of the inflow, over the run
# and over its last 5 days. Positive directions: Fram, Davis and Bering north; Barents Sea Opening east; Denmark Strait,
# Iceland–Faroe, Faroe–Shetland, Svinøy and Lofoten north-east.

using CairoMakie, Printf, Statistics

function read_monitor(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(l)]
    return Dict(String(h) => [r[k] for r in rows] for (k, h) in enumerate(header))
end

runs = isempty(ARGS) ? ["flather", "hybrid", "total_only", "clamped"] : ARGS
data = [(basename(rstrip(name, '/')), read_monitor(joinpath(name, "section_monitor.csv"))) for name in runs]
colours = Dict("flather" => :steelblue, "hybrid" => :black, "total_only" => :darkorange, "clamped" => :firebrick)
function colour(name)                 # by the scheme in the run's name, e.g. "30day-hybrid"
    for (scheme, c) in colours
        occursin(scheme, name) && return c
    end
    return :gray50
end

gateways = unique(first(split(k, "_net_")) for k in keys(data[1][2]) if occursin("_net_model_Sv", k))
order = ["denmark", "iceland_faroe", "faroe_shetland", "svinoy", "lofoten", "bso", "fram", "davis", "bering"]
gateways = sort(gateways, by = g -> something(findfirst(==(g), order), 99))
titles = Dict("fram" => "Fram Strait (north)", "bso" => "Barents Sea Opening (east)",
              "svinoy" => "Svinøy section (north-east)", "lofoten" => "Lofoten slope (north-east)",
              "davis" => "Davis Strait (north)", "bering" => "Bering Strait (north)",
              "denmark" => "Denmark Strait (north-east)", "iceland_faroe" => "Iceland–Faroe (north-east)",
              "faroe_shetland" => "Faroe–Shetland (north-east)")

rms(a) = sqrt(mean(abs2, filter(isfinite, a)))
last_days(h, n = 5) = h .≥ h[end] - 24n

println()
@printf("%-10s", "")
for g in gateways
    @printf(" %25s", first(get(titles, g, g), 25))
end
@printf("\n%-10s", "run")
for _ in gateways
    @printf(" %25s", "net rms Sv (last 5 d) | ΔT")
end
println()
for (name, d) in data
    @printf("%-10s", name)
    late = last_days(d["hour"])
    for g in gateways
        e = d["$(g)_net_model_Sv"] .- d["$(g)_net_glorys_Sv"]
        eT = d["$(g)_inflow_T_model_C"] .- d["$(g)_inflow_T_glorys_C"]
        @printf(" %8.2f (%5.2f) | %+6.2f", rms(e), rms(e[late]), mean(filter(isfinite, eT[late])))
    end
    println()
end
println("(T: mean over the last 5 days of the model minus GLORYS transport-weighted temperature of the inflow)")

n = length(gateways)
fig = Figure(size = (1500, 380 * n), fontsize = 14)
axes = []
for (r, g) in enumerate(gateways)
    axQ = Axis(fig[r, 1], title = "($('a' + 2r - 2)) $(get(titles, g, g)): net transport", xlabel = "day", ylabel = "Sv")
    axT = Axis(fig[r, 2], title = "($('a' + 2r - 1)) $(get(titles, g, g)): temperature of the inflow",
               xlabel = "day", ylabel = "°C")
    for (name, d) in data
        day = d["hour"] ./ 24
        lines!(axQ, day, d["$(g)_net_model_Sv"], color = colour(name), label = name)
        lines!(axT, day, d["$(g)_inflow_T_model_C"], color = colour(name))
    end
    d = data[1][2]
    day = d["hour"] ./ 24
    lines!(axQ, day, d["$(g)_net_glorys_Sv"], color = :gray, linestyle = :dash, linewidth = 2, label = "GLORYS")
    lines!(axT, day, d["$(g)_inflow_T_glorys_C"], color = :gray, linestyle = :dash, linewidth = 2)
    push!(axes, axQ)
end
Legend(fig[n + 1, 1:2], axes[1], orientation = :horizontal, framevisible = false)
save("section_comparison.png", fig)
println("\nsaved section_comparison.png")
