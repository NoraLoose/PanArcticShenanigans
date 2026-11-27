# Gateway transports (Davis, Fram, the BSO, Bering) of the model and of GLORYS, as monthly means, and against
# Tsubouchi et al. (2018, JPO 48, 2029-2055), the Sep 2005 - Aug 2006 mooring-based inverse estimate.
#
#   julia +1.13 --project=. analysis/validation_diagnostics.jl GLORYS --run config.toml   (once, for the GLORYS reference)
#   julia +1.13 --project=. analysis/tsubouchi_comparison.jl [config.toml]
#
# Reads the model's online transports (arctic_transports.jld2) from `output_dir` and, if present, GLORYS's from
# <validation_dir>/GLORYS/transports/. Writes two tables, printed and saved as tsubouchi_comparison.md in `output_dir` (markdown, ready to
# paste), and tsubouchi_comparison.png in `figure_dir`:
#
#   1. Calendar-month means of volume, heat and freshwater, model against GLORYS, on our reference values
#      ([transports] reference_temperature and reference_salinity).
#   2. Means over Tsubouchi et al.'s monthly blocks, model and GLORYS against Tsubouchi et al., on THEIR reference
#      values: heat and freshwater (mSv) relative to their monthly boundary-mean θ and S, block by block.
#
# Heat and freshwater per gateway depend on the reference values (the net over a closed boundary does not), so they are
# converted. With the volume Q, our heat H (relative to T₀) and our freshwater F (relative to S₀):
#     heat relative to T′       = H - ρ cₚ (T′ - T₀) Q
#     freshwater relative to S′ = Q - (S₀ / S′) (Q - F)          (since F = Q - Σ vS / S₀)
# T′ and S′ are Tsubouchi et al.'s monthly boundary-mean potential temperature and salinity (Section 2d), one per block.
#
# Tsubouchi et al.'s monthly values (validation/references/tsubouchi2018_monthly_*.csv) are read from their Figs. 5a, 7a
# and 9a: about ±0.2 Sv, ±5 TW, ±5 mSv, drawn as bands. Their "months" are 30-day blocks from 5 September 2005 (the
# tenth is 35 days); a block only partly covered by the run is compared over the covered days, which the table lists.
# Everything is positive into the Arctic.

include("../config.jl")
include("../transports.jl")
include("../validation_config.jl")   # read_glorys_transports

using CairoMakie, Dates, JLD2, Printf

start_date = Date(config["time"]["start_date"])
end_date   = Date(config["time"]["end_date"])
run_days   = Dates.value(end_date - start_date)

# (column in Tsubouchi et al.'s files, section name in sections.toml, label)
gateways = filter(g -> g[2] in first.(transport_sections(config)),
                  [("davis", "davis_strait", "Davis"), ("fram", "fram_strait", "Fram"),
                   ("bso", "barents_sea_opening", "BSO"), ("bering", "bering_strait", "Bering")])
isempty(gateways) && error("None of davis_strait, fram_strait, barents_sea_opening, bering_strait among the transport sections")

#####
##### Tsubouchi et al. (2018)
#####

# Tsubouchi et al.'s monthly boundary-mean potential temperature (°C) and salinity, Sep 2005 - Aug 2006 (their
# Section 2d), one per 30-day block. (The paper prints the 6th salinity as 36.67, a typo for 34.67.)
const T_tsubouchi = [1.26, 1.15, 1.01, 0.91, 0.94, 0.91, 0.83, 0.82, 0.88, 0.99, 1.15, 1.24]
const S_tsubouchi = [34.67, 34.64, 34.63, 34.63, 34.66, 34.67, 34.68, 34.68, 34.70, 34.68, 34.67, 34.66]
const tsubouchi_uncertainty = Dict("volume" => 0.2, "heat" => 5.0, "freshwater" => 5.0)   # reading the figures
const tsubouchi_units       = Dict("volume" => "Sv", "heat" => "TW", "freshwater" => "mSv")

function read_tsubouchi(quantity)
    file  = joinpath(@__DIR__, "..", "validation", "references", "tsubouchi2018_monthly_$quantity.csv")
    lines = filter(l -> !startswith(l, "#") && !isempty(strip(l)), readlines(file))
    header = strip.(split(lines[1], ","))
    rows   = [strip.(split(l, ",")) for l in lines[2:end]]
    values = Dict(String(header[c]) => [parse(Float64, r[c]) for r in rows] for c in 2:length(header))
    return String.(first.(rows)), values
end

block_labels, _ = read_tsubouchi("volume")
tsubouchi = Dict(q => last(read_tsubouchi(q)) for q in ("volume", "heat", "freshwater"))

# Blocks as [first day, last day + 1) in days since the run's start
block_lengths = [fill(30, 9); 35; 30; 30]
block_first   = Dates.value(Date(2005, 9, 5) - start_date) .+ [0; cumsum(block_lengths)[1:end-1]]
blocks = [(label = block_labels[b], a = block_first[b], b = block_first[b] + block_lengths[b], index = b)
          for b in eachindex(block_lengths) if block_first[b] < run_days && block_first[b] + block_lengths[b] > 0]

#####
##### Model and GLORYS
#####

# The model's averages, at the middle of their windows, in days since the start (the first output, the initial state,
# is left out); and the reference values and ρ cₚ the transports were computed with
window = config["output"]["average_hours"] / 24
(model_t, model), reference = jldopen(transport_output_file) do f
    iterations = sort(filter(!=("serialized"), keys(f["timeseries/t"])), by = s -> parse(Int, s))[2:end]
    t = [f["timeseries/t/$i"] for i in iterations] ./ 86400 .- window / 2
    values = Dict(key => [only(f["timeseries/$key/$i"]) for i in iterations] for key in keys(f["timeseries"]) if key != "t")
    reference = (T₀ = f["metadata/reference_temperature"], S₀ = f["metadata/reference_salinity"],
                 ρcₚ = f["metadata/reference_density"] * f["metadata/heat_capacity"])
    return (t, values), reference
end

# GLORYS's daily means at noon of their day, computed with the config's reference values and the same ρ cₚ
glorys_days = start_date:Day(1):end_date - Day(1)
glorys_read = read_glorys_transports(config, glorys_days)
if !isnothing(glorys_read)
    glorys_t, glorys = Dates.value.(glorys_read.dates .- start_date) .+ 0.5, glorys_read.transports
    config_reference = reference_values(config)
    (config_reference.reference_temperature, config_reference.reference_salinity) == (reference.T₀, reference.S₀) ||
        @warn "The config's reference values differ from the model run's; GLORYS's heat and freshwater may not compare"
else
    @warn "No GLORYS transports for the run's days: analysis/validation_diagnostics.jl GLORYS --run $config_file"
    glorys_t, glorys = Float64[], nothing
end

# Volume (Sv), heat (TW) and freshwater (Sv) on our reference values; heat (TW) and freshwater (mSv) on Tsubouchi's
ours(values, name, quantity) = values["$(name)_$quantity"]

function theirs(values, t, name, quantity)
    Q = values["$(name)_volume"]
    quantity == "volume" && return Q
    b = clamp.(searchsortedlast.(Ref(block_first), t), 1, length(block_first))   # Tsubouchi et al.'s block at each time
    quantity == "heat" && return values["$(name)_heat"] .- reference.ρcₚ .* (T_tsubouchi[b] .- reference.T₀) .* Q .* 1e6 ./ 1e12
    F = values["$(name)_freshwater"]
    return 1000 .* (Q .- (reference.S₀ ./ S_tsubouchi[b]) .* (Q .- F))
end

function window_mean(t, x, a, b)
    inside = findall(τ -> a ≤ τ < b, t)
    return isempty(inside) ? (NaN, 0) : (sum(x[inside]) / length(inside), length(inside))
end

fmt(x, digits) = isnan(x) ? "–" : (digits == 0 ? @sprintf("%.0f", x) : digits == 1 ? @sprintf("%.1f", x) :
                                                  digits == 2 ? @sprintf("%.2f", x) : @sprintf("%.3f", x))

#####
##### Tables
#####

io = IOBuffer()
println(io, "# Gateway transports: $experiment ($start_date to $end_date)\n")
println(io, "Positive into the Arctic. Model: online transports, ", @sprintf("%g", window * 24), " h averages. ",
            isnothing(glorys) ? "No GLORYS reference." : "GLORYS: daily means on its own grid, same sections.", "\n")

# 1. Calendar months, our reference values
months = Date[]
let m = firstdayofmonth(start_date)
    while m < end_date
        push!(months, m)
        m += Month(1)
    end
end
monthname(m) = Dates.format(m, "u yyyy")

println(io, "## Calendar-month means, model against GLORYS\n")
println(io, @sprintf("Heat relative to %g °C, freshwater relative to %g. ", reference.T₀, reference.S₀),
            "Days covered by the run in brackets when a month is not complete.\n")
units = Dict("volume" => "Sv", "heat" => "TW", "freshwater" => "Sv")
for (quantity, digits) in (("volume", 2), ("heat", 0), ("freshwater", 3))
    println(io, "**", titlecase(quantity), " ($(units[quantity]))**\n")
    head = ["Gateway"]
    for m in months
        a = max(Dates.value(m - start_date), 0)
        b = min(Dates.value(m + Month(1) - start_date), run_days)
        days = b - a
        tag = days == daysinmonth(m) ? monthname(m) : "$(monthname(m)) [$days d]"
        append!(head, ["$tag model", "GLORYS", "model − GLORYS"])
    end
    println(io, "| ", join(head, " | "), " |")
    println(io, "|", join(fill("---", length(head)), "|"), "|")
    for (_, name, label) in gateways
        row = [label]
        for m in months
            a = max(Dates.value(m - start_date), 0)
            b = min(Dates.value(m + Month(1) - start_date), run_days)
            mm, _ = window_mean(model_t, ours(model, name, quantity), a, b)
            gm = isnothing(glorys) ? NaN : first(window_mean(glorys_t, ours(glorys, name, quantity), a, b))
            append!(row, [fmt(mm, digits), fmt(gm, digits), fmt(mm - gm, digits)])
        end
        println(io, "| ", join(row, " | "), " |")
    end
    println(io)
end

# 2. Tsubouchi et al.'s blocks, their reference values
println(io, "## Against Tsubouchi et al. (2018), over their 30-day blocks\n")
reference_list = join(["$(block.label) $(T_tsubouchi[block.index]) °C / $(S_tsubouchi[block.index])" for block in blocks], ", ")
println(io, "Heat and freshwater (mSv) relative to Tsubouchi et al.'s monthly boundary-mean θ / S ($reference_list). ",
            "Tsubouchi et al.'s ",
            "values read from their figures: ±0.2 Sv, ±5 TW, ±5 mSv. Blocks start on 5 September 2005; ",
            "\"days\" is how many days of the block the run covers.\n")
for (quantity, digits) in (("volume", 2), ("heat", 0), ("freshwater", 0))
    println(io, "**", titlecase(quantity), " ($(tsubouchi_units[quantity]))**\n")
    println(io, "| Gateway | Block | Days | Model | GLORYS | Tsubouchi |")
    println(io, "|---|---|---|---|---|---|")
    for (code, name, label) in gateways, block in blocks
        a, b = max(block.a, 0), min(block.b, run_days)
        mm, _ = window_mean(model_t, theirs(model, model_t, name, quantity), a, b)
        gm = isnothing(glorys) ? NaN : first(window_mean(glorys_t, theirs(glorys, glorys_t, name, quantity), a, b))
        ts = tsubouchi[quantity][code][block.index]
        println(io, "| $label | $(block.label) | $(b - a) of $(block.b - block.a) | ", fmt(mm, digits), " | ",
                    fmt(gm, digits), " | ", fmt(ts, digits), " |")
    end
    println(io)
end

text = String(take!(io))
print(text)
table_file = joinpath(output_dir, "tsubouchi_comparison.md")
write(table_file, text)
@info "Wrote $table_file"

#####
##### Figure: daily model and GLORYS on Tsubouchi et al.'s reference values, with their blocks
#####

quantities = ("volume", "heat", "freshwater")
fig = Figure(size = (300 * length(gateways) + 100, 750))
for (row, quantity) in enumerate(quantities), (col, (code, name, label)) in enumerate(gateways)
    ax = Axis(fig[row, col]; xlabel = row == length(quantities) ? "days since $start_date" : "",
              ylabel = col == 1 ? "$quantity ($(tsubouchi_units[quantity]))" : "",
              title = row == 1 ? label : "")
    hlines!(ax, 0; color = :black, linewidth = 0.5)
    for block in blocks
        a, b = max(block.a, 0), min(block.b, run_days)
        value, δ = tsubouchi[quantity][code][block.index], tsubouchi_uncertainty[quantity]
        band!(ax, [a, b], [value - δ, value - δ], [value + δ, value + δ]; color = (:dodgerblue, 0.25))
        lines!(ax, [a, b], [value, value]; color = :dodgerblue, linewidth = 3, label = "Tsubouchi et al. (2018)")
        mm, _ = window_mean(model_t, theirs(model, model_t, name, quantity), a, b)
        lines!(ax, [a, b], [mm, mm]; color = :red, linewidth = 3, linestyle = :dash, label = "model, block mean")
    end
    lines!(ax, model_t, theirs(model, model_t, name, quantity); color = :red, label = "model")
    isnothing(glorys) || scatterlines!(ax, glorys_t, theirs(glorys, glorys_t, name, quantity); color = :black,
                                       markersize = 5, label = "GLORYS")
    row == 1 && col == 1 && axislegend(ax; position = :lb, framevisible = false, unique = true, labelsize = 10)
end
Label(fig[0, :], "$experiment: transports into the Arctic against Tsubouchi et al. (2018); " *
                 "heat and freshwater relative to their monthly boundary-mean θ and S"; fontsize = 16)

file = joinpath(figure_dir, "tsubouchi_comparison.png")
save(file, fig)
@info "Wrote $file"
