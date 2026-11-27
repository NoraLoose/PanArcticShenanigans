# The model's transports through the sections (online, arctic_transports.jld2) against GLORYS's on its own grid
# (<validation_dir>/GLORYS/transports/, from analysis/validation_diagnostics.jl): volume, heat and freshwater (rows) for each section
# (columns), positive into the Arctic. The model's averages are drawn at the middle of their averaging windows (the
# first output, the initial state, is left out); GLORYS's daily means at noon of their day. Writes
# gateway_transports.png to `figure_dir`.
#
#   julia +1.13 --project=. analysis/validation_diagnostics.jl GLORYS --run config.toml   (once, for the reference)
#   julia +1.13 --project=. plotting/gateway_transports.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../config.jl")
include("../transports.jl")
include("../validation_config.jl")   # read_glorys_transports

using CairoMakie, Dates, JLD2

start_date = config["time"]["start_date"]
section_names = first.(transport_sections(config))
quantities = ("volume", "heat", "freshwater")

# The model's averages: times (days since the start, at the middle of each window) and values, by transport name
model = jldopen(transport_output_file) do f
    iterations = sort(filter(!=("serialized"), keys(f["timeseries/t"])), by = s -> parse(Int, s))[2:end]
    t = [f["timeseries/t/$i"] for i in iterations] ./ 86400
    window = config["output"]["average_hours"] / 24
    Dict("t" => t .- window / 2,
         (key => [only(f["timeseries/$key/$i"]) for i in iterations] for key in keys(f["timeseries"]) if key != "t")...)
end

# Runs from before a section was added to sections.toml don't have it
missing_sections = filter(name -> !haskey(model, "$(name)_volume"), section_names)
isempty(missing_sections) || @info "Not in $transport_output_file (an earlier run): $(join(missing_sections, ", "))"
section_names = setdiff(section_names, missing_sections)

run_days = Date(config["time"]["start_date"]):Day(1):Date(DateTime(config["time"]["end_date"]) - Millisecond(1))
glorys   = read_glorys_transports(config, run_days)
isnothing(glorys) && @warn "No GLORYS transports for the run's days: analysis/validation_diagnostics.jl GLORYS --run $config_file"

fig = Figure(size = (300 * length(section_names) + 100, 750))
for (row, quantity) in enumerate(quantities), (col, name) in enumerate(section_names)
    key = "$(name)_$quantity"
    ax = Axis(fig[row, col]; xlabel = row == length(quantities) ? "days since $start_date" : "",
              ylabel = col == 1 ? "$quantity ($(transport_units[quantity]))" : "",
              title = row == 1 ? replace(name, "_" => " ") : "")
    hlines!(ax, 0; color = :black, linewidth = 0.5)
    lines!(ax, model["t"], model[key]; color = :red, label = "model")
    if !isnothing(glorys) && haskey(glorys.transports, key)
        t = Dates.value.(glorys.dates .- start_date) .+ 0.5
        scatterlines!(ax, t, glorys.transports[key]; color = :black, markersize = 5, label = "GLORYS")
    end
    row == 1 && col == 1 && axislegend(ax; position = :lb, framevisible = false)
end
Label(fig[0, :], "$experiment: transports into the Arctic"; fontsize = 18)

file = joinpath(figure_dir, "gateway_transports.png")
save(file, fig)
@info "Wrote $file"
