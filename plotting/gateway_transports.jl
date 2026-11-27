# The model's transports through the sections (online, arctic_transports.jld2) against GLORYS's on its own grid
# (glorys_transports.jld2, from analysis/glorys_transports.jl): volume, heat and freshwater (rows) for each section
# (columns), positive into the Arctic. The model's averages are drawn at the middle of their averaging windows (the
# first output, the initial state, is left out); GLORYS's daily means at noon of their day. Writes
# gateway_transports.png to `figure_dir`.
#
#   julia +1.13 --project=. analysis/glorys_transports.jl [config.toml]    (once, for the reference)
#   julia +1.13 --project=. plotting/gateway_transports.jl [config.toml]
#
# Uses the same config as panarctic.jl.

include("../config.jl")
include("../transports.jl")

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

glorys_file = joinpath(output_dir, "glorys_transports.jld2")
glorys = isfile(glorys_file) ? jldopen(f -> (dates = Date.(f["dates"]), transports = Dict(k => f["transports/$k"] for k in keys(f["transports"]))), glorys_file) :
                               (@warn "No $glorys_file: run analysis/glorys_transports.jl first for the GLORYS reference"; nothing)

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
        scatter!(ax, t, glorys.transports[key]; color = :black, marker = :diamond, label = "GLORYS")
    end
    row == 1 && col == 1 && axislegend(ax; position = :lb, framevisible = false)
end
Label(fig[0, :], "$experiment: transports into the Arctic"; fontsize = 18)

file = joinpath(figure_dir, "gateway_transports.png")
save(file, fig)
@info "Wrote $file"
