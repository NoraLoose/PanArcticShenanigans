# GLORYS's net transport through the model's two open edges, computed on GLORYS's own grid (its velocities, levels and
# land mask, before any interpolation), next to what the boundary data carry after interpolation onto the model grid
# (the `*_glorys_Sv` columns of a transport monitor). A difference points to the interpolation, the inpainting or the
# bathymetry rather than to GLORYS itself.
#
#   julia +1.13 --project=. plotting/glorys_section_transport.jl [config.toml] [monitor csv]
#
# Reads the GLORYS daily files already downloaded for the run (u and v for each date). Transports in Sv; Atlantic
# positive into the domain, Pacific positive out of it (as in the monitor).

include("../setup.jl")

using Printf
include("../glorys_sections.jl")

monitor_file = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(output_dir, "transport_monitor.csv")
@info @sprintf("Edges between %.1f°N and %.1f°N", extrema(edge_latitudes)...)

function read_monitor(path)
    isfile(path) || return nothing
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [parse.(Float64, split(l, ',')) for l in lines[2:end] if !isempty(l)]
    return Dict(String(h) => [r[k] for r in rows] for (k, h) in enumerate(header))
end
monitor = read_monitor(monitor_file)

# GLORYS daily means are centred at noon; the run starts at noon on dates[start_index], so snapshot n is at monitor
# hour 24 (n - start_index).
function monitor_value(column, n)
    isnothing(monitor) && return NaN
    k = findfirst(h -> h ≈ 24 * (n - start_index), monitor["hour"])
    return isnothing(k) ? NaN : monitor[column][k]
end

@printf("\n%-12s %26s %26s %26s\n", "", "Atlantic (into domain)", "Pacific (out of domain)", "net inflow")
@printf("%-12s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n", "date", "native", "model", "diff",
        "native", "model", "diff", "native", "model", "diff")
atlantic_native, pacific_native = native_edge_transports(dates)
for (n, date) in enumerate(dates)
    A, P = atlantic_native[n] / 1e6, pacific_native[n] / 1e6
    Am, Pm = monitor_value("atlantic_glorys_Sv", n), monitor_value("pacific_glorys_Sv", n)
    @printf("%-12s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", Date(date), A, Am, Am - A,
            P, Pm, Pm - P, A - P, Am - Pm, (Am - Pm) - (A - P))
end
println("\nnative: GLORYS's own grid; model: the boundary data on the model grid (from $monitor_file, NaN if absent)")
