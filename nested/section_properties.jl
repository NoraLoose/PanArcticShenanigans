# What the water crossing each gateway carries: inflow and outflow separately, and the mean section of normal velocity,
# temperature and salinity, from the daily 3-D output of a run (nested_ocean_3d.jld2, [output] ocean_3d_days).
#
#   julia --project=. section_properties.jl runs/<run>.toml
#
# Writes in the run's output_dir:
#   gateways_inout.csv     per average: volume (Sv), heat (TW) and freshwater (mSv) carried into the Arctic (inflow) and
#                          out of it (outflow) through each gateway line of gateways.jl, references as in gateways.jl
#   section_<name>.csv     the mean over the run, per face and level: distance along the line, position, depth,
#                          thickness, width, normal velocity (positive into the Arctic), T, S
# glorys_gateways_native.jl writes the same for GLORYS on its own grid (glorys_gateways_inout.csv,
# glorys_section_<name>.csv); plot_sections.jl compares the two.

ENV["PANARCTIC_ARCHITECTURE"] = "CPU"
include("setup.jl")
include("gateways.jl")
include("section_properties_common.jl")

using Oceananigans.Grids: znode

file = joinpath(output_dir, "nested_ocean_3d.jld2")
isfile(file) || error("$file not found: run with [output] ocean_3d_days")
average_days = config["output"]["ocean_3d_days"]

Tₜ, Sₜ, uₜ, vₜ = (FieldTimeSeries(file, name) for name in ("T", "S", "u", "v"))
sections = gateway_sections(Tₜ.grid)
line_start = Dict(name => points[1] for (name, points, _) in gateway_lines)

# Each wet face and level of a section: position, depth, thickness, width, normal velocity, T and S on the face
function section_records(section, u, v, T, S)
    g = v.grid
    records = []
    for (kind, i, j, sgn) in section.faces, k in 1:size(g, 3)
        wet_face(kind, i, j, k, g) || continue
        ℓ = kind === :v ? (Center(), Face(), Center()) : (Face(), Center(), Center())
        width, Δz = kind === :v ? (Δxᶜᶠᶜ(i, j, k, g), Δzᶜᶠᶜ(i, j, k, g)) : (Δyᶠᶜᶜ(i, j, k, g), Δzᶠᶜᶜ(i, j, k, g))
        velocity = sgn * (kind === :v ? v[i, j, k] : u[i, j, k])
        i₀, j₀ = kind === :v ? (i, j - 1) : (i - 1, j)
        wet₀ = !inactive_node(i₀, j₀, k, g, Center(), Center(), Center())
        wet₁ = !inactive_node(i, j, k, g, Center(), Center(), Center())
        Tᶠ = (wet₀ * T[i₀, j₀, k] + wet₁ * T[i, j, k]) / (wet₀ + wet₁)
        Sᶠ = (wet₀ * S[i₀, j₀, k] + wet₁ * S[i, j, k]) / (wet₀ + wet₁)
        push!(records, (; key = (kind, i, j, k), λ = λnode(i, j, k, g, ℓ...), φ = φnode(i, j, k, g, ℓ...),
                          z = znode(i, j, k, g, ℓ...), Δz, width, velocity, T = Tᶠ, S = Sᶠ))
    end
    return records
end

inout_file = joinpath(output_dir, "gateways_inout.csv")
open(io -> println(io, inout_header([s.name for s in sections])), inout_file, "w")
means = Dict(s.name => Dict() for s in sections)

for n in eachindex(Tₜ.times)
    Tₜ.times[n] == 0 && continue                                       # the initial state, not an average
    hour = Tₜ.times[n] / 3600 - 12 * average_days                      # centre of the averaging window
    row = Float64[hour]
    u, v, T, S = uₜ[n], vₜ[n], Tₜ[n], Sₜ[n]
    for s in sections
        records = section_records(s, u, v, T, S)
        append!(row, inflow_outflow(records, s.Tʳᵉᶠ))
        accumulate_section!(means[s.name], records)
    end
    write_gateway_row(inout_file, row)
    @info @sprintf("hour %6.1f: ", hour) *
          join((@sprintf("%s in %+.2f out %+.2f Sv", s.name, row[6g - 4], row[6g - 1]) for (g, s) in enumerate(sections)), ", ")
end

for s in sections
    write_section(joinpath(output_dir, "section_$(s.name).csv"), line_start[s.name], means[s.name])
end
@info "Wrote $inout_file and section_<name>.csv for $(join((s.name for s in sections), ", "))"
