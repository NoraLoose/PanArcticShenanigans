# The water crossing the gateways, run against GLORYS: inflow and outflow, their transport-weighted temperature and
# salinity, and the heat and freshwater errors split into a volume part (too much or too little water, carrying
# GLORYS's temperature or salinity) and a property part (the same water, warmer or saltier than GLORYS's). Then the
# mean sections of normal velocity, T and S, model beside GLORYS.
#
#   julia --project=.. plot_sections.jl <output_dir> [gateway ...] [depth=500]      (default gateways: bso fram)
#
# Needs gateways_inout.csv and section_<name>.csv (section_properties.jl) and glorys_gateways_inout.csv and
# glorys_section_<name>.csv (glorys_gateways_native.jl). Writes sections_<name>.png.

using CairoMakie, DelimitedFiles, Printf

const ρ₀cₚ = 1026 * 3991.9
const Sʳᵉᶠ = 34.8
reference_temperature(name) = name == "bering" ? -1.9 : 0.0   # as gateway_lines in gateways.jl

folder   = ARGS[1]
options  = filter(a -> occursin('=', a), ARGS[2:end])
gateways = filter(a -> !occursin('=', a), ARGS[2:end])
isempty(gateways) && (gateways = ["bso", "fram"])
maxdepth = any(startswith("depth="), options) ? parse(Float64, split(only(filter(startswith("depth="), options)), '=')[2]) : Inf

function read_csv(path)
    data, header = readdlm(path, ',', Float64; header = true)
    return Dict(String(h) => data[:, c] for (c, h) in enumerate(vec(header)))
end

model  = read_csv(joinpath(folder, "gateways_inout.csv"))
glorys = read_csv(joinpath(folder, "glorys_gateways_inout.csv"))
common = intersect(round.(Int, model["hour"]), round.(Int, glorys["hour"]))
isempty(common) && error("the model and GLORYS files have no day in common")
average(d, column) = sum(d[column][i] for i in eachindex(d["hour"]) if round(Int, d["hour"][i]) in common) / length(common)

#####
##### Inflow and outflow, and where the heat and freshwater errors come from
#####

println("\n$folder: means over the $(length(common)) days with both model and GLORYS. Positive into the Arctic.")
println("Heat relative to 0 °C (Bering −1.9 °C), freshwater to 34.8. θ, S: transport-weighted temperature and salinity.")
println("Error split: volume part = (Q − Q_G) × GLORYS's heat or freshwater per Sv; property part = the rest.\n")
@printf("%-7s %-4s %16s %14s %16s %21s %23s\n", "", "", "Q (Sv)", "θ (°C)", "S", "heat error (TW)", "freshwater error (mSv)")
@printf("%-7s %-4s %16s %14s %16s %21s %23s\n", "", "", "model GLORYS", "model GLORYS", "model GLORYS", "total vol. prop.", "total vol. prop.")
for name in gateways, direction in ("in", "out")
    column(q) = "$(name)_$(direction)_$(q)"
    Qm, Hm, Fm = (average(model,  column(q)) for q in ("volume_Sv", "heat_TW", "freshwater_mSv"))
    Qg, Hg, Fg = (average(glorys, column(q)) for q in ("volume_Sv", "heat_TW", "freshwater_mSv"))
    θ(H, Q) = reference_temperature(name) + H * 1e12 / (ρ₀cₚ * Q * 1e6)
    S(F, Q) = Sʳᵉᶠ * (1 - F * 1e3 / (Q * 1e6))
    heat_volume, freshwater_volume = (Qm - Qg) * Hg / Qg, (Qm - Qg) * Fg / Qg
    @printf("%-7s %-4s %7.2f %8.2f %6.2f %7.2f %7.3f %8.3f %7.1f %6.1f %6.1f %7.1f %7.1f %7.1f\n",
            name, direction, Qm, Qg, θ(Hm, Qm), θ(Hg, Qg), S(Fm, Qm), S(Fg, Qg),
            Hm - Hg, heat_volume, Hm - Hg - heat_volume, Fm - Fg, freshwater_volume, Fm - Fg - freshwater_volume)
end
println()

#####
##### Mean sections
#####

# A rectangle per face and level: columns between the midpoints of neighbouring positions along the line
function rectangles(section)
    s = section["s_km"]
    columns = sort(unique(s))
    gap = length(columns) > 1 ? minimum(diff(columns)) : 1.0
    edges = [columns[1] - gap / 2; (columns[1:end-1] .+ columns[2:end]) ./ 2; columns[end] + gap / 2]
    column = Dict(c => n for (n, c) in enumerate(columns))
    return [Rect2f(edges[column[s[i]]], section["z_m"][i] - section["dz_m"][i] / 2,
                   edges[column[s[i]] + 1] - edges[column[s[i]]], section["dz_m"][i]) for i in eachindex(s)]
end

percentile(x, p) = sort(x)[clamp(round(Int, p * length(x)), 1, length(x))]

for name in gateways
    sections = [read_csv(joinpath(folder, file)) for file in ("section_$name.csv", "glorys_section_$name.csv")]
    shapes = rectangles.(sections)
    deepest = min(maxdepth, maximum(maximum(-sec["z_m"] .+ sec["dz_m"] ./ 2) for sec in sections))

    fig = Figure(size = (1300, 1000))
    variables = (("velocity_ms", "normal velocity (m/s, + into the Arctic)", :balance),
                 ("T_C", "temperature (°C)", :thermal),
                 ("S", "salinity", :haline))
    for (row, (key, label, colormap)) in enumerate(variables)
        values = vcat((sec[key] for sec in sections)...)
        color_limits = key == "velocity_ms" ? (-1, 1) .* percentile(abs.(values), 0.98) :
                                       (percentile(values, 0.02), percentile(values, 0.98))
        for (col, title) in enumerate(("model", "GLORYS"))
            ax = Axis(fig[row, col]; title = "$name: $title, $label", xlabel = "distance along the line (km)",
                      ylabel = "z (m)", limits = (nothing, (-deepest, 0)))
            poly!(ax, shapes[col]; color = sections[col][key], colormap, colorrange = color_limits, strokewidth = 0)
        end
        Colorbar(fig[row, 3]; colormap, limits = color_limits, label)
    end
    path = joinpath(folder, "sections_$name.png")
    save(path, fig)
    @info "Wrote $path"
end
