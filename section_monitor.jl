# Interior diagnostics for the open-boundary comparison: transport through the main gateways, every hour, model
# against GLORYS on the model grid; and daily means of the water near the Atlantic edge. Included by panarctic.jl
# after the simulation is built, when [diagnostics] sections = true.
#
# 1. Sections are polylines of (longitude, latitude) points, so they stay at the same place when the resolution changes.
#    Each is mapped to a staircase of model faces through the cell corners nearest the line. Transport is positive to
#    the left of the direction of travel along the polyline (the `positive` entry says which way that is).
#    Hourly, to section_monitor.csv: net, positive and negative parts (summed face by face over every level), and the
#    transport-weighted temperature of the positive part. GLORYS uses its velocities interpolated onto the model's
#    faces (uᶠ, vᶠ) and the model's own wet face areas, so model and GLORYS differ only in how the water moves.
# 2. Daily means of T, S and v on the first `edge_section_rows` rows inside the Atlantic edge
#    (atlantic_edge_daily.jld2), with GLORYS on the same rows (atlantic_edge_glorys.jld2).
# 3. Daily means of the depth-integrated velocities and the free surface (barotropic_daily.jld2), so any other section
#    can be computed afterwards.

using Oceananigans.Operators: Δyᶠᶜᶜ, Δzᶠᶜᶜ
using Oceananigans.OutputWriters: AveragedTimeInterval

diagnostics = get(config, "diagnostics", Dict{String, Any}())

# (name, points, which way is positive). Ends are placed on land; the setup log reports whether they are.
gateway_sections = [
    ("fram",    [(-22.0, 78.83), (13.0, 78.83)],              "north"),   # 78°50′N mooring line, Greenland to Spitsbergen
    ("bso",     [(16.5, 76.6), (19.0, 74.45), (20.6, 70.1)],  "east"),    # Barents Sea Opening: Sørkapp, Bjørnøya, Norway
    ("svinoy",  [(-1.0, 65.0), (6.2, 62.25)],                 "north-east"), # Svinøy section: Norwegian Atlantic Current
    ("lofoten", [(9.25, 69.71), (13.31, 67.88)],              "north-east"), # Norwegian slope current off Lofoten: the old
                                                                             # jet band (v-faces i = 335–352, j = 175 at 15 km)
    ("davis",   [(-63.0, 66.6), (-52.5, 66.6)],               "north"),   # Davis Strait, Baffin Island to Greenland
    ("bering",  [(-170.3, 66.2), (-167.6, 65.55)],            "north"),   # Bering Strait, Cape Dezhnev to Cape Prince of Wales
    # Atlantic inflow to the Nordic Seas, Greenland to Scotland: how it splits between the three passages
    ("denmark", [(-34.5, 67.8), (-23.8, 65.6)],               "north-east"), # Denmark Strait, near the Látrabjarg line
    ("iceland_faroe",  [(-15.0, 64.8), (-6.9, 62.15)],        "north-east"), # Iceland–Faroe Ridge
    ("faroe_shetland", [(-6.9, 62.05), (-1.2, 60.4)],         "north-east"), # Faroe–Shetland Channel
]

# ## Geometry: corners, the staircase path and its faces

corner_λ = [λnode(i, j, 1, section_grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]
corner_φ = [φnode(i, j, 1, section_grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]

great_circle(λ₁, φ₁, λ₂, φ₂) = 2asind(sqrt(sind((φ₂ - φ₁) / 2)^2 + cosd(φ₁) * cosd(φ₂) * sind((λ₂ - λ₁) / 2)^2))

function nearest_corner(λ, φ; near = nothing, reach = 4)
    candidates = isnothing(near) ? CartesianIndices(corner_λ) :
                 CartesianIndices((max(1, near[1] - reach):min(Nx + 1, near[1] + reach),
                                   max(1, near[2] - reach):min(Ny + 1, near[2] + reach)))
    best = argmin(c -> great_circle(λ, φ, corner_λ[c], corner_φ[c]), candidates)
    return (best[1], best[2])
end

# Corners visited by the line through `points` (linear in longitude and latitude between points), joined by unit steps
function corner_path(points)
    spacing = minimum(great_circle(corner_λ[i, j], corner_φ[i, j], corner_λ[i+1, j], corner_φ[i+1, j])
                      for i in 1:Nx, j in (1, Ny ÷ 2, Ny + 1))
    path = [nearest_corner(points[1]...)]
    for p in 1:length(points) - 1
        (λa, φa), (λb, φb) = points[p], points[p + 1]
        n = ceil(Int, 8 * great_circle(λa, φa, λb, φb) / spacing)
        for s in 1:n
            α = s / n
            c = nearest_corner((1 - α) * λa + α * λb, (1 - α) * φa + α * φb; near = path[end])
            while c != path[end]                      # unit steps: i first, then j
                i, j = path[end]
                push!(path, i != c[1] ? (i + sign(c[1] - i), j) : (i, j + sign(c[2] - j)))
            end
        end
    end
    return path
end

# Faces crossed by the path, each once: (:u or :v, i, j, sign), where sign × velocity is the transport to the left.
#   corner (i, j) → (i + 1, j): v-face (i, j), +1      corner (i, j) → (i, j + 1): u-face (i, j), −1
function path_faces(path)
    signs = Dict{Tuple{Symbol, Int, Int}, Int}()
    for (a, b) in zip(path[1:end-1], path[2:end])
        Δi, Δj = b[1] - a[1], b[2] - a[2]
        face = Δj == 0 ? (:v, min(a[1], b[1]), a[2]) : (:u, a[1], min(a[2], b[2]))
        signs[face] = get(signs, face, 0) + (Δj == 0 ? Δi : -Δj)
    end
    return [(f..., s) for (f, s) in signs if s != 0]   # a face crossed back and forth drops out
end

# Wet area of every level of every face (Nfaces × Nz), the sign, and where T is averaged onto the face
function face_geometry(faces)
    A = zeros(length(faces), Nz)
    for (n, (kind, i, j, _)) in enumerate(faces), k in 1:Nz
        A[n, k] = kind === :v ?
            (wetcell(i, j, k, cpu_grid, Center(), Face(), Center()) ? Δxᶜᶠᶜ(i, j, k, cpu_grid) * Δzᶜᶠᶜ(i, j, k, cpu_grid) : 0.0) :
            (wetcell(i, j, k, cpu_grid, Face(), Center(), Center()) ? Δyᶠᶜᶜ(i, j, k, cpu_grid) * Δzᶠᶜᶜ(i, j, k, cpu_grid) : 0.0)
    end
    return A
end

struct GatewaySection
    name      :: String
    positive  :: String
    faces     :: Vector{Tuple{Symbol, Int, Int, Int}}
    area      :: Matrix{Float64}
end

gateways = map(gateway_sections) do (name, points, positive)
    path  = corner_path(points)
    faces = path_faces(path)
    area  = face_geometry(faces)
    wet   = vec(sum(area, dims = 2)) .> 0
    @info @sprintf("Section %-14s %3d faces (%d wet), corners (%d, %d) to (%d, %d), %.2f°E %.2f°N to %.2f°E %.2f°N, wet area %.2f km², positive %s",
                   name, length(faces), count(wet), path[1]..., path[end]...,
                   corner_λ[path[1]...], corner_φ[path[1]...], corner_λ[path[end]...], corner_φ[path[end]...],
                   sum(area) / 1e6, positive)
    # Both ends should be on land; a wet end face means water can pass round the end of the section
    end_faces = [path_faces(path[1:2]); path_faces(path[end-1:end])]
    any(f -> sum(face_geometry([f])) > 0, end_faces) &&
        @warn "Section $name: an end of the line is in water at this resolution; move its end point further inland"
    GatewaySection(name, positive, faces, area)
end

# ## Transport through a section from u, v and T arrays (on the CPU)

face_value(u, v, kind, i, j, k) = kind === :v ? v[i, j, k] : u[mod1(i, size(u, 1)), j, k]
function face_temperature(T, kind, i, j, k)            # mean of the two cells either side, the wet one if only one is
    a = kind === :v ? (i, j - 1) : (i - 1, j)
    T₁ = checkbounds(Bool, T, a..., k) ? T[a..., k] : NaN
    T₂ = checkbounds(Bool, T, i, j, k) ? T[i, j, k] : NaN
    return isnan(T₁) ? T₂ : isnan(T₂) ? T₁ : (T₁ + T₂) / 2
end

function gateway_transport(s::GatewaySection, u, v, T)
    net = positive = negative = heat = 0.0
    for (n, (kind, i, j, sgn)) in enumerate(s.faces), k in 1:Nz
        A = s.area[n, k]
        A == 0 && continue
        q = sgn * face_value(u, v, kind, i, j, k) * A
        net += q
        if q > 0
            positive += q
            heat += q * face_temperature(T, kind, i, j, k)
        else
            negative += q
        end
    end
    return (; net, positive, negative, temperature = positive > 0 ? heat / positive : NaN)
end

cpu_array(f) = replace(Array(interior(f)), NaN => 0.0)

# GLORYS at each daily snapshot: the T of land cells is set to NaN so it is not averaged onto wet faces
land_cells = [!wetcell(i, j, k, cpu_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
function masked_temperature(T)
    T = Array(interior(T))
    T[land_cells] .= NaN
    return T
end

glorys_gateways = [[gateway_transport(s, cpu_array(uᶠ[n]), cpu_array(vᶠ[n]), masked_temperature(Tᵉˣᵗ[n]))
                    for n in eachindex(vᶠ.times)] for s in gateways]

section_file = joinpath(output_dir, "section_monitor.csv")
open(section_file, "w") do io
    columns = ["hour"]
    for s in gateways, q in ("net", "positive", "negative"), source in ("model", "glorys")
        push!(columns, "$(s.name)_$(q)_$(source)_Sv")
    end
    for s in gateways, source in ("model", "glorys")
        push!(columns, "$(s.name)_inflow_T_$(source)_C")
    end
    println(io, join(columns, ","))
end

function section_monitor(sim)
    model = ocean_model_of(sim)
    t = model.clock.time
    u, v = cpu_array(model.velocities.u), cpu_array(model.velocities.v)
    T = masked_temperature(model.tracers.T)
    model_values = [gateway_transport(s, u, v, T) for s in gateways]
    glorys_value(n, q) = in_time([getproperty(g, q) for g in glorys_gateways[n]], t, t₀_obc, Δt_obc)
    values = Float64[(t - start_time) / 3600]
    for n in eachindex(gateways), q in (:net, :positive, :negative)
        push!(values, getproperty(model_values[n], q) / 1e6, glorys_value(n, q) / 1e6)
    end
    for n in eachindex(gateways)
        push!(values, model_values[n].temperature, glorys_value(n, :temperature))
    end
    open(section_file, "a") do io
        println(io, join((@sprintf("%.4f", x) for x in values), ","))
    end
    return nothing
end

# Callbacks also run at iteration 0: the model starts from GLORYS, so the first row should nearly agree
add_run_callback!(section_monitor, TimeInterval(1hours), :panarctic_section_monitor)
@info "Section monitor: $section_file"

# ## Daily means near the Atlantic edge, and GLORYS on the same rows

edge_rows_daily = 1:get(diagnostics, "edge_section_rows", 8)
add_run_writer!(:atlantic_edge, JLD2Writer(ocean.model, (; T = ocean.model.tracers.T, S = ocean.model.tracers.S,
                                                             v = ocean.model.velocities.v);
                                           filename = joinpath(output_dir, "atlantic_edge_daily.jld2"),
                                           indices = (:, edge_rows_daily, :),
                                           schedule = AveragedTimeInterval(1days),
                                           overwrite_files = true))

Oceananigans.OutputWriters.jldopen(joinpath(output_dir, "atlantic_edge_glorys.jld2"), "w") do file
    file["times"] = collect(vᶠ.times)          # seconds from 00:00 on dates[1]; the model's clock counts the same way
    file["start_time"] = start_time
    file["rows"] = collect(edge_rows_daily)
    file["T"] = cat((Array(interior(Tᵉˣᵗ[n]))[:, edge_rows_daily, :] for n in eachindex(Tᵉˣᵗ.times))..., dims = 4)
    file["S"] = cat((Array(interior(Sᵉˣᵗ[n]))[:, edge_rows_daily, :] for n in eachindex(Sᵉˣᵗ.times))..., dims = 4)
    file["v"] = cat((Array(interior(vᶠ[n]))[:, edge_rows_daily, :] for n in eachindex(vᶠ.times))..., dims = 4)
end

# ## Daily means of the depth-integrated velocities and the free surface

U_daily = Field(Integral(ocean.model.velocities.u, dims = 3))
V_daily = Field(Integral(ocean.model.velocities.v, dims = 3))
add_run_writer!(:barotropic, JLD2Writer(ocean.model, (; U = U_daily, V = V_daily,
                                                          η = ocean.model.free_surface.displacement);
                                        filename = joinpath(output_dir, "barotropic_daily.jld2"),
                                        schedule = AveragedTimeInterval(1days),
                                        overwrite_files = true))
@info "Daily means: atlantic_edge_daily.jld2 (rows $(edge_rows_daily)), barotropic_daily.jld2"

# ## Hourly free surface (whole domain), to follow how η changes during the first days (plotting/eta_hourly.jl)

add_run_writer!(:free_surface_hourly, JLD2Writer(ocean.model, (; η = ocean.model.free_surface.displacement);
                                                 filename = joinpath(output_dir, "free_surface_hourly.jld2"),
                                                 schedule = TimeInterval(1hours),
                                                 overwrite_files = true))
@info "Hourly free surface: free_surface_hourly.jld2"
