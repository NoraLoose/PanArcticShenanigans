# The faces of each transport section (see transports.jl) on the bathymetry of the model grid, one panel per section.
# Writes transport_sections.png to `figure_dir`. Use it to check the sections after changing the grid or the waypoints.
#
#   julia +1.13 --project=. plotting/transport_sections.jl [config.toml]
#
# Builds the grid from the config, like panarctic.jl, so run it where panarctic.jl runs.
# Faces: red where positive u or v flows in the positive direction of the section, blue where it flows against it;
# thin where closed at the surface (land). The black line is the great-circle path through the waypoints.

include("../setup.jl")
include("../transports.jl")

using CairoMakie

cpu_grid = on_architecture(CPU(), grid)
sections = gateway_sections(cpu_grid, config)
waypoints = Dict(name => (lons, lats) for (name, lons, lats) in transport_sections(config))

h_bottom = Array(interior(bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN

ccc = (Center(), Center(), Center())
x = [unit_vector(λnode(i, j, 1, cpu_grid, ccc...), φnode(i, j, 1, cpu_grid, ccc...)) for i in 1:Nx, j in 1:Ny]
nearest_cell(λ, φ) = Tuple(argmax(dot3.(x, Ref(unit_vector(λ, φ)))))

# Points along the great circles between the waypoints, as (i, j) of the nearest cell centre
function path_in_grid(lons, lats; points = 50)
    path = Point2f[]
    for n in 1:length(lons)-1
        a, b = unit_vector(lons[n], lats[n]), unit_vector(lons[n+1], lats[n+1])
        for s in range(0, 1, length = points)
            p = (1 - s) .* a .+ s .* b
            λ, φ = atand(p[2], p[1]), asind(p[3] / sqrt(dot3(p, p)))
            push!(path, Point2f(nearest_cell(λ, φ)...))
        end
    end
    return path
end

fig = Figure(size = (1100, 1000))

for (n, (name, sᵘ, sᵛ)) in enumerate(sections)
    # u face (i, j) is the west edge of cell (i, j), v face (i, j) its south edge
    u_faces, v_faces, _ = open_faces(cpu_grid, sᵘ, sᵛ)
    is_open = Set(vcat([(:u, f...) for f in u_faces], [(:v, f...) for f in v_faces]))
    segments = []
    for (kind, s) in ((:u, sᵘ), (:v, sᵛ)), I in findall(!=(0), s)
        i, j = Tuple(I)
        ends = kind == :u ? (Point2f(i - 0.5, j - 0.5), Point2f(i - 0.5, j + 0.5)) :
                            (Point2f(i - 0.5, j - 0.5), Point2f(i + 0.5, j - 0.5))
        push!(segments, (ends, s[i, j] > 0 ? :red : :blue, (kind, i, j) in is_open ? 3 : 1))
    end

    faces = vcat(findall(!=(0), sᵘ), findall(!=(0), sᵛ))
    i₁, i₂ = extrema(I[1] for I in faces)
    j₁, j₂ = extrema(I[2] for I in faces)
    pad = 8

    ax = Axis(fig[fldmod1(n, 2)...], title = name, xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
    hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-3000, 0), nan_color = :gray80)
    lines!(ax, path_in_grid(waypoints[name]...), color = :black, linewidth = 1)
    for (ends, color, width) in segments
        linesegments!(ax, [ends...], color = color, linewidth = width)
    end
    limits!(ax, i₁ - pad, i₂ + pad, j₁ - pad, j₂ + pad)
    n == 1 && Colorbar(fig[1:2, 3], hm, label = "bottom height [m]")
end

save(joinpath(figure_dir, "transport_sections.png"), fig)
@info "Wrote $(joinpath(figure_dir, "transport_sections.png"))"
