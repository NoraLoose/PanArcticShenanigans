# The transport sections (see transports.jl) on the bathymetry of the model grid: all of them on the whole domain
# (transport_sections.png), and the faces of each, zoomed in (transport_section_faces.png), written to `figure_dir`.
# Use them to check the sections after changing the grid or the waypoints.
#
#   julia +1.13 --project=. plotting/transport_sections.jl [config.toml]
#
# Builds the grid from the config, like panarctic.jl, so run it where panarctic.jl runs.
# The open boundaries, whose transports are also computed, are shown on the whole domain. In the zoomed panels, faces are
# red where positive u or v flows in the positive direction of the section (into the Arctic), blue where it flows
# against it, and thin where closed at the surface (land).

include("../setup.jl")
include("../transports.jl")

using CairoMakie

cpu_grid = on_architecture(CPU(), grid)
sections = gateway_sections(cpu_grid, config)
isempty(sections) && error("No transport sections in $config_file: add [transports.sections] (see config.toml).")
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

# The faces of a section as line segments: (ends, colour, line width). u face (i, j) is the west edge of cell (i, j),
# v face (i, j) its south edge.
function face_segments(sᵘ, sᵛ)
    u_faces, v_faces, _ = open_faces(cpu_grid, sᵘ, sᵛ)
    is_open = Set(vcat([(:u, f...) for f in u_faces], [(:v, f...) for f in v_faces]))
    segments = []
    for (kind, s) in ((:u, sᵘ), (:v, sᵛ)), I in findall(!=(0), s)
        i, j = Tuple(I)
        ends = kind == :u ? (Point2f(i - 0.5, j - 0.5), Point2f(i - 0.5, j + 0.5)) :
                            (Point2f(i - 0.5, j - 0.5), Point2f(i + 0.5, j - 0.5))
        push!(segments, (ends, s[i, j] > 0 ? :red : :blue, (kind, i, j) in is_open ? 3 : 1))
    end
    return segments
end

bathymetry!(ax) = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-3000, 0), nan_color = :gray80)

# ## All sections on the whole domain

section_colors = [:red, :darkorange, :magenta, :lime, :cyan, :yellow]   # bright, against the bathymetry colours

fig = Figure(size = (900, 800))
ax  = Axis(fig[1, 1], title = "Transport sections", xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
hm  = bathymetry!(ax)
for (n, (name, sᵘ, sᵛ)) in enumerate(sections)
    linesegments!(ax, [p for (ends, _, _) in face_segments(sᵘ, sᵛ) for p in ends],
                  color = section_colors[mod1(n, end)], linewidth = 4, label = name)
end
for (j, name, style) in ((0.5, "atlantic_boundary", :dash), (Ny + 0.5, "pacific_boundary", :dot))
    lines!(ax, [Point2f(0.5, j), Point2f(Nx + 0.5, j)], color = :black, linestyle = style, linewidth = 3, label = name)
end
Colorbar(fig[1, 2], hm, label = "bottom height [m]")
Legend(fig[2, 1:2], ax, orientation = :horizontal, nbanks = 2, framevisible = false)

save(joinpath(figure_dir, "transport_sections.png"), fig)
@info "Wrote $(joinpath(figure_dir, "transport_sections.png"))"

# ## The faces of each section

fig_faces = Figure(size = (1100, 1000))
for (n, (name, sᵘ, sᵛ)) in enumerate(sections)
    faces = vcat(findall(!=(0), sᵘ), findall(!=(0), sᵛ))
    i₁, i₂ = extrema(I[1] for I in faces)
    j₁, j₂ = extrema(I[2] for I in faces)
    pad = 8

    panel = Axis(fig_faces[fldmod1(n, 2)...], title = name, xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
    map = bathymetry!(panel)
    lines!(panel, path_in_grid(waypoints[name]...), color = :black, linewidth = 1)
    for (ends, color, width) in face_segments(sᵘ, sᵛ)
        linesegments!(panel, [ends...], color = color, linewidth = width)
    end
    limits!(panel, i₁ - pad, i₂ + pad, j₁ - pad, j₂ + pad)
    n == 1 && Colorbar(fig_faces[1:2, 3], map, label = "bottom height [m]")
end
Legend(fig_faces[3, 1:2],
       [LineElement(color = :red, linewidth = 3), LineElement(color = :blue, linewidth = 3),
        LineElement(color = :gray40, linewidth = 1), LineElement(color = :black, linewidth = 1)],
       ["face, positive u or v flows the positive way", "face, positive u or v flows the negative way",
        "thin: closed at the surface (land)", "great-circle path through the waypoints"],
       orientation = :horizontal, nbanks = 2, framevisible = false)

save(joinpath(figure_dir, "transport_section_faces.png"), fig_faces)
@info "Wrote $(joinpath(figure_dir, "transport_section_faces.png"))"
