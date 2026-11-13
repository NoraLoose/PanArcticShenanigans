# Monitors for the nested pan-Arctic runs: edge transports, domain-mean η against GLORYS's, the shelf and deep η
# changes, and the child's η against the parent's along each open edge. Included by panarctic_nested.jl.

using Oceananigans.Architectures: on_architecture
using Oceananigans.Grids: λnode, φnode, znode
using Oceananigans.Operators: Δxᶜᶜᶜ, Δxᶜᶠᶜ, Δzᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.OutputReaders: Clamp
using Oceananigans.Units: Time
using NumericalEarth.DataWrangling: bracketing_dates, default_download_directory, default_horizontal_padding
using NumericalEarth.NestedModels: band_value, boundary_location, sample_parent!
using Printf

# Net transport through the Atlantic (south, j = 1) and Pacific (north, j = Ny + 1) edges, positive into the domain (Sv).
function edge_transports(child)
    Nx, Ny, Nz = size(child.grid)
    v = on_architecture(CPU(), child.velocities.v)       # with a CPU copy of the grid in its current z★ state
    Q(j) = sum(v[i, j, k] * Δxᶜᶠᶜ(i, j, k, v.grid) * Δzᶜᶠᶜ(i, j, k, v.grid) for i in 1:Nx, k in 1:Nz)
    return (atlantic = Q(1) / 1e6, pacific = -Q(Ny + 1) / 1e6)
end

# The parent dataset's free surface over the whole child domain (2D, small), as a series whose time 0 is
# first(dates): to compare the child's domain-mean η with the parent's. The run itself holds the parent only on
# strips along the open edges.
function parent_free_surface(grid, dataset, dates; dir = default_download_directory(dataset))
    region   = BoundingBox(grid.underlying_grid; padding = default_horizontal_padding(dataset))
    native   = bracketing_dates(dataset, :free_surface, dates)
    metadata = Metadata(:free_surface; dataset, dates = native, region, dir)
    return FieldTimeSeries(metadata, Oceananigans.Architectures.architecture(grid);
                           time_indices_in_memory = length(native), time_indexing = Clamp(), start_time = first(dates))
end

# Hourly: net transport through the Atlantic (south) and Pacific (north) edges, positive into the domain; domain-mean
# η; the change in η since the first call averaged over the shelf (H < 200 m) and the deep ocean (H > 1000 m).
# With the nested model's `exchanger`: the child's η minus the parent's η that Flather sees, averaged along each
# edge (positive drains the domain). With `parent_η` (from parent_free_surface): the parent's domain-mean η on the
# child grid. Also a volume budget accumulated every time step: the volume through the open edges (from the
# barotropic transport the tracers see) and across the surface (the freshwater volume flux added to η), each as a
# change in domain-mean η, and the residual (the change in mean η minus both), which should be zero.
# Written to <output_dir>/nested_monitor.csv.
function add_nested_monitor!(simulation, child, output_dir; schedule = TimeInterval(1hours),
                             exchanger = nothing, parent_η = nothing)
    Nx, Ny, Nz = size(child.grid)
    bottom = Array(interior(bottom_height_field(child.grid)))[:, :, 1]   # in 0.117 the immersed boundary holds a bare array
    wet    = bottom .< 0
    shelf  = wet .& (bottom .> -200)
    deep   = wet .& (bottom .< -1000)

    function area_mean(a, mask, g)
        A = [Azᶜᶜᶜ(i, j, 1, g) for i in 1:Nx, j in 1:Ny]
        return sum(a[mask] .* A[mask]) / sum(A[mask])
    end

    # The seafloor depth (m, 0 on land), one row per j, for maps of the monitor's locations
    open(joinpath(output_dir, "depth.csv"), "w") do io
        foreach(j -> println(io, join(-min.(bottom[:, j], 0), ',')), 1:Ny)
    end

    file = joinpath(output_dir, "nested_monitor.csv")
    open(io -> println(io, "hour,atlantic_in_Sv,pacific_in_Sv,net_in_Sv,eta_mean_cm,shelf_change_cm,deep_change_cm,max_u,max_v,",
                       "atlantic_eta_offset_cm,pacific_eta_offset_cm,parent_eta_mean_cm,",
                       "edge_volume_cm,surface_volume_cm,budget_residual_cm,edge_inflow_from_volume_Sv,",
                       "fastest_component,fastest_i,fastest_j,fastest_k,fastest_longitude,fastest_latitude,fastest_z"), file, "w")

    # Volume budget, accumulated every time step
    arch      = Oceananigans.Architectures.architecture(child.grid)
    cpu_grid  = on_architecture(CPU(), child.grid)
    wet_area  = sum(Azᶜᶜᶜ(i, j, 1, cpu_grid) for i in 1:Nx, j in 1:Ny if wet[i, j])
    cell_area = on_architecture(arch, [wet[i, j] ? Azᶜᶜᶜ(i, j, 1, cpu_grid) : 0.0 for i in 1:Nx, j in 1:Ny])
    face_width(j) = on_architecture(arch, [Oceananigans.Grids.static_column_depthᶜᶠᵃ(i, j, cpu_grid) > 0 ?
                                           Δxᶜᶠᶜ(i, j, 1, cpu_grid) : 0.0 for i in 1:Nx])
    south_width, north_width = face_width(1), face_width(Ny + 1)
    Ṽ  = child.free_surface.filtered_state.Ṽ
    Jʷ = child.forcing.η.parameters      # freshwater volume flux, m s⁻¹, positive into the ocean
    edge_volume    = Ref(0.0)
    surface_volume = Ref(0.0)

    function accumulate_volume(sim)
        Q = sum(view(Ṽ.data, 1:Nx, 1, 1) .* south_width) - sum(view(Ṽ.data, 1:Nx, Ny + 1, 1) .* north_width)
        S = sum(view(Jʷ.data, 1:Nx, 1:Ny, 1) .* cell_area)
        edge_volume[]    += Q * sim.Δt
        surface_volume[] += S * sim.Δt
        return nothing
    end
    add_callback!(simulation, accumulate_volume, IterationInterval(1); name = :volume_budget)
    η̄₀ = Ref(NaN)
    previous = Ref((time = 0.0, η̄ = NaN, surface = 0.0))   # for the edge inflow over the last interval, from the volume


    # Child minus parent η along an edge row j (the first interior cells, where Flather reads η), Δx-weighted
    function edge_offset(η, band, j, t, g)
        isnothing(band) && return NaN
        cells = [i for i in 1:Nx if wet[i, j]]
        Δx = [Δxᶜᶜᶜ(i, j, 1, g) for i in cells]
        δη = [η[i, j] - band_value(band, i, j, 1, Time(t)) for i in cells]
        return sum(δη .* Δx) / sum(Δx)
    end
    side_band(side) = isnothing(exchanger) ? nothing : on_architecture(CPU(), exchanger.bands[side].η)
    ηᵖ = isnothing(parent_η) ? nothing : Field{Center, Center, Nothing}(child.grid)

    η₀   = Ref{Matrix{Float64}}()
    wall = Ref(time_ns())

    function nested_monitor(sim)
        v = on_architecture(CPU(), child.velocities.v)     # for the grid metrics in the current z★ state
        η = Array(interior(on_architecture(CPU(), child.free_surface.displacement)))[:, :, 1]
        isassigned(η₀) || (η₀[] = copy(η))
        δη = η .- η₀[]

        atlantic, pacific = edge_transports(child)
        η̄        = 100 * area_mean(η, wet, v.grid)
        isnan(η̄₀[]) && (η̄₀[] = η̄)
        edge_cm    = 100 * edge_volume[] / wet_area
        surface_cm = 100 * surface_volume[] / wet_area
        residual   = (η̄ - η̄₀[]) - edge_cm - surface_cm
        Δt         = time(sim) - previous[].time
        from_volume = Δt > 0 ? (wet_area * (η̄ - previous[].η̄) / 100 - (surface_volume[] - previous[].surface)) / Δt / 1e6 : NaN
        previous[] = (time = time(sim), η̄ = η̄, surface = surface_volume[])
        shelf_δη = 100 * area_mean(δη, shelf, v.grid)
        deep_δη  = 100 * area_mean(δη, deep, v.grid)
        u = on_architecture(CPU(), child.velocities.u)
        umax, iᵘ = findmax(abs, interior(u))
        vmax, iᵛ = findmax(abs, interior(v))
        # Where the faster component peaks: 1 = u (at x faces), 2 = v (at y faces)
        component, (i, j, k), location = umax ≥ vmax ? (1, Tuple(iᵘ), (Face(), Center(), Center())) :
                                                      (2, Tuple(iᵛ), (Center(), Face(), Center()))
        fastest = (component, i, j, k, λnode(i, j, k, v.grid, location...), φnode(i, j, k, v.grid, location...),
                   znode(i, j, k, v.grid, location...))
        model_hour = time(sim) / hours

        t = time(sim)
        atlantic_offset = 100 * edge_offset(η, side_band(:south), 1, t, v.grid)
        pacific_offset  = 100 * edge_offset(η, side_band(:north), Ny, t, v.grid)
        if isnothing(ηᵖ)
            η̄ᵖ = NaN
        else
            sample_parent!(ηᵖ, boundary_location(:η), parent_η, t)
            η̄ᵖ = 100 * area_mean(Array(interior(on_architecture(CPU(), ηᵖ)))[:, :, 1], wet, v.grid)
        end

        @info @sprintf("         fastest: %s = %.2f m/s at (i, j, k) = (%d, %d, %d), %.1f°E %.1f°N, z = %.0f m",
                       component == 1 ? "u" : "v", max(umax, vmax), fastest[2:4]..., fastest[5:7]...)
        @info @sprintf("%6.1f h  Atlantic in %+6.2f Sv, Pacific in %+6.2f Sv, net %+6.2f Sv; mean η %+7.2f cm (parent %+7.2f); edge η − parent: Atlantic %+5.2f, Pacific %+5.2f cm; Δη shelf %+6.2f, deep %+6.2f cm; volume: edges %+6.2f, surface %+6.2f, residual %+6.2f cm; max|u|, |v| %.2f, %.2f m/s; wall %s",
                       model_hour, atlantic, pacific, atlantic + pacific, η̄, η̄ᵖ, atlantic_offset, pacific_offset,
                       shelf_δη, deep_δη, edge_cm, surface_cm, residual, umax, vmax, prettytime(1e-9 * (time_ns() - wall[])))
        open(io -> @printf(io, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%.3f,%.3f,%.1f\n",
                           model_hour, atlantic, pacific, atlantic + pacific, η̄, shelf_δη, deep_δη, umax, vmax,
                           atlantic_offset, pacific_offset, η̄ᵖ, edge_cm, surface_cm, residual, from_volume, fastest...),
             file, "a")
        wall[] = time_ns()
        return nothing
    end

    add_callback!(simulation, nested_monitor, schedule; name = :nested_monitor)
    return file
end
