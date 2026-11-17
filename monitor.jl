# Hourly monitor of the nested ocean (panarctic.jl, `[output] monitor_hours`), written to monitor.csv in the
# experiment's output folder and read by plotting/volume_budget.jl and plotting/max_speed.jl:
#
# - net transport through the Atlantic (south, j = 1) and Pacific (north, j = Ny + 1) sides, positive into the domain;
# - domain-mean η, model and GLORYS, and the change in η since the start over the shelf (H < 200 m) and the deep ocean
#   (H > 1000 m);
# - a volume budget accumulated every time step: the volume through the open sides (from the barotropic transport the
#   tracers see) and across the surface (the freshwater volume flux added to η), each as a change in domain-mean η, and
#   the residual (the change in mean η minus both), which should be zero; and the net inflow implied by the change in
#   volume over the last interval;
# - the largest |u| and |v| and where the larger one is (grid indices, longitude, latitude, depth).
#
# Also writes depth.csv (seafloor depth, m, 0 on land, one row per j) for the maps of plotting/max_speed.jl.

using Oceananigans.Architectures: on_architecture
using Oceananigans.Grids: λnode, φnode, znode, static_column_depthᶜᶠᵃ
using Oceananigans.Operators: Δxᶜᶠᶜ, Δzᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: bottom_height_field
using NumericalEarth.NestedModels: boundary_location, sample_parent!
using Printf

include("monitor_file.jl")

# Net transport through the Atlantic and Pacific sides, positive into the domain (Sv), from the velocities at the
# boundary faces (the ones that carry T and S)
function open_side_transports(child)
    Nx, Ny, Nz = size(child.grid)
    v = on_architecture(CPU(), child.velocities.v)       # with a CPU copy of the grid in its current z★ state
    Q(j) = sum(v[i, j, k] * Δxᶜᶠᶜ(i, j, k, v.grid) * Δzᶜᶠᶜ(i, j, k, v.grid) for i in 1:Nx, k in 1:Nz)
    return (atlantic = Q(1) / 1e6, pacific = -Q(Ny + 1) / 1e6)
end

function add_monitor!(simulation, child, output_dir, parent_η; schedule = TimeInterval(1hours))
    Nx, Ny, Nz = size(child.grid)
    bottom = Array(interior(bottom_height_field(child.grid)))[:, :, 1]
    wet    = bottom .< 0
    shelf  = wet .& (bottom .> -200)
    deep   = wet .& (bottom .< -1000)

    function area_mean(a, mask, g)
        A = [Azᶜᶜᶜ(i, j, 1, g) for i in 1:Nx, j in 1:Ny]
        return sum(a[mask] .* A[mask]) / sum(A[mask])
    end

    open(depth_file(output_dir), "w") do io
        foreach(j -> println(io, join(-min.(bottom[:, j], 0), ',')), 1:Ny)
    end

    file = monitor_file(output_dir)
    open(io -> println(io, "hour,atlantic_in_Sv,pacific_in_Sv,net_in_Sv,eta_mean_cm,parent_eta_mean_cm,",
                       "shelf_change_cm,deep_change_cm,edge_volume_cm,surface_volume_cm,budget_residual_cm,",
                       "edge_inflow_from_volume_Sv,max_u,max_v,",
                       "fastest_component,fastest_i,fastest_j,fastest_k,fastest_longitude,fastest_latitude,fastest_z"), file, "w")

    # Volume budget, accumulated every time step
    arch      = Oceananigans.Architectures.architecture(child.grid)
    cpu_grid  = on_architecture(CPU(), child.grid)
    wet_area  = sum(Azᶜᶜᶜ(i, j, 1, cpu_grid) for i in 1:Nx, j in 1:Ny if wet[i, j])
    cell_area = on_architecture(arch, [wet[i, j] ? Azᶜᶜᶜ(i, j, 1, cpu_grid) : 0.0 for i in 1:Nx, j in 1:Ny])
    face_width(j) = on_architecture(arch, [static_column_depthᶜᶠᵃ(i, j, cpu_grid) > 0 ? Δxᶜᶠᶜ(i, j, 1, cpu_grid) : 0.0
                                           for i in 1:Nx])
    south_width, north_width = face_width(1), face_width(Ny + 1)
    Ṽ  = child.free_surface.filtered_state.Ṽ   # barotropic transport per unit width that the tracers see
    Jʷ = child.forcing.η.parameters            # freshwater volume flux, m s⁻¹, positive into the ocean
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

    η̄₀       = Ref(NaN)
    η₀       = Ref{Matrix{Float64}}()
    previous = Ref((time = 0.0, η̄ = NaN, surface = 0.0))   # for the net inflow over the last interval, from the volume
    ηᵖ       = Field{Center, Center, Nothing}(child.grid)

    function monitor(sim)
        v = on_architecture(CPU(), child.velocities.v)      # for the grid metrics in the current z★ state
        η = Array(interior(on_architecture(CPU(), child.free_surface.displacement)))[:, :, 1]
        isassigned(η₀) || (η₀[] = copy(η))
        δη = η .- η₀[]
        t  = time(sim)

        atlantic, pacific = open_side_transports(child)
        η̄ = 100 * area_mean(η, wet, v.grid)
        isnan(η̄₀[]) && (η̄₀[] = η̄)
        sample_parent!(ηᵖ, boundary_location(:η), parent_η, t)
        η̄ᵖ = 100 * area_mean(Array(interior(on_architecture(CPU(), ηᵖ)))[:, :, 1], wet, v.grid)

        edge_cm     = 100 * edge_volume[] / wet_area
        surface_cm  = 100 * surface_volume[] / wet_area
        residual    = (η̄ - η̄₀[]) - edge_cm - surface_cm
        Δt          = t - previous[].time
        from_volume = Δt > 0 ? (wet_area * (η̄ - previous[].η̄) / 100 - (surface_volume[] - previous[].surface)) / Δt / 1e6 : NaN
        previous[]  = (time = t, η̄ = η̄, surface = surface_volume[])

        u = on_architecture(CPU(), child.velocities.u)
        umax, iᵘ = findmax(abs, interior(u))
        vmax, iᵛ = findmax(abs, interior(v))
        # Where the faster component peaks: 1 = u (at x faces), 2 = v (at y faces)
        component, (i, j, k), location = umax ≥ vmax ? (1, Tuple(iᵘ), (Face(), Center(), Center())) :
                                                      (2, Tuple(iᵛ), (Center(), Face(), Center()))
        fastest = (component, i, j, k, λnode(i, j, k, v.grid, location...), φnode(i, j, k, v.grid, location...),
                   znode(i, j, k, v.grid, location...))

        open(io -> @printf(io, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%.3f,%.3f,%.1f\n",
                           t / hours, atlantic, pacific, atlantic + pacific, η̄, η̄ᵖ,
                           100 * area_mean(δη, shelf, v.grid), 100 * area_mean(δη, deep, v.grid),
                           edge_cm, surface_cm, residual, from_volume, umax, vmax, fastest...),
             file, "a")
        return nothing
    end

    add_callback!(simulation, monitor, schedule; name = :monitor)
    @info "Hourly monitor written to $file"
    return file
end
