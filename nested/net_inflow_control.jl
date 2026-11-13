# Prototype: hold the model's actual net inflow through the open edges at zero, as the volume constraints of NEMO
# (`ln_vol`) and ROMS (`VolCons`) do, instead of only the net of Flather's exterior transport (`zero_net_inflow`).
#
# The nested model's substeps are inside Oceananigans, so this prototype works one baroclinic step later: after each
# time step it takes the net inflow through the edges from the change in the domain's volume minus the surface
# freshwater flux, and raises the parent's sea level ηᵉˣᵗ on both edges by one uniform amount, δη = -Q / Σ √(gH) Δ,
# accumulated over the run (as the port's hybrid condition did). Flather, Uᵇ = Uᵉˣᵗ ∓ √(gH) (ηᵇ − ηᵉˣᵗ), then lets in
# √(gH) δη per unit length more on every face. Shifting Uᵉˣᵗ instead does not work: Flather turns most of it into a
# change of the edge sea level within a cell's wave-crossing time, and the shift needed reached ±140 Sv.
#
# Hourly in <output_dir>/net_inflow_control.csv: the sea-level offset (cm), the edge inflow from the volume and from the
# face transport (Sv). Needs `zero_net_inflow = true` (for the face areas). Included by panarctic_nested.jl when
# `[boundaries] net_inflow_control = true`.

using Oceananigans.Operators: Δxᶜᶠᶜ, Azᶜᶜᶜ

mutable struct NetInflowControl{W, C, E, J, H}
    face_widths  :: W        # per open side: the wet width of each boundary face (m), zero where dry, on the device
    cell_area    :: C        # horizontal area of each wet cell (m²), zero where dry, on the device
    exchanger    :: E
    freshwater   :: J        # surface freshwater volume flux added to η (m s⁻¹, positive into the ocean)
    displacement :: H        # η
    conductance  :: Float64  # Σ √(gH) Δ over the open faces (m² s⁻¹): inflow per metre of ηᵉˣᵗ
    offset       :: Float64  # added to the parent's sea level on both edges (m)
    sample_time  :: Float64  # time of the parent sample the offset was last applied to
    volume       :: Float64  # Σ η A at the end of the previous step (m³)
    edge_inflow  :: Float64  # edge inflow over the last step, from the volume (m³ s⁻¹)
end

inflow_sign(side) = side ∈ (:west, :south) ? 1 : -1

function NetInflowControl(child, exchanger)
    arch  = Oceananigans.Architectures.architecture(child.grid)
    grid  = on_architecture(CPU(), child.grid)
    Nx, Ny = size(grid)[1:2]
    areas = exchanger.net_inflow.face_areas
    sides = keys(areas)
    face_row(side) = side == :south ? 1 : Ny + 1

    widths = map(sides) do side
        j = face_row(side)
        area = Array(interior(on_architecture(CPU(), areas[side])))[:, 1, 1]
        on_architecture(arch, [area[i] > 0 ? Δxᶜᶠᶜ(i, j, 1, grid) : 0.0 for i in eachindex(area)])
    end

    bottom    = Array(interior(bottom_height_field(child.grid)))[:, :, 1]
    cell_area = on_architecture(arch, [bottom[i, j] < 0 ? Azᶜᶜᶜ(i, j, 1, grid) : 0.0 for i in 1:Nx, j in 1:Ny])

    g = 9.80665
    conductance = sum(sides) do side
        j = face_row(side)
        area = Array(interior(on_architecture(CPU(), areas[side])))[:, 1, 1]
        sum(area[i] > 0 ? sqrt(g * area[i] / Δxᶜᶠᶜ(i, j, 1, grid)) * Δxᶜᶠᶜ(i, j, 1, grid) : 0.0 for i in eachindex(area))
    end
    control = NetInflowControl(NamedTuple{sides}(widths), cell_area, exchanger, child.forcing.η.parameters,
                               child.free_surface.displacement, conductance, 0.0, Float64(exchanger.times.current),
                               0.0, 0.0)
    control.volume = domain_volume(control)
    return control
end

domain_volume(control)  = sum(view(interior(control.displacement), :, :, 1) .* control.cell_area)   # η sits at k = Nz + 1
surface_inflow(control) = sum(view(interior(control.freshwater), :, :, 1) .* control.cell_area)

# Net barotropic inflow at the open south and north faces (m³ s⁻¹), from the transport the tracers see
function face_inflow(control, child)
    Ṽ  = child.free_surface.filtered_state.Ṽ
    Ny = size(child.grid, 2)
    return sum(inflow_sign(side) * sum(view(Ṽ.data, 1:length(width), side == :south ? 1 : Ny + 1, 1) .* width)
               for (side, width) in pairs(control.face_widths))
end

function shift_exterior!(control, δη; current_only = false)
    for side_bands in control.exchanger.bands
        Base.parent(side_bands.η.current) .+= δη
        current_only || (Base.parent(side_bands.η.previous) .+= δη)
    end
    return nothing
end

function control_net_inflow!(sim, control)
    # A new parent sample arrives without the offset: give it the offset
    if control.exchanger.times.current != control.sample_time
        shift_exterior!(control, control.offset; current_only = true)
        control.sample_time = control.exchanger.times.current
    end

    volume = domain_volume(control)
    control.edge_inflow = (volume - control.volume) / sim.Δt - surface_inflow(control)
    control.volume = volume

    δη = -control.edge_inflow / control.conductance
    control.offset += δη
    shift_exterior!(control, δη)
    return nothing
end

function add_net_inflow_control!(simulation, child, exchanger, output_dir)
    control = NetInflowControl(child, exchanger)
    add_callback!(simulation, sim -> control_net_inflow!(sim, control), IterationInterval(1); name = :net_inflow_control)

    file = joinpath(output_dir, "net_inflow_control.csv")
    open(io -> println(io, "hour,eta_offset_cm,edge_inflow_from_volume_Sv,edge_inflow_from_faces_Sv"), file, "w")
    function log_control(sim)
        offset, from_volume, from_faces = 100 * control.offset, control.edge_inflow / 1e6, face_inflow(control, child) / 1e6
        @info @sprintf("         net inflow control: sea-level offset %+.2f cm; edge inflow from volume %+.2f, from faces %+.2f Sv",
                       offset, from_volume, from_faces)
        open(io -> @printf(io, "%.2f,%.4f,%.4f,%.4f\n", time(sim) / hours, offset, from_volume, from_faces), file, "a")
        return nothing
    end
    add_callback!(simulation, log_control, TimeInterval(1hours); name = :net_inflow_control_log)
    @info @sprintf("Net inflow control on: %.2f Sv per cm of edge sea level", control.conductance / 100 / 1e6)
    return control
end
