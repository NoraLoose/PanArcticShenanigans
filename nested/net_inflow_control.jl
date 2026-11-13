# Prototype: hold the model's actual net inflow through the open edges at zero, as the volume constraints of NEMO
# (`ln_vol`) and ROMS (`VolCons`) do, instead of only the net of Flather's exterior transport (`zero_net_inflow`).
#
# NEMO and ROMS remove one uniform velocity from all open faces after the boundary scheme, at every barotropic
# substep. The nested model's substeps are inside Oceananigans, so this prototype works one baroclinic step later: after
# each time step it measures the net barotropic inflow Q through the edges and shifts the Flather exterior on both edges
# by the uniform velocity -Q/A (A = Σ H Δ over the open faces), accumulated over the run. If Flather's η term changes
# slowly compared with a time step, the actual net inflow stays at zero; the exterior offset (logged hourly) shows how
# much it takes.
#
# Needs `zero_net_inflow = true` (for the face areas). Included by panarctic_nested.jl when
# `[boundaries] net_inflow_control = true`.

using Oceananigans.Operators: Δxᶜᶠᶜ

mutable struct NetInflowControl{W, E}
    face_widths :: W   # per open side: the wet width of each boundary face (m), zero where dry, on the device
    exchanger   :: E
    area        :: Float64   # Σ H Δ over the open faces (m²)
    offset      :: Float64   # depth-uniform velocity added to the exterior, positive into the domain (m s⁻¹)
    sample_time :: Float64   # time of the parent sample the offset was last applied to
end

inflow_sign(side) = side ∈ (:west, :south) ? 1 : -1

function NetInflowControl(child, exchanger)
    grid  = on_architecture(CPU(), child.grid)
    Ny    = size(grid, 2)
    areas = exchanger.net_inflow.face_areas
    sides = keys(areas)
    face_row(side) = side == :south ? 1 : Ny + 1

    widths = map(sides) do side
        j = face_row(side)
        area = Array(interior(on_architecture(CPU(), areas[side])))[:, 1, 1]
        width = [area[i] > 0 ? Δxᶜᶠᶜ(i, j, 1, grid) : 0.0 for i in eachindex(area)]
        on_architecture(Oceananigans.Architectures.architecture(child.grid), width)
    end

    A = sum(side -> sum(Base.parent(areas[side])), sides)
    return NetInflowControl(NamedTuple{sides}(widths), exchanger, Float64(A), 0.0, Float64(exchanger.times.current))
end

# Net barotropic inflow through the open south and north faces (m³ s⁻¹), from the transport the tracers see
function net_barotropic_inflow(control, child)
    Ṽ  = child.free_surface.filtered_state.Ṽ
    Ny = size(child.grid, 2)
    Q  = 0.0
    for (side, width) in pairs(control.face_widths)
        j = side == :south ? 1 : Ny + 1
        Q += inflow_sign(side) * sum(view(Ṽ.data, 1:length(width), j, 1) .* width)
    end
    return Q
end

function shift_exterior!(control, δ; current_only = false)
    for (side, depth_mean) in pairs(control.exchanger.depth_means)
        Base.parent(depth_mean.current) .+= inflow_sign(side) * δ
        current_only || (Base.parent(depth_mean.previous) .+= inflow_sign(side) * δ)
    end
    return nothing
end

function control_net_inflow!(sim, control, child)
    # A new parent sample arrives without the offset (zero_net_inflow has zeroed its net): give it the offset
    if control.exchanger.times.current != control.sample_time
        shift_exterior!(control, control.offset; current_only = true)
        control.sample_time = control.exchanger.times.current
    end

    δ = -net_barotropic_inflow(control, child) / control.area
    control.offset += δ
    shift_exterior!(control, δ)
    return nothing
end

function add_net_inflow_control!(simulation, child, exchanger)
    control = NetInflowControl(child, exchanger)
    add_callback!(simulation, sim -> control_net_inflow!(sim, control, child), IterationInterval(1);
                  name = :net_inflow_control)
    log_offset(sim) = @info @sprintf("         net inflow control: exterior offset %+.2e m/s (%+.2f Sv)", control.offset,
                                     control.offset * control.area / 1e6)
    add_callback!(simulation, log_offset, TimeInterval(1hours); name = :net_inflow_control_log)
    @info @sprintf("Net inflow control on: open-face area %.3e m²", control.area)
    return control
end
