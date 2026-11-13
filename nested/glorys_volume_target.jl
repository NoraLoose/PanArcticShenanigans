# Make the domain-mean η follow GLORYS's: before every time step, set the net inflow through the open edges (pinned by
# `net_inflow`, a `Ref`) to
#
#     T = A dη̄ᴳ/dt − F + A [(η̄ᴳ − η̄ᴳ₀) − (η̄ − η̄₀)] / τ    (the first term only with `tendency = true`),
#
# with A the wet area, η̄ᴳ GLORYS's domain-mean η (daily means, interpolated in time), F the surface freshwater volume
# flux into the ocean, η̄ the model's domain-mean η and the subscript 0 the start of the run. The first term follows
# GLORYS's volume change, the second hands the surface freshwater back through the edges, and the third relaxes any
# remaining drift over τ. GLORYS's daily means make dη̄ᴳ/dt jump once a day; without it, the relaxation alone follows
# GLORYS, smoothly but about τ behind. Included by panarctic_nested.jl when `[boundaries] follow_glorys_volume = true`.

using Oceananigans.Operators: Azᶜᶜᶜ
using NumericalEarth.NestedModels: boundary_location, sample_parent!

struct GLORYSVolumeTarget{R, C, J, H, M}
    net_inflow   :: R        # Ref read by the barotropic pin (m³ s⁻¹, positive into the domain)
    cell_area    :: C        # horizontal area of each wet cell (m²), zero where dry, on the device
    freshwater   :: J        # surface freshwater volume flux added to η (m s⁻¹, positive into the ocean)
    displacement :: H        # η
    times        :: Vector{Float64}   # times of GLORYS's daily means (s since the start of the run)
    means        :: Vector{Float64}   # GLORYS's domain-mean η at those times (m)
    area         :: Float64
    relaxation_time :: Float64
    tendency        :: Bool
    initial_means   :: M     # (model, GLORYS) domain-mean η at the start (m)
end

function GLORYSVolumeTarget(net_inflow, child, parent_η, relaxation_time, tendency)
    arch   = Oceananigans.Architectures.architecture(child.grid)
    grid   = on_architecture(CPU(), child.grid)
    Nx, Ny = size(grid)[1:2]
    bottom = Array(interior(bottom_height_field(child.grid)))[:, :, 1]
    area   = [bottom[i, j] < 0 ? Azᶜᶜᶜ(i, j, 1, grid) : 0.0 for i in 1:Nx, j in 1:Ny]
    A      = sum(area)

    ηᵖ = Field{Center, Center, Nothing}(child.grid)
    times = collect(Float64, parent_η.times)
    means = map(times) do t
        sample_parent!(ηᵖ, boundary_location(:η), parent_η, t)
        sum(Array(interior(ηᵖ))[:, :, 1] .* area) / A
    end

    target = GLORYSVolumeTarget(net_inflow, on_architecture(arch, area), child.forcing.η.parameters,
                                child.free_surface.displacement, times, means, A, Float64(relaxation_time), tendency,
                                Ref((0.0, 0.0)))
    target.initial_means[] = (model_mean(target), glorys_mean(target, child.clock.time)[1])
    return target
end

model_mean(target) = sum(view(interior(target.displacement), :, :, 1) .* target.cell_area) / target.area

# GLORYS's domain-mean η and its rate of change at time t, linear between the daily means and held before the first
function glorys_mean(target, t)
    times, means = target.times, target.means
    t ≤ times[1] && return means[1], 0.0
    n = min(searchsortedlast(times, t), length(times) - 1)
    slope = (means[n + 1] - means[n]) / (times[n + 1] - times[n])
    return means[n] + slope * (t - times[n]), slope
end

function set_net_inflow!(sim, target)
    η̄ᴳ, ∂tη̄ᴳ = glorys_mean(target, time(sim))
    η̄₀, η̄ᴳ₀  = target.initial_means[]
    F = sum(view(interior(target.freshwater), :, :, 1) .* target.cell_area)
    A = target.area
    target.net_inflow[] = A * ∂tη̄ᴳ * target.tendency - F + A * ((η̄ᴳ - η̄ᴳ₀) - (model_mean(target) - η̄₀)) / target.relaxation_time
    return nothing
end

function add_glorys_volume_target!(simulation, net_inflow, child, parent_η, relaxation_time; tendency = true)
    target = GLORYSVolumeTarget(net_inflow, child, parent_η, relaxation_time, tendency)
    add_callback!(simulation, sim -> set_net_inflow!(sim, target), IterationInterval(1); name = :glorys_volume_target)
    @info @sprintf("Net inflow follows GLORYS's domain-mean η (relaxation time %.1f days, tendency %s)",
                   relaxation_time / Oceananigans.Units.day, tendency)
    return target
end
