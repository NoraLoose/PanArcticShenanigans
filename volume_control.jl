# Volume control through the open boundaries: the net inflow through the Atlantic and Pacific sides together is pinned
# at every barotropic substep (NumericalEarth.jl#938 `net_inflow`, Oceananigans.jl#6524 `NetInflow`), as the volume
# constraints of NEMO and ROMS do; Flather still splits it between the two sides. Without it, Flather alone lets the
# domain's volume drift away from GLORYS's (the mean η gains or loses several centimetres in days).
#
# With `[boundaries] follow_glorys_volume = true` (the default), the net inflow is set before every time step to
#
#     T = A dη̄ᴳ/dt − F + A [(η̄ᴳ − η̄ᴳ₀) − (η̄ − η̄₀)] / τ    (the first term only with `follow_glorys_tendency = true`),
#
# with A the wet area, η̄ᴳ GLORYS's domain-mean η (daily means, joined by a monotone cubic so that dη̄ᴳ/dt is
# continuous), F the surface freshwater volume flux into the ocean (precipitation, evaporation, rivers, ice), η̄ the
# model's domain-mean η and the subscript 0 the start of the run. The first term follows GLORYS's volume change, the
# second hands the surface freshwater back through the open sides, and the third removes any remaining drift over
# τ = `volume_relaxation_days`. Without the first term the relaxation alone follows GLORYS, about τ behind.
# The inverse-barometer sea level has its domain mean removed, so it does not change the volume.
#
# With `follow_glorys_volume = false` and `net_inflow_Sv` set, the net inflow is that constant instead; with neither,
# Flather alone (no volume control).
#
# `net_inflow_depth_exponent` p spreads the pin's correction over the boundary columns in proportion to Hᵖ: 0 uniform
# per length, 1 a uniform velocity (NEMO, ROMS), ½ a uniform shift of the exterior η in Flather.

using Oceananigans.Operators: Azᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: bottom_height_field
using Oceananigans.OutputReaders: Clamp
using NumericalEarth.DataWrangling: bracketing_dates, default_horizontal_padding
using NumericalEarth.NestedModels: boundary_location, sample_parent!
using Printf

# Net inflow for nested_ocean_simulation (m³ s⁻¹, positive into the domain): a Ref the target below sets every time
# step, a constant, or nothing
function net_inflow_setting(boundaries)
    get(boundaries, "follow_glorys_volume", false) && return Ref(0.0)
    Sv = get(boundaries, "net_inflow_Sv", nothing)
    return isnothing(Sv) ? nothing : 1e6 * Sv
end

# GLORYS's free surface on its own grid, over the region of `grid` and bracketing `dates`, held at the ends
function parent_free_surface(grid, dataset, dates; dir_kw...)
    region   = BoundingBox(grid.underlying_grid; padding = default_horizontal_padding(dataset))
    native   = bracketing_dates(dataset, :free_surface, dates)
    metadata = Metadata(:free_surface; dataset, dates = native, region, dir_kw...)
    return FieldTimeSeries(metadata, Oceananigans.Architectures.architecture(grid);
                           time_indices_in_memory = length(native), time_indexing = Clamp(), start_time = first(dates))
end

struct GLORYSVolumeTarget{R, C, J, H, M}
    net_inflow      :: R        # Ref read by the barotropic pin (m³ s⁻¹, positive into the domain)
    cell_area       :: C        # horizontal area of each wet cell (m²), zero where dry, on the device
    freshwater      :: J        # surface freshwater volume flux added to η (m s⁻¹, positive into the ocean)
    displacement    :: H        # η
    times           :: Vector{Float64}   # times of GLORYS's daily means (s since the start of the run)
    means           :: Vector{Float64}   # GLORYS's domain-mean η at those times (m)
    slopes          :: Vector{Float64}   # dη̄ᴳ/dt at those times (m s⁻¹)
    area            :: Float64
    relaxation_time :: Float64
    tendency        :: Bool
    initial_means   :: M        # (model, GLORYS) domain-mean η at the start (m)
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
                                child.free_surface.displacement, times, means, monotone_slopes(times, means), A,
                                Float64(relaxation_time), tendency, Ref((0.0, 0.0)))
    target.initial_means[] = (model_mean(target), glorys_mean(target, child.clock.time)[1])
    return target
end

model_mean(target) = sum(view(interior(target.displacement), :, :, 1) .* target.cell_area) / target.area

# Slopes at the daily means: the harmonic mean of the neighbouring secants, zero where they differ in sign
# (Fritsch & Butland 1984), so the cubic between two means never overshoots them
function monotone_slopes(times, means)
    δ = diff(means) ./ diff(times)
    inner = [δ[n - 1] * δ[n] > 0 ? 2 / (1 / δ[n - 1] + 1 / δ[n]) : 0.0 for n in 2:length(δ)]
    return [δ[1]; inner; δ[end]]
end

# GLORYS's domain-mean η and its rate of change at time t: the cubic Hermite interpolant through the daily means and
# their slopes, held before the first
function glorys_mean(target, t)
    times, means, slopes = target.times, target.means, target.slopes
    t ≤ times[1] && return means[1], 0.0
    n = min(searchsortedlast(times, t), length(times) - 1)
    h = times[n + 1] - times[n]
    s = (t - times[n]) / h
    η̄ᴳ = (2s^3 - 3s^2 + 1) * means[n] + (s^3 - 2s^2 + s) * h * slopes[n] +
          (3s^2 - 2s^3) * means[n + 1] + (s^3 - s^2) * h * slopes[n + 1]
    ∂tη̄ᴳ = ((6s^2 - 6s) * means[n] + (6s - 6s^2) * means[n + 1]) / h +
            (3s^2 - 4s + 1) * slopes[n] + (3s^2 - 2s) * slopes[n + 1]
    return η̄ᴳ, ∂tη̄ᴳ
end

function set_net_inflow!(sim, target)
    η̄ᴳ, ∂tη̄ᴳ = glorys_mean(target, time(sim))
    η̄₀, η̄ᴳ₀  = target.initial_means[]
    F = sum(view(interior(target.freshwater), :, :, 1) .* target.cell_area)
    A = target.area
    target.net_inflow[] = A * ∂tη̄ᴳ * target.tendency - F + A * ((η̄ᴳ - η̄ᴳ₀) - (model_mean(target) - η̄₀)) / target.relaxation_time
    return nothing
end

# Hourly: the change in domain-mean η since the start, model and GLORYS, and the net inflow the pin holds
function log_volume_target(sim, target)
    η̄ᴳ, _   = glorys_mean(target, time(sim))
    η̄₀, η̄ᴳ₀ = target.initial_means[]
    @info @sprintf("    volume: Δη̄ model %+.2f cm, GLORYS %+.2f cm, net inflow %+.3f Sv",
                   100 * (model_mean(target) - η̄₀), 100 * (η̄ᴳ - η̄ᴳ₀), 1e-6 * target.net_inflow[])
    return nothing
end

function add_glorys_volume_target!(simulation, net_inflow, child, parent_η, relaxation_time; tendency = false)
    target = GLORYSVolumeTarget(net_inflow, child, parent_η, relaxation_time, tendency)
    add_callback!(simulation, sim -> set_net_inflow!(sim, target), IterationInterval(1); name = :glorys_volume_target)
    add_callback!(simulation, sim -> log_volume_target(sim, target), TimeInterval(1hours); name = :volume_log)
    @info @sprintf("Net inflow follows GLORYS's domain-mean η (relaxation time %.1f days, tendency %s)",
                   relaxation_time / days, tendency)
    return target
end
