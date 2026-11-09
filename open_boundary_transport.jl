# Barotropic transport through the open edges (Atlantic, j = 1; Pacific, j = Ny + 1). Included by panarctic.jl
# after V_obcs is built. Three parts:
#
# 1. GLORYS's net transport through each edge, from the boundary data, at each daily snapshot.
# 2. With `barotropic = "hybrid"` in [boundaries]: the hybrid condition replaces plain Flather on the barotropic
#    transport. On each edge, Flather about ηᵉˣᵗ + C, with C chosen every barotropic substep so that the edge's net
#    transport equals GLORYS's (a uniform η shift; the transport changes in proportion to √(gH)). Then one uniform
#    velocity w on both edges sets the total net inflow, chosen by `total_net_inflow_Sv`:
#      "glorys_volume"  the model's volume changes at the rate GLORYS's does over the same wet domain:
#                       total = A d⟨ηᴳ⟩/dt − (model's surface freshwater flux). Only the rate of change of GLORYS's mean
#                       free surface enters, so a datum offset drops out; the model's own surface fluxes (ice growth
#                       and melt, rivers, P − E) are accounted for every time step.
#      a number (Sv)    a fixed total
#      "glorys"         GLORYS's own total (w = 0). It need not close on the model grid.
#    A datum difference between GLORYS's and the model's free surface then moves no water, and the throughflow
#    follows GLORYS's. Tested in panarctic-tools/sandbox/obc-vortex.
# 3. An hourly monitor, transport_monitor.csv in the output directory: model and GLORYS net transport per edge, the
#    hybrid's total, the surface freshwater volume flux (positive into the ocean), the area-weighted mean free surface
#    and the residual of the volume budget.

using Oceananigans.Grids: static_column_depthᶜᶠᵃ, static_column_depthᶜᶜᵃ
using Oceananigans.Operators: Δxᶜᶠᶜ, Azᶜᶜᶜ

const g_obc = Oceananigans.defaults.gravitational_acceleration

# ## 1. Edge geometry and GLORYS's net transport per edge

edge_rows = (1, Ny + 1)
Δx_edge, H_edge = CUDA.@allowscalar begin
    ([Δxᶜᶠᶜ(i, j, 1, grid) for i in 1:Nx, j in edge_rows],
     [max(static_column_depthᶜᶠᵃ(i, j, grid), 0.0) for i in 1:Nx, j in edge_rows])
end

function glorys_edge_transport(n, e)
    t = Oceananigans.Units.Time(vᶠ.times[n])
    Q = 0.0
    CUDA.@allowscalar for i in 1:Nx
        Q += vertical_integral(i, edge_rows[e], grid, vᶠ, t, Δzᶜᶠᶜ, Center(), Face(), Center()) * Δx_edge[i, e]
    end
    return Q
end

atlantic_glorys = [glorys_edge_transport(n, 1) for n in eachindex(vᶠ.times)]   # into the domain
pacific_glorys  = [glorys_edge_transport(n, 2) for n in eachindex(vᶠ.times)]   # northward: out of the domain
glorys_net      = atlantic_glorys .- pacific_glorys                              # net inflow

t₀_obc = vᶠ.times[1]
Δt_obc = vᶠ.times[2] - vᶠ.times[1]
all(diff(vᶠ.times) .≈ Δt_obc) || error("the boundary data must be evenly spaced in time")

@info @sprintf("GLORYS at the start: Atlantic %.2f Sv in, Pacific %.2f Sv out, net %.3f Sv",
               atlantic_glorys[start_index] / 1e6, pacific_glorys[start_index] / 1e6, glorys_net[start_index] / 1e6)

# Wet columns, their areas, and the rate of change of GLORYS's volume over the same wet domain
wet_column, cell_area = CUDA.@allowscalar begin
    ([static_column_depthᶜᶜᵃ(i, j, grid) > 0 for i in 1:Nx, j in 1:Ny], [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny])
end
wet_area    = cell_area .* wet_column
domain_area = sum(wet_area)
all(ηᵉˣᵗ.times .≈ vᶠ.times) || error("GLORYS η and v must share their snapshot times")
glorys_mean_η = [sum(Array(interior(ηᵉˣᵗ[n]))[:, :, 1] .* wet_area) / domain_area for n in eachindex(ηᵉˣᵗ.times)]
glorys_volume_rate = domain_area .* diff(glorys_mean_η) ./ Δt_obc    # between snapshots n and n + 1

# Linear interpolation of a snapshot series to time t (clamped to the first and last snapshots)
@inline function in_time(series, t, t₀, Δt)
    N = length(series)
    s = clamp((t - t₀) / Δt, zero(t), convert(typeof(t), N - 1))
    n = min(floor(Int, s) + 1, N - 1)
    α = s - (n - 1)
    return @inbounds (1 - α) * series[n] + α * series[n + 1]
end

# The value of a series defined between snapshots (length N - 1) at time t
@inline between_snapshots(series, t, t₀, Δt) = @inbounds series[clamp(floor(Int, (t - t₀) / Δt) + 1, 1, length(series))]

# ## 2. The hybrid condition

# C for one edge: the √(gH)-weighted mean of ηᵇ − ηᵉˣᵗ along it. jᶠ: the boundary face row; jᵢ, jₕ: the interior and
# halo rows of η either side of it; jₑ: the row of ηᵉˣᵗ used by atlantic_V_obc / pacific_V_obc.
@inline function edge_offset(jᶠ, jᵢ, jₕ, jₑ, grid, clock, fields, ηᵉˣᵗ)
    t  = Oceananigans.Units.Time(clock.time)
    η  = fields.η
    kη = grid.Nz + 1
    num = zero(grid)
    den = zero(grid)
    for i′ in 1:grid.Nx
        H  = static_column_depthᶜᶠᵃ(i′, jᶠ, grid)
        c  = sqrt(g_obc * max(H, zero(H)))
        Δx = Δxᶜᶠᶜ(i′, jᶠ, 1, grid)
        ηᵇ = @inbounds (η[i′, jᵢ, kη] + η[i′, jₕ, kη]) / 2
        num += ifelse(H > 0, c * (ηᵇ - @inbounds(ηᵉˣᵗ[i′, jₑ, 1, t])) * Δx, zero(grid))
        den += c * Δx
    end
    return ifelse(den > 0, num / den, zero(grid))
end

@inline function hybrid_velocity(t, p)
    target = ifelse(p.follow_volume, between_snapshots(p.rate, t, p.t₀, p.Δt) - @inbounds(p.freshwater[1]), p.target)
    return ifelse(p.use_total, (in_time(p.net, t, p.t₀, p.Δt) - target) / p.area, zero(t))
end

@inline function atlantic_V_hybrid(i, k, grid, clock, fields, p)
    V, ηᵉ = atlantic_V_obc(i, k, grid, clock, fields, p)
    C = edge_offset(1, 1, 0, 1, grid, clock, fields, p.η)
    H = max(static_column_depthᶜᶠᵃ(i, 1, grid), zero(grid))
    return (V - H * hybrid_velocity(clock.time, p), ηᵉ + C)
end

@inline function pacific_V_hybrid(i, k, grid, clock, fields, p)
    V, ηᵉ = pacific_V_obc(i, k, grid, clock, fields, p)
    C = edge_offset(grid.Ny + 1, grid.Ny, grid.Ny + 1, grid.Ny, grid, clock, fields, p.η)
    H = max(static_column_depthᶜᶠᵃ(i, grid.Ny + 1, grid), zero(grid))
    return (V + H * hybrid_velocity(clock.time, p), ηᵉ + C)
end

barotropic_scheme = get(bc, "barotropic", "flather")
total_setting = get(bc, "total_net_inflow_Sv", "glorys_volume")
follow_volume = total_setting == "glorys_volume"
use_total     = total_setting != "glorys"
target_inflow = total_setting isa Number ? total_setting * 1e6 : 0.0
open_area = sum(H_edge .* Δx_edge)

# The model's surface freshwater volume flux over the wet domain (positive into the ocean), updated every time step
surface_freshwater = on_architecture(arch, [0.0])
wet_area_device = on_architecture(arch, wet_area)
surface_freshwater_flux(ocean_model) = NumericalEarth.Oceans.extract_freshwater_flux(ocean_model.tracers.S.boundary_conditions.top.condition)

function update_surface_freshwater!(sim)
    Jʷ = surface_freshwater_flux(sim.model.ocean.model)
    fill!(surface_freshwater, sum(interior(Jʷ) .* wet_area_device))
    return nothing
end

if barotropic_scheme == "hybrid"
    p_hybrid = (v = vᶠ, η = ηᵉˣᵗ, net = on_architecture(arch, glorys_net), rate = on_architecture(arch, glorys_volume_rate),
                freshwater = surface_freshwater, t₀ = t₀_obc, Δt = Δt_obc, area = open_area, target = target_inflow,
                follow_volume, use_total)
    V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
        south = GravityWaveRadiationBoundaryCondition(atlantic_V_hybrid, discrete_form = true, parameters = p_hybrid),
        north = GravityWaveRadiationBoundaryCondition(pacific_V_hybrid,  discrete_form = true, parameters = p_hybrid))
    @info "Barotropic transport: hybrid (each edge's net to GLORYS's; total: " *
          (follow_volume ? "the model's volume follows GLORYS's)" :
           use_total ? @sprintf("%.3f Sv)", target_inflow / 1e6) : "GLORYS's own)")
elseif barotropic_scheme == "flather"
    @info "Barotropic transport: Flather"
else
    error("Unknown barotropic scheme \"$barotropic_scheme\"; use \"flather\" or \"hybrid\"")
end

# ## 3. Hourly transport monitor

monitor_file = joinpath(output_dir, "transport_monitor.csv")
open(monitor_file, "w") do io
    println(io, "hour,atlantic_model_Sv,atlantic_glorys_Sv,pacific_model_Sv,pacific_glorys_Sv,net_model_Sv,",
                "net_glorys_Sv,hybrid_total_Sv,glorys_volume_rate_Sv,surface_freshwater_Sv,eta_mean_m,budget_residual_Sv")
end
previous_monitor = Ref((t = NaN, η = NaN))

function transport_monitor(sim)
    model = sim.model.ocean.model
    t = model.clock.time
    V = Array(interior(model.free_surface.barotropic_velocities.V))[:, :, 1]
    η = Array(interior(model.free_surface.displacement))[:, :, 1]
    atlantic = sum(V[:, 1] .* Δx_edge[:, 1])
    pacific  = sum(V[:, Ny + 1] .* Δx_edge[:, 2])
    freshwater = sum(Array(interior(surface_freshwater_flux(model)))[:, :, 1] .* wet_area)
    mean_η = sum(η .* wet_area) / domain_area
    rate = between_snapshots(glorys_volume_rate, t, t₀_obc, Δt_obc)
    total = follow_volume ? rate - freshwater : use_total ? target_inflow : in_time(glorys_net, t, t₀_obc, Δt_obc)
    # A d⟨η⟩/dt over the last hour minus what the boundaries and the surface supplied now (zero if the budget closes)
    residual = domain_area * (mean_η - previous_monitor[].η) / (t - previous_monitor[].t) - (atlantic - pacific + freshwater)
    previous_monitor[] = (; t, η = mean_η)
    open(monitor_file, "a") do io
        @printf(io, "%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.4f\n", (t - start_time) / 3600,
                atlantic / 1e6, in_time(atlantic_glorys, t, t₀_obc, Δt_obc) / 1e6,
                pacific / 1e6, in_time(pacific_glorys, t, t₀_obc, Δt_obc) / 1e6,
                (atlantic - pacific) / 1e6, in_time(glorys_net, t, t₀_obc, Δt_obc) / 1e6,
                total / 1e6, rate / 1e6, freshwater / 1e6, mean_η, residual / 1e6)
    end
    return nothing
end
