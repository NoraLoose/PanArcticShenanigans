# Barotropic transport through the open edges (Atlantic, j = 1; Pacific, j = Ny + 1). Included by panarctic.jl
# after V_obcs is built. Three parts:
#
# 1. GLORYS's net transport through each edge, from the boundary data, at each daily snapshot.
# 2. With `barotropic = "hybrid"` in [boundaries] (or "total_only", "clamped" for comparison; see section 2): the
#    hybrid condition replaces plain Flather on the barotropic transport. On each edge, Flather about ηᵉˣᵗ + C, with C
#    chosen every barotropic substep so that the edge's net transport equals its target (a uniform η shift; the
#    transport changes in proportion to √(gH)). The target is GLORYS's transport through the same line on GLORYS's own
#    grid (`per_edge_target = "native"`), or after interpolation onto the model grid ("interpolated"). Then one uniform
#    velocity w on both edges sets the total net inflow, chosen by `total_net_inflow_Sv`:
#      "glorys_transport"  GLORYS's net transport (the sum of the per-edge targets) plus a slow nudge,
#                          A (Δ⟨ηᴳ⟩ − Δ⟨η⟩) / τ with τ = `volume_nudging_days`, that keeps the change of the model's
#                          domain-mean free surface close to GLORYS's over weeks (changes from the start only, so a
#                          datum drops out). Transports are trusted over days, sea level over weeks.
#      "glorys_volume"     the model's volume changes at GLORYS's rate every step: A d⟨ηᴳ⟩/dt − surface freshwater flux
#      a number (Sv)       a fixed total
#      "glorys"            no total correction (w = 0)
#    A datum difference between GLORYS's and the model's free surface then moves no water. Tested in
#    panarctic-tools/sandbox/obc-vortex.
# 3. An hourly monitor, transport_monitor.csv in the output directory: model and GLORYS net transport per edge (on the
#    model grid and on GLORYS's own grid), the target total and its nudge, the surface freshwater volume flux (positive
#    into the ocean), the area-weighted mean free surface of the model and of GLORYS, and the volume budget residual.

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

# The same transports on GLORYS's own grid (reads the daily u, v files once)
include("glorys_sections.jl")
atlantic_native, pacific_native = native_edge_transports(dates)
length(atlantic_native) == length(vᶠ.times) || error("one GLORYS file per boundary snapshot expected")

@info @sprintf("GLORYS at the start, model grid: Atlantic %.2f Sv in, Pacific %.2f Sv out; own grid: %.2f, %.2f Sv",
               atlantic_glorys[start_index] / 1e6, pacific_glorys[start_index] / 1e6,
               atlantic_native[start_index] / 1e6, pacific_native[start_index] / 1e6)

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

# ## 2. Constrained conditions on the barotropic transport
#
# Three methods share one boundary function, all with Flather underneath (V = Vᵉˣᵗ ∓ c(ηᵇ − ηᵉˣᵗ), c = √(gH)):
#   hybrid      each edge's net set to GLORYS's by a uniform η shift C (∝ √(gH)), then the total by w
#   total_only  Flather on each edge, then the total by w (NEMO `ln_vol`, ROMS `VolCons`)
#   clamped     V = Vᵉˣᵗ (ηᵉˣᵗ replaced by ηᵇ, so the Flather term vanishes), then the total by w (MITgcm prescribed
#               velocities with `useOBCSbalance`)
# w is one uniform velocity on every open face (∝ H), removed from the inflow, that brings the total net inflow to the
# target set by `total_net_inflow_Sv`.

# Sums along one edge: Σ c (ηᵇ − ηᵉˣᵗ) Δx and Σ c Δx. jᶠ: the boundary face row; jᵢ, jₕ: the interior and halo rows of
# η either side of it; jₑ: the row of ηᵉˣᵗ used by atlantic_V_obc / pacific_V_obc.
@inline function edge_mismatch(jᶠ, jᵢ, jₕ, jₑ, grid, clock, fields, ηᵉˣᵗ)
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
    return num, den
end

# Per-edge corrections and the total. Per edge: hybrid shifts ηᵉˣᵗ by C (√(gH)-weighted); clamped subtracts a uniform
# velocity u (∝ H), used only to move from the interpolated to the native target. Then w (∝ H) on both edges sets the
# total. Returns (C_atlantic, C_pacific, u_atlantic, u_pacific, w).
@inline function constraints(grid, clock, fields, p)
    t = clock.time
    nA, dA = edge_mismatch(1, 1, 0, 1, grid, clock, fields, p.η)
    nP, dP = edge_mismatch(grid.Ny + 1, grid.Ny, grid.Ny + 1, grid.Ny, grid, clock, fields, p.η)
    VA = in_time(p.atlantic, t, p.t₀, p.Δt)                     # Σ Vᵉˣᵗ Δx on the model grid: into the domain
    VP = in_time(p.pacific,  t, p.t₀, p.Δt)                     # out of the domain
    TA = in_time(p.atlantic_target, t, p.t₀, p.Δt)              # the per-edge targets
    TP = in_time(p.pacific_target,  t, p.t₀, p.Δt)
    CA = ifelse(p.edge_shift & (dA > 0), (nA - (VA - TA)) / dA, zero(grid))
    CP = ifelse(p.edge_shift & (dP > 0), (nP + (VP - TP)) / dP, zero(grid))
    uA = ifelse(p.clamped, (VA - TA) / p.atlantic_area, zero(grid))
    uP = ifelse(p.clamped, (VP - TP) / p.pacific_area,  zero(grid))
    flather = !p.clamped
    inflow = (VA - ifelse(flather, nA, zero(grid)) + CA * dA - uA * p.atlantic_area) -
             (VP + ifelse(flather, nP, zero(grid)) - CP * dP - uP * p.pacific_area)
    nudge  = p.domain_area * p.nudging_rate * (in_time(p.glorys_Δη, t, p.t₀, p.Δt) - @inbounds(p.model_Δη[1]))
    target = ifelse(p.follow_volume, between_snapshots(p.rate, t, p.t₀, p.Δt) - @inbounds(p.freshwater[1]),
             ifelse(p.follow_transport, TA - TP + nudge, p.target))
    w = ifelse(p.use_total, (inflow - target) / p.area, zero(grid))
    return CA, CP, uA, uP, w
end

@inline function atlantic_V_constrained(i, k, grid, clock, fields, p)
    V, ηᵉ = atlantic_V_obc(i, k, grid, clock, fields, p)
    CA, CP, uA, uP, w = constraints(grid, clock, fields, p)
    H  = max(static_column_depthᶜᶠᵃ(i, 1, grid), zero(grid))
    ηᵇ = @inbounds (fields.η[i, 0, grid.Nz + 1] + fields.η[i, 1, grid.Nz + 1]) / 2
    return (V - H * (uA + w), ifelse(p.clamped, ηᵇ, ηᵉ + CA))
end

@inline function pacific_V_constrained(i, k, grid, clock, fields, p)
    V, ηᵉ = pacific_V_obc(i, k, grid, clock, fields, p)
    CA, CP, uA, uP, w = constraints(grid, clock, fields, p)
    H  = max(static_column_depthᶜᶠᵃ(i, grid.Ny + 1, grid), zero(grid))
    ηᵇ = @inbounds (fields.η[i, grid.Ny, grid.Nz + 1] + fields.η[i, grid.Ny + 1, grid.Nz + 1]) / 2
    return (V - H * (uP - w), ifelse(p.clamped, ηᵇ, ηᵉ + CP))
end

barotropic_scheme = get(bc, "barotropic", "flather")
total_setting    = get(bc, "total_net_inflow_Sv", "glorys_transport")
follow_transport = total_setting == "glorys_transport"
follow_volume    = total_setting == "glorys_volume"
use_total        = total_setting != "glorys"
target_inflow    = total_setting isa Number ? total_setting * 1e6 : 0.0
nudging_days     = get(bc, "volume_nudging_days", 14)
nudging_rate     = nudging_days > 0 ? 1 / (nudging_days * days) : 0.0
edge_target      = get(bc, "per_edge_target", "native")
edge_target in ("native", "interpolated") || error("per_edge_target must be \"native\" or \"interpolated\"")
atlantic_target  = edge_target == "native" ? atlantic_native : atlantic_glorys
pacific_target   = edge_target == "native" ? pacific_native  : pacific_glorys
atlantic_area, pacific_area = sum(H_edge[:, 1] .* Δx_edge[:, 1]), sum(H_edge[:, 2] .* Δx_edge[:, 2])
open_area = atlantic_area + pacific_area
glorys_Δη = glorys_mean_η .- in_time(glorys_mean_η, start_time, t₀_obc, Δt_obc)    # change since the start of the run

# The model's surface freshwater volume flux over the wet domain (positive into the ocean) and the change of its
# area-weighted mean free surface since the start, updated every time step
surface_freshwater = on_architecture(arch, [0.0])
model_Δη = on_architecture(arch, [0.0])
model_mean_η₀ = Ref(NaN)
wet_area_device = on_architecture(arch, wet_area)
surface_freshwater_flux(ocean_model) = NumericalEarth.Oceans.extract_freshwater_flux(ocean_model.tracers.S.boundary_conditions.top.condition)

function update_surface_freshwater!(sim)
    ocean_model = sim.model.ocean.model
    fill!(surface_freshwater, sum(interior(surface_freshwater_flux(ocean_model)) .* wet_area_device))
    mean_η = sum(interior(ocean_model.free_surface.displacement) .* wet_area_device) / domain_area
    isnan(model_mean_η₀[]) && (model_mean_η₀[] = mean_η)
    fill!(model_Δη, mean_η - model_mean_η₀[])
    return nothing
end

constrained_schemes = ("hybrid", "total_only", "clamped")
if barotropic_scheme in constrained_schemes
    p_constrained = (v = vᶠ, η = ηᵉˣᵗ, atlantic = on_architecture(arch, atlantic_glorys),
                     pacific = on_architecture(arch, pacific_glorys),
                     atlantic_target = on_architecture(arch, atlantic_target),
                     pacific_target = on_architecture(arch, pacific_target),
                     rate = on_architecture(arch, glorys_volume_rate), glorys_Δη = on_architecture(arch, glorys_Δη),
                     freshwater = surface_freshwater, model_Δη, t₀ = t₀_obc, Δt = Δt_obc, area = open_area,
                     atlantic_area, pacific_area, domain_area, nudging_rate, target = target_inflow,
                     follow_transport, follow_volume, use_total, edge_shift = barotropic_scheme == "hybrid",
                     clamped = barotropic_scheme == "clamped")
    V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
        south = GravityWaveRadiationBoundaryCondition(atlantic_V_constrained, discrete_form = true, parameters = p_constrained),
        north = GravityWaveRadiationBoundaryCondition(pacific_V_constrained,  discrete_form = true, parameters = p_constrained))
    @info "Barotropic transport: $barotropic_scheme (per-edge target: GLORYS $edge_target; total: " *
          (follow_transport ? "GLORYS's transport, volume nudged over $nudging_days days)" :
           follow_volume ? "the model's volume follows GLORYS's every step)" :
           use_total ? @sprintf("%.3f Sv)", target_inflow / 1e6) : "no correction)")
elseif barotropic_scheme == "flather"
    @info "Barotropic transport: Flather"
else
    error("Unknown barotropic scheme \"$barotropic_scheme\"; use \"flather\", \"hybrid\", \"total_only\" or \"clamped\"")
end

# ## 3. Hourly transport monitor

monitor_file = joinpath(output_dir, "transport_monitor.csv")
open(monitor_file, "w") do io
    println(io, "hour,atlantic_model_Sv,atlantic_glorys_Sv,pacific_model_Sv,pacific_glorys_Sv,net_model_Sv,",
                "net_glorys_Sv,target_total_Sv,glorys_volume_rate_Sv,surface_freshwater_Sv,eta_mean_m,glorys_eta_mean_m,",
                "budget_residual_Sv,atlantic_native_Sv,pacific_native_Sv,nudge_Sv")
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
    TA, TP = in_time(atlantic_target, t, t₀_obc, Δt_obc), in_time(pacific_target, t, t₀_obc, Δt_obc)
    nudge = isnan(model_mean_η₀[]) ? 0.0 :
            domain_area * nudging_rate * (in_time(glorys_Δη, t, t₀_obc, Δt_obc) - (mean_η - model_mean_η₀[]))
    total = barotropic_scheme == "flather" ? NaN :
            follow_transport ? TA - TP + nudge : follow_volume ? rate - freshwater :
            use_total ? target_inflow : NaN
    # A d⟨η⟩/dt over the last hour minus what the boundaries and the surface supplied now (zero if the budget closes)
    residual = domain_area * (mean_η - previous_monitor[].η) / (t - previous_monitor[].t) - (atlantic - pacific + freshwater)
    previous_monitor[] = (; t, η = mean_η)
    open(monitor_file, "a") do io
        @printf(io, "%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.4f,%.4f,%.4f,%.4f\n", (t - start_time) / 3600,
                atlantic / 1e6, in_time(atlantic_glorys, t, t₀_obc, Δt_obc) / 1e6,
                pacific / 1e6, in_time(pacific_glorys, t, t₀_obc, Δt_obc) / 1e6,
                (atlantic - pacific) / 1e6, in_time(glorys_net, t, t₀_obc, Δt_obc) / 1e6,
                total / 1e6, rate / 1e6, freshwater / 1e6, mean_η, in_time(glorys_mean_η, t, t₀_obc, Δt_obc),
                residual / 1e6, in_time(atlantic_native, t, t₀_obc, Δt_obc) / 1e6,
                in_time(pacific_native, t, t₀_obc, Δt_obc) / 1e6, follow_transport ? nudge / 1e6 : NaN)
    end
    return nothing
end
