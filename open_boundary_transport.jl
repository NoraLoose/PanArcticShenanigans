# Barotropic transport through the open edges (Atlantic, j = 1; Pacific, j = Ny + 1). Included by panarctic.jl
# after V_obcs is built. Three parts:
#
# 1. GLORYS's net transport through each edge, from the boundary data, at each daily snapshot.
# 2. With `barotropic = "hybrid"` in [boundaries]: the hybrid condition replaces plain Flather on the barotropic
#    transport. On each edge, Flather about ηᵉˣᵗ + C, with C chosen every barotropic substep so that the edge's net
#    transport equals GLORYS's (a uniform η shift; the transport changes in proportion to √(gH)). Then, if
#    `total_net_inflow_Sv` is a number, one uniform velocity w on both edges brings the total net inflow to it; with
#    "glorys" the total is GLORYS's own (w = 0). A datum difference between GLORYS's and the model's free surface then
#    moves no water, and the throughflow follows GLORYS's. Tested in panarctic-tools/sandbox/obc-vortex.
# 3. An hourly monitor, transport_monitor.csv in the output directory: model and GLORYS net transport per edge, the
#    domain-integrated surface freshwater volume flux and the domain-mean free surface (to close the volume budget
#    before tying the total to the surface fluxes).

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

# Linear interpolation of a snapshot series to time t (clamped to the first and last snapshots)
@inline function in_time(series, t, t₀, Δt)
    N = length(series)
    s = clamp((t - t₀) / Δt, zero(t), convert(typeof(t), N - 1))
    n = min(floor(Int, s) + 1, N - 1)
    α = s - (n - 1)
    return @inbounds (1 - α) * series[n] + α * series[n + 1]
end

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

@inline hybrid_velocity(t, p) = (in_time(p.net, t, p.t₀, p.Δt) - p.target) / p.area

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
total_setting = get(bc, "total_net_inflow_Sv", "glorys")
total_from_glorys = total_setting == "glorys"
target_inflow = total_from_glorys ? 0.0 : total_setting * 1e6
open_area = sum(H_edge .* Δx_edge)
corrected_net = total_from_glorys ? zero(glorys_net) : glorys_net     # w = 0 when the total is GLORYS's

if barotropic_scheme == "hybrid"
    p_hybrid = (v = vᶠ, η = ηᵉˣᵗ, net = on_architecture(arch, corrected_net), t₀ = t₀_obc, Δt = Δt_obc,
                area = open_area, target = target_inflow)
    V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
        south = GravityWaveRadiationBoundaryCondition(atlantic_V_hybrid, discrete_form = true, parameters = p_hybrid),
        north = GravityWaveRadiationBoundaryCondition(pacific_V_hybrid,  discrete_form = true, parameters = p_hybrid))
    @info "Barotropic transport: hybrid (each edge's net to GLORYS's; total: " *
          (total_from_glorys ? "GLORYS's)" : @sprintf("%.3f Sv)", target_inflow / 1e6))
elseif barotropic_scheme == "flather"
    @info "Barotropic transport: Flather"
else
    error("Unknown barotropic scheme \"$barotropic_scheme\"; use \"flather\" or \"hybrid\"")
end

# ## 3. Hourly transport monitor

monitor_file = joinpath(output_dir, "transport_monitor.csv")
wet_column, cell_area = CUDA.@allowscalar begin
    ([static_column_depthᶜᶜᵃ(i, j, grid) > 0 for i in 1:Nx, j in 1:Ny], [Azᶜᶜᶜ(i, j, 1, grid) for i in 1:Nx, j in 1:Ny])
end
open(monitor_file, "w") do io
    println(io, "hour,atlantic_model_Sv,atlantic_glorys_Sv,pacific_model_Sv,pacific_glorys_Sv,net_model_Sv,",
                "net_glorys_Sv,hybrid_w_mm_per_s,surface_freshwater_Sv,eta_mean_m")
end

function transport_monitor(sim)
    model = sim.model.ocean.model
    t = model.clock.time
    V = Array(interior(model.free_surface.barotropic_velocities.V))[:, :, 1]
    η = Array(interior(model.free_surface.displacement))[:, :, 1]
    atlantic = sum(V[:, 1] .* Δx_edge[:, 1])
    pacific  = sum(V[:, Ny + 1] .* Δx_edge[:, 2])
    Jʷ = Array(interior(NumericalEarth.Oceans.extract_freshwater_flux(model.tracers.S.boundary_conditions.top.condition)))
    freshwater = sum(Jʷ[:, :, 1][wet_column] .* cell_area[wet_column])   # sign as NumericalEarth's carrying flux
    w = barotropic_scheme == "hybrid" ? (in_time(corrected_net, t, t₀_obc, Δt_obc) - target_inflow) / open_area : 0.0
    open(monitor_file, "a") do io
        @printf(io, "%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f\n", (t - start_time) / 3600,
                atlantic / 1e6, in_time(atlantic_glorys, t, t₀_obc, Δt_obc) / 1e6,
                pacific / 1e6, in_time(pacific_glorys, t, t₀_obc, Δt_obc) / 1e6,
                (atlantic - pacific) / 1e6, in_time(glorys_net, t, t₀_obc, Δt_obc) / 1e6,
                1e3 * w, freshwater / 1e6, sum(η[wet_column]) / count(wet_column))
    end
    return nothing
end
