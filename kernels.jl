# Functions called inside GPU kernels: the open boundary conditions, the hybrid constraint, and the rotation of GLORYS
# velocities into the grid frame. panarctic.jl includes this file only when it has not been loaded in this Julia
# session or its contents changed. Redefining a function, even with identical text, makes every kernel that uses it
# stale and the whole time step compiles again (~15 min), so keep these definitions here, keep them free of global
# variables (other than constants), and put everything else in the other files.

using Oceananigans.Operators: Δzᶜᶠᶜ, Δxᶜᶠᶜ, intrinsic_vector
using Oceananigans.Grids: static_column_depthᶜᶠᵃ
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node

const g_obc = Oceananigans.defaults.gravitational_acceleration

# ## Rotation of GLORYS velocities from east/north into the grid frame

@inline grid_frame_u(i, j, k, grid, uₑ, vₑ) = intrinsic_vector(i, j, k, grid, uₑ, vₑ)[1]
@inline grid_frame_v(i, j, k, grid, uₑ, vₑ) = intrinsic_vector(i, j, k, grid, uₑ, vₑ)[2]

# ## GLORYS values on the open edges (discrete-form boundary conditions)

@inline atlantic_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1,          k, Oceananigans.Units.Time(clock.time)]
@inline  pacific_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny,    k, Oceananigans.Units.Time(clock.time)]
@inline pacific_face_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny+1, k, Oceananigans.Units.Time(clock.time)]

# ## Depth integrals on the edges, and the Flather data (V, η) for the barotropic transport

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

@inline function vertical_integral(i, j, grid, u, t, Δz, ℓx, ℓy, ℓz)
    U = zero(eltype(grid))
    @inbounds for k in 1:grid.Nz
        wet = wetcell(i, j, k, grid, ℓx, ℓy, ℓz)
        U += ifelse(wet, u[i, j, k, t] * Δz(i, j, k, grid), zero(U))
    end
    return U
end

@inline function atlantic_V_obc(i, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time)
    V = vertical_integral(i, 1, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, 1, 1, t])
end

@inline function pacific_V_obc(i, k, grid, clock, fields, p)
    j = grid.Ny+1
    t = Oceananigans.Units.Time(clock.time)
    V = vertical_integral(i, j, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, grid.Ny, 1, t])
end

# ## Snapshot series in time

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

# ## Constrained conditions on the barotropic transport (see open_boundary_transport.jl, section 2)

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
