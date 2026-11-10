# Balanced barotropic initial state, by two elliptic solves on the model grid. Included by panarctic.jl after the
# initial conditions are set, when [initialization] balance = "projection" or "full".
#
# GLORYS's daily-mean velocities, interpolated onto our grid, give a depth-integrated transport U whose divergence is
# ~50 m/day at the grid scale (initial_vorticity.jl): an unbalanced part the model removes by gravity waves in its
# first hours, and whose large-scale remainder can excite basin modes. Both solves use the operator
# L φ = ∇·(H ∇φ) on wet columns, with the face weights of the model's own transport (H = wet depth at the face).
# Coasts and the two open edges are no-flux for the correction, so the edge transports are left to the boundary
# conditions.
#
# 1. Projection: L φ = ∇·U − d, then u ← u − ∇φ at every level of every interior face. d is the area mean of ∇·U
#    (−dη/dt of the domain mean from the edge transports), so only the divergence that varies in space is removed.
#    The correction is uniform in depth: the vertical shear and the vorticity of U/H are unchanged.
# 2. ("full") Balanced η: g L η = ∇·(Gᵇᶜ + Gᶜᵒʳ), Gᵇᶜ = −Σ ∂pHY′ Δz (from the initial T and S), Gᶜᵒʳ = −Σ f × u Δz
#    (the corrected u), so the depth-integrated momentum tendency has no divergence at t = 0: no gravity waves are
#    launched. η keeps GLORYS's area mean (its datum at the edges). By default advection, friction and wind are left
#    out; [initialization] balance_terms = ["advection", "drag"] adds the depth-integrated momentum advection (the
#    model's WENO vector-invariant operator, with its numerical dissipation) and the quadratic bottom drag (Cd and u_bg
#    as in [ocean]) to the right-hand side. Wind is zero in the no-wind runs.
#
# Prints how much each step changes, writes balanced_initialization.jld2/.png (with the barotropic streamfunction of
# the corrected transport) in output_dir. With [diagnostics] initial_vorticity on, that diagnostic then measures the
# balanced state.

using Oceananigans.Operators: Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δxᶜᶠᶜ, Δyᶠᶜᶜ, Azᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ, ∂xᶠᶜᶜ, ∂yᶜᶠᶜ
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Advection: U_dot_∇u, U_dot_∇v, U_dot_∇u_hydrostatic_metric, U_dot_∇v_hydrostatic_metric
using Oceananigans.Grids: static_column_depthᶜᶜᵃ, λnode, φnode
using Oceananigans.OutputWriters: jldopen
using CairoMakie, Statistics

# SparseArrays is a dependency of Oceananigans, not of this project: load it by its identity
sparse_arrays = Base.require(Base.PkgId(Base.UUID("2f01184e-e22b-5df5-ae63-d93ebab69eaf"), "SparseArrays"))

balance_mode = get(get(config, "initialization", Dict()), "balance", false)
balance_mode in (false, "projection", "full") ||
    error("[initialization] balance must be \"projection\" or \"full\" (or absent)")
balance_terms = String.(get(get(config, "initialization", Dict()), "balance_terms", String[]))
all(in(("advection", "drag")), balance_terms) ||
    error("[initialization] balance_terms may contain \"advection\" and \"drag\"")
balance_drag = let ocean_settings = get(config, "ocean", Dict())
    (Cd  = Float64(get(ocean_settings, "bottom_drag_coefficient", 0.003)),
     ubg = Float64(get(ocean_settings, "bottom_drag_background_velocity", 0.0)))
end

balance_mode != false && let model = ocean.model
    g     = model.free_surface.gravitational_acceleration
    cgrid = on_architecture(CPU(), model.grid)
    Nzc   = size(cgrid, 3)
    k₀    = Nzc   # level for horizontal metrics (they do not depend on k)

    # ## Geometry: wet depth on faces, wet columns, weights, components

    Hᵘ = zeros(Nx + 1, Ny)   # u-faces; i = 1 and Nx + 1 are land (east/west edges)
    Hᵛ = zeros(Nx, Ny + 1)   # v-faces; j = 1 and Ny + 1 are the open edges
    for j in 1:Ny, i in 2:Nx, k in 1:Nzc
        wetcell(i, j, k, cgrid, Face(), Center(), Center()) && (Hᵘ[i, j] += Δzᶠᶜᶜ(i, j, k, cgrid))
    end
    for j in 1:Ny+1, i in 1:Nx, k in 1:Nzc
        wetcell(i, j, k, cgrid, Center(), Face(), Center()) && (Hᵛ[i, j] += Δzᶜᶠᶜ(i, j, k, cgrid))
    end

    id = zeros(Int, Nx, Ny)
    n = 0
    for j in 1:Ny, i in 1:Nx
        static_column_depthᶜᶜᵃ(i, j, cgrid) > 0 && (id[i, j] = (n += 1))
    end
    area = zeros(n)
    for j in 1:Ny, i in 1:Nx
        id[i, j] > 0 && (area[id[i, j]] = Azᶜᶜᶜ(i, j, k₀, cgrid))
    end

    # Interior faces joining two wet columns: (cell a, cell b, weight Δ⊥ H / Δ∥); b is on the positive side
    links = Tuple{Int, Int, Float64}[]
    for j in 1:Ny, i in 2:Nx
        Hᵘ[i, j] > 0 && id[i-1, j] > 0 && id[i, j] > 0 &&
            push!(links, (id[i-1, j], id[i, j], Δyᶠᶜᶜ(i, j, k₀, cgrid) * Hᵘ[i, j] / Δxᶠᶜᶜ(i, j, k₀, cgrid)))
    end
    for j in 2:Ny, i in 1:Nx
        Hᵛ[i, j] > 0 && id[i, j-1] > 0 && id[i, j] > 0 &&
            push!(links, (id[i, j-1], id[i, j], Δxᶜᶠᶜ(i, j, k₀, cgrid) * Hᵛ[i, j] / Δyᶜᶠᶜ(i, j, k₀, cgrid)))
    end

    # Connected components (union–find); each needs one pinned cell, as the problem is Neumann
    parent_of = collect(1:n)
    root(c) = (while parent_of[c] != c; parent_of[c] = parent_of[parent_of[c]]; c = parent_of[c]; end; c)
    for (a, b, _) in links
        ra, rb = root(a), root(b)
        ra != rb && (parent_of[ra] = rb)
    end
    component = [root(c) for c in 1:n]
    pinned = Set(unique(component))

    # Σ_neighbours w (φₙ − φ_c) = rhs, with one row per component replaced by φ = 0
    I, J, V = Int[], Int[], Float64[]
    for (a, b, w) in links
        for (r, c, x) in ((a, a, -w), (b, b, -w), (a, b, w), (b, a, w))
            r ∈ pinned || (push!(I, r); push!(J, c); push!(V, x))
        end
    end
    for p in pinned
        push!(I, p); push!(J, p); push!(V, 1.0)
    end
    L = sparse_arrays.sparse(I, J, V, n, n)

    function solve(rhs)   # rhs: Az × (source) per cell; made compatible per component, then solved
        total = zeros(n); weight = zeros(n)
        for c in 1:n
            total[component[c]] += rhs[c]; weight[component[c]] += area[c]
        end
        b = [rhs[c] - area[c] * total[component[c]] / weight[component[c]] for c in 1:n]
        for p in pinned
            b[p] = 0
        end
        return L \ b
    end
    field_of(x) = [id[i, j] > 0 ? x[id[i, j]] : NaN for i in 1:Nx, j in 1:Ny]

    # ## Depth-integrated transport and its divergence (Az ∇·U per wet column)

    function transports(u, v)
        Ux, Vy = zeros(Nx + 1, Ny), zeros(Nx, Ny + 1)
        for j in 1:Ny, i in 2:Nx, k in 1:Nzc
            wetcell(i, j, k, cgrid, Face(), Center(), Center()) && (Ux[i, j] += u[i, j, k] * Δzᶠᶜᶜ(i, j, k, cgrid))
        end
        for j in 1:Ny+1, i in 1:Nx, k in 1:Nzc
            wetcell(i, j, k, cgrid, Center(), Face(), Center()) && (Vy[i, j] += v[i, j, k] * Δzᶜᶠᶜ(i, j, k, cgrid))
        end
        return Ux, Vy
    end
    function area_divergence(Fx, Fy)
        d = zeros(n)
        for j in 1:Ny, i in 1:Nx
            c = id[i, j]
            c > 0 || continue
            d[c] = Δyᶠᶜᶜ(i+1, j, k₀, cgrid) * Fx[i+1, j] - Δyᶠᶜᶜ(i, j, k₀, cgrid) * Fx[i, j] +
                   Δxᶜᶠᶜ(i, j+1, k₀, cgrid) * Fy[i, j+1] - Δxᶜᶠᶜ(i, j, k₀, cgrid) * Fy[i, j]
        end
        return d
    end
    cm_per_day(d) = 100 * 86400 .* d ./ area   # −dη/dt from Az ∇·U

    function smooth(a; m = 3)
        s = fill(NaN, size(a))
        for j in axes(a, 2), i in axes(a, 1)
            isfinite(a[i, j]) || continue
            s[i, j] = mean(filter(isfinite, a[max(1, i-m):min(end, i+m), max(1, j-m):min(end, j+m)]))
        end
        return s
    end
    rms(a) = sqrt(mean(abs2, filter(isfinite, a)))

    # ## 1. Projection (as a function: it is repeated after η changes, since on the z-star grid the layer
    # thicknesses, and so the transports, scale with 1 + η/H)

    function project!()
        local u, v, uₐ, vₐ, D₀, D₁, φ, largest, Ux, Vy
        u = on_architecture(CPU(), model.velocities.u)
        v = on_architecture(CPU(), model.velocities.v)
        uₐ, vₐ = Array(interior(u)), Array(interior(v))
        size(uₐ) == (Nx + 1, Ny, Nzc) && size(vₐ) == (Nx, Ny + 1, Nzc) ||
            error("unexpected velocity sizes $(size(uₐ)), $(size(vₐ))")
        D₀ = area_divergence(transports(uₐ, vₐ)...)
        φ  = solve(D₀)
        largest = (0.0, :u, 0, 0)
        for j in 1:Ny, i in 2:Nx
            Hᵘ[i, j] > 0 && id[i-1, j] > 0 && id[i, j] > 0 || continue
            Δu = -(φ[id[i, j]] - φ[id[i-1, j]]) / Δxᶠᶜᶜ(i, j, k₀, cgrid)
            abs(Δu) > largest[1] && (largest = (abs(Δu), :u, i, j))
            for k in 1:Nzc
                wetcell(i, j, k, cgrid, Face(), Center(), Center()) && (uₐ[i, j, k] += Δu)
            end
        end
        for j in 2:Ny, i in 1:Nx
            Hᵛ[i, j] > 0 && id[i, j-1] > 0 && id[i, j] > 0 || continue
            Δv = -(φ[id[i, j]] - φ[id[i, j-1]]) / Δyᶜᶠᶜ(i, j, k₀, cgrid)
            abs(Δv) > largest[1] && (largest = (abs(Δv), :v, i, j))
            for k in 1:Nzc
                wetcell(i, j, k, cgrid, Center(), Face(), Center()) && (vₐ[i, j, k] += Δv)
            end
        end
        Ux, Vy = transports(uₐ, vₐ)
        D₁ = area_divergence(Ux, Vy)
        set!(model; u = uₐ, v = vₐ, intrinsic_velocities = true)
        return D₀, D₁, Vy, largest
    end
    function report_largest(largest)
        Δ, kind, i, j = largest
        H = kind === :u ? Hᵘ[i, j] : Hᵛ[i, j]
        φ_ = φnode(i, j, 1, cgrid, kind === :u ? Face() : Center(), kind === :u ? Center() : Face(), Center())
        λ_ = λnode(i, j, 1, cgrid, kind === :u ? Face() : Center(), kind === :u ? Center() : Face(), Center())
        @printf("   largest velocity change %.2f cm/s, at the %s-face (%d, %d), %.2f°E %.2f°N, wet depth %.0f m\n",
                100Δ, kind, i, j, λ_, φ_, H)
    end

    D₀, D₁, Vy, largest = project!()
    before, after = field_of(cm_per_day(D₀)), field_of(cm_per_day(D₁))
    println("\nBalanced initialization, step 1 (projection):")
    report_largest(largest)
    nordic(a) = a[200:380, 60:240]
    @printf("   initial dη/dt = −∇·U, rms: raw %.1f → %.2f cm/day, 7×7 mean %.1f → %.2f cm/day; Nordic Seas 31×31 mean %.2f → %.3f cm/day\n",
            rms(before), rms(after), rms(smooth(before)), rms(smooth(after)), rms(smooth(nordic(before); m = 15)),
            rms(smooth(nordic(after); m = 15)))
    @printf("   domain-mean dη/dt from the edge transports: %.2f cm/day (kept)\n",
            100 * 86400 * sum(D₁) / sum(area) * -1)

    # ## 2. Balanced η

    ηᴳ = Array(interior(model.free_surface.displacement))[:, :, 1]
    η_new = copy(ηᴳ)
    if balance_mode == "full"
        u = on_architecture(CPU(), model.velocities.u)   # corrected, with halos filled by set!
        v = on_architecture(CPU(), model.velocities.v)
        w = on_architecture(CPU(), model.velocities.w)
        pHY = on_architecture(CPU(), model.pressure.pHY′)
        coriolis = model.coriolis
        advection = model.advection.momentum
        U₃ = (; u, v, w)
        with_advection, with_drag = "advection" in balance_terms, "drag" in balance_terms
        Cd, ubg = balance_drag.Cd, balance_drag.ubg
        Gx, Gy = zeros(Nx + 1, Ny), zeros(Nx, Ny + 1)   # interior faces only; the edges are the boundary's
        for j in 1:Ny, i in 2:Nx
            kb = 0
            for k in 1:Nzc
                wetcell(i, j, k, cgrid, Face(), Center(), Center()) || continue
                kb == 0 && (kb = k)
                G = ∂xᶠᶜᶜ(i, j, k, cgrid, pHY) + x_f_cross_U(i, j, k, cgrid, coriolis, (u, v))
                with_advection && (G += U_dot_∇u(i, j, k, cgrid, advection, U₃) +
                                        U_dot_∇u_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃))
                Gx[i, j] -= G * Δzᶠᶜᶜ(i, j, k, cgrid)
            end
            if with_drag && kb > 0
                ub = u[i, j, kb]
                vb = (v[i-1, j, kb] + v[i, j, kb] + v[i-1, j+1, kb] + v[i, j+1, kb]) / 4
                Gx[i, j] -= Cd * sqrt(ub^2 + vb^2 + ubg^2) * ub
            end
        end
        for j in 2:Ny, i in 1:Nx
            kb = 0
            for k in 1:Nzc
                wetcell(i, j, k, cgrid, Center(), Face(), Center()) || continue
                kb == 0 && (kb = k)
                G = ∂yᶜᶠᶜ(i, j, k, cgrid, pHY) + y_f_cross_U(i, j, k, cgrid, coriolis, (u, v))
                with_advection && (G += U_dot_∇v(i, j, k, cgrid, advection, U₃) +
                                        U_dot_∇v_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃))
                Gy[i, j] -= G * Δzᶜᶠᶜ(i, j, k, cgrid)
            end
            if with_drag && kb > 0
                vb = v[i, j, kb]
                ub = (u[i, j-1, kb] + u[i+1, j-1, kb] + u[i, j, kb] + u[i+1, j, kb]) / 4
                Gy[i, j] -= Cd * sqrt(ub^2 + vb^2 + ubg^2) * vb
            end
        end
        isempty(balance_terms) || @info "Balanced η includes: pressure, Coriolis, $(join(balance_terms, ", "))"
        η̃ = solve(area_divergence(Gx, Gy) ./ g)

        # Each component keeps GLORYS's area mean
        shift_sum = Dict{Int, Float64}(); weight = Dict{Int, Float64}()
        for j in 1:Ny, i in 1:Nx
            c = id[i, j]
            c > 0 || continue
            r = component[c]
            shift_sum[r] = get(shift_sum, r, 0.0) + area[c] * (ηᴳ[i, j] - η̃[c])
            weight[r] = get(weight, r, 0.0) + area[c]
        end
        for j in 1:Ny, i in 1:Nx
            c = id[i, j]
            c > 0 && (η_new[i, j] = η̃[c] + shift_sum[component[c]] / weight[component[c]])
        end

        # Check: divergence of the depth-integrated momentum tendency, before and after
        function tendency_divergence(η)
            Fx, Fy = copy(Gx), copy(Gy)
            for j in 1:Ny, i in 2:Nx
                Hᵘ[i, j] > 0 && id[i-1, j] > 0 && id[i, j] > 0 &&
                    (Fx[i, j] -= g * Hᵘ[i, j] * (η[i, j] - η[i-1, j]) / Δxᶠᶜᶜ(i, j, k₀, cgrid))
            end
            for j in 2:Ny, i in 1:Nx
                Hᵛ[i, j] > 0 && id[i, j-1] > 0 && id[i, j] > 0 &&
                    (Fy[i, j] -= g * Hᵛ[i, j] * (η[i, j] - η[i, j-1]) / Δyᶜᶠᶜ(i, j, k₀, cgrid))
            end
            return field_of(area_divergence(Fx, Fy) ./ area)
        end
        t₀, t₁ = tendency_divergence(ηᴳ), tendency_divergence(η_new)
        change = [id[i, j] > 0 ? η_new[i, j] - ηᴳ[i, j] : NaN for i in 1:Nx, j in 1:Ny]
        set!(model; η = reshape(η_new, Nx, Ny, 1))

        println("Balanced initialization, step 2 (balanced η):")
        @printf("   η − GLORYS η: rms %.1f cm, 7×7 mean %.1f cm, largest %.1f cm\n", 100 * rms(change),
                100 * rms(smooth(change)), 100 * maximum(abs, filter(isfinite, change)))
        @printf("   divergence of the depth-integrated momentum tendency, rms: %.2e → %.2e m s⁻² (7×7 mean %.2e → %.2e)\n",
                rms(t₀), rms(t₁), rms(smooth(t₀)), rms(smooth(t₁)))

        # The new η changes the z-star layer thicknesses, and so the transports: project once more, on a fresh copy
        # of the grid (its thicknesses now follow the new η; the solver's weights keep the old ones, a small effect)
        cgrid = on_architecture(CPU(), model.grid)
        D₀′, D₁′, Vy, largest′ = project!()
        after = field_of(cm_per_day(D₁′))
        println("Balanced initialization, step 1 again (the layer thicknesses changed with η):")
        report_largest(largest′)
        @printf("   initial dη/dt = −∇·U, rms: raw %.1f → %.2f cm/day, 7×7 mean %.1f → %.2f cm/day\n",
                rms(field_of(cm_per_day(D₀′))), rms(after), rms(smooth(field_of(cm_per_day(D₀′)))), rms(smooth(after)))
    end

    # Barotropic streamfunction of the final transport (Sv), at cell corners, from the western land eastward:
    # ψ(i + 1, j) = ψ(i, j) + Δx V. Along j = 1 it is the cumulative transport through the Atlantic edge.
    ψ = zeros(Nx + 1, Ny + 1)
    for j in 1:Ny+1, i in 1:Nx
        ψ[i+1, j] = ψ[i, j] + Δxᶜᶠᶜ(i, j, k₀, cgrid) * Vy[i, j] / 1e6
    end

    # ## Output

    jldopen(joinpath(output_dir, "balanced_initialization.jld2"), "w") do file
        file["mode"] = balance_mode
        file["dηdt_before_cm_per_day"], file["dηdt_after_cm_per_day"] = -before, -after
        file["η_glorys"], file["η_balanced"] = ηᴳ, η_new
        file["ψ_Sv"] = ψ
    end

    irange, jrange = 200:380, 60:240
    window(a) = a[irange, jrange]
    wet2 = id .> 0
    anomaly(a) = (b = ifelse.(wet2, a, NaN); b .- mean(filter(isfinite, window(b))))
    fig = Figure(size = (1700, 1050), fontsize = 13)
    kwη = (colormap = :balance, colorrange = (-40, 40), nan_color = :gray80)
    panels = [("(a) GLORYS η′ (cm)", 100 .* anomaly(ηᴳ), kwη),
              ("(b) balanced η′ (cm)", 100 .* anomaly(η_new), kwη),
              ("(c) balanced − GLORYS η (cm)", 100 .* ifelse.(wet2, η_new .- ηᴳ, NaN),
               (colormap = :balance, colorrange = (-15, 15), nan_color = :gray80)),
              ("(e) dη/dt before, 7×7 mean (cm/day)", -smooth(before),
               (colormap = :balance, colorrange = (-500, 500), nan_color = :gray80)),
              ("(f) dη/dt after, 7×7 mean (cm/day)", -smooth(after),
               (colormap = :balance, colorrange = (-5, 5), nan_color = :gray80))]
    for (m, (title, a, kw)) in enumerate(panels)
        r, c = m ≤ 3 ? (1, 2m - 1) : (2, 2(m - 3) + 1)   # row 2: (d) is drawn at column 1 below
        ax = Axis(fig[r, c], title = title, aspect = DataAspect())
        hm = heatmap!(ax, irange, jrange, window(a); kw...)
        Colorbar(fig[r, c + 1], hm)
    end
    ax = Axis(fig[2, 1], title = "(d) barotropic streamfunction ψ (Sv), contours every 2 Sv", aspect = DataAspect())
    ψw = ψ[irange, jrange]
    ψw = ψw .- ψw[1, 1]
    hm = heatmap!(ax, irange, jrange, ψw; colormap = :balance, colorrange = (-20, 20))
    contour!(ax, irange, jrange, ψw; levels = -60:2:60, color = :black, linewidth = 0.5)
    Colorbar(fig[2, 2], hm)
    save(joinpath(output_dir, "balanced_initialization.png"), fig)
    @info "Wrote balanced_initialization.jld2 and balanced_initialization.png in $output_dir"
end
