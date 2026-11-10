# Which term of the depth-integrated momentum equation moves water onto or off the shelves. Included by panarctic.jl
# when [diagnostics] tendencies = true; a callback at 0, ½, 1, 2, 3, 4, 6, 9, 12, 18, 24, 36 and 48 hours.
#
# The free surface changes only through the depth-integrated volume equation, ∂η/∂t = −∇·U + Jʷ. Over a region R
# (the shelves, H < 200 m; the slopes, 200–1000 m; the deep ocean, H > 1000 m) the mean η therefore accelerates as
#
#     d²⟨η⟩_R/dt² = −(1/A_R) ∮_R (∂U/∂t)·n dl  (+ dJʷ/dt),
#
# and ∂U/∂t is the sum of the depth-integrated momentum terms, each computed here with the model's own discretization
# and current fields (as in its tendency kernels):
#
#     advection             −Σₖ (U_dot_∇u + its metric term) Δz   (WENO vector invariant; its numerical dissipation is
#                                                                   the model's only horizontal friction: there is no
#                                                                   explicit horizontal viscosity)
#     coriolis              −Σₖ (f × u) Δz
#     baroclinic_pressure   −Σₖ ∂x pHY′ Δz                         (density)
#     barotropic_pressure   −g D ∂x η                              (D the column thickness at the face; in the
#                                                                   barotropic substeps)
#     bottom_drag           −Cd √(u_b² + v_b² + u_bg²) u_b          (Cd and u_bg as in [ocean], default 0.003 and 0;
#                                                                   applied implicitly by the model, estimated here
#                                                                   from the bottom velocity)
#     surface_stress        −(top flux of u)                        (wind and ice stress; zero without wind)
#
# Vertical mixing (CATKE) only moves momentum within the column, so its depth integral is the two stresses above.
# Each term's contribution is printed in cm h⁻² of the region's mean η, with their sum and the region's current
# dη/dt in cm h⁻¹. At t = 0 of a balanced run the pressure and Coriolis terms cancel in divergence by construction, so
# what is left is what the balance did not include; later on, the term that grows is the one breaking the balance.
# Also written to tendency_budget.csv in the output directory.

using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δxᶜᶠᶜ, Δyᶠᶜᶜ, Azᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Advection: U_dot_∇u, U_dot_∇v, U_dot_∇u_hydrostatic_metric, U_dot_∇v_hydrostatic_metric
using Oceananigans.Grids: static_column_depthᶜᶜᵃ

tendency_terms = (:advection, :coriolis, :baroclinic_pressure, :barotropic_pressure, :bottom_drag, :surface_stress)
tendency_file  = joinpath(output_dir, "tendency_budget.csv")
open(tendency_file, "w") do io
    println(io, "hour,region,deta_dt_cm_per_h,", join(("$(t)_cm_per_h2" for t in tendency_terms), ","), ",sum_cm_per_h2")
end

let ocean_settings = get(config, "ocean", Dict())
    global tendency_drag = (Cd  = Float64(get(ocean_settings, "bottom_drag_coefficient", 0.003)),
                            ubg = Float64(get(ocean_settings, "bottom_drag_background_velocity", 0.0)))
end

# Regions by static depth, away from the open edges (whose faces are set by the boundary conditions)
tendency_regions = let cgrid = on_architecture(CPU(), grid)
    H = [static_column_depthᶜᶜᵃ(i, j, cgrid) for i in 1:Nx, j in 1:Ny]
    inside = [3 ≤ j ≤ Ny - 2 for i in 1:Nx, j in 1:Ny]
    (shelf = inside .& (0 .< H .< 200), slope = inside .& (200 .≤ H .< 1000), deep = inside .& (H .≥ 1000))
end

top_flux(field) = (c = field.boundary_conditions.top.condition; c isa Oceananigans.Fields.AbstractField ? Array(interior(c))[:, :, 1] : nothing)

function tendency_budget(sim)
    model = ocean_model_of(sim)
    hour  = (model.clock.time - start_time) / 3600
    cgrid = on_architecture(CPU(), model.grid)          # the current z-star thicknesses
    u     = on_architecture(CPU(), model.velocities.u)
    v     = on_architecture(CPU(), model.velocities.v)
    w     = on_architecture(CPU(), model.velocities.w)
    pHY   = on_architecture(CPU(), model.pressure.pHY′)
    η     = Array(interior(model.free_surface.displacement))[:, :, 1]
    U₃    = (; u, v, w)
    advection = model.advection.momentum
    coriolis  = model.coriolis
    g     = model.free_surface.gravitational_acceleration
    Nzc   = size(cgrid, 3)
    τx, τy = top_flux(model.velocities.u), top_flux(model.velocities.v)
    Cd, ubg = tendency_drag.Cd, tendency_drag.ubg

    Gx = Dict(t => zeros(Nx + 1, Ny) for t in tendency_terms)
    Gy = Dict(t => zeros(Nx, Ny + 1) for t in tendency_terms)
    Ux, Vy = zeros(Nx + 1, Ny), zeros(Nx, Ny + 1)

    for j in 1:Ny, i in 2:Nx
        D, kb = 0.0, 0
        for k in 1:Nzc
            wetcell(i, j, k, cgrid, Face(), Center(), Center()) || continue
            kb == 0 && (kb = k)
            Δz = Δzᶠᶜᶜ(i, j, k, cgrid)
            D += Δz
            Gx[:advection][i, j] -= (U_dot_∇u(i, j, k, cgrid, advection, U₃) +
                                     U_dot_∇u_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃)) * Δz
            Gx[:coriolis][i, j]  -= x_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
            Gx[:baroclinic_pressure][i, j] -= ∂xᶠᶜᶜ(i, j, k, cgrid, pHY) * Δz
            Ux[i, j] += u[i, j, k] * Δz
        end
        D > 0 || continue
        Gx[:barotropic_pressure][i, j] = -g * D * (η[i, j] - η[i-1, j]) / Δxᶠᶜᶜ(i, j, Nzc, cgrid)
        ub = u[i, j, kb]
        vb = (v[i-1, j, kb] + v[i, j, kb] + v[i-1, j+1, kb] + v[i, j+1, kb]) / 4
        Gx[:bottom_drag][i, j] = -Cd * sqrt(ub^2 + vb^2 + ubg^2) * ub
        isnothing(τx) || (Gx[:surface_stress][i, j] = -τx[i, j])
    end
    for j in 3:Ny-1, i in 1:Nx       # interior faces only; the regions stay two rows away from the edges
        D, kb = 0.0, 0
        for k in 1:Nzc
            wetcell(i, j, k, cgrid, Center(), Face(), Center()) || continue
            kb == 0 && (kb = k)
            Δz = Δzᶜᶠᶜ(i, j, k, cgrid)
            D += Δz
            Gy[:advection][i, j] -= (U_dot_∇v(i, j, k, cgrid, advection, U₃) +
                                     U_dot_∇v_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃)) * Δz
            Gy[:coriolis][i, j]  -= y_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
            Gy[:baroclinic_pressure][i, j] -= ∂yᶜᶠᶜ(i, j, k, cgrid, pHY) * Δz
            Vy[i, j] += v[i, j, k] * Δz
        end
        D > 0 || continue
        Gy[:barotropic_pressure][i, j] = -g * D * (η[i, j] - η[i, j-1]) / Δyᶜᶠᶜ(i, j, Nzc, cgrid)
        vb = v[i, j, kb]
        ub = (u[i, j-1, kb] + u[i+1, j-1, kb] + u[i, j, kb] + u[i+1, j, kb]) / 4
        Gy[:bottom_drag][i, j] = -Cd * sqrt(ub^2 + vb^2 + ubg^2) * vb
        isnothing(τy) || (Gy[:surface_stress][i, j] = -τy[i, j])
    end

    # Net inflow of a face field into region R per unit area of R, in m s⁻¹ (for U) or m s⁻² (for ∂U/∂t)
    function region_inflow(Fx, Fy, mask)
        total, area = 0.0, 0.0
        for j in 1:Ny, i in 1:Nx
            mask[i, j] || continue
            total -= Δyᶠᶜᶜ(i+1, j, Nzc, cgrid) * Fx[i+1, j] - Δyᶠᶜᶜ(i, j, Nzc, cgrid) * Fx[i, j] +
                     Δxᶜᶠᶜ(i, j+1, Nzc, cgrid) * Fy[i, j+1] - Δxᶜᶠᶜ(i, j, Nzc, cgrid) * Fy[i, j]
            area  += Azᶜᶜᶜ(i, j, Nzc, cgrid)
        end
        return total / area
    end

    acceleration = 100 * 3600^2      # m s⁻² → cm h⁻²
    @printf("\nDepth-integrated momentum budget at hour %.1f: each term's d²⟨η⟩/dt² (cm/h²), and dη/dt (cm/h)\n", hour)
    @printf("%-7s %9s %10s %9s %11s %11s %9s %9s %9s\n", "region", "dη/dt", "advection", "Coriolis",
            "baroclinic", "barotropic", "drag", "surface", "sum")
    open(tendency_file, "a") do io
        for (name, mask) in pairs(tendency_regions)
            rate = 100 * 3600 * region_inflow(Ux, Vy, mask)
            terms = [acceleration * region_inflow(Gx[t], Gy[t], mask) for t in tendency_terms]
            @printf("%-7s %9.3f %10.3f %9.3f %11.3f %11.3f %9.3f %9.3f %9.3f\n", name, rate, terms..., sum(terms))
            println(io, join((@sprintf("%.3f", hour), name, @sprintf("%.5f", rate),
                              (@sprintf("%.5f", x) for x in terms)..., @sprintf("%.5f", sum(terms))), ","))
        end
    end
    return nothing
end

tendency_hours = filter(h -> h ≤ config["time"]["stop_hours"], [0, 0.5, 1, 2, 3, 4, 6, 9, 12, 18, 24, 36, 48])
add_callback!(simulation, tendency_budget, SpecifiedTimes((start_time .+ tendency_hours .* hours)...))
@info "Momentum budget by depth class at hours $(tendency_hours): $tendency_file"
