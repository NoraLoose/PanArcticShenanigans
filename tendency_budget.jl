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
#
# The baroclinic term is also split (second table, tendency_pressure_split.csv). On the z-star grid the model takes the
# pressure gradient at constant z with a chain-rule correction, level by level:
#
#     ∂x pHY′|_z = δx pHY′ / Δx  −  ∂x_z · ℑxz(∂z pHY′)
#                  along_levels     level_slope
#
# with ∂x_z the slope of the level (from znode), zero when η = 0. pHY′ is integrated on the underlying grid through
# immersed cells, whose T and S are masked to zero, and ∂z pHY′ is not masked, so the level_slope part of the bottom cell
# sees that water. The table gives both parts, and the baroclinic term and its level_slope part in the bottom wet cell
# and in the top cell of each face alone (∂z pHY′ at the surface face comes from pHY′'s halo).
#
# The advection term is split too (third table, tendency_advection_split.csv): the model's WENO vector-invariant total,
# the same velocities through a centred, non-dissipative VectorInvariant() (enstrophy-conserving vorticity, energy-
# conserving vertical advection and kinetic-energy gradient), their difference (WENO's numerical dissipation, our only
# lateral friction), and the WENO total's three parts: vorticity flux (ζ + f is not in it; only ζ × u), kinetic-energy
# gradient (with the divergence part of the upwinding) and vertical advection.
#
# Closure check: one model step after each budget time the change in each region's dη/dt is measured and printed
# beside the predicted sum ("measured" is (dη/dt(t + Δt) − dη/dt(t))/Δt).

using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δxᶜᶠᶜ, Δyᶠᶜᶜ, Azᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.Operators: δxᶠᶜᶜ, δyᶜᶠᶜ, Δx⁻¹ᶠᶜᶜ, Δy⁻¹ᶜᶠᶜ
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Advection: U_dot_∇u, U_dot_∇v, U_dot_∇u_hydrostatic_metric, U_dot_∇v_hydrostatic_metric
using Oceananigans.Advection: VectorInvariant, horizontal_advection_U, horizontal_advection_V,
                              vertical_advection_U, vertical_advection_V, bernoulli_head_U, bernoulli_head_V
using Oceananigans.Grids: static_column_depthᶜᶜᵃ

tendency_terms = (:advection, :coriolis, :baroclinic_pressure, :barotropic_pressure, :bottom_drag, :surface_stress)
tendency_file  = joinpath(output_dir, "tendency_budget.csv")
open(tendency_file, "w") do io
    println(io, "hour,region,deta_dt_cm_per_h,", join(("$(t)_cm_per_h2" for t in tendency_terms), ","), ",sum_cm_per_h2")
end

pressure_split_terms = (:along_levels, :level_slope, :bottom_cell, :bottom_cell_level_slope, :top_cell, :top_cell_level_slope)
pressure_split_file  = joinpath(output_dir, "tendency_pressure_split.csv")
open(pressure_split_file, "w") do io
    println(io, "hour,region,", join(("$(t)_cm_per_h2" for t in pressure_split_terms), ","))
end

advection_split_terms = (:weno, :centered, :weno_dissipation, :vorticity_flux, :kinetic_energy_gradient, :vertical_advection)
advection_split_file  = joinpath(output_dir, "tendency_advection_split.csv")
open(advection_split_file, "w") do io
    println(io, "hour,region,", join(("$(t)_cm_per_h2" for t in advection_split_terms), ","))
end
centered_advection = VectorInvariant()

closure_file = joinpath(output_dir, "tendency_closure.csv")
open(closure_file, "w") do io
    println(io, "hour,region,predicted_cm_per_h2,measured_cm_per_h2")
end
predicted_acceleration = Dict{Symbol, Float64}()
rate_before_step       = Dict{Symbol, Float64}()
last_budget_time       = Ref(NaN)

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

# The surface momentum flux J. NumericalEarth's coupled ocean uses an IMEX condition, J = Jₑ + λ u_top, with the
# ice–ocean (and air–ocean) stress split into an explicit flux Jₑ and an implicit coefficient λ; reading only a plain
# Field (as before) gave zero surface stress in coupled runs.
function top_flux(field, top_values)
    c = field.boundary_conditions.top.condition
    c isa Oceananigans.Fields.AbstractField && return Array(interior(c))[:, :, 1]
    if c isa Oceananigans.BoundaryConditions.IMEXFlux
        Jₑ, λ = c.explicit_flux, c.implicit_coefficient
        if Jₑ isa Oceananigans.Fields.AbstractField && λ isa Oceananigans.Fields.AbstractField
            return Array(interior(Jₑ))[:, :, 1] .+ Array(interior(λ))[:, :, 1] .* top_values
        end
        @warn "Surface momentum flux of an unexpected form ($(typeof(Jₑ)), $(typeof(λ))); surface stress left out" maxlog = 1
    end
    return nothing
end

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
    τx = top_flux(model.velocities.u, Array(interior(u))[:, :, size(cgrid, 3)])
    τy = top_flux(model.velocities.v, Array(interior(v))[:, :, size(cgrid, 3)])
    Cd, ubg = tendency_drag.Cd, tendency_drag.ubg

    Gx = Dict(t => zeros(Nx + 1, Ny) for t in tendency_terms)
    Gy = Dict(t => zeros(Nx, Ny + 1) for t in tendency_terms)
    Px = Dict(t => zeros(Nx + 1, Ny) for t in pressure_split_terms)
    Py = Dict(t => zeros(Nx, Ny + 1) for t in pressure_split_terms)
    Ax = Dict(t => zeros(Nx + 1, Ny) for t in advection_split_terms)
    Ay = Dict(t => zeros(Nx, Ny + 1) for t in advection_split_terms)
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
            Ax[:weno][i, j]     -= (U_dot_∇u(i, j, k, cgrid, advection, U₃) +
                                    U_dot_∇u_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃)) * Δz
            Ax[:centered][i, j] -= (U_dot_∇u(i, j, k, cgrid, centered_advection, U₃) +
                                    U_dot_∇u_hydrostatic_metric(i, j, k, cgrid, centered_advection, U₃, U₃)) * Δz
            Ax[:vorticity_flux][i, j]          -= horizontal_advection_U(i, j, k, cgrid, advection, u, v) * Δz
            Ax[:kinetic_energy_gradient][i, j] -= bernoulli_head_U(i, j, k, cgrid, advection, u, v) * Δz
            Ax[:vertical_advection][i, j]      -= vertical_advection_U(i, j, k, cgrid, advection, U₃) * Δz
            Gx[:coriolis][i, j]  -= x_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
            ∂xp = ∂xᶠᶜᶜ(i, j, k, cgrid, pHY)
            ∂xpᴸ = δxᶠᶜᶜ(i, j, k, cgrid, pHY) * Δx⁻¹ᶠᶜᶜ(i, j, k, cgrid)
            Gx[:baroclinic_pressure][i, j] -= ∂xp * Δz
            Px[:along_levels][i, j] -= ∂xpᴸ * Δz
            Px[:level_slope][i, j]  -= (∂xp - ∂xpᴸ) * Δz
            if k == kb
                Px[:bottom_cell][i, j] -= ∂xp * Δz
                Px[:bottom_cell_level_slope][i, j] -= (∂xp - ∂xpᴸ) * Δz
            end
            if k == Nzc
                Px[:top_cell][i, j] -= ∂xp * Δz
                Px[:top_cell_level_slope][i, j] -= (∂xp - ∂xpᴸ) * Δz
            end
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
            Ay[:weno][i, j]     -= (U_dot_∇v(i, j, k, cgrid, advection, U₃) +
                                    U_dot_∇v_hydrostatic_metric(i, j, k, cgrid, advection, U₃, U₃)) * Δz
            Ay[:centered][i, j] -= (U_dot_∇v(i, j, k, cgrid, centered_advection, U₃) +
                                    U_dot_∇v_hydrostatic_metric(i, j, k, cgrid, centered_advection, U₃, U₃)) * Δz
            Ay[:vorticity_flux][i, j]          -= horizontal_advection_V(i, j, k, cgrid, advection, u, v) * Δz
            Ay[:kinetic_energy_gradient][i, j] -= bernoulli_head_V(i, j, k, cgrid, advection, u, v) * Δz
            Ay[:vertical_advection][i, j]      -= vertical_advection_V(i, j, k, cgrid, advection, U₃) * Δz
            Gy[:coriolis][i, j]  -= y_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
            ∂yp = ∂yᶜᶠᶜ(i, j, k, cgrid, pHY)
            ∂ypᴸ = δyᶜᶠᶜ(i, j, k, cgrid, pHY) * Δy⁻¹ᶜᶠᶜ(i, j, k, cgrid)
            Gy[:baroclinic_pressure][i, j] -= ∂yp * Δz
            Py[:along_levels][i, j] -= ∂ypᴸ * Δz
            Py[:level_slope][i, j]  -= (∂yp - ∂ypᴸ) * Δz
            if k == kb
                Py[:bottom_cell][i, j] -= ∂yp * Δz
                Py[:bottom_cell_level_slope][i, j] -= (∂yp - ∂ypᴸ) * Δz
            end
            if k == Nzc
                Py[:top_cell][i, j] -= ∂yp * Δz
                Py[:top_cell_level_slope][i, j] -= (∂yp - ∂ypᴸ) * Δz
            end
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
            predicted_acceleration[name] = sum(terms)
            rate_before_step[name] = rate
            last_budget_time[] = model.clock.time
            println(io, join((@sprintf("%.3f", hour), name, @sprintf("%.5f", rate),
                              (@sprintf("%.5f", x) for x in terms)..., @sprintf("%.5f", sum(terms))), ","))
        end
    end

    @printf("Baroclinic term split (cm/h²): along_levels + level_slope = baroclinic; bottom and top cells alone\n")
    @printf("%-7s %13s %12s %12s %19s %9s %16s\n", "region", "along_levels", "level_slope", "bottom_cell",
            "bottom_level_slope", "top_cell", "top_level_slope")
    open(pressure_split_file, "a") do io
        for (name, mask) in pairs(tendency_regions)
            parts = [acceleration * region_inflow(Px[t], Py[t], mask) for t in pressure_split_terms]
            @printf("%-7s %13.3f %12.3f %12.3f %19.3f %9.3f %16.3f\n", name, parts...)
            println(io, join((@sprintf("%.3f", hour), name, (@sprintf("%.5f", x) for x in parts)...), ","))
        end
    end

    Ax[:weno_dissipation] .= Ax[:weno] .- Ax[:centered]
    Ay[:weno_dissipation] .= Ay[:weno] .- Ay[:centered]
    @printf("Advection split (cm/h²): weno = centered + weno_dissipation; weno = vorticity_flux + ke_gradient + vertical\n")
    @printf("%-7s %9s %9s %17s %15s %12s %9s\n", "region", "weno", "centered", "weno_dissipation", "vorticity_flux",
            "ke_gradient", "vertical")
    open(advection_split_file, "a") do io
        for (name, mask) in pairs(tendency_regions)
            parts = [acceleration * region_inflow(Ax[t], Ay[t], mask) for t in advection_split_terms]
            @printf("%-7s %9.3f %9.3f %17.3f %15.3f %12.3f %9.3f\n", name, parts...)
            println(io, join((@sprintf("%.3f", hour), name, (@sprintf("%.5f", x) for x in parts)...), ","))
        end
    end
    return nothing
end

# One step later: the measured change of each region's dη/dt, against the sum predicted at the budget time
function tendency_closure(sim)
    model = ocean_model_of(sim)
    # Oceananigans calls every callback once at start-up whatever its schedule: act only one step after a budget
    isapprox(model.clock.time, last_budget_time[] + Δt; atol = 1) || return nothing
    hour  = (model.clock.time - start_time) / 3600
    cgrid = on_architecture(CPU(), model.grid)
    u     = on_architecture(CPU(), model.velocities.u)
    v     = on_architecture(CPU(), model.velocities.v)
    Nzc   = size(cgrid, 3)
    Ux, Vy = zeros(Nx + 1, Ny), zeros(Nx, Ny + 1)
    for j in 1:Ny, i in 2:Nx, k in 1:Nzc
        wetcell(i, j, k, cgrid, Face(), Center(), Center()) && (Ux[i, j] += u[i, j, k] * Δzᶠᶜᶜ(i, j, k, cgrid))
    end
    for j in 3:Ny-1, i in 1:Nx, k in 1:Nzc
        wetcell(i, j, k, cgrid, Center(), Face(), Center()) && (Vy[i, j] += v[i, j, k] * Δzᶜᶠᶜ(i, j, k, cgrid))
    end
    @printf("Closure over one step (%.0f s), hour %.2f: predicted sum vs measured d(dη/dt)/dt (cm/h²)\n", Δt, hour)
    open(closure_file, "a") do io
        for (name, mask) in pairs(tendency_regions)
            total, area = 0.0, 0.0
            for j in 1:Ny, i in 1:Nx
                mask[i, j] || continue
                total -= Δyᶠᶜᶜ(i+1, j, Nzc, cgrid) * Ux[i+1, j] - Δyᶠᶜᶜ(i, j, Nzc, cgrid) * Ux[i, j] +
                         Δxᶜᶠᶜ(i, j+1, Nzc, cgrid) * Vy[i, j+1] - Δxᶜᶠᶜ(i, j, Nzc, cgrid) * Vy[i, j]
                area  += Azᶜᶜᶜ(i, j, Nzc, cgrid)
            end
            rate = 100 * 3600 * total / area
            measured = (rate - rate_before_step[name]) / (Δt / 3600)
            @printf("   %-7s predicted %9.3f   measured %9.3f\n", name, predicted_acceleration[name], measured)
            println(io, join((@sprintf("%.3f", hour), name, @sprintf("%.5f", predicted_acceleration[name]),
                              @sprintf("%.5f", measured)), ","))
        end
    end
    return nothing
end

# [diagnostics] tendency_hours overrides the default times, e.g. every 10 minutes in the first hour to check that the
# terms integrate to the change in dη/dt
default_tendency_hours = [0, 0.5, 1, 2, 3, 4, 6, 9, 12, 18, 24, 36, 48]
tendency_hours = filter(h -> h ≤ config["time"]["stop_hours"],
                        Float64.(get(config["diagnostics"], "tendency_hours", default_tendency_hours)))
# Budget times rounded to whole steps, so that SpecifiedTimes does not shorten a step to reach them. The closure callback
# is added first: when a budget time falls one step after another, the closure must read the earlier dη/dt first.
tendency_times = sort(unique(round.(tendency_hours .* hours ./ Δt) .* Δt))
add_callback!(simulation, tendency_closure, SpecifiedTimes((start_time .+ tendency_times .+ Δt)...))
add_callback!(simulation, tendency_budget,  SpecifiedTimes((start_time .+ tendency_times)...))
@info "Momentum budget by depth class at hours $(tendency_times ./ hours): $tendency_file"
