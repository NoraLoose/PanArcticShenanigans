# Mass-consistent momentum advection at bottom-cell thickness jumps: Oceananigans.jl PR #6533 (issue #6532), as a
# run-time override until it is released (it is on Oceananigans 0.118, and NumericalEarth is on 0.117).
#
# With partial bottom cells, WENOVectorInvariant lets a jet grow in a thin bottom cell next to a thicker one (deep
# currents of ~1 m/s, in the interior and at the open sides): its kinetic-energy gradient and vertical term do not use
# the transports that carry the u-cell's mass. Where the wet thickness of the bottom cells around a velocity face
# changes by more than 1 %, momentum is advected in flux form (with the scheme's own kinetic-energy-gradient scheme
# and the curvature terms); everywhere else the vector-invariant form is unchanged. Immersed cells are ignored, so a
# whole-cell grid (and every face away from a partial bottom cell) is exactly as before. Only partial-cell grids
# (PCBIBG) are affected.
#
# [grid] partial_cell_fix = true (default) includes this file in panarctic.jl, before the model is built.
# Remove this file once the Oceananigans that NumericalEarth uses contains PR #6533.

@eval Oceananigans.Advection begin
    using Oceananigans.Grids: inactive_node, topology
    using Oceananigans.Operators: Δrᶜᶜᶜ
    using Oceananigans.ImmersedBoundaries: PCBIBG

    @inline mass_consistent_wet_thickness(i, j, k, grid) =
        ifelse(inactive_node(i, j, k, grid, Center(), Center(), Center()), zero(grid), Δrᶜᶜᶜ(i, j, k, grid))

    # immersed cells (thickness 0) are ignored
    @inline function mass_consistent_wet_thickness_changes(a, b, c, d)
        hi = max(a, b, c, d)
        lo = min(ifelse(a == 0, hi, a), ifelse(b == 0, hi, b), ifelse(c == 0, hi, c), ifelse(d == 0, hi, d))
        return 100 * (hi - lo) > hi
    end

    # T-cells i-2 … i+1 enter the transports of the flux-form tendency at u-face i
    @inline mass_consistent_jump_x(i, j, k, grid) =
        mass_consistent_wet_thickness_changes(mass_consistent_wet_thickness(i-2, j, k, grid), mass_consistent_wet_thickness(i-1, j, k, grid),
                                              mass_consistent_wet_thickness(i,   j, k, grid), mass_consistent_wet_thickness(i+1, j, k, grid))

    @inline mass_consistent_jump_y(i, j, k, grid) =
        mass_consistent_wet_thickness_changes(mass_consistent_wet_thickness(i, j-2, k, grid), mass_consistent_wet_thickness(i, j-1, k, grid),
                                              mass_consistent_wet_thickness(i, j,   k, grid), mass_consistent_wet_thickness(i, j+1, k, grid))

    # vector-invariant schemes whose kinetic-energy gradient uses a reconstruction scheme that also has momentum fluxes
    const VectorInvariantWithFluxScheme = VectorInvariant{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any,
                                                          <:Union{AbstractCenteredAdvectionScheme, AbstractUpwindBiasedAdvectionScheme}}

    # `if`, not `ifelse`: only one of the two advection forms is computed at each face
    @inline function U_dot_∇u(i, j, k, grid::PCBIBG, scheme::VectorInvariantWithFluxScheme, U)
        if topology(grid, 1) !== Flat && mass_consistent_jump_x(i, j, k, grid)
            flux_form = scheme.kinetic_energy_gradient_scheme
            return div_𝐯u(i, j, k, grid, flux_form, U, U.u) + U_dot_∇u_hydrostatic_metric(i, j, k, grid, flux_form, U, U)
        else
            return horizontal_advection_U(i, j, k, grid, scheme, U.u, U.v) +
                   vertical_advection_U(i, j, k, grid, scheme, U) +
                   bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v)
        end
    end

    @inline function U_dot_∇v(i, j, k, grid::PCBIBG, scheme::VectorInvariantWithFluxScheme, U)
        if topology(grid, 2) !== Flat && mass_consistent_jump_y(i, j, k, grid)
            flux_form = scheme.kinetic_energy_gradient_scheme
            return div_𝐯v(i, j, k, grid, flux_form, U, U.v) + U_dot_∇v_hydrostatic_metric(i, j, k, grid, flux_form, U, U)
        else
            return horizontal_advection_V(i, j, k, grid, scheme, U.u, U.v) +
                   vertical_advection_V(i, j, k, grid, scheme, U) +
                   bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v)
        end
    end
end

@info "Partial-cell fix on: flux-form momentum advection at bottom-cell thickness jumps (Oceananigans.jl#6533)"
