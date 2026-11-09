# Barotropic vorticity tendency of the initial state: what spins the depth-integrated flow up or down at t = 0.
# Included by panarctic.jl right after the initial conditions are set, when [diagnostics] initial_vorticity = true
# (or "only", which stops panarctic.jl afterwards, before the atmosphere is loaded and before anything compiles for
# the time step).
#
# The depth-integrated momentum tendency, using the model's own discretization and fields at t = 0 (no wind, no
# advection, no friction):
#
#     ∂U/∂t = Gᵇᶜ + Gᵑ + Gᶜᵒʳ,   Gᵇᶜ = −Σₖ ∂x pHY′ Δz,   Gᵑ = −g H ∂x η,   Gᶜᵒʳ = −Σₖ (f × u) Δz
#
# pHY′ is the model's hydrostatic pressure anomaly from the initial T and S, H the column depth at the face. The curl
# of Gᵇᶜ + Gᵑ is the bottom pressure torque J(p_b, H) (its baroclinic part is JEBAR); the curl of Gᶜᵒʳ is −βV − f∇·U.
# A state in balance with our topography has a small total; a state that will drive barotropic flow does not.
# Also: the initial η tendency −∇·U, from the GLORYS velocities on our grid and bathymetry.
#
# Writes initial_vorticity.jld2 (whole domain) and initial_vorticity.png (Nordic Seas) in output_dir and prints rms
# values; compare runs with different bathymetry or levels by these numbers and maps.

using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δxᶜᶠᶜ, Δyᶠᶜᶜ, Azᶠᶠᶜ, Azᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.Coriolis: x_f_cross_U, y_f_cross_U
using Oceananigans.Grids: static_column_depthᶜᶜᵃ
using Oceananigans.OutputWriters: jldopen
using CairoMakie, Statistics

let model = ocean.model
    g        = model.free_surface.gravitational_acceleration
    cgrid    = on_architecture(CPU(), model.grid)
    pHY      = on_architecture(CPU(), model.pressure.pHY′)
    u        = on_architecture(CPU(), model.velocities.u)
    v        = on_architecture(CPU(), model.velocities.v)
    η        = Array(interior(model.free_surface.displacement))[:, :, 1]
    coriolis = model.coriolis   # parameters only, no arrays
    U₃       = (u, v)
    Nzc      = size(cgrid, 3)

    # Depth-integrated forces at u-faces (i, j) and v-faces (i, j), and the depth-integrated flow
    Gxᵇᶜ, Gxᵑ, Gxᶜᵒʳ, Ux = (zeros(Nx + 1, Ny) for _ in 1:4)
    Gyᵇᶜ, Gyᵑ, Gyᶜᵒʳ, Vy = (zeros(Nx, Ny + 1) for _ in 1:4)
    for j in 1:Ny, i in 2:Nx
        H = 0.0
        for k in 1:Nzc
            wetcell(i, j, k, cgrid, Face(), Center(), Center()) || continue
            Δz = Δzᶠᶜᶜ(i, j, k, cgrid)
            H += Δz
            Gxᵇᶜ[i, j]  -= ∂xᶠᶜᶜ(i, j, k, cgrid, pHY) * Δz
            Gxᶜᵒʳ[i, j] -= x_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
            Ux[i, j]    += u[i, j, k] * Δz
        end
        H > 0 && (Gxᵑ[i, j] = -g * H * (η[i, j] - η[i-1, j]) / Δxᶠᶜᶜ(i, j, Nzc, cgrid))
    end
    for j in 1:Ny+1, i in 1:Nx
        interior_face = 2 ≤ j ≤ Ny   # forces on the open edges (j = 1, Ny + 1) are set by the boundary condition
        H = 0.0
        for k in 1:Nzc
            wetcell(i, j, k, cgrid, Center(), Face(), Center()) || continue
            Δz = Δzᶜᶠᶜ(i, j, k, cgrid)
            H += Δz
            Vy[i, j] += v[i, j, k] * Δz
            interior_face || continue
            Gyᵇᶜ[i, j]  -= ∂yᶜᶠᶜ(i, j, k, cgrid, pHY) * Δz
            Gyᶜᵒʳ[i, j] -= y_f_cross_U(i, j, k, cgrid, coriolis, U₃) * Δz
        end
        interior_face && H > 0 && (Gyᵑ[i, j] = -g * H * (η[i, j] - η[i, j-1]) / Δyᶜᶠᶜ(i, j, Nzc, cgrid))
    end

    # Curl at corners (i, j), as Oceananigans' ζ₃ᶠᶠᶜ; NaN unless all four faces around the corner are wet
    xwet = [i > 1 && i ≤ Nx && wetcell(i, j, Nzc, cgrid, Face(), Center(), Center()) for i in 1:Nx+1, j in 1:Ny]
    ywet = [j > 1 && j ≤ Ny && wetcell(i, j, Nzc, cgrid, Center(), Face(), Center()) for i in 1:Nx, j in 1:Ny+1]
    function curl(Gx, Gy)
        c = fill(NaN, Nx, Ny)
        for j in 2:Ny, i in 2:Nx
            xwet[i, j] && xwet[i, j-1] && ywet[i, j] && ywet[i-1, j] || continue
            c[i, j] = (Δyᶜᶠᶜ(i, j, Nzc, cgrid) * Gy[i, j] - Δyᶜᶠᶜ(i-1, j, Nzc, cgrid) * Gy[i-1, j] -
                       Δxᶠᶜᶜ(i, j, Nzc, cgrid) * Gx[i, j] + Δxᶠᶜᶜ(i, j-1, Nzc, cgrid) * Gx[i, j-1]) / Azᶠᶠᶜ(i, j, Nzc, cgrid)
        end
        return c
    end
    jebar   = curl(Gxᵇᶜ, Gyᵇᶜ)                       # baroclinic bottom pressure torque
    ηtorque = curl(Gxᵑ, Gyᵑ)                         # g J(η, H)
    bpt     = jebar .+ ηtorque                       # total bottom pressure torque
    cor     = curl(Gxᶜᵒʳ, Gyᶜᵒʳ)                     # −βV − f∇·U
    total   = bpt .+ cor

    # Initial η tendency −∇·U (m/s), at cell centres
    wetc = [static_column_depthᶜᶜᵃ(i, j, cgrid) > 0 for i in 1:Nx, j in 1:Ny]
    dηdt = fill(NaN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        wetc[i, j] || continue
        dηdt[i, j] = -(Δyᶠᶜᶜ(i+1, j, Nzc, cgrid) * Ux[i+1, j] - Δyᶠᶜᶜ(i, j, Nzc, cgrid) * Ux[i, j] +
                       Δxᶜᶠᶜ(i, j+1, Nzc, cgrid) * Vy[i, j+1] - Δxᶜᶠᶜ(i, j, Nzc, cgrid) * Vy[i, j]) / Azᶜᶜᶜ(i, j, Nzc, cgrid)
    end
    depth = [wetc[i, j] ? static_column_depthᶜᶜᵃ(i, j, cgrid) : NaN for i in 1:Nx, j in 1:Ny]

    # Box average over (2n + 1)² points, ignoring NaN: the torques are noisy at the grid scale, and what drives the
    # basin-scale response is their large-scale part
    function smooth(a; n = 3)
        s = fill(NaN, size(a))
        for j in axes(a, 2), i in axes(a, 1)
            isfinite(a[i, j]) || continue
            block = filter(isfinite, a[max(1, i-n):min(end, i+n), max(1, j-n):min(end, j+n)])
            s[i, j] = mean(block)
        end
        return s
    end
    fields = (; jebar, ηtorque, bpt, cor, total)
    smoothed = map(smooth, fields)

    irange, jrange = 200:380, 60:240
    rms(a) = sqrt(mean(abs2, filter(isfinite, a)))
    window(a) = a[irange, jrange]
    scale = 1e-9
    println("\nInitial barotropic vorticity tendency [1e-9 m s⁻²], rms (raw / 7×7 box mean), bathymetry \"$bathymetry_source\", vertical \"$vertical_levels\" ($Nzc levels):")
    @printf("%-28s %20s %20s\n", "", "whole domain", "Nordic Seas window")
    rows = ("baroclinic BPT (JEBAR)" => :jebar, "g J(η, H)" => :ηtorque, "bottom pressure torque" => :bpt,
             "Coriolis: −βV − f∇·U" => :cor, "total" => :total)
    for (label, key) in rows
        a, s = getproperty(fields, key), getproperty(smoothed, key)
        @printf("%-28s %9.2f / %8.2f %9.2f / %8.2f\n", label, rms(a) / scale, rms(s) / scale,
                rms(window(a)) / scale, rms(window(s)) / scale)
    end
    @printf("%-28s %9.2f / %8.2f %9.2f / %8.2f   (smaller is closer to balance)\n", "total / bottom pressure torque",
            rms(total) / rms(bpt), rms(smoothed.total) / rms(smoothed.bpt),
            rms(window(total)) / rms(window(bpt)), rms(window(smoothed.total)) / rms(window(smoothed.bpt)))
    @printf("Initial η tendency −∇·U: rms %.1f cm/day (whole domain), %.1f cm/day (Nordic Seas); 7×7 mean %.1f / %.1f\n",
            100 * 86400 * rms(dηdt), 100 * 86400 * rms(window(dηdt)),
            100 * 86400 * rms(smooth(dηdt)), 100 * 86400 * rms(window(smooth(dηdt))))

    jldopen(joinpath(output_dir, "initial_vorticity.jld2"), "w") do file
        for (name, a) in pairs(fields)
            file[String(name)] = a
        end
        file["dηdt"], file["depth"] = dηdt, depth
        file["bathymetry"], file["vertical"] = bathymetry_source, vertical_levels
    end

    # Figure: Nordic Seas window
    fig = Figure(size = (1700, 1050), fontsize = 13)
    limit = 3 * rms(window(smoothed.total)) / scale
    torque_kw = (colormap = :balance, colorrange = (-limit, limit), nan_color = :gray80)
    panels = [("(a) bottom pressure torque", smoothed.bpt), ("(b) Coriolis: −βV − f∇·U", smoothed.cor),
              ("(c) total = (a) + (b)", smoothed.total), ("(d) baroclinic part of (a) (JEBAR)", smoothed.jebar)]
    for (n, (title, a)) in enumerate(panels)
        ax = Axis(fig[(n - 1) ÷ 2 + 1, 2 * ((n - 1) % 2) + 1], title = title * ", 7×7 mean (1e-9 m s⁻²)",
                  aspect = DataAspect())
        hm = heatmap!(ax, irange, jrange, window(a) ./ scale; torque_kw...)
        contour!(ax, irange, jrange, window(depth); levels = [500, 1000, 2000, 3000], color = :black, linewidth = 0.4)
        Colorbar(fig[(n - 1) ÷ 2 + 1, 2 * ((n - 1) % 2) + 2], hm)
    end
    ax = Axis(fig[1, 5], title = "(e) initial dη/dt = −∇·U (cm/day)", aspect = DataAspect())
    hm = heatmap!(ax, irange, jrange, 100 * 86400 .* window(dηdt); colormap = :balance, colorrange = (-20, 20),
                  nan_color = :gray80)
    Colorbar(fig[1, 6], hm)
    ax = Axis(fig[2, 5], title = "(f) depth (m), \"$bathymetry_source\", $Nzc levels", aspect = DataAspect())
    hm = heatmap!(ax, irange, jrange, window(depth); colormap = :deep, colorrange = (0, 4000), nan_color = :gray80)
    Colorbar(fig[2, 6], hm)
    save(joinpath(output_dir, "initial_vorticity.png"), fig)
    @info "Wrote initial_vorticity.jld2 and initial_vorticity.png in $output_dir"
end

if get(get(config, "diagnostics", Dict()), "initial_vorticity", false) == "only"
    error("Stopped on purpose after the initial vorticity diagnostic ([diagnostics] initial_vorticity = \"only\"); " *
          "this is not a failure.")
end
