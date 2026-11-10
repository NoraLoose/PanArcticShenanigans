# Pan-Arctic coupled ocean–sea ice simulation (test configuration).
# Settings are read from a TOML file: `julia +1.13 --project=. panarctic.jl [config.toml]`. See README.
# Figures and movies: scripts in plotting/ (one per figure).

# Wall-clock time of each phase; rerun.jl runs again from the same initial state without rebuilding (see README)
include("timing.jl")

include("setup.jl")
finish_phase!("setup: config, grid and bathymetry")

using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition, NormalFlowBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ, intrinsic_vector
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: DatasetRestoring, UVelocity, VVelocity
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Printf
using ClimaSeaIce: IncrementalRemapping

# Functions used inside GPU kernels: loaded once per Julia session (again only if kernels.jl changed), so that a second
# run in the same REPL reuses the compiled time step instead of compiling it again
kernels_hash = hash(read(joinpath(@__DIR__, "kernels.jl")))
if !@isdefined(loaded_kernels_hash) || loaded_kernels_hash != kernels_hash
    include("kernels.jl")
    loaded_kernels_hash = kernels_hash
    @info "Loaded kernels.jl: the time step compiles at the first step"
else
    @info "kernels.jl already loaded in this session: reusing the compiled kernels"
end

# ## GLORYS boundary data

# T and S need full inpainting: capped inpainting leaves gaps in deep basins (Baffin Bay overturning, #2).
Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature;  dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity;     dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
ηᵉˣᵗ = glorys_free_surface()

# Datum test ([initialization] eta_offset = "mean" or a number in m): subtract a constant from GLORYS η everywhere and
# at all times, in the initial state and the boundary data. The physics does not depend on the datum (the hybrid and
# clamped edges are immune to it, and the volume target uses only changes in η), but z-star levels do: z = z̃(1 + η/H) + η.
# Any change in the run is then a numerical effect of the absolute η.
eta_offset = get(get(config, "initialization", Dict()), "eta_offset", 0)
if eta_offset != 0
    η_shift = if eta_offset == "mean"
        let offset_grid = on_architecture(CPU(), grid)
            column_area = [Oceananigans.Grids.static_column_depthᶜᶜᵃ(i, j, offset_grid) > 0 ?
                           Oceananigans.Operators.Azᶜᶜᶜ(i, j, 1, offset_grid) : 0.0 for i in 1:Nx, j in 1:Ny]
            sum(Array(interior(ηᵉˣᵗ[start_index]))[:, :, 1] .* column_area) / sum(column_area)
        end
    else
        Float64(eta_offset)
    end
    parent(ηᵉˣᵗ.data) .-= η_shift
    @warn @sprintf("Datum test: GLORYS η shifted by %+.3f m (initial state and boundary data)", -η_shift)
end

# Rotate GLORYS velocities from east/north into the grid frame (grid_frame_u, grid_frame_v in kernels.jl).

for n in eachindex(uᵉˣᵗ.times)
    uₙ, vₙ = uᵉˣᵗ[n], vᵉˣᵗ[n]
    uᵍ = compute!(Field(KernelFunctionOperation{Center, Center, Center}(grid_frame_u, grid, uₙ, vₙ)))
    vᵍ = compute!(Field(KernelFunctionOperation{Center, Center, Center}(grid_frame_v, grid, uₙ, vₙ)))
    interior(uₙ) .= interior(uᵍ)
    interior(vₙ) .= interior(vᵍ)
    fill_halo_regions!((uₙ, vₙ))
end

# Interpolate the rotated velocities to the faces where the model's u and v live.

uᶠ = FieldTimeSeries{Face, Center, Center}(grid, uᵉˣᵗ.times; time_indexing = uᵉˣᵗ.time_indexing)
vᶠ = FieldTimeSeries{Center, Face, Center}(grid, vᵉˣᵗ.times; time_indexing = vᵉˣᵗ.time_indexing)

for n in eachindex(uᵉˣᵗ.times)
    set!(uᶠ[n], compute!(Field(@at((Face, Center, Center), uᵉˣᵗ[n]))))
    set!(vᶠ[n], compute!(Field(@at((Center, Face, Center), vᵉˣᵗ[n]))))
end
@info "Rotated GLORYS velocities into the grid frame and interpolated them to faces"
finish_phase!("GLORYS boundary data: load, rotate, interpolate")

# ## Rest-state test (optional): [initialization] rest = "uniform", "homogeneous" or "glorys"
#
# A check for motion made by the numerics alone (pressure gradients over sloping topography, partial cells, z-star in
# thin columns). The ocean starts at rest (u = v = 0, η = rest_eta, default 0), the edges are held at zero velocity (the
# GLORYS boundary velocities are set to zero and η to rest_eta, so use barotropic = "clamped" with a total of 0 Sv), and
# the run is ocean only: no atmosphere and no sea-ice coupling, so no surface fluxes of heat, salt, freshwater or
# momentum. "uniform" replaces T and S by their mean on each level; "homogeneous" by one value each (the mean over wet
# cells), so the exact solution is rest for any uniform η, and anything that moves comes from the z-star levels
# (z = z̃(1 + η/H) + η tilts them over slopes when η ≠ 0) or the partial cells; "glorys" keeps GLORYS's T and S, so
# whatever moves is the adjustment to its density field.
rest_state = get(get(config, "initialization", Dict()), "rest", "none")
rest_state in ("none", "uniform", "homogeneous", "glorys") ||
    error("[initialization] rest must be \"uniform\", \"homogeneous\" or \"glorys\" (or absent)")
rest_test  = rest_state != "none"
# [initialization] ocean_only = true runs the ocean alone (no atmosphere, no sea-ice coupling) without a rest test:
# GLORYS initial and boundary data as usual, but no surface stress, heat or freshwater fluxes
ocean_only = rest_test || get(get(config, "initialization", Dict()), "ocean_only", false)
ocean_only && !rest_test && @warn "Ocean only: no atmosphere and no sea-ice coupling ([initialization] ocean_only = true)"
if rest_test
    get(config["initialization"], "balance", false) == false || error("[initialization] rest and balance cannot be combined")
    eta_offset == 0 || error("[initialization] rest and eta_offset cannot be combined; use rest_eta")
    rest_eta = Float64(get(config["initialization"], "rest_eta", 0.0))
    for series in (uᵉˣᵗ, vᵉˣᵗ, uᶠ, vᶠ)
        fill!(parent(series.data), 0)
    end
    fill!(parent(ηᵉˣᵗ.data), rest_eta)
    rest_eta == 0 || @warn @sprintf("Rest-state test: uniform η = %+.3f m at the start and in the boundary data", rest_eta)
    if rest_state == "uniform"
        rest_grid = on_architecture(CPU(), grid)
        rest_wet = [!immersed_inactive_node(i, j, k, rest_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
        for (name, series) in (("T", Tᵉˣᵗ), ("S", Sᵉˣᵗ))
            c₀ = Array(interior(series[start_index]))
            profile = [count(rest_wet[:, :, k]) > 0 ? sum(c₀[:, :, k][rest_wet[:, :, k]]) / count(rest_wet[:, :, k]) : 0.0
                       for k in 1:Nz]
            uniform = repeat(reshape(profile, 1, 1, Nz), Nx, Ny, 1)
            for n in eachindex(series.times)
                set!(series[n], uniform)
                fill_halo_regions!(series[n])
            end
            @info @sprintf("Rest-state test: %s uniform on each level, from %.3f (bottom) to %.3f (top)", name, profile[1], profile[end])
        end
    end
    if rest_state == "homogeneous"
        rest_grid = on_architecture(CPU(), grid)
        rest_wet = [!immersed_inactive_node(i, j, k, rest_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
        for (name, series) in (("T", Tᵉˣᵗ), ("S", Sᵉˣᵗ))
            c₀ = Array(interior(series[start_index]))
            value = sum(c₀[rest_wet]) / count(rest_wet)
            for n in eachindex(series.times)
                set!(series[n], value)
                fill_halo_regions!(series[n])
            end
            @info @sprintf("Rest-state test: %s = %.3f everywhere", name, value)
        end
    end
    @warn "Rest-state test ($rest_state): starting at rest, edges at zero velocity, ocean only (no atmosphere or sea ice)"
end

# ## Open boundary conditions
#
# Velocities: oblique radiation (Raymond & Kuo 1984), nudged toward GLORYS on inflow, on both the normal (v)
# and the tangential (u) component. Tracers: reservoirs (as in MOM6) that relax toward GLORYS over a
# length scale on inflow. Barotropic transport: Flather. Timescale and length scales are in the config.

bc = config["boundaries"]
velocity_scheme()  = ObliqueRadiation(inflow_timescale = bc["velocity_inflow_days"] * days)
atlantic_tracer_scheme() = TracerReservoir(inflow_length_scale = bc["tracer_inflow_km_atlantic"] * kilometers)
pacific_tracer_scheme()  = TracerReservoir(inflow_length_scale = bc["tracer_inflow_km_pacific"]  * kilometers)

u_obcs = FieldBoundaryConditions(
    south = ValueBoundaryCondition(atlantic_obc, discrete_form = true, parameters = uᶠ, scheme = velocity_scheme()),
    north = ValueBoundaryCondition(pacific_obc,  discrete_form = true, parameters = uᶠ, scheme = velocity_scheme()))

v_obcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(atlantic_obc,    discrete_form = true, parameters = vᶠ, scheme = velocity_scheme()),
    north = NormalFlowBoundaryCondition(pacific_face_obc, discrete_form = true, parameters = vᶠ, scheme = velocity_scheme()))

T_obcs = FieldBoundaryConditions(
    south = ValueBoundaryCondition(atlantic_obc, discrete_form = true, parameters = Tᵉˣᵗ, scheme = atlantic_tracer_scheme()),
    north = ValueBoundaryCondition(pacific_obc,  discrete_form = true, parameters = Tᵉˣᵗ, scheme = pacific_tracer_scheme()))

S_obcs = FieldBoundaryConditions(
    south = ValueBoundaryCondition(atlantic_obc, discrete_form = true, parameters = Sᵉˣᵗ, scheme = atlantic_tracer_scheme()),
    north = ValueBoundaryCondition(pacific_obc,  discrete_form = true, parameters = Sᵉˣᵗ, scheme = pacific_tracer_scheme()))

V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
    south = GravityWaveRadiationBoundaryCondition(atlantic_V_obc, discrete_form = true, parameters = (v = vᶠ, η = ηᵉˣᵗ)),
    north = GravityWaveRadiationBoundaryCondition(pacific_V_obc,  discrete_form = true, parameters = (v = vᶠ, η = ηᵉˣᵗ)))

# GLORYS's net transport per edge, the hybrid condition (replaces V_obcs when barotropic = "hybrid"), and the hourly
# transport monitor
include("open_boundary_transport.jl")
finish_phase!("open-boundary transports (GLORYS on its own grid)")

# ## Sponge (optional, off by default)
#
# Gaussian rim on the open edges. Tracers relax over 1 day, velocities over 20 minutes.

if bc["sponge"]
    sponge_width = 4   # cells
    sponge_mask  = zeros(Nx, Ny, Nz)
    for j in 1:Ny
        d = min(j - 1, Ny - j)   # cells to the nearest open edge
        sponge_mask[:, j, :] .= exp(-(d / sponge_width)^2)
    end

    FT = DatasetRestoring(Metadata(:temperature; dates, dataset, region, dir_kw...), grid; rate = 1/1days, mask = sponge_mask)
    FS = DatasetRestoring(Metadata(:salinity;    dates, dataset, region, dir_kw...), grid; rate = 1/1days, mask = sponge_mask)
    Fu = DatasetRestoring(uᵉˣᵗ, nothing, on_architecture(arch, sponge_mask), UVelocity(), 1/20minutes)
    Fv = DatasetRestoring(vᵉˣᵗ, nothing, on_architecture(arch, sponge_mask), VVelocity(), 1/20minutes)
    forcing = (T = FT, S = FS, u = Fu, v = Fv)
else
    forcing = NamedTuple()
end

# ## Ocean and sea ice

Δt      = config["time"]["dt_minutes"] * minutes
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)
ocean_settings = get(config, "ocean", Dict())

# Optional [ocean] horizontal_biharmonic_viscosity (m⁴/s): without it, momentum is dissipated only by the WENO
# advection scheme. GLORYS12 uses 2e10 m⁴/s at 1/12° (Lellouche et al. 2018), about 9e10 at 15 km if scaled with Δx³.
# Adding it changes the model type (one compilation); different values then reuse the compiled time step.
biharmonic_viscosity = Float64(get(ocean_settings, "horizontal_biharmonic_viscosity", 0))
closure = if biharmonic_viscosity > 0
    @info "Horizontal biharmonic viscosity $(biharmonic_viscosity) m⁴/s"
    (CATKEVerticalDiffusivity(), HorizontalScalarBiharmonicDiffusivity(ν = biharmonic_viscosity))
else
    CATKEVerticalDiffusivity()
end

# Optional [ocean] bottom_drag_coefficient (NumericalEarth's default is 0.003) and bottom_drag_background_velocity
# (m/s, added in quadrature to the resolved speed; default 0). NEMO's defaults are 0.001 and 0.05 (rn_Cd0 = 1e-3,
# rn_ke0 = 2.5e-3 m² s⁻²).
drag_kw = (; (Symbol(k) => Float64(ocean_settings[k])
              for k in ("bottom_drag_coefficient", "bottom_drag_background_velocity") if haskey(ocean_settings, k))...)
isempty(drag_kw) || @info "Bottom drag settings: $drag_kw"

ocean = ocean_simulation(grid; drag_kw...,
                         clock = Clock(time = start_time),
                         free_surface = SplitExplicitFreeSurface(grid; cfl = 0.7, fixed_Δt = Δt),
                         momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
                         tracer_advection = WENO(; order=7, time_discretization, minimum_buffer_upwind_order=1),
                         closure,
                         forcing,
                         boundary_conditions = (u = u_obcs, v = v_obcs,
                                                T = T_obcs, S = S_obcs,
                                                V = V_obcs))

# EVP dynamics (default). Incremental remapping (Lipscomb & Hunke 2004) transports ice area and volume together, so
# the thickness h = volume / area stays bounded where the concentration is small; it needs forward-Euler steps
# (the NumericalEarth default) and an ice Courant number below one. WENO advected h and ℵ separately and inflated
# h to hundreds of metres in nearly ice-free cells (#29).
sea_ice = sea_ice_simulation(grid, ocean; clock = Clock(time = start_time),
                             advection = IncrementalRemapping(), timestepper = :ForwardEuler)
finish_phase!("ocean and sea-ice models")

# ## Initial conditions: GLORYS (ocean and sea ice)

# uᶠ, vᶠ are already rotated into the grid frame (above), so tell set! not to rotate them again:
# since Oceananigans 0.116 (#6316) it treats velocities as east/north on immersed curvilinear grids too.
n₀ = start_index   # start_date in dates
set!(ocean.model, T = Tᵉˣᵗ[n₀], S = Sᵉˣᵗ[n₀], u = uᶠ[n₀], v = vᶠ[n₀], η = ηᵉˣᵗ[n₀], intrinsic_velocities = true)
set!(sea_ice.model, h = Metadatum(:sea_ice_thickness;     date = start_date, dataset, region, dir_kw...),
                    ℵ = Metadatum(:sea_ice_concentration; date = start_date, dataset, region, dir_kw...))

# GLORYS sea-ice fields are NaN on land and open water; after interpolation make them consistent:
# concentration in [0, 1], and thickness and concentration either both zero or both positive.
h = Array(interior(sea_ice.model.ice_thickness))
ℵ = Array(interior(sea_ice.model.ice_concentration))
@info @sprintf("GLORYS ice: %d / %d NaN cells; max h %.2f m, max ℵ %.3f",
               count(isnan, h), count(isnan, ℵ), maximum(filter(!isnan, h)), maximum(filter(!isnan, ℵ)))
replace!(h, NaN => 0)
replace!(ℵ, NaN => 0)
clamp!(ℵ, 0, 1)
h .= max.(h, 0)
no_ice = (h .== 0) .| (ℵ .== 0)
h[no_ice] .= 0
ℵ[no_ice] .= 0
set!(sea_ice.model, h = h, ℵ = ℵ)

# set! from an array fills only the interior; the halos still hold the NaN from the GLORYS fields.
fill_halo_regions!(sea_ice.model.ice_thickness)
fill_halo_regions!(sea_ice.model.ice_concentration)

# Warn about unstable stratification (N² < -1e-6 s⁻² between wet cells) in the initial state.

cpu_grid = Oceananigans.Architectures.on_architecture(CPU(), grid)
wet = [!immersed_inactive_node(i, j, k, cpu_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
b = Field(Oceananigans.Models.buoyancy_operation(ocean.model))
compute!(b)
N² = Array(interior(compute!(Field(∂z(b)))))
unstable = [(i, j, k) for i in 1:Nx, j in 1:Ny, k in 2:Nz if wet[i, j, k-1] && wet[i, j, k] && N²[i, j, k] < -1e-6]

if isempty(unstable)
    @info "Initial state is gravitationally stable (no N² < -1e-6 s⁻² between wet cells)"
else
    @warn "$(length(unstable)) faces with N² < -1e-6 s⁻² in the initial state; the five most unstable (i, j, k):"
    for (i, j, k) in first(sort(unstable, by = c -> N²[c...]), 5)
        @warn @sprintf("   (%d, %d, %d)   N² = %.2e s⁻²", i, j, k, N²[i, j, k])
    end
end

# Balanced barotropic initial state ([initialization] balance = "projection" or "full")
haskey(config, "initialization") && include("balanced_initialization.jl")

# Barotropic vorticity tendency of the initial state ([diagnostics] initial_vorticity = true, or "only" to stop here)
get(get(config, "diagnostics", Dict()), "initial_vorticity", false) != false && include("initial_vorticity.jl")
finish_phase!("initial state: GLORYS, stability check, balance")

# ## Atmosphere and coupled model

# JRA55 for the run's own dates, so the winds match the GLORYS initial and boundary data ("multi_year"),
# or the 1990–91 repeat year whatever the start date ("repeat_year"). Same dates as GLORYS (including the margin
# days), so its time also counts from 00:00 on dates[1].

if !ocean_only
    jra55 = if config["atmosphere"]["jra55"] == "multi_year"
        (; dataset = MultiYearJRA55(), start_date = dates[1], end_date = dates[end])
    elseif config["atmosphere"]["jra55"] == "repeat_year"
        (; dataset = RepeatYearJRA55())
    else
        error("Unknown JRA55 forcing \"$(config["atmosphere"]["jra55"])\"; use \"multi_year\" or \"repeat_year\"")
    end

    atmosphere    = JRA55PrescribedAtmosphere(arch; jra55..., dir_kw...)

    # Test of the initial adjustment: with `winds = false` in [atmosphere], the JRA55 winds are set to zero (heat and
    # freshwater fluxes still use the other JRA55 fields, with no wind in the bulk formulae). The whole record is held in
    # memory so that nothing is reloaded from file.
    #
    # With `pressure = false`, the sea-level pressure is set to 101325 Pa everywhere. NumericalEarth applies it to the ocean
    # as −∇(p/ρ₀) on every level (an inverse-barometer response of about −1 cm per hPa), so `winds = false` alone still
    # forces the ocean with the January 1993 pressure systems.
    no_winds    = !get(config["atmosphere"], "winds", true)
    no_pressure = !get(config["atmosphere"], "pressure", true)
    if no_winds || no_pressure
        atmosphere = JRA55PrescribedAtmosphere(arch; jra55..., dir_kw...,
                                               time_indices_in_memory = length(atmosphere.velocities.u.times))
    end
    if no_winds
        for u in atmosphere.velocities
            fill!(parent(u.data), 0)
        end
        @warn "JRA55 winds set to zero ([atmosphere] winds = false)"
    end
    if no_pressure
        fill!(parent(atmosphere.pressure.data), 101325)
        @warn "JRA55 sea-level pressure set to 101325 Pa everywhere ([atmosphere] pressure = false)"
    end
    radiation     = JRA55PrescribedRadiation(arch; jra55..., dir_kw...)
    land          = JRA55PrescribedLand(grid; jra55..., dir_kw...)

    for component in (atmosphere, radiation, land)
        component.clock.time = start_time
    end

    coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation, clock = Clock(time = start_time))
end

# ## Simulation

stop_time = start_time + config["time"]["stop_hours"] * hours
simulation = if ocean_only    # the ocean simulation itself, with no coupling
    ocean.Δt = Δt
    ocean.stop_time = stop_time
    ocean
else
    Simulation(coupled_model; Δt, stop_time)
end

finish_phase!("atmosphere and simulation")

# Progress, timing, transport and section monitors, budgets and output writers (also used by rerun.jl)
include("run_outputs.jl")
finish_phase!("monitors and output writers")

# A host copy of the initial state (ocean only), so that rerun.jl can start again from it in this session.
# prognostic_state copies fields to the host; it returns the clock itself (written to file when checkpointing), so the
# clock is copied here. On a CPU the host copy would share memory with the model, hence the deepcopy.
if ocean_only
    initial_state = Oceananigans.prognostic_state(ocean.model)
    arch isa CPU && (initial_state = deepcopy(initial_state))
    let clock = ocean.model.clock
        global initial_state = merge(initial_state, (; clock = (; clock.time, clock.iteration, clock.last_Δt,
                                                                   clock.last_stage_Δt, clock.stage)))
    end
    finish_phase!("copy of the initial state for rerun.jl")
end

start_run_timing!()
run!(simulation)
report_run_timing(simulation)
