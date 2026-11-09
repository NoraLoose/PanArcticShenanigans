# Pan-Arctic coupled ocean–sea ice simulation (test configuration).
# Settings are read from a TOML file: `julia +1.13 --project=. panarctic.jl [config.toml]`. See README.
# Figures and movies: scripts in plotting/ (one per figure).

include("setup.jl")

using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition, NormalFlowBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ, intrinsic_vector
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: DatasetRestoring, UVelocity, VVelocity
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Printf
using ClimaSeaIce: IncrementalRemapping

# ## GLORYS boundary data

# T and S need full inpainting: capped inpainting leaves gaps in deep basins (Baffin Bay overturning, #2).
Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature;  dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity;     dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
ηᵉˣᵗ = glorys_free_surface()

# Rotate GLORYS velocities from east/north into the grid frame.

@inline grid_frame_u(i, j, k, grid, uₑ, vₑ) = intrinsic_vector(i, j, k, grid, uₑ, vₑ)[1]
@inline grid_frame_v(i, j, k, grid, uₑ, vₑ) = intrinsic_vector(i, j, k, grid, uₑ, vₑ)[2]

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

# ## Open boundary conditions
#
# Velocities: oblique radiation (Raymond & Kuo 1984), nudged toward GLORYS on inflow, on both the normal (v)
# and the tangential (u) component. Tracers: reservoirs (as in MOM6) that relax toward GLORYS over a
# length scale on inflow. Barotropic transport: Flather. Timescale and length scales are in the config.

bc = config["boundaries"]
velocity_scheme()  = ObliqueRadiation(inflow_timescale = bc["velocity_inflow_days"] * days)
atlantic_tracer_scheme() = TracerReservoir(inflow_length_scale = bc["tracer_inflow_km_atlantic"] * kilometers)
pacific_tracer_scheme()  = TracerReservoir(inflow_length_scale = bc["tracer_inflow_km_pacific"]  * kilometers)

@inline atlantic_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1,          k, Oceananigans.Units.Time(clock.time)]
@inline  pacific_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny,    k, Oceananigans.Units.Time(clock.time)]
@inline pacific_face_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny+1, k, Oceananigans.Units.Time(clock.time)]

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

V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
    south = GravityWaveRadiationBoundaryCondition(atlantic_V_obc, discrete_form = true, parameters = (v = vᶠ, η = ηᵉˣᵗ)),
    north = GravityWaveRadiationBoundaryCondition(pacific_V_obc,  discrete_form = true, parameters = (v = vᶠ, η = ηᵉˣᵗ)))

# GLORYS's net transport per edge, the hybrid condition (replaces V_obcs when barotropic = "hybrid"), and the hourly
# transport monitor
include("open_boundary_transport.jl")

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
closure = CATKEVerticalDiffusivity()
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

ocean = ocean_simulation(grid;
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

# ## Atmosphere and coupled model

# JRA55 for the run's own dates, so the winds match the GLORYS initial and boundary data ("multi_year"),
# or the 1990–91 repeat year whatever the start date ("repeat_year"). Same dates as GLORYS (including the margin
# days), so its time also counts from 00:00 on dates[1].

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
if !get(config["atmosphere"], "winds", true)
    atmosphere = JRA55PrescribedAtmosphere(arch; jra55..., dir_kw...,
                                           time_indices_in_memory = length(atmosphere.velocities.u.times))
    for u in atmosphere.velocities
        fill!(parent(u.data), 0)
    end
    @warn "JRA55 winds set to zero ([atmosphere] winds = false)"
end
radiation     = JRA55PrescribedRadiation(arch; jra55..., dir_kw...)
land          = JRA55PrescribedLand(grid; jra55..., dir_kw...)

for component in (atmosphere, radiation, land)
    component.clock.time = start_time
end

coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation, clock = Clock(time = start_time))

# ## Simulation

stop_time = start_time + config["time"]["stop_hours"] * hours
simulation = Simulation(coupled_model; Δt, stop_time)

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    sea_ice = sim.model.sea_ice
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    h = sea_ice.model.ice_thickness
    msg = @sprintf("time: %s, iter: %d, extrema(T, S): (%.1f, %.1f) °C (%.1f, %.1f) psu, extrema(u) (%.2e, %.2e, %.2e) max(h): %.2f m, wall: %s",
                   prettytime(sim), iteration(sim),
                   extrema(T)..., extrema(S)..., maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(h),
                   prettytime(1e-9 * (time_ns() - wall_time[])))
    @info msg
    wall_time[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, TimeInterval(1hours))

function ice_progress(sim)
    uᵢ, vᵢ = sim.model.sea_ice.model.velocities
    @info @sprintf("    ice: max |u| = %.2f m/s, max |v| = %.2f m/s", maximum(abs, uᵢ), maximum(abs, vᵢ))
    return nothing
end
add_callback!(simulation, ice_progress, TimeInterval(1hours))
add_callback!(simulation, transport_monitor, TimeInterval(1hours))
add_callback!(simulation, update_surface_freshwater!, IterationInterval(1))

# Gateway transports (Greenland–Scotland passages, Svinøy, Lofoten slope, Barents Sea Opening, Fram, Davis and Bering
# Straits), daily means near the Atlantic edge and of the barotropic flow
get(get(config, "diagnostics", Dict()), "sections", false) && include("section_monitor.jl")

# ## Output: hourly fields near the surface

u, v, w = ocean.model.velocities
h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration
𝒱 = @at((Center, Center, Center), sqrt(u^2 + v^2))
he = h * ℵ
ocean_outputs = merge(ocean.model.tracers, (; 𝒱))

sea_ice_outputs = (; he, uᵢ = sea_ice.model.velocities.u, vᵢ = sea_ice.model.velocities.v)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            filename = ocean_surface_file,
                                            indices = (:, :, grid.Nz-2),
                                            schedule = TimeInterval(1hours),
                                            overwrite_files = true)

sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                              filename = sea_ice_surface_file,
                                              schedule = TimeInterval(1hours),
                                              overwrite_files = true)

run!(simulation)
