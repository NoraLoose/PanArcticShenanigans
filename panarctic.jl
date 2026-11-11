# Pan-Arctic coupled ocean–sea ice simulation (test configuration).
# Settings are read from a TOML file: `julia +1.13 --project=. panarctic.jl [config.toml]`. See README.
# Figures and movies: scripts in plotting/ (one per figure).

include("setup.jl")

using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Printf
using ClimaSeaIce: IncrementalRemapping

# ## Ocean, nested in GLORYS
#
# nested_ocean_simulation (NumericalEarth.jl#676) reads GLORYS on strips along the open sides, the Atlantic
# (j = 1) and Pacific (j = Ny) edges, and starts the ocean from GLORYS on start_date. On the open sides:
# velocities: oblique radiation (Raymond & Kuo 1984), nudged toward GLORYS on inflow, on both the normal and
# the tangential component. Tracers: reservoirs (as in MOM6) that relax toward GLORYS over a length scale on
# inflow. Barotropic transport: Flather. Near the open sides, the bathymetry is blended toward GLORYS's
# seafloor, and T, S and the baroclinic velocities can relax toward GLORYS in a band (relaxation_hours).
# The child grid is a copy of `grid` with Open sides and the blended bathymetry: use child.grid below.

bc = config["boundaries"]
Δt = config["time"]["dt_minutes"] * minutes
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

ocean = nested_ocean_simulation(grid, dataset;
                                dates = (start_date, end_date),
                                dir_kw...,
                                Δt,
                                open_sides = (:south, :north),
                                normal_velocity_scheme = ObliqueRadiation(inflow_timescale = bc["velocity_inflow_days"] * days),
                                tracer_scheme = TracerReservoir(inflow_length_scale = bc["tracer_inflow_km"] * kilometers),
                                relaxation_rate = bc["relaxation_hours"] > 0 ? 1 / (bc["relaxation_hours"] * hours) : nothing,
                                free_surface = SplitExplicitFreeSurface(grid; cfl = 0.7, fixed_Δt = Δt),
                                momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
                                tracer_advection = WENO(; order=7, time_discretization, minimum_buffer_upwind_order=1),
                                closure = CATKEVerticalDiffusivity())

child = ocean.model.child   # the regional ocean model; ocean.model also holds the GLORYS parent

# ## Sea ice
#
# Ice enters and leaves through the open sides, with GLORYS thickness and concentration on inflow, and starts
# from GLORYS on start_date.
# EVP dynamics (default). Incremental remapping (Lipscomb & Hunke 2004) transports ice area and volume together, so
# the thickness h = volume / area stays bounded where the concentration is small; it needs forward-Euler steps
# (the NumericalEarth default) and an ice Courant number below one. WENO advected h and ℵ separately and inflated
# h to hundreds of metres in nearly ice-free cells (#29).
sea_ice = nested_sea_ice_simulation(ocean, dataset; dates = (start_date, end_date), dir_kw...,
                                    advection = IncrementalRemapping(), timestepper = :ForwardEuler)

# nested_sea_ice_simulation sets the GLORYS NaN (land, open water) to zero. Also make the concentration lie in
# [0, 1], and thickness and concentration either both zero or both positive.
h = Array(interior(sea_ice.model.ice_thickness))
ℵ = Array(interior(sea_ice.model.ice_concentration))
@info @sprintf("GLORYS ice: max h %.2f m, max ℵ %.3f", maximum(h), maximum(ℵ))
clamp!(ℵ, 0, 1)
h .= max.(h, 0)
no_ice = (h .== 0) .| (ℵ .== 0)
h[no_ice] .= 0
ℵ[no_ice] .= 0
set!(sea_ice.model, h = h, ℵ = ℵ)

# Warn about unstable stratification (N² < -1e-6 s⁻² between wet cells) in the initial state.

cpu_grid = on_architecture(CPU(), child.grid)
wet = [!immersed_inactive_node(i, j, k, cpu_grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
b = Field(Oceananigans.Models.buoyancy_operation(child))
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
# or the 1990–91 repeat year whatever the start date ("repeat_year"). Its time, like the model's, counts
# from start_date; one extra day at the end covers the last time step.

jra55 = if config["atmosphere"]["jra55"] == "multi_year"
    (; dataset = MultiYearJRA55(), start_date, end_date = end_date + Day(1))
elseif config["atmosphere"]["jra55"] == "repeat_year"
    (; dataset = RepeatYearJRA55())
else
    error("Unknown JRA55 forcing \"$(config["atmosphere"]["jra55"])\"; use \"multi_year\" or \"repeat_year\"")
end

atmosphere    = JRA55PrescribedAtmosphere(arch; jra55..., dir_kw...)
radiation     = JRA55PrescribedRadiation(arch; jra55..., dir_kw...)
land          = JRA55PrescribedLand(child.grid; jra55..., dir_kw...)

coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation)

# ## Simulation

simulation = Simulation(coupled_model; Δt, stop_time)

wall_time = Ref(time_ns())

function progress(sim)
    T = child.tracers.T
    S = child.tracers.S
    u, v, w = child.velocities
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

# ## Output: time averages and checkpoints
#
# Averages over `average_hours` (daily by default) of the 3D ocean state, the free surface and the sea ice,
# in single precision. Checkpoints of the full model state every `checkpoint_days`, for restarts.

output  = config["output"]
restart = output["restart"]

(; u, v) = child.velocities
η = child.free_surface.displacement
ocean_outputs   = (; T = child.tracers.T, S = child.tracers.S, u, v, η)
sea_ice_outputs = (; h = sea_ice.model.ice_thickness, ℵ = sea_ice.model.ice_concentration,
                     uᵢ = sea_ice.model.velocities.u, vᵢ = sea_ice.model.velocities.v)

# Attached to the coupled simulation, whose run loop accumulates the averages every time step.
# On a restart, append to the existing files instead of overwriting them.
simulation.output_writers[:ocean] = JLD2Writer(child, ocean_outputs;
                                         filename = ocean_output_file,
                                         schedule = AveragedTimeInterval(output["average_hours"] * hours),
                                         overwrite_files = !restart)

simulation.output_writers[:sea_ice] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                                 filename = sea_ice_output_file,
                                                 schedule = AveragedTimeInterval(output["average_hours"] * hours),
                                                 overwrite_files = !restart)

simulation.output_writers[:checkpointer] = Checkpointer(coupled_model;
                                                        schedule = TimeInterval(output["checkpoint_days"] * days),
                                                        dir = output_dir,
                                                        prefix = checkpoint_prefix)

# Continue from the latest checkpoint: it replaces the initial conditions and sets all clocks.
# @allowscalar: Oceananigans restores the open-boundary radiation buffers with an element-wise copy into
# GPU arrays (small boundary rows, so this is cheap).
restart && CUDA.@allowscalar set!(simulation; checkpoint = :latest)

run!(simulation)
