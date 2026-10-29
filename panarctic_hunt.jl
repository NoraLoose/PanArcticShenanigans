# Pan-Arctic coupled ocean–sea ice simulation (test configuration).
# Settings are read from a TOML file: `julia +1.13 --project=. panarctic.jl [config.toml]`. See README.

using Pkg
Pkg.activate(".")

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir, GravityWaveRadiationBoundaryCondition, NormalFlowBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ, intrinsic_vector
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: on_architecture
using NumericalEarth.DataWrangling: DatasetRestoring, UVelocity, VVelocity
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Oceananigans.OrthogonalSphericalShellGrids: LambertConformalConicGrid
using Dates, CUDA, Printf
using CopernicusMarine   # enables the GLORYS download extension
using TOML

config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

arch = config["architecture"] == "GPU" ? GPU() : CPU()

# ## Grid
#
# Polar stereographic (LambertConformalConicGrid, λ₀ = −40°). Open boundaries on the Atlantic
# (j = 1) and Pacific (j = Ny) edges; east/west edges are land. `spacing_km` is the largest
# spacing (at the pole; cells shrink ~15% toward the Atlantic edge).

Nz      = config["grid"]["Nz"]
depth   = config["grid"]["depth"] * meters
spacing = config["grid"]["spacing_km"] * kilometers
z       = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

x₁, x₂  = -3280e3, 3616e3   # west–east projected extent (m)
y₁, y₂  = -4100e3, 2700e3   # south–north projected extent (m)
underlying_grid = LambertConformalConicGrid(arch;
                                            size               = (round(Int, (x₂ - x₁) / spacing), round(Int, (y₂ - y₁) / spacing), Nz),
                                            x                  = (x₁, x₂),
                                            y                  = (y₁, y₂),
                                            central_longitude  = -40,   # λ₀: mid-Atlantic
                                            standard_parallel  = 90,
                                            latitude_of_origin = 90,
                                            z,
                                            halo = (7, 7, 7))
φ_south = 43   # southern limit of GLORYS data (grid reaches 43.6°N)

Nx, Ny = size(underlying_grid)[1:2]
@info "Grid: $Nx × $Ny × $Nz, largest spacing $(spacing / 1e3) km"

data_dir   = config["paths"]["data_dir"]
dir_kw     = isempty(data_dir) ? (;) : (; dir = data_dir)
output_dir = mkpath(config["paths"]["output_dir"])
make_figures = config["output"]["figures"]   # bathymetry and grid-spacing PNGs
make_movies  = config["output"]["movies"]    # GLORYS free-surface and results movies

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...)
bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

using CairoMakie

h_bottom = Array(interior(bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN

if make_figures
    fig = Figure(size = (700, 600))
    ax = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
              title = "Arctic bathymetry", aspect = DataAspect())
    hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-depth, 0), nan_color = :gray80)
    Colorbar(fig[1, 2], hm, label = "bottom height [m]")

    save(joinpath(output_dir, "bathymetry.png"), fig)
end

if make_figures
    Δx = Array(interior(compute!(Field(xspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
    Δy = Array(interior(compute!(Field(yspacings(underlying_grid, Center(), Center(), Center()))), :, :, 1)) ./ 1e3   # km
    dry = isnan.(h_bottom)
    Δx[dry] .= NaN
    Δy[dry] .= NaN

    fig = Figure(size = (1000, 450))
    ax1 = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j", title = "Δx [km]", aspect = DataAspect())
    hm1 = heatmap!(ax1, Δx, colormap = :viridis)
    Colorbar(fig[1, 2], hm1)
    ax2 = Axis(fig[1, 3], xlabel = "grid i", title = "Δy [km]", aspect = DataAspect())
    hm2 = heatmap!(ax2, Δy, colormap = :viridis)
    Colorbar(fig[1, 4], hm2)

    save(joinpath(output_dir, "grid_spacing.png"), fig)
end

# ## GLORYS boundary data
#
# Daily means (centred at noon) from the start date to the end of the run.

start_date = DateTime(config["time"]["start_date"])
run_days   = ceil(Int, config["time"]["stop_hours"] / 24)
dates      = start_date : Day(1) : start_date + Day(max(run_days, 1))
dataset = GLORYSDaily()
region  = BoundingBox(longitude = (-180, 180), latitude = (φ_south, 90))

# T and S need full inpainting: capped inpainting leaves gaps in deep basins (Baffin Bay overturning, #2).
Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature;  dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity;     dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
ηᵉˣᵗ = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))

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

if make_movies
    iter = Observable(1)
    ηv = @lift begin
        η = Array(interior(ηᵉˣᵗ[$iter], :, :, 1))
        η[isnan.(h_bottom)] .= NaN
        η
    end
    fig = Figure(size = (700, 600))
    ax  = Axis(fig[1, 1], xlabel = "grid i", ylabel = "grid j",
               title = "GLORYS free surface", aspect = DataAspect())
    hm  = heatmap!(ax, ηv, colormap = :balance, colorrange = (-1.5, 0.2), nan_color = :gray80)
    Colorbar(fig[1, 2], hm, label = "Free surface [m]")
    CairoMakie.record(fig, joinpath(output_dir, "arctic_free_surface.mp4"), 1:length(dates), framerate = 4) do i
        iter[] = i
    end
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

closure = (CATKEVerticalDiffusivity(minimum_tke=1e-7))
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

ocean = ocean_simulation(grid;
                         free_surface = SplitExplicitFreeSurface(grid; substeps=30),
                         momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
                         tracer_advection = WENO(; order=7, time_discretization, minimum_buffer_upwind_order=1),
                         closure,
                         forcing,
                         boundary_conditions = (u = u_obcs, v = v_obcs,
                                                T = T_obcs, S = S_obcs,
                                                V = V_obcs))

sea_ice = sea_ice_simulation(grid, ocean; advection = WENO(order = 5))   # EVP dynamics (default) and ice advection

# ## Initial conditions: GLORYS (ocean and sea ice)

set!(ocean.model, T = Tᵉˣᵗ[1], S = Sᵉˣᵗ[1], u = uᶠ[1], v = vᶠ[1])
set!(sea_ice.model, h = Metadatum(:sea_ice_thickness;     date = dates[1], dataset, region, dir_kw...),
                    ℵ = Metadatum(:sea_ice_concentration; date = dates[1], dataset, region, dir_kw...))

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

# Sea-ice fields are not gap-filled (zeros are real), so coastal cells can pick up NaN from land.
for f in (sea_ice.model.ice_thickness, sea_ice.model.ice_concentration)
    a = Array(interior(f))
    @info @sprintf("%s: %d NaN cells set to 0; max %.2f", summary(f), count(isnan, a), maximum(filter(!isnan, a)))
    replace!(a, NaN => 0)
    set!(f, a)
end

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
# or the 1990–91 repeat year whatever the start date ("repeat_year", e.g. for spin-ups).

jra55 = if config["atmosphere"]["jra55"] == "multi_year"
    (; dataset = MultiYearJRA55(), start_date = dates[1], end_date = dates[end])
elseif config["atmosphere"]["jra55"] == "repeat_year"
    (; dataset = RepeatYearJRA55())
else
    error("Unknown JRA55 forcing \"$(config["atmosphere"]["jra55"])\"; use \"multi_year\" or \"repeat_year\"")
end

atmosphere    = JRA55PrescribedAtmosphere(arch; jra55..., dir_kw...)
radiation     = JRA55PrescribedRadiation(arch; jra55..., dir_kw...)
land          = JRA55PrescribedLand(grid; jra55..., dir_kw...)
coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation)

# ## Simulation

Δt        = config["time"]["dt_minutes"] * minutes
stop_time = config["time"]["stop_hours"] * hours
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

# ## Output: hourly fields near the surface

u, v, w = ocean.model.velocities
h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration
𝒱 = @at((Center, Center, Center), sqrt(u^2 + v^2))
he = h * ℵ
ocean_outputs = merge(ocean.model.tracers, (; 𝒱))

sea_ice_outputs = (; he, uᵢ = sea_ice.model.velocities.u, vᵢ = sea_ice.model.velocities.v)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            filename = joinpath(output_dir, "arctic_ocean_surface.jld2"),
                                            indices = (:, :, grid.Nz-2),
                                            schedule = TimeInterval(1hours),
                                            overwrite_files = true)

sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                              filename = joinpath(output_dir, "arctic_sea_ice_surface.jld2"),
                                              schedule = TimeInterval(1hours),
                                              overwrite_files = true)

# Local debugging, not for committing. Placed just before `run!(simulation)`.
# Takes the first 6 time steps one at a time and reports, after each, which fields hold NaN and where the
# first one is; then stops the script, so the full run does not start.

function report_nans(sim)
    ocean = sim.model.ocean.model
    ice   = sim.model.sea_ice.model
    fields = (T  = ocean.tracers.T,           S  = ocean.tracers.S,
              u  = ocean.velocities.u,        v  = ocean.velocities.v,
              η  = ocean.free_surface.displacement,
              h  = ice.ice_thickness,         ℵ  = ice.ice_concentration,
              uᵢ = ice.velocities.u,          vᵢ = ice.velocities.v)
    ℓ = (Center(), Center(), Center())
    @info "After iteration $(iteration(sim)):"
    for (name, f) in pairs(fields)
        a = Array(interior(f))
        n = count(isnan, a)
        if n == 0
            @info "    $name: no NaN"
        else
            I = Tuple(findfirst(isnan, a))
            i, j = min(I[1], Nx), min(I[2], Ny)
            @info @sprintf("    %s: %d NaN, first at %s ≈ %.1f°E %.1f°N", name, n, string(I),
                           Oceananigans.Grids.λnode(i, j, 1, cpu_grid, ℓ...), Oceananigans.Grids.φnode(i, j, 1, cpu_grid, ℓ...))
        end
    end
    return nothing
end

report_nans(simulation)
for n in 1:6
    simulation.stop_iteration = n
    run!(simulation)
    report_nans(simulation)
end
exit()

run!(simulation)

# ## Movie of SST, SSS, surface speed and ice volume

if make_movies
    To = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "T")
    So = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "S")
    Uo = FieldTimeSeries(joinpath(output_dir, "arctic_ocean_surface.jld2"),   "𝒱")
    hi = FieldTimeSeries(joinpath(output_dir, "arctic_sea_ice_surface.jld2"), "he")

    times = To.times
    n = Observable(length(times))

    title = @lift "Arctic — hour " * string(round(Int, times[$n] / hours))

    Tₙ = @lift(To[$n])
    Sₙ = @lift(So[$n])
    Uₙ = @lift(Uo[$n])
    hₙ = @lift(hi[$n])

    fig = Figure(size = (1100, 1000))
    fig[0, 1:4] = Label(fig, title, fontsize = 20, tellwidth = false)

    ax = Axis(fig[1, 1], ylabel = "grid j", aspect = DataAspect())
    hm_T = heatmap!(ax, Tₙ, colormap = :thermal, colorrange = (-2, 8), nan_color = :gray80)
    ax = Axis(fig[1, 3], aspect = DataAspect())
    hm_S = heatmap!(ax, Sₙ, colormap = :haline, colorrange = (32.5, 35.5), nan_color = :gray80)
    ax = Axis(fig[2, 1], xlabel = "grid i", ylabel = "grid j", aspect = DataAspect())
    hm_U = heatmap!(ax, Uₙ, colormap = Reverse(:solar), colorrange = (0, 0.5), nan_color = :gray80)
    ax = Axis(fig[2, 3], xlabel = "grid i", aspect = DataAspect())
    hm_h = heatmap!(ax, hₙ, colormap = Reverse(:blues), colorrange = (0.01, 1.5), lowclip = :gray80)
    Colorbar(fig[1, 2], hm_T, label = "SST [°C]")
    Colorbar(fig[1, 4], hm_S, label = "SSS [psu]")
    Colorbar(fig[2, 2], hm_U, label = "Surface speed [ms⁻¹]")
    Colorbar(fig[2, 4], hm_h, label = "ice volume per area [m]")

    CairoMakie.record(fig, joinpath(output_dir, "arctic_sea.mp4"), 1:length(times), framerate = 4) do i
        n[] = i
    end
end
