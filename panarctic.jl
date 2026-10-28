# # The Arctic Ocean: a pan-Arctic coupled ocean–sea ice simulation
#
# *A short regional test of a pan-Arctic domain on a rotated grid.*
#
# This script adapts the [Barents Sea regional tutorial](02_barents_sea_regional.jl) to the **whole Arctic
# Ocean**: an Atlantic open boundary reaching down to ~50°N and a domain that also connects to the Pacific
# through the Bering side. The science target is the same — Atlantic/Pacific inflow, an ice edge, the polar
# front — but the geometry forces one structural change: the domain arcs across the high Arctic near the North
# Pole, where a plain latitude–longitude grid is singular.
#
# !!! warning "This is a *test* configuration"
#     The default config runs 0.15° for 10 days (see `config.toml`). The goal is to
#     shake out the grid orientation and the open-boundary plumbing — *not* to produce a validated simulation. Several
#     choices below (the rotation, the box bounds, the velocity boundary conditions) are flagged as things to
#     verify and tune.
#
# !!! warning "Hardware and data requirements"
#     GPU tutorial. The first run downloads GLORYS12 and JRA55 over the **whole Arctic cap** (45–90°N, all
#     longitudes) — substantially more than the Barents box, and possibly not pre-staged in the workshop
#     cache. Keep the date range short (below) on a first run.

using Pkg
Pkg.activate(".")

# Upload packages

DEPOT_PATH

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.BoundaryConditions: NormalRadiation, GravityWaveRadiationBoundaryCondition, NormalFlowBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Dates, CUDA, Printf, Base64
using CopernicusMarine   # enables the GLORYS download extension
using TOML

# ## Settings
#
# Run settings (grid size, dates, run length, paths) live in a TOML file, so nobody has to edit this script to
# change them. Pass the file on the command line, `julia +1.13 --project=. panarctic.jl my_config.toml`;
# without an argument the script reads `config.toml` next to it. Credentials for the downloads stay in
# environment variables (see the README), never in the config file.

config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

arch = config["architecture"] == "GPU" ? GPU() : CPU()

# A small helper to embed an mp4 inline in the notebook (self-contained, plays in JupyterLab):

mp4_html(path) = HTML(string("<video autoplay loop muted playsinline controls ",
                             "src=\"data:video/mp4;base64,", base64encode(read(path)),
                             "\" style=\"max-width:100%\"></video>"))

# ## A rotated pan-Arctic grid
#
# A plain `LatitudeLongitudeGrid` **cannot contain the North Pole**: 90°N is a coordinate singularity and cell
# areas collapse toward it (usable only to ~84°N). We instead use a `RotatedLatitudeLongitudeGrid`, which is
# the *same* structured, `Bounded` lat–lon machinery (so the open-boundary, Flather and sponge code below
# still applies) but with the singular grid pole **moved off the domain**.
#
# We put the grid's north pole on the geographic equator (`north_pole = (λ₀, 0)`): keeping the latitude at 0
# centres the box on the geographic North Pole, and the longitude `λ₀` simply **spins the cap about the pole**.
# We orient it so the two ocean gateways land on the *meridional* edges — the **south** edge on the Atlantic
# (limited to ~50°N) and the **north** edge across the Bering Sea — while the **east/west** edges fall on land.
# `longitude`/`latitude` below are in this **rotated frame**.
#
# !!! tip "Tune the rotation from the diagnostic"
#     With the pole on the equator (`φ_p = 0`) the cap is centred on the North Pole and an edge midpoint sits
#     near geographic latitude `90 − |φ|`. So `λ₀` spins the cap, while `φ₁`/`φ₂` set how far the south/north
#     edges reach: lower `φ₂` to pull the **north** edge in from the open Pacific (e.g. `30` → ~60°N, across
#     the Bering Sea), raise `φ₁` to pull the **south** edge back from the deep Atlantic (e.g. `-30` → ~60°N).
#     Run the grid + bathymetry + edge-coordinate cells below (cheap — no GLORYS) to check after each change.

const λ₁, λ₂ = 151, 212   # rotated longitude (east ↔ west edges → continents)
const φ₁, φ₂ = -34.8, 27  # rotated latitude: south = Atlantic (~50°N) ↔ north (cuts across Bering Sea)

res = config["grid"]["resolution"]   # degrees; 0.15 ≈ TOPAZ2, 0.075 ≈ TOPAZ5

Nx  = round(Int, (λ₂ - λ₁) / res)
Ny  = round(Int, (φ₂ - φ₁) / res)
Nz  = config["grid"]["Nz"]

depth = 4000meters
z = ExponentialDiscretization(Nz, -depth, 0; scale = depth/4, mutable = true)

underlying_grid = RotatedLatitudeLongitudeGrid(arch;
                                               size       = (Nx, Ny, Nz),
                                               longitude  = (λ₁, λ₂),
                                               latitude   = (φ₁, φ₂),
                                               north_pole = (145, 0),   # λ₀: spin the cap (south=Atlantic, north=Bering)
                                               z,
                                               halo = (7, 7, 7))

# Downloads land in `data_dir` from the config, or in each product's default cache when that is empty.
# Output goes to `output_dir`:

data_dir   = config["paths"]["data_dir"]
dir_kw     = isempty(data_dir) ? (;) : (; dir = data_dir)
output_dir = mkpath(config["paths"]["output_dir"])

bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...)
bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)

using CairoMakie, SixelTerm

h_bottom = Array(interior(bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN

fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1], xlabel = "grid i (rotated)", ylabel = "grid j (rotated)",
          title = "Arctic bathymetry (rotated grid)", aspect = DataAspect())
hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-depth, 0), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "bottom height [m]")
fig

# The grid spacing varies across the rotated cap. `xspacings`/`yspacings` at cell centres return the metric
# (in metres) as 2D arrays — heatmap them in km to see where the grid is finest/coarsest and how the rotation
# distorts it:

# (`xspacings`/`yspacings` return a lazy `KernelFunctionOperation`; wrap in a `Field` + `compute!` so we get a
# concrete array we can pull off the GPU — `Array(::KernelFunctionOperation)` would scalar-index on the device.)

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
fig

# ## Open boundary conditions from GLORYS12
#
# Same recipe as the Barents tutorial — Flather (1976) for the barotropic mode, Orlanski/Marchesiello
# radiation for baroclinic velocities and tracers — opened on the two meridional gateways: the **south** edge
# (the Atlantic, down to ~50°N) and the **north** edge (across the Bering Sea). The **east** and **west** edges
# fall on land (the Bering Strait archipelago / continental coasts) and stay closed walls. Both open edges are
# fed by GLORYS12; the data is cropped to the cap.
#
# !!! warning "Velocity components on a rotated grid"
#     GLORYS `u`/`v` are geographic **east/north** components, while the normal-flow and Flather conditions
#     expect velocities along the **rotated grid** axes. On this rotated grid the two frames differ, so the
#     velocity boundary conditions may be mis-oriented until the components are rotated into the grid frame.
#     The tracer (`T`, `S`) and free-surface (`η`) conditions are unaffected (scalars). The north/south normal
#     velocity is the grid-`v`, fed from GLORYS `v_velocity` — treat this as the main thing to validate.

dates   = DateTime(config["time"]["start_date"]) : Day(1) : DateTime(config["time"]["end_date"])
dataset = GLORYSDaily()
region  = BoundingBox(longitude = (-180, 180), latitude = (45, 90))

# T and S use NumericalEarth's default inpainting, which fills every gap in GLORYS. A capped number of passes
# (e.g. `inpainting=100`) leaves zeros below GLORYS's sea floor in enclosed deep basins such as Baffin Bay;
# partial bottom cells are sampled at the centre of the full cell, below the sea floor, and pick those zeros
# up as spuriously fresh, light bottom water that overturns the whole column.
Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature;  dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity;     dates, dataset, region, dir_kw...), grid, time_indices_in_memory=length(dates))
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity;   dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
ηᵉˣᵗ = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kw...), grid, inpainting=100, time_indices_in_memory=length(dates))
nothing #hide

# A quick look at the external free surface (again in grid-index space):

iter = Observable(1)
ηv = @lift begin
    η = Array(interior(ηᵉˣᵗ[$iter], :, :, 1))
    η[isnan.(h_bottom)] .= NaN
    η
end
fig = Figure(size = (700, 600))
ax  = Axis(fig[1, 1], xlabel = "grid i (rotated)", ylabel = "grid j (rotated)",
           title = "GLORYS free surface", aspect = DataAspect())
hm  = heatmap!(ax, ηv, colormap = :balance, colorrange = (-1.5, 0.2), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "Free surface [m]")
CairoMakie.record(fig, joinpath(output_dir, "arctic_free_surface.mp4"), 1:length(dates), framerate = 4) do i
    iter[] = i
end
mp4_html(joinpath(output_dir, "arctic_free_surface.mp4"))

# Discrete boundary functions evaluate the external `FieldTimeSeries` at the boundary index and the current
# clock time. We need the south and north edges; `v` is the normal velocity there:

@inline  south_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1,       k, Oceananigans.Units.Time(clock.time)]
@inline  north_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny, k, Oceananigans.Units.Time(clock.time)]
@inline north_v_obc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny+1, k, Oceananigans.Units.Time(clock.time)]
nothing #hide

# NormalRadiation timescales à la Marchesiello: 1 day on inflow, infinite on outflow. The normal velocity (`v`) and
# the tracers are open on south and north; `u` (normal to the closed east/west walls) keeps its default
# condition:

v_obcs = FieldBoundaryConditions(
    south = NormalFlowBoundaryCondition(south_obc,   discrete_form = true, parameters = vᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)),
    north = NormalFlowBoundaryCondition(north_v_obc, discrete_form = true, parameters = vᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)))

T_obcs = FieldBoundaryConditions(
    south = ValueBoundaryCondition(south_obc, discrete_form = true, parameters = Tᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)),
    north = ValueBoundaryCondition(north_obc, discrete_form = true, parameters = Tᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)))

S_obcs = FieldBoundaryConditions(
    south = ValueBoundaryCondition(south_obc, discrete_form = true, parameters = Sᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)),
    north = ValueBoundaryCondition(north_obc, discrete_form = true, parameters = Sᵉˣᵗ, scheme = NormalRadiation(inflow_timescale = 1days)))

# The Flather condition acts on the barotropic transport `V` with external state `(∫vᵉˣᵗ dz, ηᵉˣᵗ)`,
# `immersed_peripheral_node` skipping solid cells so the integral is the wet-column transport. Only south and
# north (the `V` edges) are open:

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

@inline function south_V_obc(i, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time)
    V = vertical_integral(i, 1, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, 1, 1, t])
end

@inline function north_V_obc(i, k, grid, clock, fields, p)
    j = grid.Ny+1
    t = Oceananigans.Units.Time(clock.time)
    V = vertical_integral(i, j, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, grid.Ny, 1, t])
end

V_obcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
    south = GravityWaveRadiationBoundaryCondition(south_V_obc, discrete_form = true, parameters = (v = vᵉˣᵗ, η = ηᵉˣᵗ)),
    north = GravityWaveRadiationBoundaryCondition(north_V_obc, discrete_form = true, parameters = (v = vᵉˣᵗ, η = ηᵉˣᵗ)))

# ## ... and a sponge behind them
#
# The sponge restores toward GLORYS in a thin rim just inside the open edges. On a rotated grid the geographic
# `(λ, φ)` are no longer aligned with the box edges, so instead of an analytic mask of `(λ, φ)` we build the
# mask in **grid-index space**: a Gaussian rim that is 1 on the open **south/north** edges and decays a few
# cells inward (no rim on the closed east/west walls). `DatasetRestoring` accepts a plain array mask (indexed
# by `i, j, k`):

sponge_width = 4   # cells
sponge_mask  = zeros(Nx, Ny, Nz)
for j in 1:Ny
    d = min(j - 1, Ny - j)                          # distance (in cells) to the nearest open (south/north) edge
    sponge_mask[:, j, :] .= exp(-(d / sponge_width)^2)
end

# Tracers relax on the gentle 1-day timescale; velocities get a much stronger ~20-minute edge nudge so the
# near-boundary interior stays matched to the prescribed boundary:

FT = DatasetRestoring(Metadata(:temperature; dates, dataset, region, dir_kw...), grid; rate = 1/1days,     mask = sponge_mask)
Fu = DatasetRestoring(Metadata(:u_velocity;  dates, dataset, region, dir_kw...), grid; rate = 1/20minutes, mask = sponge_mask, inpainting=100)
Fv = DatasetRestoring(Metadata(:v_velocity;  dates, dataset, region, dir_kw...), grid; rate = 1/20minutes, mask = sponge_mask, inpainting=100)
FS = DatasetRestoring(Metadata(:salinity;    dates, dataset, region, dir_kw...), grid; rate = 1/1days,     mask = sponge_mask)

# ## The ocean component
#
# Assembled by `ocean_simulation` with realistic defaults plus our open boundaries and sponge forcings. Fewer
# free-surface substeps than the Barents run, since the coarse grid has a gentler barotropic CFL:

closure = (CATKEVerticalDiffusivity(minimum_tke=1e-7))
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

ocean = ocean_simulation(grid;
                         free_surface = SplitExplicitFreeSurface(grid; substeps=30),
                         momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
                         tracer_advection = WENO(; order=7, time_discretization, minimum_buffer_upwind_order=1),
                         closure,
                         forcing = (T = FT, S = FS, u = Fu, v = Fv),
                         boundary_conditions = (v = v_obcs,
                                                T = T_obcs, S = S_obcs,
                                                V = V_obcs))

# ## The sea-ice component
#
# Thermodynamics only (`dynamics = nothing`), wired to the ocean below:

sea_ice = sea_ice_simulation(grid, ocean; dynamics=nothing)

# ## Initial conditions
#
# Mid-winter 1993 from GLORYS12 for the ocean, and ECCO4 for the ice thickness and concentration:

set!(ocean.model, T = Tᵉˣᵗ[1], S = Sᵉˣᵗ[1])
set!(sea_ice.model, h = Metadatum(:sea_ice_thickness;     date=dates[1], dataset=ECCO4Monthly(), dir_kw...),
                    ℵ = Metadatum(:sea_ice_concentration; date=dates[1], dataset=ECCO4Monthly(), dir_kw...))

# Check that the initial state is gravitationally stable: N² = ∂b/∂z between vertically adjacent wet cells
# (face k lies between cells k-1 and k). Real data has tiny inversions; anything below -1e-6 s⁻² is suspect.

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

# ## The atmosphere and the coupled model

atmosphere    = JRA55PrescribedAtmosphere(arch; dir_kw...)
radiation     = JRA55PrescribedRadiation(arch; dir_kw...)
land          = JRA55PrescribedLand(grid; dir_kw...)
coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation)

# A short run — a handful of hours — just to see whether the boundaries hold and the fields stay finite.
# The progress line prints the velocity extrema every couple of steps: if the open boundaries are wired
# sanely these stay bounded; a blow-up at the edge is the signature of a misconfigured OBC.

# Time step and run length come from the config (`stop_hours`; e.g. 3 for a quick check):

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

# ## Output
#
# Surface fields from both components, on a short schedule so the few-step run still yields a frame or two:

u, v, w = ocean.model.velocities
h = sea_ice.model.ice_thickness
ℵ = sea_ice.model.ice_concentration
𝒱 = @at((Center, Center, Center), sqrt(u^2 + v^2))
he = h * ℵ
ocean_outputs = merge(ocean.model.tracers, (; 𝒱))

sea_ice_outputs = (; he)

ocean.output_writers[:surface] = JLD2Writer(ocean.model, ocean_outputs;
                                            filename = joinpath(output_dir, "arctic_ocean_surface.jld2"),
                                            indices = (:, :, grid.Nz-2),
                                            schedule = TimeInterval(1hours),
                                            overwrite_files = true)

sea_ice.output_writers[:surface] = JLD2Writer(sea_ice.model, sea_ice_outputs;
                                              filename = joinpath(output_dir, "arctic_sea_ice_surface.jld2"),
                                              schedule = TimeInterval(1hours),
                                              overwrite_files = true)

# The big red button:

run!(simulation)

# ## Results
#
# Four surface fields — SST, SSS, surface speed, and ice volume — in grid-index space. With only a few hours
# of integration this is a sanity check, not a climate story: look for finite, bounded fields and sensible
# behaviour at the open edges.

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
mp4_html(joinpath(output_dir, "arctic_sea.mp4"))
