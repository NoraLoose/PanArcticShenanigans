# GLORYS's daily-mean volume, heat and freshwater transport through the four gateways of gateways.jl, computed on
# GLORYS's own 1/12° grid: no interpolation onto the model grid, no inpainting, GLORYS's own coastline and sea floor.
# This is the reference for the model's gateways.csv (glorys_gateways_run.jl point-samples GLORYS onto the model grid,
# which does not conserve transport at 1°; its numbers are not GLORYS's).
#
#   julia --project=. glorys_gateways_native.jl runs/glorys-sqrtH-4d.toml
#
# Writes, in the run's output_dir,
#   glorys_gateways.csv            through each gateway line (gateway_lines in gateways.jl): the reference;
#                                  plot_gateways.jl draws it as black diamonds
#   glorys_gateways_modelline.csv  through the model's staircase for that line (its cell corners at this resolution),
#                                  still on GLORYS's grid: the difference from the first file is what the model's
#                                  placement of the line alone changes
# and prints a table per day plus the sum of the four (the net inflow to the Arctic, which only a change in its volume,
# runoff and Fury and Hecla Strait can balance: it should be a few tenths of a Sv, not several Sv).
#
# Downloads, once, small boxes around each gateway: daily uo, vo, thetao, so (all depths) and the static deptho.
# Builds the model grid on the CPU only to find the model's staircase; needs no GPU.
#
# Method. GLORYS's daily product is on a regular longitude-latitude grid with every variable at cell centres. Each line
# is turned into a staircase through the corners of GLORYS's cells (as gateways.jl does on the model grid), so every
# step is the face between two GLORYS cells. The velocity on a face is the mean of the two cells' (u on a meridional
# step, v on a zonal one), at every level where both are wet; a face with land on either side carries nothing. Level
# thicknesses come from GLORYS's levels, with the deepest one cut at the shallower of the two columns' deptho (partial
# bottom cells). Heat relative to the gateway's reference temperature, freshwater relative to 34.8 (Shu et al. 2023).

ENV["PANARCTIC_ARCHITECTURE"] = "CPU"
include("setup.jl")
include("gateways.jl")

using NumericalEarth.DataWrangling: BoundingBox, Metadatum, default_download_directory, NCDatasets, Downloads
using NumericalEarth.DataWrangling.GLORYS: GLORYSDaily, GLORYSStatic
using .NCDatasets: Dataset, dimnames

const R_earth = 6.371e6

glorys_dir = isempty(data_dir) ? default_download_directory(GLORYSDaily()) : data_dir

# ## GLORYS in a box around one line

# Box around the points, 0.5° wider on every side (the download adds two more GLORYS cells)
function line_box(points)
    λ, φ = first.(points), last.(points)
    return BoundingBox(longitude = (minimum(λ) - 0.5, maximum(λ) + 0.5), latitude = (minimum(φ) - 0.5, maximum(φ) + 0.5))
end

function read_variable(name, dataset, region; date = nothing)
    meta = isnothing(date) ? Metadatum(name; dataset, region, dir = glorys_dir) :
                             Metadatum(name; dataset, date, region, dir = glorys_dir)
    path = Downloads.download(meta)
    ds   = Dataset(path)
    var  = ds[NumericalEarth.DataWrangling.GLORYS.GLORYS_dataset_variable_names[name]]
    lon  = Float64.(ds["longitude"][:])
    lat  = Float64.(ds["latitude"][:])
    dims = dimnames(var)
    if "depth" in dims
        dims[1:3] == ("longitude", "latitude", "depth") || error("unexpected dimensions $dims in $path")
        data  = Float64.(coalesce.(var[:, :, :, 1], NaN))
        depth = Float64.(ds["depth"][:])
    else
        dims[1:2] == ("longitude", "latitude") || error("unexpected dimensions $dims in $path")
        data  = Float64.(coalesce.(ndims(var) == 3 ? var[:, :, 1] : var[:, :], NaN))
        depth = nothing
    end
    close(ds)
    return (; lon, lat, depth, data)
end

# Thicknesses of GLORYS's levels from its level depths (interfaces halfway between, the top one at the surface)
function level_interfaces(depth)
    zf = vcat(0.0, (depth[1:end-1] .+ depth[2:end]) ./ 2)
    return vcat(zf, depth[end] + (depth[end] - zf[end]))
end

# ## Staircase through GLORYS's cell corners
#
# Corner (I, J) sits at (lon[1] + (I - 3/2) Δλ, lat[1] + (J - 3/2) Δφ): cell (i, j) has corners I = i, i + 1 and
# J = j, j + 1. Steps follow gateways.jl: corner (I, J) → (I + 1, J) crosses the v-face between cells (I, J - 1) and
# (I, J), positive northward; (I, J) → (I, J + 1) crosses the u-face between cells (I - 1, J) and (I, J), with sign −1
# for an eastward u, so a positive transport is always to the left of the direction of travel.

struct LonLatCorners
    λ₀ :: Float64
    φ₀ :: Float64
    Δλ :: Float64
    Δφ :: Float64
    NI :: Int
    NJ :: Int
end

function LonLatCorners(lon, lat)
    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    all(d -> isapprox(d, Δλ; rtol = 1e-3), diff(lon)) && all(d -> isapprox(d, Δφ; rtol = 1e-3), diff(lat)) ||
        error("GLORYS's grid is expected to be regular in longitude and latitude")
    return LonLatCorners(lon[1] - Δλ / 2, lat[1] - Δφ / 2, Δλ, Δφ, length(lon) + 1, length(lat) + 1)
end

corner_lonlat(c::LonLatCorners, I, J) = (c.λ₀ + (I - 1) * c.Δλ, c.φ₀ + (J - 1) * c.Δφ)

function nearest_corner(c::LonLatCorners, λ, φ)
    I = 1 + round(Int, (mod(λ - c.λ₀ + 180, 360) - 180) / c.Δλ)
    J = 1 + round(Int, (φ - c.φ₀) / c.Δφ)
    (1 ≤ I ≤ c.NI && 1 ≤ J ≤ c.NJ) || error(@sprintf("point %.2f°E %.2f°N is outside the downloaded box", λ, φ))
    return (I, J)
end

# Corners visited by the line through `points` (straight in longitude-latitude between points), in unit steps
function lonlat_corner_path(c::LonLatCorners, points)
    path = [nearest_corner(c, points[1]...)]
    for p in 1:length(points) - 1
        (λa, φa), (λb, φb) = points[p], points[p + 1]
        Δλ = mod(λb - λa + 180, 360) - 180
        n  = 4 * ceil(Int, max(abs(Δλ) / c.Δλ, abs(φb - φa) / c.Δφ)) + 1
        for s in 1:n
            α = s / n
            target = nearest_corner(c, λa + α * Δλ, φa + α * (φb - φa))
            while target != path[end]
                I, J = path[end]
                # step along the axis with the larger remaining distance, so the staircase hugs the line
                if abs(target[1] - I) ≥ abs(target[2] - J)
                    push!(path, (I + sign(target[1] - I), J))
                else
                    push!(path, (I, J + sign(target[2] - J)))
                end
            end
        end
    end
    return path
end

# Faces as (kind, a, b, sign): the two cells on either side, `a` on the side the transport counts from
function lonlat_faces(path)
    signs = Dict{Tuple{Symbol, Int, Int}, Int}()
    for (a, b) in zip(path[1:end-1], path[2:end])
        face = b[2] == a[2] ? (:v, min(a[1], b[1]), a[2]) : (:u, a[1], min(a[2], b[2]))
        signs[face] = get(signs, face, 0) + (b[2] == a[2] ? b[1] - a[1] : a[2] - b[2])
    end
    faces = [(f..., s) for (f, s) in signs if s != 0]
    sort!(faces, by = f -> (f[2], f[3]))
    return faces
end

# ## Transport through the faces of one day

struct GatewayBox
    corners :: LonLatCorners
    lat     :: Vector{Float64}
    zf      :: Vector{Float64}     # level interfaces, m, positive down
    deptho  :: Matrix{Float64}     # sea floor depth, m (NaN on land)
end

# (cells on either side, the face's length) for a face; cells are (i, j) of GLORYS's arrays
function face_geometry(box::GatewayBox, kind, I, J)
    c = box.corners
    if kind === :v
        cells = ((I, J - 1), (I, J))
        φ = c.φ₀ + (J - 1) * c.Δφ
        return cells, R_earth * cosd(φ) * deg2rad(c.Δλ)
    else
        cells = ((I - 1, J), (I, J))
        return cells, R_earth * deg2rad(c.Δφ)
    end
end

inside(data, (i, j)) = 1 ≤ i ≤ size(data, 1) && 1 ≤ j ≤ size(data, 2)
wet_floor(h) = isfinite(h) && h > 0   # deptho is missing or 0 on land

function native_transport(box::GatewayBox, faces, Tʳᵉᶠ, u, v, T, S)
    volume = heat = freshwater = 0.0
    Nk = length(box.zf) - 1
    for (kind, I, J, sgn) in faces
        (a, b), L = face_geometry(box, kind, I, J)
        (inside(u, a) && inside(u, b)) || continue
        floor = min(box.deptho[a...], box.deptho[b...])
        wet_floor(floor) || continue
        velocity = kind === :v ? v : u
        for k in 1:Nk
            box.zf[k] ≥ floor && break
            Δz = min(box.zf[k + 1], floor) - box.zf[k]
            va, vb = velocity[a..., k], velocity[b..., k]
            (isfinite(va) && isfinite(vb)) || continue
            q = sgn * (va + vb) / 2 * L * Δz
            volume += q
            Tᶠ = (T[a..., k] + T[b..., k]) / 2
            Sᶠ = (S[a..., k] + S[b..., k]) / 2
            if isfinite(Tᶠ) && isfinite(Sᶠ)
                heat       += ρ₀cₚ * q * (Tᶠ - Tʳᵉᶠ)
                freshwater += q * (1 - Sᶠ / Sʳᵉᶠ)
            end
        end
    end
    return volume, heat, freshwater
end

# ## The model's staircase for each line, as (longitude, latitude) of its corners

model_grid = on_architecture(CPU(), grid)
model_sections = gateway_sections(model_grid)
model_corner_λ = [λnode(i, j, 1, model_grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]
model_corner_φ = [φnode(i, j, 1, model_grid, Face(), Face(), Center()) for i in 1:Nx+1, j in 1:Ny+1]
model_lines = map(gateway_lines) do (name, points, Tʳᵉᶠ)
    path = corner_path(points, model_corner_λ, model_corner_φ)
    [(model_corner_λ[c...], model_corner_φ[c...]) for c in path]
end

# ## Boxes, staircases and checks

run_dates = [start_date + Day(n) for n in 0:ceil(Int, run_days) - 1]

boxes = map(enumerate(gateway_lines)) do (g, (name, points, Tʳᵉᶠ))
    all_points = vcat(points, model_lines[g])
    region = line_box(all_points)
    static = read_variable(:depth, GLORYSStatic(), region)
    sample = read_variable(:u_velocity, GLORYSDaily(), region; date = run_dates[1])
    (static.lon ≈ sample.lon && static.lat ≈ sample.lat) || error("$name: deptho and uo are on different grids")
    box = GatewayBox(LonLatCorners(sample.lon, sample.lat), sample.lat, level_interfaces(sample.depth), static.data)

    lines = map((("line", points), ("model line", model_lines[g]))) do (label, pts)
        path  = lonlat_corner_path(box.corners, pts)
        faces = lonlat_faces(path)
        face_floor(f) = ((a, b), _ = face_geometry(box, f[1], f[2], f[3]); min(box.deptho[a...], box.deptho[b...]))
        wet   = count(f -> wet_floor(face_floor(f)), faces)
        width = sum(f -> face_geometry(box, f[1], f[2], f[3])[2], faces)
        @info @sprintf("GLORYS %-6s %-10s %5d faces (%d wet), %.0f km of staircase", name, label, length(faces), wet, width / 1e3)
        ends = (lonlat_faces(path[1:2]), lonlat_faces(path[end-1:end]))
        for f in Iterators.flatten(ends)
            wet_floor(face_floor(f)) &&
                @warn "GLORYS $name ($label): an end of the line is in water on GLORYS's grid"
        end
        faces
    end
    (; name, region, box, Tʳᵉᶠ, faces = lines[1], model_faces = lines[2])
end

# ## Transports, day by day

files = (joinpath(output_dir, "glorys_gateways.csv"), joinpath(output_dir, "glorys_gateways_modelline.csv"))
for file in files
    open(io -> println(io, gateway_header(model_sections)), file, "w")
end

@printf("\n%-12s", "date")
for b in boxes
    @printf(" %15s", b.name)
end
@printf(" %15s\n", "sum")
@printf("%-12s", "")
for _ in 1:length(boxes) + 1
    @printf(" %7s %7s", "line", "model")
end
println()

for (n, date) in enumerate(run_dates)
    rows = (Float64[24 * (n - 1) + 12], Float64[24 * (n - 1) + 12])   # GLORYS daily means are centred at noon
    for b in boxes
        uᴳ = read_variable(:u_velocity,  GLORYSDaily(), b.region; date).data
        vᴳ = read_variable(:v_velocity,  GLORYSDaily(), b.region; date).data
        Tᴳ = read_variable(:temperature, GLORYSDaily(), b.region; date).data
        Sᴳ = read_variable(:salinity,    GLORYSDaily(), b.region; date).data
        for (row, faces) in zip(rows, (b.faces, b.model_faces))
            Q, H, F = native_transport(b.box, faces, b.Tʳᵉᶠ, uᴳ, vᴳ, Tᴳ, Sᴳ)
            append!(row, (Q / 1e6, H / 1e12, F / 1e3))
        end
    end
    foreach(write_gateway_row, files, rows)
    @printf("%-12s", Date(date))
    for g in 1:length(boxes)
        @printf(" %+7.2f %+7.2f", rows[1][3g - 1], rows[2][3g - 1])
    end
    @printf(" %+7.2f %+7.2f\n", sum(rows[1][3g - 1] for g in 1:length(boxes)), sum(rows[2][3g - 1] for g in 1:length(boxes)))
end
println("\nvolume transport, Sv, positive into the Arctic; line = the gateway line, model = the model's staircase for it")
@info "Wrote $(files[1]) and $(files[2])"
