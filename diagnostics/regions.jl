# Ocean regions of the model domain, for regional averages and integrals (heat content, profiles, ...), in sets:
#
#     include("diagnostics/regions.jl")
#     regions = region_masks(grid; set = "PAN")   # (; set, codes, names, masks): masks[code] is a Bool matrix (Nx, Ny),
#                                                 # true in the region's ocean columns; codes in a fixed order
#
#   "PAN"   15 sub-basins that cover the whole domain without overlapping, adapted from the PAN regions of
#           hycom-validation (regions_pan.py); see below.
#   "TSAR"  the 4 boxes of Muilwijk et al. (2022) where the observed profiles of TS_Arctic_Regions_Observed_1970_2017.nc
#           are (Polyakov et al., 2020), as in hycom-validation (regions_tsar.py): parts of the domain only.
#
# The regions are defined by longitude, latitude and bottom depth only, not by grid indices, so they work for any grid
# (any spacing, and GLORYS's or WOA's own grids). Boxes are (latitude range, longitude range), with longitude in
# (-180, 180]: lat₀ ≤ φ < lat₁ and lon₀ < λ ≤ lon₁.

using Oceananigans.Grids: λnode, φnode, topology, Periodic
using Oceananigans.Architectures: on_architecture
using TOML

# The sections of the transports (sections.toml) bound some of the regions, so regions and transports always agree
region_section(name) = TOML.parsefile(joinpath(@__DIR__, "..", "sections.toml"))[name]

# ## PAN
#
# Assigned in the order of `pan_codes`: the Arctic basins first (see below), then each region claims the ocean columns
# in its area that no earlier region has claimed, so no two regions overlap. The areas leave no part of the domain out,
# so together the regions cover every ocean column (`pan_masks` checks both).

const pan_names = Dict(
    "EUBA" => "Eurasian Basin",
    "AMBA" => "Amerasian Basin",
    "BRTS" => "Barents Sea",
    "KARA" => "Kara Sea",
    "LPES" => "Laptev + East Siberian Seas",
    "CHSH" => "Chukchi + Beaufort Shelf",
    "CACH" => "Canadian Arctic Archipelago",
    "BFFB" => "Baffin Bay",
    "HUDS" => "Hudson Bay + Strait",
    "NRDS" => "Nordic Seas",
    "NRTS" => "North Sea + English Channel",
    "BALT" => "Baltic Sea",
    "LBRD" => "Labrador Sea",
    "BRGS" => "Bering Sea",
    "ESPNA" => "Subpolar North Atlantic",
)

# Order of assignment: Arctic basins, Arctic shelf seas, Nordic Seas, European marginal seas, sub-Arctic.
const pan_codes = ["EUBA", "AMBA", "BRTS", "KARA", "LPES", "CHSH", "CACH", "BFFB", "HUDS",
                   "NRDS", "NRTS", "BALT", "LBRD", "BRGS", "ESPNA"]

in_box(λ, φ, ((lat₀, lat₁), (lon₀, lon₁))) = (lat₀ ≤ φ < lat₁) & (lon₀ < λ ≤ lon₁)
in_boxes(λ, φ, boxes) = any(in_box(λ, φ, box) for box in boxes)

# ## Arctic basins
#
# EUBA and AMBA are the columns deeper than deep_basin_depth inside `arctic_ocean_boxes` that are connected (through such
# columns) to the North Pole, with any shallower patches they enclose (rises on the ridges), and everything north of
# Fram Strait (`fram_strait_box`, north of the latitude of the fram_strait section in sections.toml), so that the Fram
# Strait section is the boundary of the Arctic basins against the Nordic Seas. South of 76°N the boxes cut the troughs into the Barents and
# Kara Seas (the St. Anna Trough reaches ~77°N); shallow sills elsewhere (Nares, Barrow and Bering Straits) keep out
# Baffin Bay, the Archipelago and the Pacific. The Lomonosov Ridge (~140°E / ~25°W, as in regions_pan.py) splits them:
# EUBA is the part in `eurasian_basin_boxes`.

const deep_basin_depth = 500   # minimum bottom depth (m) of the Arctic basins EUBA and AMBA
const fram_strait_latitude = only(unique(region_section("fram_strait")["lat"]))   # 78°50′N
const fram_strait_box  = ((fram_strait_latitude, 90), (-45, 16))
const arctic_ocean_boxes = [fram_strait_box,
                            ((76, 90), (16, 100)),     # Barents and Kara Seas, incl. the troughs into them
                            ((70, 90), (100, 180)),
                            ((70, 90), (-180, -45))]
const eurasian_basin_boxes = [((-90, 90), (-25, 140))]

# ## Greenland-Scotland Ridge
#
# The southern boundary of the Nordic Seas: the sections north_sea, iceland_scotland and denmark_strait of sections.toml
# (from Norway north of Shetland to Scotland, along the Wyville Thomson and Iceland-Faroe Ridges to Iceland, across Iceland
# and Denmark Strait to Greenland), so the Nordic Seas are what these transports flow into. Waypoints joined by great
# circles as in the transports, extended over land along the latitude of their ends to 45°W in Greenland and 16°E in
# Norway. `ridge_latitude` needs the waypoints in order of longitude along the whole ridge.

const greenland_scotland_ridge = let
    ridge = unique(sort([(λ, φ) for name in ("denmark_strait", "iceland_scotland", "north_sea")
                                for (λ, φ) in zip(region_section(name)["lon"], region_section(name)["lat"])]))
    all(diff(first.(ridge)) .> 0) || error("The Greenland-Scotland Ridge sections in sections.toml must follow each other " *
                                           "in longitude, without two waypoints at the same longitude")
    ridge
end
const ridge_longitudes = (-45, 16)   # the ends of the ridge beyond the sections, on land

# Latitude at longitude λ of the great circle through (λ₀, φ₀) and (λ₁, φ₁)
great_circle_latitude(λ, (λ₀, φ₀), (λ₁, φ₁)) =
    atand((tand(φ₀) * sind(λ₁ - λ) + tand(φ₁) * sind(λ - λ₀)) / sind(λ₁ - λ₀))

# Along the sections a great circle between waypoints; beyond their ends, over land, the latitude of the end (a great
# circle there would bulge north of the North Sea and Skagerrak, which then belong to no region)
function ridge_latitude(λ)
    (west, east) = (first(greenland_scotland_ridge), last(greenland_scotland_ridge))
    ridge_longitudes[1] ≤ λ < first(west) && return last(west)
    first(east) < λ ≤ ridge_longitudes[2] && return last(east)
    for (a, b) in zip(greenland_scotland_ridge[1:end-1], greenland_scotland_ridge[2:end])
        first(a) ≤ λ ≤ first(b) && return great_circle_latitude(λ, a, b)
    end
    return NaN
end

# ## The other regions: the area of each, as a function of longitude and latitude

const region_areas = Dict(
    # Shelf seas, and the shallow rim of the Arctic basins next to them
    "BRTS"  => (λ, φ) -> in_boxes(λ, φ, [((70, 90), (16, 32)),           # western Barents Sea, north of the Norwegian coast
                                         ((-90, 90), (32, 58))]),        # eastern Barents Sea + White Sea
    "KARA"  => (λ, φ) -> in_box(λ, φ, ((-90, 90), (58, 100))),           # incl. the Ob and Yenisei Gulfs
    "LPES"  => (λ, φ) -> in_box(λ, φ, ((66, 90), (100, 180))),
    "CHSH"  => (λ, φ) -> in_boxes(λ, φ, [((66, 90), (-180, -158)),       # Chukchi Sea, north of Bering Strait (66°N)
                                         ((66, 90), (-158, -123))]),     # Beaufort Sea
    "CACH"  => (λ, φ) -> in_boxes(λ, φ, [((70, 90), (-123, -75)),        # Archipelago
                                         ((83, 90), (-75, -45)),         # Lincoln Sea
                                         ((66, 70), (-123, -65)),        # Fury and Hecla Strait, Gulf of Boothia
                                         ((63, 66), (-95, -65))]),       # Foxe Basin
    "BFFB"  => (λ, φ) -> in_box(λ, φ, ((67, 83), (-75, -45))),           # south to Davis Strait (~67°N), incl. Nares Strait
    "HUDS"  => (λ, φ) -> in_box(λ, φ, ((45, 66), (-100, -65))),          # east to the mouth of Hudson Strait (~65°W)
    # Nordic Seas: between the Greenland-Scotland Ridge and Fram Strait, and the Norwegian coast from Tromsø to Vardø
    "NRDS"  => (λ, φ) -> (-45 < λ ≤ 16 && ridge_latitude(λ) ≤ φ < fram_strait_latitude) || in_box(λ, φ, ((66, 70), (16, 32))),
    # European marginal seas
    "NRTS"  => (λ, φ) -> in_boxes(λ, φ, [((58, 62), (-3.5, 9)),          # North Sea + English Channel: east of Orkney,
                                         ((55, 58), (-4.3, 9)),          # of Scotland's east coast (Forth, Moray Firth),
                                         ((51, 55), (-2.5, 9)),          # of England's (the Irish Sea is west of ~3°W),
                                         ((48.4, 51), (-5.2, 9))]),      # and the Channel to Land's End and Ushant
    "BALT"  => (λ, φ) -> in_box(λ, φ, ((-90, 66), (9, 32))),             # Baltic Sea + Kattegat + Skagerrak east of 9°E
    # Sub-Arctic
    "LBRD"  => (λ, φ) -> in_box(λ, φ, ((52, 67), (-67, -43))),           # north to Davis Strait
    "BRGS"  => (λ, φ) -> in_boxes(λ, φ, [((-90, 66), (100, 180)),        # Bering Sea, Asian side
                                         ((-90, 66), (-180, -100))]),    # Bering Sea, American side
    # The rest of the Atlantic sector, south of the Greenland-Scotland Ridge: Irminger Sea, Iceland Basin, Rockall
    # Trough, Bay of Biscay
    "ESPNA" => (λ, φ) -> in_box(λ, φ, ((-90, 67), (-65, 0))),          # to 67°N: the Denmark Strait section bulges north of 66°N
)

# ## Connectivity

const neighbours = CartesianIndex.(((1, 0), (-1, 0), (0, 1), (0, -1)))   # across cell faces

# The neighbour of cell I across a face, or nothing beyond the domain edge. With `periodic`, the first index wraps around
# (longitude on a latitude-longitude grid).
function neighbour(cells, I, Δ, periodic)
    J = I + Δ
    periodic && (J = CartesianIndex(mod1(J[1], size(cells, 1)), J[2]))
    return checkbounds(Bool, cells, J) ? J : nothing
end

# The cells of `cells` connected through each other to `start`
function connected_cells(cells, start; periodic = false)
    connected = falses(size(cells))
    connected[start] = true
    queue = [start]
    while !isempty(queue)
        I = pop!(queue)
        for Δ in neighbours
            J = neighbour(cells, I, Δ, periodic)
            if !isnothing(J) && cells[J] && !connected[J]
                connected[J] = true
                push!(queue, J)
            end
        end
    end
    return connected
end

# The cells of `candidates` connected through each other to the northernmost of them
connected_to_pole(φ, candidates; periodic = false) =
    connected_cells(candidates, argmax(ifelse.(candidates, φ, -Inf)); periodic)

# `basin` with the patches of other ocean columns it encloses: those that touch neither land nor the domain edge
function fill_enclosed(basin, ocean; periodic = false)
    filled = copy(basin)
    seen   = copy(basin)
    for I in findall(ocean .& .!basin)
        seen[I] && continue
        patch = connected_cells(ocean .& .!basin, I; periodic)
        seen .|= patch
        enclosed = all(findall(patch)) do P
            all(neighbours) do Δ
                J = neighbour(ocean, P, Δ, periodic)
                !isnothing(J) && ocean[J]
            end
        end
        enclosed && (filled .|= patch)
    end
    return filled
end

# ## PAN masks: a Dict of region code => BitArray (λ in (-180, 180])

function pan_masks(λ, φ, bottom_depth; periodic_longitude = false)
    ocean = bottom_depth .> 0

    deep    = (bottom_depth .> deep_basin_depth) .& in_boxes.(λ, φ, Ref(arctic_ocean_boxes))
    basins  = fill_enclosed(connected_to_pole(φ, deep; periodic = periodic_longitude), ocean; periodic = periodic_longitude)
    arctic  = basins .| (ocean .& in_box.(λ, φ, Ref(fram_strait_box)))
    eurasian = in_boxes.(λ, φ, Ref(eurasian_basin_boxes))
    masks = Dict{String, BitArray}("EUBA" => arctic .& eurasian, "AMBA" => arctic .& .!eurasian)

    assigned = copy(arctic)
    for code in pan_codes[3:end]
        masks[code] = ocean .& .!assigned .& region_areas[code].(λ, φ)
        assigned .|= masks[code]
    end
    assigned == ocean || error("$(count(ocean .& .!assigned)) ocean columns are in no region: extend the areas in regions.jl.")
    return masks
end

# ## TSAR
#
# The boxes of Muilwijk et al. (2022), from the Lat_bnds and Lon_bnds of TS_Arctic_Regions_Observed_1970_2017.nc: every
# ocean column inside, shelf or basin. They do not overlap, and leave most of the domain out.

const tsar_names = Dict(
    "WEUR" => "Western Eurasian Basin",
    "EEUR" => "Eastern Eurasian Basin",
    "CHUK" => "Chukchi Sea",
    "BEAU" => "Beaufort Gyre",
)

const tsar_codes = ["WEUR", "EEUR", "CHUK", "BEAU"]

const tsar_boxes = Dict(
    "WEUR" => [((80, 90), (0, 55))],
    "EEUR" => [((76, 84), (100, 140))],
    "CHUK" => [((72, 81), (170, 180)), ((72, 81), (-180, -155))],   # across the dateline
    "BEAU" => [((70, 80), (-155, -125))],
)

tsar_masks(λ, φ, bottom_depth) =
    Dict{String, BitArray}(code => (bottom_depth .> 0) .& in_boxes.(λ, φ, Ref(tsar_boxes[code])) for code in tsar_codes)

# ## Masks of a set

const region_sets = Dict("PAN" => (codes = pan_codes, names = pan_names), "TSAR" => (codes = tsar_codes, names = tsar_names))

"""
    region_masks(λ, φ, bottom_depth; set = "PAN", periodic_longitude = false)

The regions of `set` ("PAN" or "TSAR") from cell-centre longitudes `λ`, latitudes `φ` and bottom depths (positive down;
≤ 0 on land), all of the same size: `(; set, codes, names, masks)`, with `masks` a `Dict` of region code => `BitArray`
and `codes` and `names` in the set's order. `periodic_longitude = true` connects the first and last columns (the first
index) across the dateline, as on a global latitude-longitude grid.
"""
function region_masks(λ::AbstractArray, φ::AbstractArray, bottom_depth::AbstractArray; set = "PAN", periodic_longitude = false)
    haskey(region_sets, set) || error("Unknown region set \"$set\": one of $(join(keys(region_sets), ", "))")
    λ = @. mod(λ + 180, 360) - 180   # to [-180, 180) ...
    λ[λ .== -180] .= 180              # ... and (-180, 180]
    φ = min.(φ, prevfloat(90.0))      # the boxes exclude their upper latitude: keep the pole itself inside
    masks = set == "PAN" ? pan_masks(λ, φ, bottom_depth; periodic_longitude) : tsar_masks(λ, φ, bottom_depth)
    (; codes, names) = region_sets[set]
    return (; set, codes, names = [names[code] for code in codes], masks)
end

"""
    region_masks(grid; set = "PAN", within = nothing)

The regions of `set` on the ocean columns of an `ImmersedBoundaryGrid` (the model's, e.g. the grid saved with the
output, or one of GLORYS's or WOA's cells): `(; set, codes, names, masks)`, the masks `BitMatrix`es of size (Nx, Ny).
`within`, a Bool matrix of the same size, leaves the other columns out of every region (e.g. GLORYS's columns outside
the model domain). A grid that is periodic in longitude wraps around across the dateline.
"""
function region_masks(grid; set = "PAN", within = nothing)
    grid = on_architecture(CPU(), grid)
    Nx, Ny = size(grid)[1:2]
    λ = [λnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    φ = [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
    bottom_depth = -grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]   # works for a Field or an offset array
    isnothing(within) || (bottom_depth[.!within] .= 0)
    return region_masks(λ, φ, bottom_depth; set, periodic_longitude = topology(grid, 1) === Periodic)
end
