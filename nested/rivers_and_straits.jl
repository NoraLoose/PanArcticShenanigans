# Fresh water into the Norwegian Coastal Current: is the Baltic connected to the North Sea on the model grid, and how much
# river water reaches the model ocean in each region compared with what JRA55 discharges there?
#
#   julia --project=. rivers_and_straits.jl runs/<run>.toml
#
# CPU only, no run needed. Prints
#   1. whether the Baltic's surface water connects to the Skagerrak, and a map of the Danish straits (column depth);
#   2. JRA55 river discharge by region at its mouths, against what the model's river routing (JRA55PrescribedLand, the
#      same call as panarctic_nested.jl) delivers to ocean cells in that region, averaged over the first days of the run.

ENV["PANARCTIC_ARCHITECTURE"] = "CPU"
include("setup.jl")
include("jra55_read_fix.jl")

using Oceananigans.Grids: inactive_node, λnode, φnode
using Oceananigans.Operators: Azᶜᶜᶜ, Δzᶜᶜᶜ

Nx, Ny, Nz = size(grid)
C = Center()
cell_λ = [mod(λnode(i, j, 1, grid, C, C, C) + 180, 360) - 180 for i in 1:Nx, j in 1:Ny]
cell_φ = [φnode(i, j, 1, grid, C, C, C) for i in 1:Nx, j in 1:Ny]
wet    = [!inactive_node(i, j, Nz, grid, C, C, C) for i in 1:Nx, j in 1:Ny]
depth  = [sum(Δzᶜᶜᶜ(i, j, k, grid) for k in 1:Nz if !inactive_node(i, j, k, grid, C, C, C); init = 0.0) for i in 1:Nx, j in 1:Ny]

function nearest_wet(λ, φ)
    d = [wet[i, j] ? ((cell_λ[i, j] - λ) * cosd(φ))^2 + (cell_φ[i, j] - φ)^2 : Inf for i in 1:Nx, j in 1:Ny]
    return Tuple(argmin(d))
end

# ## 1. The Danish straits

baltic, skagerrak = nearest_wet(19.5, 57.5), nearest_wet(9.5, 58.0)    # Gotland Basin, Skagerrak
reached = falses(Nx, Ny)
frontier = [baltic]
reached[baltic...] = true
while !isempty(frontier)
    i, j = pop!(frontier)
    for (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))
        (1 ≤ a ≤ Nx && 1 ≤ b ≤ Ny && wet[a, b] && !reached[a, b]) || continue
        reached[a, b] = true
        push!(frontier, (a, b))
    end
end
@printf("\nBaltic (cell %s) and Skagerrak (cell %s): %s\n", baltic, skagerrak,
        reached[skagerrak...] ? "CONNECTED at the surface" : "NOT CONNECTED: the Baltic is a closed lake on this grid")

in_straits = [(i, j) for i in 1:Nx, j in 1:Ny if 8.5 ≤ cell_λ[i, j] ≤ 14.5 && 54.0 ≤ cell_φ[i, j] ≤ 58.5]
irange = extrema(first.(in_straits))
jrange = extrema(last.(in_straits))
symbol(i, j) = !wet[i, j] ? '.' : depth[i, j] < 20 ? '1' : depth[i, j] < 50 ? '2' : depth[i, j] < 100 ? '3' :
                                  depth[i, j] < 200 ? '4' : '5'
println("\nDanish straits, cells i = $(irange[1])–$(irange[2]) (columns), j = $(jrange[2]) down to $(jrange[1]) (rows):")
println(". land   1 < 20 m   2 < 50 m   3 < 100 m   4 < 200 m   5 deeper   B reached from the Baltic, S Skagerrak cell")
for j in jrange[2]:-1:jrange[1]
    row = map(irange[1]:irange[2]) do i
        (i, j) == skagerrak ? 'S' : (wet[i, j] && reached[i, j] && cell_λ[i, j] > 12.5) ? 'B' : symbol(i, j)
    end
    @printf("%4d %s   %.1f°E %.1f°N\n", j, String(row), cell_λ[irange[1], j], cell_φ[irange[1], j])
end

# ## 2. Rivers

regions = [("Baltic",          (λ, φ) -> 13.5 ≤ λ ≤ 31 && 53 ≤ φ ≤ 66),
           ("Danish straits",  (λ, φ) -> 8 ≤ λ < 13.5 && 54.5 ≤ φ ≤ 59.5),
           ("North Sea",       (λ, φ) -> -4 ≤ λ < 8 && 50 ≤ φ ≤ 59.5),
           ("Norwegian coast", (λ, φ) -> (4 ≤ λ ≤ 31 && 59.5 < φ ≤ 72) && !(13.5 ≤ λ ≤ 31 && φ ≤ 66)),
           ("Kara (Ob, Yenisei)", (λ, φ) -> 55 ≤ λ ≤ 100 && 64 ≤ φ ≤ 78),
           ("Laptev (Lena)",   (λ, φ) -> 100 < λ ≤ 145 && 68 ≤ φ ≤ 78),
           ("Mackenzie",       (λ, φ) -> -142 ≤ λ ≤ -125 && 66 ≤ φ ≤ 72)]
region_of(λ, φ) = something(findfirst(r -> r[2](mod(λ + 180, 360) - 180, φ), regions), 0)

jra55 = (; dataset = MultiYearJRA55(), start_date = dates[1], end_date = dates[2], dir_kw...)
land = JRA55PrescribedLand(grid; jra55...)
rivers = land.freshwater_flux.rivers
routing = land.river_routing.rivers
source_grid = on_architecture(CPU(), rivers.grid)

records = 1:min(10, length(rivers.times))
flux = sum(Array(interior(rivers[n]))[:, :, 1] for n in records) / length(records)   # kg m⁻² s⁻¹

source = zeros(length(regions))
for I in findall(flux .> 0)
    λ, φ = λnode(I[1], I[2], 1, source_grid, C, C, C), φnode(I[1], I[2], 1, source_grid, C, C, C)
    r = region_of(λ, φ)
    r > 0 && (source[r] += flux[I] * Azᶜᶜᶜ(I[1], I[2], 1, source_grid) / 1000)
end

delivered = zeros(length(regions))
target_i, target_j, offsets = Array(routing.target_i), Array(routing.target_j), Array(routing.offsets)
outlet_i, outlet_j, weight = Array(routing.contribution_outlet_i), Array(routing.contribution_outlet_j), Array(routing.contribution_weight)
for c in eachindex(target_i)
    i, j = target_i[c], target_j[c]
    r = region_of(cell_λ[i, j], cell_φ[i, j])
    r == 0 && continue
    for m in offsets[c]:offsets[c + 1] - 1
        delivered[r] += weight[m] * flux[outlet_i[m], outlet_j[m]] * Azᶜᶜᶜ(i, j, Nz, grid) / 1000
    end
end

@printf("\nRiver discharge, mean of the first %d JRA55 records from %s (m³/s; 1 mSv = 1000 m³/s)\n", length(records), dates[1])
@printf("%-20s %14s %18s\n", "region", "JRA55 mouths", "delivered to ocean")
for (r, (name, _)) in enumerate(regions)
    @printf("%-20s %14.0f %18.0f\n", name, source[r], delivered[r])
end
println("""
Mouths are counted where JRA55 puts them, delivery where the model's routing puts the water (nearest wet cell within
about 3°, spread over 1.2°), so water can be counted in a neighbouring region. Large shortfalls summed over neighbouring
regions are water that never reaches the ocean.""")
