# GLORYS's own sea floor and vertical levels, so the model can start from GLORYS on (nearly) GLORYS's geometry.
# Included by setup.jl when [grid] bathymetry = "glorys" or vertical = "glorys".
#
# glorys_vertical_faces(): GLORYS's 50 levels, with the same interfaces NumericalEarth uses when it reads GLORYS
#   (midway between level centres, the surface at 0, the bottom one spacing below the deepest midpoint).
# glorys_bottom_height(grid): GLORYS's `deptho` (static file) regridded onto `grid` with the same steps as
#   NumericalEarth's regrid_bathymetry (interpolation passes, minimum depth, removal of minor basins). NumericalEarth
#   cannot read `deptho` itself (it treats it as 3-D and has no value for land), so the file is read here. Land is set
#   to `land_height` before interpolating, as `height_above_water` does in regrid_bathymetry.

using NumericalEarth.DataWrangling: NCDatasets, metadata_path
using NumericalEarth.DataWrangling.GLORYS: GLORYSStatic
using .NCDatasets: Dataset, dimnames
using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Utils: launch!
using Downloads, Printf

function glorys_vertical_faces()
    metadatum = Metadatum(:temperature; date = first(dates), dataset, region, dir_kw...)
    Downloads.download(metadatum)                                  # skips a file already there
    ds = Dataset(metadata_path(metadatum))
    zc = -reverse(Float64.(ds["depth"][:]))                        # level centres, deepest first
    close(ds)
    zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
    return vcat(zf[1] - (zc[2] - zc[1]), zf, 0.0)
end

function glorys_bottom_height(target_grid; minimum_depth = 15, interpolation_passes = 2, major_basins = 1,
                              land_height = 50)
    metadatum = Metadatum(:depth; dataset = GLORYSStatic(), region, dir_kw...)
    Downloads.download(metadatum)
    path = metadata_path(metadatum)
    ds = Dataset(path)
    var = ds["deptho"]
    names = dimnames(var)
    lon = Float64.(ds["longitude"][:])
    lat = Float64.(ds["latitude"][:])
    depth = Float64.(coalesce.(var[:, :], NaN))                    # positive down; missing on land
    close(ds)
    names[1] == "longitude" || (depth = permutedims(depth))
    size(depth) == (length(lon), length(lat)) || error("unexpected shape $(size(depth)) of deptho in $path")

    z = ifelse.(isnan.(depth) .| (depth .≤ 0), Float64(land_height), -depth)

    Δλ, Δφ = lon[2] - lon[1], lat[2] - lat[1]
    TX = isapprox(length(lon) * Δλ, 360; atol = Δλ / 2) ? Periodic : Bounded
    native_grid = LatitudeLongitudeGrid(CPU(); size = (length(lon), length(lat), 1),
                                        longitude = (lon[1] - Δλ / 2, lon[end] + Δλ / 2),
                                        latitude  = (lat[1] - Δφ / 2, min(90, lat[end] + Δφ / 2)),
                                        z = (0, 1), topology = (TX, Bounded, Bounded), halo = (10, 10, 1))
    native_z = Field{Center, Center, Nothing}(native_grid)
    set!(native_z, z)
    fill_halo_regions!(native_z)

    cpu_target = on_architecture(CPU(), target_grid)
    target_z = NumericalEarth.Bathymetry.interpolate_bathymetry_in_passes(native_z, cpu_target;
                                                                          passes = interpolation_passes)
    minimum_depth > 0 &&
        launch!(CPU(), cpu_target, :xy, NumericalEarth.Bathymetry._enforce_minimum_depth!, target_z, minimum_depth)
    NumericalEarth.Bathymetry.remove_minor_basins!(target_z, major_basins)
    fill_halo_regions!(target_z)

    @info @sprintf("GLORYS bathymetry from %s: %d of %d columns wet, deepest %.0f m", path,
                   count(<(0), interior(target_z)), length(interior(target_z)), -minimum(interior(target_z)))
    return target_z   # on the CPU
end
