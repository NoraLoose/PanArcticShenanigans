# GLORYS's own sea floor (`deptho`, static file) regridded onto the model grid, for [grid] bathymetry = "GLORYS", so the
# gateways and sills have the depths GLORYS's transports were computed over. Same steps as NumericalEarth's
# regrid_bathymetry (interpolation passes, minimum depth, removal of minor basins). NumericalEarth cannot read `deptho`
# as a bathymetry (it has no value on land), so the file is read here and land is set to `land_height` first, as
# `height_above_water` does in regrid_bathymetry. Adapted from the port (fjp/hybrid-obc).

using NumericalEarth.DataWrangling: NCDatasets, metadata_path, BoundingBox, Metadatum
using NumericalEarth.DataWrangling: default_horizontal_padding
using NumericalEarth.DataWrangling.GLORYS: GLORYSStatic, GLORYSDaily
using .NCDatasets: Dataset, dimnames
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Utils: launch!
using Downloads

function glorys_bottom_height(target_grid; minimum_depth = 15, interpolation_passes = 2, major_basins = 1, land_height = 50)
    region = BoundingBox(target_grid; padding = default_horizontal_padding(GLORYSDaily()))
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
    periodic = isapprox(length(lon) * Δλ, 360; atol = Δλ / 2)
    TX = periodic ? Periodic : Bounded
    west = lon[1] - Δλ / 2
    east = periodic ? west + 360 : lon[end] + Δλ / 2
    native_grid = LatitudeLongitudeGrid(CPU(); size = (length(lon), length(lat), 1),
                                        longitude = (west, east),
                                        latitude  = (lat[1] - Δφ / 2, min(90, lat[end] + Δφ / 2)),
                                        z = (0, 1), topology = (TX, Bounded, Bounded), halo = (10, 10, 1))
    native_z = Field{Center, Center, Nothing}(native_grid)
    set!(native_z, z)
    fill_halo_regions!(native_z)

    target_z = NumericalEarth.Bathymetry.interpolate_bathymetry_in_passes(native_z, target_grid; passes = interpolation_passes)
    minimum_depth > 0 &&
        launch!(CPU(), target_grid, :xy, NumericalEarth.Bathymetry._enforce_minimum_depth!, target_z, minimum_depth)
    NumericalEarth.Bathymetry.remove_minor_basins!(target_z, major_basins)
    fill_halo_regions!(target_z)

    @info @sprintf("GLORYS bathymetry from %s: %d of %d columns wet, deepest %.0f m", path,
                   count(<(0), interior(target_z)), length(interior(target_z)), -minimum(interior(target_z)))
    return target_z
end
