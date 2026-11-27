# Observed annual-mean temperature and salinity profiles of the 4 TSAR regions (regions.jl "TSAR"), 1970–2017, a
# validation dataset: TS_Arctic_Regions_Observed_1970_2017.nc (Muilwijk et al., 2022, J. Climate,
# doi:10.1175/JCLI-D-22-0349.1; observations of Polyakov et al., 2020, doi:10.3389/fmars.2020.00491; data: Muilwijk and
# Polyakov, 2022, NSF Arctic Data Center, doi:10.18739/A20P0WS23, CC0): potential temperature and practical salinity,
# already averaged over each region, on 1 m levels from 1 to 1000 m. The averages mix observations unevenly in space
# and time (see Polyakov et al., 2020). The file is looked for in <validation_dir>/TSAR/raw and downloaded from the
# Arctic Data Center if missing (3 MB); a copy you already have can be linked there instead.
#
# The file's two deepest levels (999 and 1000 m) are left out: their temperatures are 26–30 °C in most regions and
# years (all 165 values at 1000 m, 16 at 999 m), which is not physical; everything above is within −2 to 5 °C.
#
# As profiles that are already regional means, they give the temperature and salinity of diagnostics/ts_profiles.jl but
# no contents (no volumes). tsar_ts_profiles! keeps them in the same form as the other datasets' profiles.
# Needs validation_config.jl and diagnostics/ts_profiles.jl.

using NumericalEarth.DataWrangling: NCDatasets
using .NCDatasets: Dataset
using Downloads: Downloads

const tsar_url  = "https://cn.dataone.org/cn/v2/resolve/urn:uuid:2553d855-f866-4a78-b232-e3e4c53544d8"
const tsar_size = 3082908   # bytes (MD5 5f21e636b648ac5d496b2f4205af7f25)

# The observations file in <validation_dir>/TSAR/raw, downloaded from the Arctic Data Center if it is not there
function tsar_file()
    file = joinpath(raw_dir("TSAR"), "TS_Arctic_Regions_Observed_1970_2017.nc")
    if !isfile(file)
        @info "Downloading the TSAR observations from the NSF Arctic Data Center"
        Downloads.download(tsar_url, file * ".part")
        filesize(file * ".part") == tsar_size || error("The TSAR download is not the expected file ($tsar_url)")
        mv(file * ".part", file)
    end
    return file
end

"""
    tsar_ts_profiles!(options; set_names)

The observed profiles as the TSAR region set's diagnostics: <validation_dir>/TSAR/ts_profiles/tsar/1970-2017.jld2 (or
the file's years), with "temperature" and "salinity" (level, region, year), the levels from the bottom up, and "dates"
the years. Only for the TSAR region set, if it is in `set_names`.
"""
function tsar_ts_profiles!(options; set_names)
    "TSAR" in set_names || return @warn "The TSAR observations are for the TSAR region set, which is not in region_sets"

    path = tsar_file()
    ds = Dataset(path)
    names = [strip(String(collect(ds["region"][r, :])), ['\0', ' ']) for r in 1:size(ds["region"], 1)]
    years = round.(Int, ds["year"][:])
    depth = Float64.(ds["depth"][:])
    θ = Float64.(coalesce.(ds["Observed_temp"][:, :, :], NaN))   # (region, year, depth)
    S = Float64.(coalesce.(ds["Observed_sal"][:, :, :], NaN))
    close(ds)

    regions = region_masks(zeros(1, 1), zeros(1, 1), zeros(1, 1); set = "TSAR")   # codes and names only
    order = [findfirst(==(name), names) for name in regions.names]
    any(isnothing, order) && error("The regions of $path ($(join(names, ", "))) are not those of the TSAR set")

    deepest = findlast(<(999), depth)   # the levels at 999 and 1000 m are left out (see above)
    levels  = deepest:-1:1               # from the bottom up, as the other datasets' levels
    results = [(; temperature = permutedims(θ[order, y, levels]), salinity = permutedims(S[order, y, levels]))
               for y in eachindex(years)]

    save_dataset_regions("TSAR", "ts_profiles", regions, nothing; level_depths = depth[levels])
    save_dataset_period("TSAR", "ts_profiles", "TSAR", "$(first(years))-$(last(years))", results; dates = string.(years))
end
