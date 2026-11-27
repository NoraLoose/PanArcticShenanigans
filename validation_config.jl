# Reads the validation settings (validation.toml) and sets the paths of the validation datasets and their diagnostics,
# in validation_dir (paths.toml, or your paths_local.toml). Included by the scripts that compute diagnostics of
# validation datasets (analysis/validation_diagnostics.jl) and by those that compute a run's diagnostics or compare it
# with them; the settings are the same for every run, so that all are computed alike.
#
# <validation_dir>/<dataset>/raw                              downloads and copies of the dataset
# <validation_dir>/<dataset>/<family>/<set>/regions.jld2      the dataset's regions of a region set, and the settings
# <validation_dir>/<dataset>/<family>/<set>/<period>.jld2     its diagnostics over one period (a day, a decade, ...)
# with <family> a family of diagnostics (diagnostics/): ocean_content, ts_profiles.

using TOML

isdefined(Main, :paths) || include("paths.jl")   # paths, data_dir (already there if config.jl was included)

# The value after `--name` on the command line, or `default`
function command_line_option(name, default = nothing)
    i = findfirst(==("--$name"), ARGS)
    return isnothing(i) ? default : ARGS[i + 1]
end

validation_file = joinpath(@__DIR__, "validation.toml")
validation      = TOML.parsefile(validation_file)
haskey(paths, "validation_dir") || @warn "No validation_dir in $paths_file: using validation_data (see paths.toml)"
validation_dir  = abspath(get(paths, "validation_dir", "validation_data"))
@info "Validation settings from $validation_file; datasets and their diagnostics in $validation_dir"

# Heat, salt and freshwater content (diagnostics/ocean_content.jl): the same settings for runs and datasets
function content_settings(validation)
    settings = validation["ocean_content"]
    return (; reference_temperature = Float64(settings["reference_temperature"]),
              reference_salinity    = Float64(settings["reference_salinity"]),
              depth_bounds          = vcat(Float64.(settings["depth_bounds"]), Inf))   # Inf: the whole water column
end

# The region sets (diagnostics/regions.jl) of every regional diagnostic
region_set_names(validation) = String.(validation["region_sets"])

dataset_names = sort(collect(keys(validation["datasets"])))
dataset_options(name) = haskey(validation["datasets"], name) ? validation["datasets"][name] :
                        error("No dataset $name in $validation_file: one of $(join(dataset_names, ", "))")

raw_dir(dataset) = mkpath(joinpath(validation_dir, dataset, "raw"))
diagnostics_dir(dataset, family, set)       = mkpath(joinpath(validation_dir, dataset, family, lowercase(set)))
regions_file(dataset, family, set)          = joinpath(diagnostics_dir(dataset, family, set), "regions.jld2")
period_file(dataset, family, set, period)   = joinpath(diagnostics_dir(dataset, family, set), "$period.jld2")
period_files(dataset, family, set) =
    filter(f -> endswith(f, ".jld2") && f != "regions.jld2", readdir(diagnostics_dir(dataset, family, set)))

# ## GLORYS transports (validation/glorys.jl), for the sections of a run's config: need transports.jl
#
# <validation_dir>/GLORYS/transports/<key>/<day>.jld2, a folder per set of sections and reference values ([transports]
# of the config, which the run uses online too), named by a checksum of them and described in its sections.toml, so
# that runs with the same sections share the days.

using SHA: sha1
using JLD2: jldopen
using Dates: Date

# Raised when the way the transports are computed changes (2: each face counted once where two segments meet, #52), so
# that GLORYS's transports computed before are not reused
const transport_method = 2

function glorys_transport_dir(config)
    sections   = transport_sections(config)
    references = reference_values(config)
    key = bytes2hex(sha1(repr((sections, references, transport_method))))[1:12]
    dir = mkpath(joinpath(validation_dir, "GLORYS", "transports", key))
    description = joinpath(dir, "sections.toml")
    isfile(description) || open(description, "w") do io
        TOML.print(io, Dict("reference_temperature" => references.reference_temperature,
                            "reference_salinity"    => references.reference_salinity,
                            "sections" => Dict(name => Dict("lon" => lons, "lat" => lats) for (name, lons, lats) in sections)))
    end
    return dir
end

"""
    read_glorys_transports(config, days)

GLORYS's transports through the sections of `config` on those of `days` that have been computed
(analysis/validation_diagnostics.jl GLORYS --run): `(; dates, transports)`, transports a `Dict` of
"<section>_<quantity>" => one value per date; `nothing` if none has.
"""
function read_glorys_transports(config, days)
    dir   = glorys_transport_dir(config)
    dates = filter(date -> isfile(joinpath(dir, "$date.jld2")), Date.(days))
    isempty(dates) && return nothing
    transports = Dict{String, Vector{Float64}}()
    for date in dates
        jldopen(joinpath(dir, "$date.jld2")) do file
            for key in keys(file["transports"])
                push!(get!(transports, key, Float64[]), file["transports/$key"])
            end
        end
    end
    return (; dates, transports)
end
