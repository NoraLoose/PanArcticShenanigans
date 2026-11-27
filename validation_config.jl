# Reads the settings shared by every run for the regional diagnostics (validation.toml), so that all runs are computed
# alike. Included by the scripts that compute or plot a run's regional diagnostics.

using TOML

validation_file = joinpath(@__DIR__, "validation.toml")
validation      = TOML.parsefile(validation_file)
@info "Shared settings from $validation_file"

# Heat, salt and freshwater content (ocean_content.jl): the same settings for every run
function content_settings(validation)
    settings = validation["ocean_content"]
    return (; reference_temperature = Float64(settings["reference_temperature"]),
              reference_salinity    = Float64(settings["reference_salinity"]),
              depth_bounds          = vcat(Float64.(settings["depth_bounds"]), Inf))   # Inf: the whole water column
end

# The region sets (regions.jl) of every regional diagnostic
region_set_names(validation) = String.(validation["region_sets"])
