# Reads the settings shared by every run for the regional diagnostics (validation.toml), so that all runs are computed
# alike. Included by the scripts that compute a run's regional diagnostics.

using TOML

validation_file = joinpath(@__DIR__, "validation.toml")
validation      = TOML.parsefile(validation_file)
@info "Shared settings from $validation_file"

# Heat, salt and freshwater content and profiles: the same settings for every run
function regional_settings(validation)
    settings = validation["regional_diagnostics"]
    return (; reference_temperature = Float64(settings["reference_temperature"]),
              reference_salinity    = Float64(settings["reference_salinity"]),
              depth_bounds          = vcat(Float64.(settings["depth_bounds"]), Inf))   # Inf: the whole water column
end

region_set_names(validation) = String.(validation["regional_diagnostics"]["region_sets"])
