# Reads the config file (first command-line argument; default config.toml) and the folders (paths.toml), and sets the
# paths. Output and figures go to a folder named after the experiment, inside output_dir and figure_dir.
# Included by setup.jl and by plotting scripts that don't need the grid.

using Pkg
Pkg.activate(@__DIR__)

using TOML

config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

include("paths.jl")   # paths (paths.toml or paths_local.toml), data_dir, dir_kw
haskey(config, "paths") &&
    @warn "$config_file has a [paths] section, which is no longer read: the folders come from $paths_file"

experiment = get(config, "experiment", "")   # "" (or no setting): directly in output_dir and figure_dir
output_dir = mkpath(joinpath(paths["output_dir"], experiment))   # model output
figure_dir = mkpath(joinpath(paths["figure_dir"], experiment))   # figures and movies
run_config_file = joinpath(output_dir, "config.toml")   # copy of the config of the latest run, written by panarctic.jl
run_paths_file  = joinpath(output_dir, "paths.toml")    # and of its folders

# Output files, written by panarctic.jl and read by plotting/surface_movie.jl
ocean_output_file   = joinpath(output_dir, "arctic_ocean.jld2")
sea_ice_output_file = joinpath(output_dir, "arctic_sea_ice.jld2")
checkpoint_prefix   = "arctic_checkpoint"   # checkpoints in output_dir, written by panarctic.jl
transport_output_file         = joinpath(output_dir, "arctic_transports.jld2")           # online, by panarctic.jl
offline_transport_output_file = joinpath(output_dir, "arctic_transports_offline.jld2")   # by analysis/transports.jl

# Derived diagnostics of the run, small and worth keeping: in diagnostics_dir/<experiment> if set, else in output_dir
run_diagnostics_dir = isempty(get(paths, "diagnostics_dir", "")) ? output_dir :
                      mkpath(joinpath(paths["diagnostics_dir"], experiment))

# by analysis/ocean_diagnostics.jl, for a region set (regions.jl)
ocean_content_file(set) = joinpath(run_diagnostics_dir, "ocean_content_$(lowercase(set)).jld2")
ts_profiles_file(set)   = joinpath(run_diagnostics_dir, "ts_profiles_$(lowercase(set)).jld2")
