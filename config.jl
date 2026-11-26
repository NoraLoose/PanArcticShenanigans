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

# by analysis/regional_diagnostics.jl, for a region set (regions.jl)
regional_output_file(set) = joinpath(run_diagnostics_dir, "arctic_regional_diagnostics_$(lowercase(set)).jld2")
