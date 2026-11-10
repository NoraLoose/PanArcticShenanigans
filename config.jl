# Reads the config file (first command-line argument; default config.toml) and sets the paths. Output and figures go
# to a folder named after the experiment, inside output_dir and figure_dir.
# Included by setup.jl and by plotting scripts that don't need the grid.

using Pkg
Pkg.activate(@__DIR__)

using TOML

config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

data_dir   = config["paths"]["data_dir"]
dir_kw     = isempty(data_dir) ? (;) : (; dir = data_dir)
experiment = get(config, "experiment", "")   # "" (or no setting): directly in output_dir and figure_dir
output_dir = mkpath(joinpath(config["paths"]["output_dir"], experiment))   # model output
figure_dir = mkpath(joinpath(config["paths"]["figure_dir"], experiment))   # figures and movies
run_config_file = joinpath(output_dir, "config.toml")   # copy of the config of the latest run, written by panarctic.jl

# Output files, written by panarctic.jl and read by plotting/surface_movie.jl
ocean_output_file   = joinpath(output_dir, "arctic_ocean.jld2")
sea_ice_output_file = joinpath(output_dir, "arctic_sea_ice.jld2")
checkpoint_prefix   = "arctic_checkpoint"   # checkpoints in output_dir, written by panarctic.jl
