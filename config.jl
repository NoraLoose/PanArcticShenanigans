# Reads the config file (first command-line argument; default config.toml) and sets the paths.
# Included by setup.jl and by plotting scripts that don't need the grid.

using Pkg
Pkg.activate(@__DIR__)

using TOML

config_file = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : ARGS[1]
config = TOML.parsefile(config_file)
@info "Settings from $config_file" config

data_dir   = config["paths"]["data_dir"]
dir_kw     = isempty(data_dir) ? (;) : (; dir = data_dir)
output_dir = mkpath(config["paths"]["output_dir"])   # model output
figure_dir = mkpath(config["paths"]["figure_dir"])   # figures and movies

# Output files, written by panarctic.jl and read by plotting/surface_movie.jl
ocean_output_file   = joinpath(output_dir, "arctic_ocean_daily.jld2")
sea_ice_output_file = joinpath(output_dir, "arctic_sea_ice_daily.jld2")
checkpoint_prefix   = "arctic_checkpoint"   # checkpoints in output_dir, written by panarctic.jl
