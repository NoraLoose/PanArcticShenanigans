# The data, output, figure and validation folders, the same for every config file: your own in paths_local.toml (a copy
# of paths.toml, not in the repo), otherwise the defaults in paths.toml. Included by config.jl, and by the validation
# scripts, which need no config file (validation_config.jl).

using TOML

paths_file = joinpath(@__DIR__, isfile(joinpath(@__DIR__, "paths_local.toml")) ? "paths_local.toml" : "paths.toml")
paths = TOML.parsefile(paths_file)["paths"]
@info "Folders from $paths_file" paths

data_dir = paths["data_dir"]
dir_kw   = isempty(data_dir) ? (;) : (; dir = data_dir)
