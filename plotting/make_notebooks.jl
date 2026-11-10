# Turns each plotting script into a Jupyter notebook in plotting/notebooks/ with Literate.jl, to look at the figures.
#
#   julia +1.13 plotting/make_notebooks.jl [config.toml]
#
# Comment lines become text cells. Each notebook starts with a cell that runs from the repository root (so the includes
# and the config's relative paths resolve as when running the script) and builds on the CPU, since the Jupyter kernel
# runs there. The config file is a setting in that cell: the argument (relative to the repository root), or config.toml.
# Each figure is shown after it is saved, a movie at the end. Rerun this after changing a plotting script.
# Literate is installed in a temporary environment, so the project's dependencies stay as they are.

using Pkg
Pkg.activate(temp = true)
Pkg.add("Literate"; io = devnull)
using Literate

repository  = dirname(@__DIR__)
config_file = isempty(ARGS) ? "config.toml" : ARGS[1]
isfile(joinpath(repository, config_file)) || error("Config file $config_file not found (relative to $repository).")

const preamble = """
# ## Settings
#
# The config and the architecture. The cell runs from the repository root, as the script does, so the includes and
# relative paths in the config resolve the same way.

isfile("panarctic.jl") || cd(joinpath("..", ".."))   # from plotting/notebooks to the repository root
repository = pwd()

config_file  = "$config_file"   # settings, as for panarctic.jl
architecture = "CPU"           # the Jupyter kernel runs on the CPU

using TOML
settings = TOML.parsefile(config_file)
settings["architecture"] = architecture
notebook_config = joinpath(mktempdir(), basename(config_file))
open(io -> TOML.print(io, settings), notebook_config, "w")
empty!(ARGS); push!(ARGS, notebook_config)   # read by config.jl, as the command-line argument of the script

# ## The script

"""

# Show the movie the script wrote, embedded in the notebook
movie_cell(file) = """

#-
using Base64
HTML("<video controls width=700 src='data:video/mp4;base64,\$(base64encode(read(joinpath(figure_dir, "$file"))))'></video>")
"""

function preprocess(script)
    script = replace(script, r"include\(\"\.\./(.+?)\"\)" => s"include(joinpath(repository, \"\1\"))")
    script = replace(script, r"^([ \t]+)# "m => s"\1## ")              # indented comments stay in the code (Literate: ## → #)

    # The script's opening comments (its description) stay first; the settings cell goes after them
    lines  = split(script, '\n')
    n      = something(findfirst(l -> !(startswith(l, "#") || isempty(strip(l))), lines), length(lines) + 1)
    header = join(lines[1:n-1], '\n')
    body   = join(lines[n:end], '\n')

    # Show each figure in a cell of its own after it is saved, and the movie at the end
    body  = replace(body, r"^(save\(joinpath\(figure_dir, \".+?\"\), (\w+)\))$"m => s"\1\n\n#-\n\2\n\n#-")
    movie = match(r"record\(fig, joinpath\(figure_dir, \"(.+?\.mp4)\"\)", script)
    movie_display = isnothing(movie) ? "" : movie_cell(movie[1])

    return header * "\n\n" * preamble * body * movie_display
end

plotting   = @__DIR__
notebooks  = mkpath(joinpath(plotting, "notebooks"))
scripts    = filter(f -> endswith(f, ".jl") && f != basename(@__FILE__), readdir(plotting))

for script in scripts
    Literate.notebook(joinpath(plotting, script), notebooks; preprocess, execute = false, credit = false)
end
