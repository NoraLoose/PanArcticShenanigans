# GLORYS's daily-mean volume, heat and freshwater transport through the four gateways of gateways.jl, on the grid of a
# run file, for the days of that run: <output_dir>/glorys_gateways.csv, drawn by plot_gateways.jl as black dots.
# Builds only the grid, on the CPU, so it needs no model and no GPU.
#
#   julia --project=. glorys_gateways_run.jl runs/glorys-sqrtH-4d.toml

ENV["PANARCTIC_ARCHITECTURE"] = "CPU"
include("setup.jl")
include("gateways.jl")

glorys_gateways(on_architecture(CPU(), grid), GLORYSDaily(), dates, output_dir; days = run_days, dir_kw...)
