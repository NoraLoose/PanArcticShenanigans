# PanArcticShenanigans

Panarctic configuration, adapted from the Barents Sea regional example [here](https://github.com/NumericalEarth/BjerknesWorkshop/blob/main/tutorials/day2/notebooks/02_barents_sea_regional.ipynb).

https://github.com/user-attachments/assets/6f6c94c5-4a65-4a45-ad01-e083ba5aa039


## Requirements

- **Julia 1.13.** `Manifest.toml` is written by Julia 1.13 and does not work with 1.12 (`Zstd_jll` is part of Julia from 1.13 on, so 1.12 reports it as not installed). With [juliaup](https://github.com/JuliaLang/juliaup): `juliaup add 1.13`, then run with `julia +1.13`.
- **Data credentials** for the downloads, as environment variables (kept out of the config file so passwords are never committed):
  - GLORYS12 (ocean initial and boundary data): a [Copernicus Marine](https://marine.copernicus.eu) account, set as `COPERNICUSMARINE_SERVICE_USERNAME` and `COPERNICUSMARINE_SERVICE_PASSWORD`.
  - ECCO4 (initial sea ice): NASA Earthdata WebDAV credentials, set as `ECCO_USERNAME` and `ECCO_WEBDAV_PASSWORD`. See NumericalEarth's [ECCO README](https://github.com/NumericalEarth/NumericalEarth.jl/blob/main/src/DataWrangling/ECCO/README.md) for where to find them (they are not your Earthdata login password).
  - JRA55 (atmosphere, rivers) and GEBCO (bathymetry, ~7.5 GB) need no account. Regridding GEBCO onto a new grid needs ~55 GB of memory (e.g. `--mem=80G`); the result is cached.

## Running

```bash
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia +1.13 --project=. panarctic.jl
```

Settings (grid, dates, time step, run length, data and output folders) are in [`config.toml`](config.toml). To use your own, copy it and pass it as an argument:

```bash
cp config.toml config_myname.toml   # edit as needed
julia +1.13 --project=. panarctic.jl config_myname.toml
```

The first run downloads tens of GB (GLORYS, JRA55), so set `data_dir` to a disk with space. For a quick check, set `stop_hours = 3`.

The first run is slow: it downloads the data and fills gaps in GLORYS (cached afterwards). Later runs start in a few minutes.

Output (`arctic_*.jld2`) is written to `output_dir`.

### Timing

The log reports how long each set-up phase takes (grid and bathymetry, GLORYS data, open-boundary transports, model
construction, initial state, monitors), then, inside `run!`, the initialization, the first time step (where the time
step is compiled when the model type is new in the Julia session) and the seconds per model hour of the later steps.
A summary table is printed when the run ends.

### Running again in the same session (ocean only)

Building the model and compiling its first time step can take 15 minutes or more. For ocean-only runs
(`[initialization] ocean_only = true`, or a rest test), `rerun.jl` starts again from the same initial state in the
same Julia session, without loading data, balancing or compiling:

```julia
julia> empty!(ARGS); push!(ARGS, "first_run.toml"); include("panarctic.jl")   # builds the model, runs, keeps the initial state
julia> empty!(ARGS); push!(ARGS, "second_run.toml"); include("rerun.jl")      # same model, new output_dir, from the start again
```

The second config may change `[paths]` (each run needs its own `output_dir`), `stop_hours` (up to the last GLORYS day
the first run loaded) and `[diagnostics]`. Anything else is part of the model or its initial state; `rerun.jl` stops
with a message if it differs, and a fresh session is needed. A run stopped with Ctrl-C can be rerun the same way.

## Figures and movies

Plotting is separate from the model run. There is one script per figure or movie; each takes the same config file as `panarctic.jl` and writes to `output_dir`:

```bash
julia +1.13 --project=. plotting/bathymetry.jl          [config.toml]   # bathymetry.png
julia +1.13 --project=. plotting/grid_spacing.jl        [config.toml]   # grid_spacing.png
julia +1.13 --project=. plotting/glorys_free_surface.jl [config.toml]   # arctic_free_surface.mp4 (GLORYS boundary data)
julia +1.13 --project=. plotting/surface_movie.jl       [config.toml]   # arctic_sea.mp4 (SST, SSS, surface speed, ice volume) from the model output
```

`surface_movie.jl` only reads the output files and doesn't need a GPU. The other three build the grid (and `glorys_free_surface.jl` reads the GLORYS data), so run them where `panarctic.jl` runs.
