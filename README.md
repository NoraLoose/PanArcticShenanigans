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

Settings (grid, dates, time step, output) are in [`config.toml`](config.toml). To use your own, copy it and pass it as an argument:

```bash
cp config.toml config_myname.toml   # edit as needed
julia +1.13 --project=. panarctic.jl config_myname.toml
```

The data, output and figure folders are in [`paths.toml`](paths.toml), the same for every config file. To use your own, copy it to `paths_local.toml` (not in the repo) and edit that; it is then used instead of `paths.toml`:

```bash
cp paths.toml paths_local.toml   # edit as needed
```

The first run downloads tens of GB (GLORYS, JRA55), so set `data_dir` to a disk with space. For a quick check, set `end_date = 1993-01-01T03:00` (3 hours).

The first run is slow: it downloads the data and fills gaps in GLORYS (cached afterwards). Later runs start in a few minutes.

The ocean and sea ice are nested in GLORYS with NumericalEarth's `nested_ocean_simulation` and `nested_sea_ice_simulation` ([NumericalEarth.jl#676](https://github.com/NumericalEarth/NumericalEarth.jl/pull/676), not yet released; `Project.toml` pins its branch). They read GLORYS along the open Atlantic and Pacific edges, set the open boundary conditions, and start the run from GLORYS on `start_date`.

Each run has an experiment name (`experiment` in the config). Its output goes to `output_dir/<experiment>/` (default `output/test/`, relative to where you run), together with copies of the config file that ran it and of its folders (`config.toml`, `paths.toml`, replaced by each run of the experiment): time averages (daily by default) of the 3D ocean (T, S, u, v), the free surface and the sea ice (`arctic_ocean.jld2`, `arctic_sea_ice.jld2`), and checkpoints of the full model state (`arctic_checkpoint*.jld2`). The averaging and checkpoint intervals are set in `[output]` in the config.

To continue a run from its latest checkpoint, set `restart = true` and run again with the same experiment name. The run still ends at `end_date`.

## Transports

Volume, heat, salt and freshwater transports through Fram Strait, the Barents Sea Opening, Bering Strait and Davis Strait (Sv, TW, kg/s, Sv; positive into the Arctic), and through the two open boundaries (`atlantic_boundary`, `pacific_boundary`; positive into the domain). The sections are defined in `[transports]` in the config as longitude/latitude waypoints, so they don't depend on the grid: the faces they cross are found on whatever grid the run uses. Both ends of a section should be on land; a warning says when one isn't. Heat and freshwater are relative to `reference_temperature` and `reference_salinity`.

- **Online** (`online = true`, the default): computed every time step and averaged like the rest of the output, written to `arctic_transports.jld2` in `output_dir/<experiment>/`. These include the covariance of velocity and T, S within each average.
- **Offline**, from the time-averaged output of a run, on the CPU:

  ```bash
  julia +1.13 --project=. analysis/transports.jl [config.toml]   # arctic_transports_offline.jld2
  ```

  Uses the grid saved with the output. Being computed from averaged u, v, T, S, these miss the covariance within each average.

To check the sections on the grid (e.g. after changing the grid or a section), plot their faces on the bathymetry with `plotting/transport_sections.jl`.

## Validation

Comparisons with GLORYS, the World Ocean Atlas 2023, observed profiles (TSAR) and gateway transports (Tsubouchi et al., 2018): heat, salt and freshwater content and temperature and salinity profiles by region, and transports. One command computes what is missing and makes the figures:

```bash
julia +1.13 --project=. analysis/compare_with_validation.jl [config.toml]
```

See [validation/README.md](validation/README.md) for the datasets, regions and figures. The datasets go to `validation_dir` (in `paths.toml`); the shared settings are in [`validation.toml`](validation.toml).

## Figures and movies

Plotting is separate from the model run. There is one script per figure or movie; each takes the same config file as `panarctic.jl` and writes to `figure_dir/<experiment>/` (default `figures/test/`), separate from the model output:

```bash
julia +1.13 --project=. plotting/bathymetry.jl          [config.toml]   # bathymetry.png
julia +1.13 --project=. plotting/grid_spacing.jl        [config.toml]   # grid_spacing.png
julia +1.13 --project=. plotting/glorys_free_surface.jl [config.toml]   # arctic_free_surface.mp4 (GLORYS boundary data)
julia +1.13 --project=. plotting/surface_movie.jl       [config.toml]   # arctic_sea.mp4 (SST, SSS, surface speed, ice volume) from the time-averaged output
julia +1.13 --project=. plotting/transport_sections.jl  [config.toml]   # transport_sections.png, transport_section_faces.png (the transport sections on the bathymetry)
```

`surface_movie.jl` only reads the output files and doesn't need a GPU. The others build the grid (and `glorys_free_surface.jl` reads the GLORYS data), so run them where `panarctic.jl` runs.
