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

The first run downloads tens of GB (GLORYS, JRA55), so set `data_dir` to a disk with space. For a quick check, set `end_date = 1993-01-01T15:00` (3 hours).

The first run is slow: it downloads the data and fills gaps in GLORYS (cached afterwards). Later runs start in a few minutes.

Output is written to `output_dir`: daily means of the 3D ocean (T, S, u, v), the free surface and the sea ice (`arctic_ocean_daily.jld2`, `arctic_sea_ice_daily.jld2`). The averaging interval is set in `[output]` in the config.

## Figures and movies

Plotting is separate from the model run. There is one script per figure or movie; each takes the same config file as `panarctic.jl` and writes to `output_dir`:

```bash
julia +1.13 --project=. plotting/bathymetry.jl          [config.toml]   # bathymetry.png
julia +1.13 --project=. plotting/grid_spacing.jl        [config.toml]   # grid_spacing.png
julia +1.13 --project=. plotting/glorys_free_surface.jl [config.toml]   # arctic_free_surface.mp4 (GLORYS boundary data)
julia +1.13 --project=. plotting/surface_movie.jl       [config.toml]   # arctic_sea.mp4 (SST, SSS, surface speed, ice volume) from the daily means
```

`surface_movie.jl` only reads the output files and doesn't need a GPU. The other three build the grid (and `glorys_free_surface.jl` reads the GLORYS data), so run them where `panarctic.jl` runs.
