# PanArcticShenanigans

Panarctic configuration, adapted from the Barents Sea regional example [here](https://github.com/NumericalEarth/BjerknesWorkshop/blob/main/tutorials/day2/notebooks/02_barents_sea_regional.ipynb).

https://github.com/user-attachments/assets/6f6c94c5-4a65-4a45-ad01-e083ba5aa039


## Requirements

- **Julia 1.13.** `Manifest.toml` is written by Julia 1.13 and does not work with 1.12 (`Zstd_jll` is part of Julia from 1.13 on, so 1.12 reports it as not installed). With [juliaup](https://github.com/JuliaLang/juliaup): `juliaup add 1.13`, then run with `julia +1.13`.
- **An NVIDIA GPU.** The script uses `arch = GPU()`. The 0.15° × 20-level setup needs about 4 GB of GPU memory.
- **Data credentials** for the downloads:
  - GLORYS12 (ocean initial and boundary data): a [Copernicus Marine](https://marine.copernicus.eu) account, set as `COPERNICUSMARINE_SERVICE_USERNAME` and `COPERNICUSMARINE_SERVICE_PASSWORD`.
  - ECCO4 (initial sea ice): NASA Earthdata WebDAV credentials, set as `ECCO_USERNAME` and `ECCO_WEBDAV_PASSWORD`. See NumericalEarth's [ECCO README](https://github.com/NumericalEarth/NumericalEarth.jl/blob/main/src/DataWrangling/ECCO/README.md) for where to find them (they are not your Earthdata login password).
  - JRA55 (atmosphere, rivers) and ETOPO (bathymetry) need no account.

## Running

```bash
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia +1.13 --project=. panarctic.jl
```

Environment variables the script reads:

| Variable | Default | Purpose |
|---|---|---|
| `DATA_DIR` | each dataset's own cache | Where downloads are stored. The first run downloads tens of GB (GLORYS, JRA55), so point this at a disk with space. |
| `STOP_HOURS` | `240` (10 days) | Run length in simulated hours. Use e.g. `STOP_HOURS=3` for a quick check. |

The first run is slow: it downloads the data and fills gaps in GLORYS (cached afterwards). Later runs start in a few minutes.

**What a healthy run looks like** (0.15°, one GH200): a warning listing ~29,000 slightly unstable faces near the surface (expected, see #3), then hourly progress lines with max |u|, |v| around 0.4–1.5 m/s, max |w| around 1 cm/s, and no NaNs. The model runs at about 0.24 s per simulated hour.

Output (`arctic_*.jld2`, `*.mp4`) is written to the current directory.
