# Validation

Compare a run with observations and reanalyses: heat, salt and freshwater content and mean temperature and salinity
profiles by region, and transports through the gateways.

## Quickstart

```bash
julia +1.13 --project=. analysis/compare_with_validation.jl config_myname.toml
```

This computes whatever is missing (the datasets' diagnostics, the run's own) and writes the figures to `figure_dir`.
Several runs go in the same figures, told apart by their experiment names:

```bash
julia +1.13 --project=. analysis/compare_with_validation.jl config_a.toml config_b.toml
```

Their figures go to a folder `compare_<a>_<b>` in the figure folder. `--colors red,royalblue` sets the runs' colours, in
the order of the config files. The
first time, it downloads WOA23 (about 75 GB transferred, 7 GB kept) and computes GLORYS for the run's days; after
that, a new run of the same period only needs its own diagnostics. Set `validation_dir` (where datasets and their
diagnostics go) in your `paths_local.toml` (see the main [README](../README.md)).

## Datasets

| Dataset | What | Period | Grid | Gives | Settings |
|---|---|---|---|---|---|
| GLORYS | Copernicus GLORYS12 reanalysis, daily means | 1993–present: the run's days | 1/12°, 50 levels | contents, profiles, transports | [validation.toml](../validation.toml) |
| WOA23 | NOAA World Ocean Atlas 2023, climatologies (Reagan et al., 2024) | 1985–94, 1995–2004, 2005–14, 2015–22 | 0.25°, 102 levels | contents, profiles | [validation.toml](../validation.toml) |
| Tsubouchi2018 | mooring-based gateway transports (Tsubouchi et al., 2018) | Sep 2005 – Aug 2006, monthly | 4 gateways | transports | [references/](references) |

The gridded datasets (GLORYS, WOA23) go through the model's own code on their own grids, over the part inside the model
domain: their regions follow their own coasts and sea floor, so their contents differ from the model's by a small
offset as well. WOA23 at 0.25° has no monthly climatologies for 1985–94 and 1995–2004; those use the seasonal ones.

Every dataset is downloaded by the code from where it is published (GLORYS from Copernicus Marine with your login; WOA23
from NOAA NCEI), into `<validation_dir>/<dataset>/raw`. A copy you already have can be linked there instead.

## Regions

| Region set | Regions | Covers |
|---|---|---|
| PAN | 15 sub-basins, adapted from hycom-validation | the whole domain, without overlaps |
| TSAR | 4 boxes: Western and Eastern Eurasian Basin, Chukchi Sea, Beaufort Gyre (Muilwijk et al., 2022) | where observed profiles are averaged |

<img src="docs/regions_pan.png" width="49%"> <img src="docs/regions_tsar.png" width="49%">

Defined by longitude, latitude and bottom depth only (`regions.jl`), so they work for any grid spacing and on each
dataset's grid. `julia +1.13 --project=. plotting/regions.jl config.toml` draws them on a run's grid.

## Figures

| File (in `figure_dir`) | Shows |
|---|---|
| `regional_{heat,freshwater,salt}_content_<set>.png` | content by region over the run: the run(s), GLORYS, each WOA23 decade |
| `regional_{temperature,salinity}_profiles_<set>.png` | mean profiles by region: model and GLORYS on the first and last day (several runs: each run's last day), WOA23 per decade |
| `gateway_transports.png` | volume, heat and freshwater through each section: model and GLORYS |
| `tsubouchi_comparison.png`, `tsubouchi_comparison.md` (in `output_dir`) | monthly gateway transports against Tsubouchi et al. (2018), for runs in Sep 2005 – Aug 2006 |

## Settings and files

- [`validation.toml`](../validation.toml): the settings, the same for every run so that all are computed alike: region
  sets, reference temperature and salinity, depth bounds, and the datasets with their options. The sections and the
  transports' reference values stay in each run's `[transports]`, as the model uses them online too.
- [`references/`](references): small reference tables that exist nowhere else in usable form (e.g. values read from a
  paper's figures), kept in the repo with their source and uncertainty in the header. Anything large or downloadable
  goes in `validation_dir` instead.
- `<validation_dir>/<dataset>/raw/`: downloads and copies of the dataset.
- `<validation_dir>/<dataset>/diagnostics/<set>/`: its regional diagnostics, `regions.jld2` and one file per period (a
  day for GLORYS, a decade for WOA23), shared by all runs. They are computed again when the settings change.
- `<validation_dir>/GLORYS/transports/<key>/`: GLORYS's transports, a file per day, for one set of sections and reference
  values (described in its `sections.toml`).

## Running the steps alone

```bash
julia +1.13 --project=. analysis/validation_diagnostics.jl WOA23
julia +1.13 --project=. analysis/validation_diagnostics.jl GLORYS --run config.toml   # or --from 1993-01-01 --to 1993-02-01
julia +1.13 --project=. analysis/regional_diagnostics.jl config.toml
julia +1.13 --project=. plotting/regional_diagnostics.jl config.toml
julia +1.13 --project=. plotting/gateway_transports.jl config.toml
julia +1.13 --project=. analysis/tsubouchi_comparison.jl config.toml
```

## Adding a dataset

1. A reader in `validation/<dataset>.jl`. For a gridded dataset, a function that builds an Oceananigans grid of its
   cells and sets T and S for each period, then passes them to `regional_diagnostics` (see `woa23.jl`).
2. A `[datasets.<dataset>]` section in `validation.toml`, with its options.
3. A line for it in `analysis/validation_diagnostics.jl` (and in `analysis/compare_with_validation.jl`), and a row in
   the table above.
4. Its lines in `plotting/regional_diagnostics.jl`, read like GLORYS's or WOA23's.
