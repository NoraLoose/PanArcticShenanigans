# Nested pan-Arctic model

The pan-Arctic ocean, and optionally sea ice and the JRA55 atmosphere, nested in GLORYS with NumericalEarth's
`nested_ocean_model` ([PR #676](https://github.com/NumericalEarth/NumericalEarth.jl/pull/676)), on the same grid as
`panarctic.jl`. It has its own environment: the nested model needs Oceananigans 0.117, while `panarctic.jl` is on
0.116.

NumericalEarth is pinned in `Project.toml` to the branch of
[PR #923](https://github.com/NumericalEarth/NumericalEarth.jl/pull/923), which is #676 plus the option
`zero_net_inflow`. Once #923 is merged the pin moves to #676, and to a release once #676 is merged.

## Running

From this folder:

```bash
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia +1.13 --project=. panarctic_nested.jl                 # settings in config.toml
julia +1.13 --project=. panarctic_nested.jl my_config.toml  # or your own copy
```

The default is a quick check: 100 km (≈1°), ocean only, 4 days; under half an hour on an A100 once the
data are downloaded. For the production grid set `spacing_km = 15` and `dt_minutes = 10`. Credentials for GLORYS are the same
environment variables as for `panarctic.jl` (see the main README).

## What it does

- **Boundaries** (from `nested_ocean_model`): Atlantic (j = 1) and Pacific (j = Ny) edges open. Oblique radiation on
  u, v (9 h inflow nudging); tracer reservoirs (20 km); Flather toward GLORYS's depth-mean velocity times the model's
  depth, with Chapman on η; a sponge on T, S and the baroclinic velocity; the seafloor blended to GLORYS's near the
  open edges. The initial state is GLORYS interpolated onto the grid.
- **`zero_net_inflow`** (`[boundaries]`): removes one depth-uniform velocity from the Flather exterior on both edges,
  so that it has no net inflow. At 1°, ocean only, it reduces the error in domain-mean η after one day from −2.0 cm to
  −0.2 cm relative to GLORYS (−7.8 → −5.6 cm after four days).
- **Coupled** (`[coupling] coupled = true`): sea ice nested in GLORYS, JRA55 atmosphere, radiation and rivers.
  `pressure = false` (default) sets the JRA55 sea-level pressure to a uniform 101325 Pa, because GLORYS's η has no
  inverse-barometer response.

## Output

In `output_dir`:

- `nested_monitor.csv`, hourly, also printed to the log:
  - net transport through the Atlantic and Pacific edges, positive into the domain (Sv);
  - domain-mean η, next to GLORYS's on the model grid;
  - the change in η since the start over the shelf (H < 200 m) and the deep ocean (H > 1000 m);
  - the model's η minus GLORYS's along each edge.
- `nested_surface.jld2` (u, v, T, S at the surface, every 3 hours) and `nested_free_surface.jld2` (η, hourly).
