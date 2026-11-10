# Transports through the sections in [transports] and the open boundaries (see transports.jl), computed offline from
# the time-averaged output of panarctic.jl.
#
#   julia +1.13 --project=. analysis/transports.jl [config.toml]
#
# Reads arctic_ocean.jld2 from `output_dir` and writes arctic_transports_offline.jld2 there. Runs on the CPU and reads one
# average at a time. Uses the grid saved with the output, so the faces are those of the run, even if the config's grid
# has changed since. The cell thicknesses follow the averaged free surface.
#
# The transports are computed from the averaged u, v, T and S, so they miss the covariance of velocity and tracers
# within each average; the online transports (arctic_transports.jld2) include it.

include("../config.jl")
include("../transports.jl")

using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Models.HydrostaticFreeSurfaceModels: update_zstar_scaling!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using JLD2

saved = Dict(name => FieldTimeSeries(ocean_output_file, name; backend = OnDisk()) for name in ("u", "v", "T", "S", "η"))
grid  = saved["u"].grid
times = saved["u"].times

u, v = XFaceField(grid), YFaceField(grid)
T, S = CenterField(grid), CenterField(grid)

sections   = vcat(gateway_sections(grid, config), boundary_sections(grid))
transports = section_transports(grid, u, v, T, S, sections; reference_values(config)...)
values     = Dict(name => zeros(length(times)) for name in keys(transports))

for n in eachindex(times)
    for (field, name) in ((u, "u"), (v, "v"), (T, "T"), (S, "S"))
        interior(field) .= interior(saved[name][n])
    end
    # T and S on the open-boundary faces average the edge cell with the halo: fill it with the edge value. (Not u, v:
    # their default boundary conditions would zero the velocities on the open boundaries.)
    fill_halo_regions!((T, S))
    update_zstar_scaling!(grid, saved["η"][n])

    for (name, transport) in transports
        values[name][n] = compute!(transport)[1, 1, 1]
    end
    @info @sprintf("%d / %d: %s", n, length(times), prettytime(times[n]))
end

jldopen(offline_transport_output_file, "w") do file
    save_transport_metadata!(file, config, sections;
                             reference_density = default_reference_density, heat_capacity = default_heat_capacity)
    file["times"] = times
    for (name, value) in values
        file["transports/$name"] = value
    end
end

@info "Wrote $offline_transport_output_file"
