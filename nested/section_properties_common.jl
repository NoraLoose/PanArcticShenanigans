# Shared by section_properties.jl (model) and glorys_gateways_native.jl (GLORYS): inflow and outflow through a section,
# and its mean section, from per-face-and-level records (key, λ, φ, z, Δz, width, velocity, T, S; key[4] is the
# level). Needs gateways.jl.

# On a staircase, each face measures the flow normal to itself, not to the line: flow along the line goes in through
# one face and out through the next. So transports are first added up over stretches of `bin_km` along the line, level
# by level, and only then split into inflow and outflow. Bins start at the line's first point, the same for the model
# and for GLORYS. Each bin keeps its transport, heat and freshwater, area-mean T and S, and mean normal velocity.
function binned(records, start, Tʳᵉᶠ; bin_km = 30)
    bins = Dict()
    for r in records
        q = r.velocity * r.width * r.Δz
        isfinite(q) && isfinite(r.T) && isfinite(r.S) || continue
        s = 6371 * deg2rad(great_circle(start..., r.λ, r.φ))
        key = (floor(Int, s / bin_km), r.key[4])
        A = r.width * r.Δz
        b = get(bins, key, (; λ = r.λ, φ = r.φ, z = 0.0, A = 0.0, width = 0.0, q = 0.0, heat = 0.0, freshwater = 0.0, TA = 0.0, SA = 0.0))
        bins[key] = merge(b, (; z = b.z + r.z * A, A = b.A + A, width = b.width + r.width, q = b.q + q,
                               heat = b.heat + ρ₀cₚ * q * (r.T - Tʳᵉᶠ), freshwater = b.freshwater + q * (1 - r.S / Sʳᵉᶠ),
                               TA = b.TA + r.T * A, SA = b.SA + r.S * A))
    end
    return [(; key, s = (key[1] + 1/2) * bin_km, b.λ, b.φ, z = b.z / b.A, Δz = b.A / b.width, b.width, velocity = b.q / b.A,
               T = b.TA / b.A, S = b.SA / b.A, b.q, b.heat, b.freshwater) for (key, b) in bins]
end

# Volume (Sv), heat (TW) and freshwater (mSv) carried by the inflow (bins with net flow into the Arctic) and the outflow
function inflow_outflow(bins)
    sums = zeros(6)
    for b in bins
        n = b.q > 0 ? 0 : 3
        sums[n + 1] += b.q / 1e6
        sums[n + 2] += b.heat / 1e12
        sums[n + 3] += b.freshwater / 1e3
    end
    return sums
end

inout_header(names) = join(["hour"; [n * "_" * d * "_" * q for n in names for d in ("in", "out")
                                     for q in ("volume_Sv", "heat_TW", "freshwater_mSv")]], ',')

# Add one day's bins to the running sums per bin and level
function accumulate_section!(means, records)
    for r in records
        isfinite(r.velocity) && isfinite(r.T) && isfinite(r.S) || continue
        m = get!(means, r.key, (; r.s, r.λ, r.φ, r.z, r.Δz, r.width, velocity = 0.0, T = 0.0, S = 0.0, n = 0))
        means[r.key] = merge(m, (; velocity = m.velocity + r.velocity, T = m.T + r.T, S = m.S + r.S, n = m.n + 1))
    end
    return means
end

# The means of the bins (distance along the line at the bin's centre)
function write_section(path, means)
    open(path, "w") do io
        println(io, "s_km,lon,lat,z_m,dz_m,width_m,velocity_ms,T_C,S")
        for m in values(means)
            println(io, join((@sprintf("%.6g", x) for x in (m.s, m.λ, m.φ, m.z, m.Δz, m.width,
                                                            m.velocity / m.n, m.T / m.n, m.S / m.n)), ','))
        end
    end
end
