# Shared by section_properties.jl (model) and glorys_gateways_native.jl (GLORYS): inflow and outflow through a section,
# and its mean section, from per-face-and-level records (λ, φ, z, Δz, width, velocity, T, S). Needs gateways.jl.

# Volume (Sv), heat (TW) and freshwater (mSv) carried by the inflow (velocity > 0, into the Arctic) and by the outflow
function inflow_outflow(records, Tʳᵉᶠ)
    sums = zeros(6)
    for r in records
        q = r.velocity * r.width * r.Δz
        isfinite(q) && isfinite(r.T) && isfinite(r.S) || continue
        n = q > 0 ? 0 : 3
        sums[n + 1] += q / 1e6
        sums[n + 2] += ρ₀cₚ * q * (r.T - Tʳᵉᶠ) / 1e12
        sums[n + 3] += q * (1 - r.S / Sʳᵉᶠ) / 1e3
    end
    return sums
end

inout_header(names) = join(["hour"; [n * "_" * d * "_" * q for n in names for d in ("in", "out")
                                     for q in ("volume_Sv", "heat_TW", "freshwater_mSv")]], ',')

# Add one day's records to the running sums per face and level
function accumulate_section!(means, records)
    for r in records
        isfinite(r.velocity) && isfinite(r.T) && isfinite(r.S) || continue
        m = get!(means, r.key, (; r.λ, r.φ, r.z, r.Δz, r.width, velocity = 0.0, T = 0.0, S = 0.0, n = 0))
        means[r.key] = merge(m, (; velocity = m.velocity + r.velocity, T = m.T + r.T, S = m.S + r.S, n = m.n + 1))
    end
    return means
end

# The means, with the distance along the line from its first point
function write_section(path, start, means)
    open(path, "w") do io
        println(io, "s_km,lon,lat,z_m,dz_m,width_m,velocity_ms,T_C,S")
        for m in values(means)
            s = 6371 * deg2rad(great_circle(start..., m.λ, m.φ))
            println(io, join((@sprintf("%.6g", x) for x in (s, m.λ, m.φ, m.z, m.Δz, m.width,
                                                            m.velocity / m.n, m.T / m.n, m.S / m.n)), ','))
        end
    end
end
