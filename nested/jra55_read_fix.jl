# Workaround for NumericalEarth's `jra55_read_data` (DataWrangling/JRA55/JRA55_field_time_series.jl): with
# `Cyclical()` time indexing, the records held in memory wrap from the end of the run back to its start, so the
# indices into the year's JRA55 file are not sorted. NumericalEarth splits them where the index equals 1, which
# exists only when the run starts on the first record of the file (1 January). Any other start date (here
# 1 September 2005) fails with `MethodError: no method matching -(::Nothing, ::Int64)` while building the river
# routing, or later when a reload wraps. Read the indices in sorted order and put them back in the requested order.
# To be fixed upstream; remove this file then.

@eval NumericalEarth.DataWrangling.JRA55 function jra55_read_data(ds, name, i, j, nn)
    issorted(nn) && return ds[name][i, j, nn]
    order = sortperm(nn)
    data  = ds[name][i, j, nn[order]]
    return data[:, :, invperm(order)]
end
