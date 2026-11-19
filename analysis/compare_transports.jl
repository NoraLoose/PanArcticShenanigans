# Compare the online transports (arctic_transports.jld2) of two experiments: for every section and quantity, the
# largest difference between their averages at the times both have, and the largest value, for scale.
#
#   julia +1.13 --project=. analysis/compare_transports.jl output/<experiment_a> output/<experiment_b>
#
# Used to check that a change leaves the model's answer alone (smoke test before and after), or to see what it moved.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JLD2, Printf

length(ARGS) == 2 || error("Usage: analysis/compare_transports.jl <output folder a> <output folder b>")
files = map(dir -> joinpath(dir, "arctic_transports.jld2"), ARGS)

jldopen(files[1]) do a
    jldopen(files[2]) do b
        units = a["metadata/units"]
        names = sort(intersect(keys(a["timeseries"]), keys(b["timeseries"])))
        @printf("%-40s %12s %12s   %s\n", "transport", "max |a - b|", "max |a|", "averages compared")
        for name in names
            name == "t" && continue
            iterations = sort(intersect(filter(!=("serialized"), keys(a["timeseries/$name"])),
                                        filter(!=("serialized"), keys(b["timeseries/$name"]))), by = s -> parse(Int, s))
            isempty(iterations) && continue
            va = [only(a["timeseries/$name/$i"]) for i in iterations]
            vb = [only(b["timeseries/$name/$i"]) for i in iterations]
            unit = units[last(split(name, "_"))]
            @printf("%-40s %12.3e %12.3e   %d   %s\n", name, maximum(abs, va .- vb), maximum(abs, va), length(iterations), unit)
        end
    end
end
