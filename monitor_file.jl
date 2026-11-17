# Files of the monitor (monitor.jl), and reading them back for the plotting scripts without loading the model.

monitor_file(output_dir) = joinpath(output_dir, "monitor.csv")
depth_file(output_dir)   = joinpath(output_dir, "depth.csv")

# The header of monitor.csv and a function returning a column by name
function read_monitor(output_dir)
    text   = readlines(monitor_file(output_dir))
    header = split(text[1], ',')
    rows   = [parse.(Float64, split(line, ',')) for line in text[2:end] if !isempty(line)]
    column(name) = [row[findfirst(==(name), header)] for row in rows]
    return header, column
end
