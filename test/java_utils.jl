module JavaTestUtils

using ArchimedLight

export locate_java_jar, with_java_case, read_semicolon_table, read_table, parse_float, filter_rows, sum_columns, compare_java_julia

"""Locate the Archimed 2018 reference JAR used as oracle for integration tests."""
function locate_java_jar()
    if haskey(ENV, "ARCHIMED_JAR")
        jar = ENV["ARCHIMED_JAR"]
        isfile(jar) && return jar
    end

    candidates = String[
        # Typical location when ArchimedLight and archimed-2018-source live side-by-side
        normpath(joinpath(@__DIR__, "..", "ARCHIMED", "archimed-2018-source", "archimed-lib-2018", "target", "archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar")),
        # Fallback when archimed-2018-source is a sibling of the repo root
        normpath(joinpath(@__DIR__, "..", "..", "ARCHIMED", "archimed-2018-source", "archimed-lib-2018", "target", "archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar")),
        normpath(joinpath(@__DIR__, "..", "..", "archimed-2018-source", "archimed-lib-2018", "target", "archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar"))
    ]
    for jar in candidates
        isfile(jar) && return jar
    end
    return nothing
end

function clean_output_dirs!(dir::AbstractString)
    for entry in readdir(dir)
        path = joinpath(dir, entry)
        if isdir(path) && startswith(lowercase(entry), "output")
            rm(path; recursive=true, force=true)
        end
    end
end

function output_directory_from_yaml(path::AbstractString)
    base = dirname(path)
    for line in readlines(path)
        stripped = strip(line)
        startswith(stripped, "#") && continue
        if startswith(stripped, "output_directory")
            parts = split(stripped, ":", limit=2)
            if length(parts) == 2
                val = strip(parts[2])
                val = split(val, "#", limit=2)[1]
                val = strip(val)
                val = replace(val, '"' => "")
                val = replace(val, "'" => "")
                isempty(val) && return joinpath(base, "output")
                return joinpath(base, val)
            end
        end
    end
    return joinpath(base, "output")
end

"""
    with_java_case(cfgdir, config; jar, extra_args=String[], f)

Copy the Java integration test directory to a temporary workspace, execute the ARCHIMED
reference jar with the provided configuration, and pass execution details to `f`.
`f` receives `(success::Bool, log_text::String, run_dirs::Vector{String})` and should return
whatever data the caller needs. Temporary artifacts are deleted once `f` returns.
"""
function with_java_case(f::Function, cfgdir::AbstractString, config::AbstractString; jar::AbstractString, extra_args::Vector{String}=String[])
    mktempdir() do tmp
        dest = joinpath(tmp, "case")
        cp(cfgdir, dest; force=true)
        clean_output_dirs!(dest)
        config_path = joinpath(dest, config)
        base_outdir = output_directory_from_yaml(config_path)
        rm(base_outdir; recursive=true, force=true)
        mkpath(base_outdir)
        log_path = joinpath(base_outdir, "log")

        args = vcat([config], extra_args)
        command = Cmd(vcat(["java", "-Xmx1024m", "-jar", jar], args))
        success = false
        Base.Filesystem.cd(dest) do
            open(log_path, "w") do logio
                pipeline_cmd = pipeline(command; stdout=logio, stderr=logio)
                try
                    run(pipeline_cmd)
                    success = true
                catch e
                    if isa(e, Base.ProcessFailedException)
                        success = false
                    else
                        rethrow(e)
                    end
                end
            end
        end

        log_text = isfile(log_path) ? read(log_path, String) : ""
        run_dirs = sort([joinpath(base_outdir, name) for name in readdir(base_outdir) if isdir(joinpath(base_outdir, name))])
        return f(success, log_text, run_dirs)
    end
end

function read_table(path::AbstractString; delim::AbstractString=";")
    lines = readlines(path)
    isempty(lines) && return String[], Dict{String,String}[]
    header = split(strip(lines[1]), delim)
    rows = Dict{String,String}[]
    for line in lines[2:end]
        stripped = strip(line)
        isempty(stripped) && continue
        values = split(stripped, delim)
        row = Dict{String,String}()
        for (col, val) in zip(header, values)
            row[col] = val
        end
        push!(rows, row)
    end
    return header, rows
end

"""Read a semi-colon delimited file into a vector of dictionaries keyed by header."""
read_semicolon_table(path::AbstractString) = read_table(path; delim=";")

parse_float(row::Dict{String,String}, key::String) = parse(Float64, row[key])

filter_rows(rows::Vector{Dict{String,String}}, key::String, value::AbstractString) = [row for row in rows if get(row, key, nothing) == value]

function sum_columns(rows::Vector{Dict{String,String}}, cols::Vector{String})
    totals = Dict{String,Float64}()
    for col in cols
        accumulator = 0.0
        for row in rows
            if haskey(row, col)
                try
                    accumulator += parse(Float64, row[col])
                catch
                    # ignore unparsable values
                end
            end
        end
        totals[col] = accumulator
    end
    return totals
end

const DEFAULT_METRICS = [
    "Ri_PAR_0_f", "Ri_NIR_0_f",
    "Ra_PAR_0_f", "Ra_NIR_0_f",
    "Ri_PAR_0_q", "Ri_NIR_0_q",
    "Ra_PAR_0_q", "Ra_NIR_0_q"
]

function compare_java_julia(cfgdir::AbstractString, config::AbstractString; jar::AbstractString, metrics::Vector{String}=DEFAULT_METRICS, java_args::Vector{String}=String[], julia_outdir::Union{Nothing,String}=nothing)
    java_success = Ref(false)
    java_rows = Ref(Vector{Dict{String,String}}())
    with_java_case(cfgdir, config; jar=jar, extra_args=java_args) do success, _log, run_dirs
        java_success[] = success
        if success && !isempty(run_dirs)
            csv_path = joinpath(run_dirs[1], "component_values.csv")
            _, rows = read_table(csv_path; delim=";")
            java_rows[] = rows
        else
            java_rows[] = Dict{String,String}[]
        end
        return nothing
    end

    java_success[] || error("Java simulation failed for $config")
    rows_java = java_rows[]
    java_totals = sum_columns(rows_java, metrics)

    cfgpath = joinpath(cfgdir, config)
    outdir = isnothing(julia_outdir) ? mktempdir() : julia_outdir
    res = ArchimedLight.Runner.run_from_config(cfgpath; outdir=outdir)
    julia_csv = joinpath(outdir, "components.csv")
    _, rows_julia = read_table(julia_csv; delim=",")
    julia_totals = sum_columns(rows_julia, metrics)

    diffs = Dict{String,Float64}()
    rel = Dict{String,Float64}()
    for col in metrics
        jv = get(java_totals, col, 0.0)
        lv = get(julia_totals, col, 0.0)
        diff = lv - jv
        diffs[col] = diff
        denom = max(abs(jv), 1e-9)
        rel[col] = abs(diff) / denom
    end

    return (java=java_totals, julia=julia_totals, diff=diffs, rel=rel)
end

end # module
