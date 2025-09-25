module JavaTestUtils

using ArchimedLight

export locate_java_jar, with_java_case, read_semicolon_table, read_table, parse_float, filter_rows, sum_columns, compare_java_julia,
       reference_tag, load_reference_case, load_reference_table

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

const REFERENCE_ROOT = joinpath(@__DIR__, "java_reference")

struct ReferenceCase
    case::String
    config::String
    tag::String
    path::String
    meta::Dict{String,Any}
    log::String
    run_paths::Vector{String}
end

reference_tag(config::AbstractString, extra_args::Vector{String}=String[]) = begin
    base = replace(config, r"[^A-Za-z0-9]+" => "_")
    if isempty(extra_args)
        return base
    end
    extra = join(extra_args, "_")
    extra = replace(extra, r"[^A-Za-z0-9=.-]+" => "_")
    extra = replace(extra, "--" => "")
    extra = replace(extra, "__" => "_")
    extra = strip(extra, '_')
    isempty(extra) ? base : string(base, "__", extra)
end

function _read_metadata(path::AbstractString)
    meta = Dict{String,Any}()
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        parts = split(stripped, "=", limit=2)
        length(parts) == 2 || continue
        key = strip(parts[1])
        valstr = strip(parts[2])
        if startswith(valstr, "\"") && endswith(valstr, "\"")
            meta[key] = replace(valstr[2:end-1], "\\\"" => "\"")
        elseif valstr == "true" || valstr == "false"
            meta[key] = valstr == "true"
        elseif startswith(valstr, "[") && endswith(valstr, "]")
            inner = strip(valstr[2:end-1])
            if isempty(inner)
                meta[key] = String[]
            else
                items = split(inner, ",")
                arr = String[]
                for item in items
                    item = strip(item)
                    if startswith(item, "\"") && endswith(item, "\"")
                        push!(arr, replace(item[2:end-1], "\\\"" => "\""))
                    else
                        push!(arr, item)
                    end
                end
                meta[key] = arr
            end
        else
            meta[key] = valstr
        end
    end
    return meta
end

function load_reference_case(case::AbstractString, config::AbstractString; extra_args::Vector{String}=String[])
    tag = reference_tag(config, extra_args)
    case_dir = joinpath(REFERENCE_ROOT, case, tag)
    isdir(case_dir) || error("Reference artefacts missing for $case/$tag. Run scripts/update_java_references.jl")
    meta_path = joinpath(case_dir, "metadata.toml")
    isfile(meta_path) || error("Missing metadata for $case/$tag")
    meta = _read_metadata(meta_path)
    runs = String[]
    for entry in get(meta, "runs", String[])
        push!(runs, joinpath(case_dir, String(entry)))
    end
    log_path = joinpath(case_dir, "log.txt")
    log_text = isfile(log_path) ? read(log_path, String) : ""
    return ReferenceCase(String(case), String(config), tag, case_dir, meta, log_text, runs)
end

function load_reference_table(case::AbstractString, config::AbstractString, filename::AbstractString; extra_args::Vector{String}=String[], delim::AbstractString=";")
    ref = load_reference_case(case, config; extra_args=extra_args)
    isempty(ref.run_paths) && error("No run outputs stored for $(ref.case)/$(ref.tag)")
    table_path = joinpath(ref.run_paths[1], filename)
    isfile(table_path) || error("File $(filename) missing in reference run $(ref.case)/$(ref.tag)")
    header, rows = read_table(table_path; delim=delim)
    return header, rows, ref
end

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

function compare_java_julia(cfgdir::AbstractString, config::AbstractString; metrics::Vector{String}=DEFAULT_METRICS, java_args::Vector{String}=String[], julia_outdir::Union{Nothing,String}=nothing, jar::Union{Nothing,AbstractString}=nothing)
    case = basename(cfgdir)
    _, rows_java, ref = load_reference_table(case, config, "component_values.csv"; extra_args=java_args)
    get(ref.meta, "success", true) || error("Reference run for $(case)/$(ref.tag) is marked as failure")
    java_totals = sum_columns(rows_java, metrics)

    cfgpath = joinpath(cfgdir, config)
    outdir = isnothing(julia_outdir) ? mktempdir() : julia_outdir
    ArchimedLight.Runner.run_from_config(cfgpath; outdir=outdir)
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

"""
    compare_and_report(cfgdir, config; metrics=DEFAULT_METRICS, tol=1e-3, strict=false)

Uses stored Java reference outputs and the Julia runner to compute aggregate metrics.
Returns the same tuple as `compare_java_julia`.
"""
function compare_and_report(cfgdir::AbstractString, config::AbstractString; metrics::Vector{String}=DEFAULT_METRICS, tol::Float64=1e-3, strict::Bool=false, java_args::Vector{String}=String[], julia_outdir::Union{Nothing,String}=nothing, jar::Union{Nothing,AbstractString}=nothing)
    cmp = compare_java_julia(cfgdir, config; metrics=metrics, java_args=java_args, julia_outdir=julia_outdir, jar=jar)
    @info "Java vs Julia comparison" cmp
    for col in metrics
        jv = get(cmp.java, col, 0.0)
        lv = get(cmp.julia, col, 0.0)
        dv = get(cmp.diff, col, 0.0)
        rv = get(cmp.rel, col, 0.0)
        @info "metric" col "java" jv "julia" lv "diff" dv "rel" rv
        if strict && rv > tol
            error("Relative difference for metric $col is $rv > $tol (java=$jv julia=$lv diff=$dv)")
        end
    end
    return cmp
end

end # module
