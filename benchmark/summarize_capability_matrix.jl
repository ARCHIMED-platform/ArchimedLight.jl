using Printf

function _usage()
    return """
    Usage:
      julia --project=test/gpu benchmark/summarize_capability_matrix.jl result1.tsv [result2.tsv ...]

    Summarizes TSV files written by benchmark/local_realistic_backends.jl.
    Rows may come from different machines and may include either old or new
    metadata columns. The summary highlights backend fallback status, reducer
    capabilities, stack-profile copy requirements, hit-density metrics, and
    stage-split preparation/public API timings when present.
    """
end

function _read_tsv(path::AbstractString)
    rows = Vector{Dict{String,String}}()
    open(path, "r") do io
        header_line = readline(io)
        columns = split(chomp(header_line), '\t'; keepempty=true)
        for line in eachline(io)
            isempty(strip(line)) && continue
            cells = split(chomp(line), '\t'; keepempty=true)
            row = Dict{String,String}("__source" => path)
            for (i, column) in pairs(columns)
                row[column] = i <= length(cells) ? cells[i] : ""
            end
            push!(rows, row)
        end
    end
    return rows
end

_cell(row, key, default="") = get(row, key, default)

function _cell_any(row, keys::Tuple, default="")
    for key in keys
        value = _cell(row, key)
        isempty(value) || return value
    end
    return default
end

function _float_cell(row, key)
    raw = _cell(row, key)
    isempty(raw) && return 0.0
    return try
        parse(Float64, raw)
    catch
        0.0
    end
end

function _float_cell_any(row, keys::Tuple)
    raw = _cell_any(row, keys)
    isempty(raw) && return 0.0
    return try
        parse(Float64, raw)
    catch
        0.0
    end
end

function _int_cell(row, key)
    raw = _cell(row, key)
    isempty(raw) && return 0
    return try
        parse(Int, raw)
    catch
        0
    end
end

function _int_cell_any(row, keys::Tuple)
    raw = _cell_any(row, keys)
    isempty(raw) && return 0
    return try
        parse(Int, raw)
    catch
        0
    end
end

function _backend_key(row)
    return (
        _cell(row, "run_label"),
        _cell(row, "workload", _cell(row, "selftest", "unknown")),
        _cell(row, "backend"),
        _cell(row, "resolved_backend"),
        _cell(row, "fallback_reason"),
        _cell(row, "reduction_capabilities"),
        _cell(row, "edge_accumulation"),
        _cell(row, "julia_version"),
        _cell(row, "os"),
        _cell(row, "arch"),
    )
end

function _row_status(row)
    fallback = _cell(row, "fallback_reason")
    resolved = _cell(row, "resolved_backend")
    backend = _cell(row, "backend")
    reductions = _cell(row, "reduction_capabilities")
    traced = _int_cell_any(row, ("traced_dirs", "profile_traced_dirs"))
    copy_required = _int_cell_any(row, ("copy_required_dirs", "profile_copy_required_dirs"))
    reduced = _int_cell_any(row, ("reduced_dirs", "profile_reduced_dirs"))
    fallback != "" && fallback != "none" && return "fallback"
    isempty(reductions) && return isempty(resolved) ? "metadata" : "no-raycore-data"
    traced > 0 && copy_required == 0 && reduced > 0 && return "readback-free-candidate"
    traced > 0 && copy_required > 0 && return "host-readback-required"
    backend == resolved || isempty(resolved) || return "resolved-different"
    return "raycore"
end

function _shape_label(row)
    instances = _int_cell(row, "raycore_tlas_instances")
    geometries = _int_cell(row, "raycore_tlas_geometries")
    instances == 0 && geometries == 0 && return ""
    return string(instances, "/", geometries)
end

function _face_label(row)
    faces = _int_cell(row, "raycore_expanded_face_count")
    upper = _int_cell(row, "raycore_expanded_face_instance_upper_bound")
    faces == 0 && upper == 0 && return ""
    return string(faces, "/", upper)
end

function _scratch_label(row)
    dense_pairs = _int_cell(row, "raycore_dense_edge_pairs")
    edge_capacity = _int_cell(row, "raycore_edge_key_capacity")
    dense_pairs == 0 && edge_capacity == 0 && return ""
    return string(dense_pairs, "/", edge_capacity)
end

function _print_summary(rows)
    isempty(rows) && (println("No rows found."); return rows)
    sort!(
        rows;
        by=row -> (
            _cell(row, "workload", _cell(row, "selftest", "unknown")),
            _cell(row, "backend"),
            _cell(row, "edge_accumulation"),
            _cell(row, "run_label"),
            _cell(row, "resolved_backend"),
            _cell(row, "fallback_reason"),
        ),
    )

    println("Capability matrix summary")
    @printf(
        "%-16s %-22s %-18s %-18s %-24s %-18s %-24s %-11s %-16s %-18s %-34s %-45s %9s %8s %8s %8s %8s %9s %9s %9s %9s %8s %9s %9s %-24s\n",
        "label",
        "workload",
        "backend",
        "edge_mode",
        "resolved",
        "mode",
        "ref",
        "tlas",
        "faces",
        "scratch",
        "fallback",
        "reductions",
        "turtles",
        "dirs",
        "copyreq",
        "hits",
        "occup%",
        "trace",
        "count",
        "edge",
        "copy",
        "max",
        "prepare",
        "public",
        "status",
    )
    for row in rows
        input_turtles = _int_cell(row, "profile_input_turtle_count")
        profile_turtles = _int_cell(row, "profile_unique_turtle_count")
        turtle_label = input_turtles == 0 && profile_turtles == 0 ? "" : string(input_turtles, "/", profile_turtles)
        @printf(
            "%-16s %-22s %-18s %-18s %-24s %-18s %-24s %-11s %-16s %-18s %-34s %-45s %9s %8d %8d %8d %8.2f %9.3f %9.3f %9.3f %9.3f %8d %9.3f %9.3f %-24s\n",
            _cell(row, "run_label"),
            _cell(row, "workload", _cell(row, "selftest", "unknown")),
            _cell(row, "backend"),
            _cell(row, "edge_accumulation"),
            _cell(row, "resolved_backend"),
            _cell(row, "geometry_mode"),
            _cell(row, "ref_instancing_status"),
            _shape_label(row),
            _face_label(row),
            _scratch_label(row),
            _cell(row, "fallback_reason"),
            _cell(row, "reduction_capabilities"),
            turtle_label,
            _int_cell_any(row, ("traced_dirs", "profile_traced_dirs")),
            _int_cell_any(row, ("copy_required_dirs", "profile_copy_required_dirs")),
            _int_cell_any(row, ("total_hits", "profile_total_hits", "trace_total_hits")),
            _float_cell_any(row, ("occupied", "profile_occupied", "trace_occupied")),
            _float_cell_any(row, ("trace_ms", "profile_trace_ms")),
            _float_cell_any(row, ("count_area_ms", "profile_count_area_ms")),
            _float_cell_any(row, ("edge_ms", "profile_edge_ms")),
            _float_cell_any(row, ("copy_ms", "profile_copy_ms")),
            _int_cell_any(row, ("max_seen", "profile_max_seen", "trace_max_seen")),
            _float_cell_any(row, ("prepare_seconds", "prepare_median_seconds")),
            _float_cell_any(row, ("public_output_seconds", "reuse_median_seconds", "cached_median_seconds")),
            _row_status(row),
        )
    end

    groups = Dict{Tuple,Int}()
    for row in rows
        key = _backend_key(row)
        groups[key] = get(groups, key, 0) + 1
    end
    println()
    println("Grouped rows: ", length(groups), " unique run/workload/backend states from ", length(rows), " rows.")
    fallbacks = count(row -> _row_status(row) == "fallback", rows)
    readback = count(row -> _row_status(row) == "host-readback-required", rows)
    candidates = count(row -> _row_status(row) == "readback-free-candidate", rows)
    println("Fallback rows: ", fallbacks)
    println("Host-readback-required rows: ", readback)
    println("Readback-free candidate rows: ", candidates)
    return rows
end

function main(args=ARGS)
    if isempty(args) || any(arg -> arg in ("-h", "--help"), args)
        print(_usage())
        return nothing
    end
    rows = Dict{String,String}[]
    for path in args
        append!(rows, _read_tsv(path))
    end
    return _print_summary(rows)
end

main()
