include(joinpath(@__DIR__, "backend_comparison.jl"))

function _local_realistic_configs()
    raw = get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_CONFIGS", "")
    if isempty(strip(raw))
        return [joinpath(REPO_ROOT, "example_2", "config.yml")]
    end
    return [abspath(strip(path)) for path in split(raw, ',') if !isempty(strip(path))]
end

function local_realistic_main()
    output = get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT", BENCH_OUTPUT)
    rows = benchmark_configs(_local_realistic_configs(); output=output)
    println("ArchimedLight local realistic backend benchmark")
    for row in rows
        @printf(
            "%-24s %-24s median=%8.4f s alloc=%8.1f MiB PAR=%12.3f NIR=%12.3f\n",
            row.config,
            row.backend,
            row.median_seconds,
            row.median_alloc_mib,
            row.incident_par,
            row.incident_nir,
        )
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    local_realistic_main()
end
