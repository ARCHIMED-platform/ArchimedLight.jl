include(joinpath(@__DIR__, "backend_comparison.jl"))

function _complexity_directions()
    raw = get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_DIRECTIONS", "6,16,46")
    return Union{Nothing,Int}[parse(Int, strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

function _complexity_configs()
    raw = get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_CONFIGS", "")
    if isempty(strip(raw))
        return [joinpath(REPO_ROOT, "example_2", "config.yml")]
    end
    return [abspath(strip(path)) for path in split(raw, ',') if !isempty(strip(path))]
end

function scene_complexity_main()
    output = get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_OUTPUT", BENCH_OUTPUT)
    rows = benchmark_configs(
        _complexity_configs();
        directions=_complexity_directions(),
        output=output,
    )
    println("ArchimedLight RasterGPU scene-complexity sweep")
    @printf("%-24s %-22s %8s %10s %10s\n", "config", "backend", "dirs", "median s", "alloc MiB")
    for row in rows
        @printf(
            "%-24s %-22s %8d %10.4f %10.1f\n",
            row.config,
            row.backend,
            row.directions,
            row.median_seconds,
            row.median_alloc_mib,
        )
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    scene_complexity_main()
end
