using ArchimedLight
using KernelAbstractions
using Printf
using Statistics

const REPO_ROOT = dirname(@__DIR__)
const BENCH_CONFIG = get(ENV, "ARCHIMEDLIGHT_BENCH_CONFIG", joinpath(REPO_ROOT, "example_2", "config.yml"))
const BENCH_BACKENDS = [
    strip(lowercase(value))
    for value in split(
        get(ENV, "ARCHIMEDLIGHT_BENCH_BACKENDS", "normal_cpu,rastergpu_ka_cpu,rastergpu_metal_gpu"),
        ',',
    ) if !isempty(strip(value))
]
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_SAMPLES", "3"))
const BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WARMUPS", "1"))
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "128"))
const BENCH_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", "256"))
const BENCH_TILE_SIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_SIZE", "1"))
const BENCH_TILE_FACE_CAPACITY =
    parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_FACE_CAPACITY", "128"))
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto"))
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_BENCH_OUTPUT", "")

function _optional_metal_backend()
    "rastergpu_metal_gpu" in BENCH_BACKENDS || "all" in BENCH_BACKENDS || return nothing
    Sys.isapple() && Sys.ARCH == :aarch64 || return nothing
    Base.find_package("Metal") === nothing && return nothing
    metal = Base.require(Main, :Metal)
    getproperty(metal, :functional)() || return nothing
    array_type = getproperty(metal, :MtlArray)
    return KernelAbstractions.get_backend(array_type(zeros(Float32, 1)))
end

function _rastergpu_case(name::String, backend)
    interception = ArchimedLight.RasterGPUBackend(
        backend=backend,
        workgroupsize=BENCH_WORKGROUPSIZE,
        max_hits_per_pixel=BENCH_MAX_HITS,
        tile_size=BENCH_TILE_SIZE,
        tile_face_capacity=BENCH_TILE_FACE_CAPACITY,
        edge_accumulation=BENCH_EDGE_ACCUMULATION,
    )
    return (
        name=name,
        interception_backend=interception,
        scattering_backend=ArchimedLight.RasterGPUScatteringBackend(interception),
    )
end

function _backend_cases()
    allowed = Set(["all", "normal_cpu", "rastergpu_ka_cpu", "rastergpu_metal_gpu"])
    unknown = setdiff(Set(BENCH_BACKENDS), allowed)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_BENCH_BACKENDS value(s): $(join(sort!(collect(unknown)), ", ")). " *
        "Use normal_cpu, rastergpu_ka_cpu, rastergpu_metal_gpu, or all.",
    )

    selected(name) = "all" in BENCH_BACKENDS || name in BENCH_BACKENDS
    cases = NamedTuple[]
    selected("normal_cpu") && push!(
        cases,
        (name="normal_cpu", interception_backend=:raster_cpu, scattering_backend=nothing),
    )
    selected("rastergpu_ka_cpu") && push!(
        cases,
        _rastergpu_case("rastergpu_ka_cpu", KernelAbstractions.CPU()),
    )
    if selected("rastergpu_metal_gpu")
        metal = _optional_metal_backend()
        if metal === nothing
            @warn "RasterGPU Metal benchmark requested but Metal is unavailable; skipping it."
        else
            push!(cases, _rastergpu_case("rastergpu_metal_gpu", metal))
        end
    end
    isempty(cases) && error("No benchmark backend is available.")
    return cases
end

function _step_totals(step)
    totals = ArchimedLight._light_step_energy_totals(step)
    return (
        incident_par=Float64(totals.incident_par_total),
        incident_nir=Float64(totals.incident_nir_total),
        absorbed_par=Float64(totals.absorbed_par_total),
        absorbed_nir=Float64(totals.absorbed_nir_total),
    )
end

function _measure_case(scene, models, row, options, case)
    for _ in 1:BENCH_WARMUPS
        sim = ArchimedLight.LightSimulation(
            scene,
            models;
            options=options,
            interception_backend=case.interception_backend,
            scattering_backend=case.scattering_backend,
        )
        ArchimedLight.run_light(sim, row)
    end

    seconds = Float64[]
    allocations = Float64[]
    totals = nothing
    resolved_backend = ""
    for _ in 1:BENCH_SAMPLES
        GC.gc()
        sim = ArchimedLight.LightSimulation(
            scene,
            models;
            options=options,
            interception_backend=case.interception_backend,
            scattering_backend=case.scattering_backend,
        )
        timed = @timed ArchimedLight.run_light(sim, row)
        push!(seconds, timed.time)
        push!(allocations, timed.bytes / 1024.0^2)
        totals = _step_totals(timed.value)
        resolved_backend = sim.cache === nothing ? "" : string(nameof(typeof(sim.cache.resolved_interception_backend)))
    end
    return (
        backend=case.name,
        resolved_backend=resolved_backend,
        median_seconds=median(seconds),
        min_seconds=minimum(seconds),
        median_alloc_mib=median(allocations),
        totals...,
    )
end

function _write_rows(rows, path::AbstractString)
    isempty(path) && return rows
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), '\t'))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), '\t'))
        end
    end
    return rows
end

function benchmark_configs(
    config_paths::Vector{String};
    directions::Vector{Union{Nothing,Int}}=Union{Nothing,Int}[nothing],
    output::AbstractString=BENCH_OUTPUT,
)
    BENCH_SAMPLES > 0 || error("ARCHIMEDLIGHT_BENCH_SAMPLES must be positive.")
    cases = _backend_cases()
    rows = NamedTuple[]
    for config_path in config_paths
        isfile(config_path) || error("Benchmark config does not exist: $config_path")
        base_options, scene, meteo, models = ArchimedLight.read_config(config_path)
        for direction_count in directions
            options = direction_count === nothing ? base_options :
                      ArchimedLight.LightOptions(base_options; turtle_sectors=direction_count)
            prepared_meteo = ArchimedLight.prepare_meteo(meteo, options)
            row = first(prepared_meteo.rows)
            for case in cases
                measurement = _measure_case(scene, models, row, options, case)
                push!(
                    rows,
                    merge(
                        (
                            config=basename(dirname(config_path)),
                            directions=options.turtle_sectors,
                            pixel_size=options.pixel_size,
                            scattering=options.scattering,
                        ),
                        measurement,
                    ),
                )
            end
        end
    end
    _write_rows(rows, output)
    return rows
end

function backend_comparison_main()
    rows = benchmark_configs([BENCH_CONFIG])
    println("ArchimedLight backend comparison")
    @printf("%-24s %-24s %10s %10s %10s\n", "backend", "resolved", "median s", "min s", "alloc MiB")
    for row in rows
        @printf(
            "%-24s %-24s %10.4f %10.4f %10.1f\n",
            row.backend,
            row.resolved_backend,
            row.median_seconds,
            row.min_seconds,
            row.median_alloc_mib,
        )
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    backend_comparison_main()
end
