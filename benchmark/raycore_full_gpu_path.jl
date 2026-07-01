const _FULL_BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_METAL", ""))
const _FULL_BENCH_WANTS_METAL = _FULL_BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _FULL_BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using KernelAbstractions

const FULL_GPU_BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_SAMPLES", "5"))
const FULL_GPU_BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_WARMUPS", "2"))
const FULL_GPU_BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_OUTPUT", "")

function _required_metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        error("Metal full GPU path benchmark requires Apple Silicon macOS.")
    end
    if Base.find_package("Metal") === nothing
        error("Metal.jl is not available. Run this benchmark with an environment that includes Metal.")
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _full_gpu_workload()
    config = get(
        ENV,
        "ARCHIMEDLIGHT_FULL_GPU_BENCH_CONFIG",
        joinpath(@__DIR__, "..", "test", "fast_fixtures", "simpleplant_16_toric", "input", "config.yml"),
    )
    options, scene, _meteo, models = ArchimedLight.read_config(config)
    options = ArchimedLight.LightOptions(
        options;
        cache_radiation=true,
        scattering=true,
        pixel_size=parse(Float64, get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_PIXEL_SIZE", "0.40")),
        turtle_sectors=parse(Int, get(ENV, "ARCHIMEDLIGHT_FULL_GPU_BENCH_SECTORS", "6")),
    )
    sky = ArchimedLight.SkyState(180.0, 55.0, 420.0, 220.0, 0.70, 0.30)
    return (; config, scene, models, options, sky)
end

function _cache_tlas_source(cache)
    data = cache.raycore_data
    data === nothing && return :none
    data.backend isa KernelAbstractions.CPU && return :native_cpu
    return data.tlas.backend isa KernelAbstractions.CPU ? :cpu_built_static : :native_device
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

function _measure_full_gpu_path(workload, backend)
    interception = ArchimedLight.RaycoreInterceptionBackend(
        backend=backend,
        max_hits_per_pixel=64,
        validate=true,
    )
    GC.gc()
    started = time_ns()
    sim = ArchimedLight.LightSimulation(
        workload.scene,
        workload.models;
        options=workload.options,
        interception_backend=interception,
        scattering_backend=ArchimedLight.RaycoreScatteringBackend(interception),
    )
    step = ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=900.0)
    elapsed_ms = (time_ns() - started) / 1e6
    source = _cache_tlas_source(sim.cache)
    return (
        full_run_ms=elapsed_ms,
        cache_backend=string(typeof(sim.cache.resolved_interception_backend)),
        cache_tlas_source=source,
        geometry_mode=sim.cache.raycore_data.geometry_mode,
        chunked_tlas=sim.cache.raycore_data.chunked_tlas,
        totals=_step_totals(step),
    )
end

function _print_rows(io, rows)
    fields = (
        :backend,
        :sample,
        :full_run_ms,
        :cache_backend,
        :cache_tlas_source,
        :geometry_mode,
        :chunked_tlas,
        :incident_par,
        :incident_nir,
        :absorbed_par,
        :absorbed_nir,
    )
    println(io, join(string.(fields), '\t'))
    for row in rows
        println(io, join((string(getproperty(row, field)) for field in fields), '\t'))
    end
end

function main()
    workload = _full_gpu_workload()
    backend = _required_metal_backend()
    @info "Raycore full GPU path benchmark" config=workload.config samples=FULL_GPU_BENCH_SAMPLES warmups=FULL_GPU_BENCH_WARMUPS
    for _ in 1:FULL_GPU_BENCH_WARMUPS
        _measure_full_gpu_path(workload, backend)
    end

    rows = NamedTuple[]
    for sample in 1:FULL_GPU_BENCH_SAMPLES
        result = _measure_full_gpu_path(workload, backend)
        flat = merge(
            (
                backend=:raycore_metal_gpu,
                sample=sample,
                full_run_ms=result.full_run_ms,
                cache_backend=result.cache_backend,
                cache_tlas_source=result.cache_tlas_source,
                geometry_mode=result.geometry_mode,
                chunked_tlas=result.chunked_tlas,
            ),
            result.totals,
        )
        push!(rows, flat)
        @info "Raycore full GPU path sample" flat
    end

    if isempty(FULL_GPU_BENCH_OUTPUT)
        _print_rows(stdout, rows)
    else
        open(FULL_GPU_BENCH_OUTPUT, "w") do io
            _print_rows(io, rows)
        end
        @info "Wrote Raycore full GPU path benchmark results" path=FULL_GPU_BENCH_OUTPUT rows=length(rows)
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
