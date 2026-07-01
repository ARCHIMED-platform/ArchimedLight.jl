const _CACHE_BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_METAL", ""))
const _CACHE_BENCH_WANTS_METAL = _CACHE_BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _CACHE_BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using KernelAbstractions

const CACHE_BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_SAMPLES", "3"))
const CACHE_BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_WARMUPS", "1"))
const CACHE_BENCH_STEPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_STEPS", "4"))
const CACHE_BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_OUTPUT", "")

function _cache_required_metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        error("Metal static-cache benchmark requires Apple Silicon macOS.")
    end
    if Base.find_package("Metal") === nothing
        error("Metal.jl is not available. Run this benchmark with an environment that includes Metal.")
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _cache_workload()
    config = get(
        ENV,
        "ARCHIMEDLIGHT_CACHE_BENCH_CONFIG",
        joinpath(@__DIR__, "..", "test", "fast_fixtures", "simpleplant_16_toric", "input", "config.yml"),
    )
    paving_flag = get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_PAVING", "")
    options, scene, meteo, models =
        isempty(paving_flag) ?
        ArchimedLight.read_config(config) :
        ArchimedLight.read_config(config; plot_paving_override=parse(Int, paving_flag))
    options = ArchimedLight.LightOptions(
        options;
        cache_radiation=true,
        scattering=true,
        meteo_range=nothing,
        pixel_size=parse(Float64, get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_PIXEL_SIZE", "0.40")),
        turtle_sectors=parse(Int, get(ENV, "ARCHIMEDLIGHT_CACHE_BENCH_SECTORS", "6")),
    )
    rows = ArchimedLight.prepare_meteo(meteo, options).rows
    if isempty(rows)
        sky = ArchimedLight.SkyState(180.0, 55.0, 420.0, 220.0, 0.70, 0.30)
        meteo_rows = NamedTuple[]
    else
        meteo_rows = collect(Iterators.take(Iterators.cycle(rows), CACHE_BENCH_STEPS))
        sky = ArchimedLight.compute_sky(first(meteo_rows), options)
    end
    return (; config, scene, models, options, sky, meteo_rows)
end

function _cache_backend(backend)
    interception = ArchimedLight.RaycoreInterceptionBackend(
        backend=backend,
        max_hits_per_pixel=64,
        validate=true,
    )
    return interception, ArchimedLight.RaycoreScatteringBackend(interception)
end

function _new_sim(workload, backend)
    interception, scattering = _cache_backend(backend)
    return ArchimedLight.LightSimulation(
        workload.scene,
        workload.models;
        options=workload.options,
        interception_backend=interception,
        scattering_backend=scattering,
    )
end

function _raycore_id(sim)
    sim.cache === nothing && return UInt(0)
    sim.cache.raycore_data === nothing && return UInt(0)
    return objectid(sim.cache.raycore_data)
end

function _run_series_rows!(sim, rows)
    out = Vector{Any}(undef, length(rows))
    for i in eachindex(rows)
        out[i] = ArchimedLight.run_light(sim, rows[i])
    end
    return out
end

function _measure_cache_lifecycle(workload, backend)
    sim = _new_sim(workload, backend)

    GC.gc()
    started = time_ns()
    cache = ArchimedLight._ensure_light_cache!(sim)
    cold_prepare_ms = (time_ns() - started) / 1e6
    prepared_raycore_id = _raycore_id(sim)
    entries_after_prepare = length(cache.entries)

    started = time_ns()
    ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=900.0)
    warm_single_populate_ms = (time_ns() - started) / 1e6
    id_after_single_populate = _raycore_id(sim)
    entries_after_single_populate = length(cache.entries)

    started = time_ns()
    ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=900.0)
    warm_single_replay_ms = (time_ns() - started) / 1e6
    id_after_single_replay = _raycore_id(sim)
    entries_after_single_replay = length(cache.entries)

    multi_populate_ms =
        if isempty(workload.meteo_rows)
            0.0
        else
            started = time_ns()
            _run_series_rows!(sim, workload.meteo_rows)
            (time_ns() - started) / 1e6
        end
    id_after_multi_populate = _raycore_id(sim)
    entries_after_multi_populate = length(cache.entries)

    multi_replay_ms =
        if isempty(workload.meteo_rows)
            0.0
        else
            started = time_ns()
            _run_series_rows!(sim, workload.meteo_rows)
            (time_ns() - started) / 1e6
        end
    id_after_multi_replay = _raycore_id(sim)
    entries_after_multi_replay = length(cache.entries)

    rebuild_sim = _new_sim(workload, backend)
    started = time_ns()
    ArchimedLight._ensure_light_cache!(rebuild_sim)
    rebuild_ms = (time_ns() - started) / 1e6

    same_raycore_single = prepared_raycore_id == id_after_single_populate == id_after_single_replay
    same_raycore_multi = prepared_raycore_id == id_after_multi_populate == id_after_multi_replay
    return (
        backend=:raycore_metal_gpu,
        cold_prepare_ms=cold_prepare_ms,
        warm_single_populate_ms=warm_single_populate_ms,
        warm_single_replay_ms=warm_single_replay_ms,
        multi_populate_ms=multi_populate_ms,
        multi_replay_ms=multi_replay_ms,
        rebuild_ms=rebuild_ms,
        steps=length(workload.meteo_rows),
        warm_single_replay_per_step_ms=warm_single_replay_ms,
        multi_populate_per_step_ms=isempty(workload.meteo_rows) ? 0.0 : multi_populate_ms / length(workload.meteo_rows),
        multi_replay_per_step_ms=isempty(workload.meteo_rows) ? 0.0 : multi_replay_ms / length(workload.meteo_rows),
        entries_after_prepare=entries_after_prepare,
        entries_after_single_populate=entries_after_single_populate,
        entries_after_single_replay=entries_after_single_replay,
        entries_after_multi_populate=entries_after_multi_populate,
        entries_after_multi_replay=entries_after_multi_replay,
        same_raycore_single=same_raycore_single,
        same_raycore_multi=same_raycore_multi,
        raycore_tlas_source=cache.raycore_data === nothing ? :none :
                            cache.raycore_data.tlas.backend isa KernelAbstractions.CPU ? :native_cpu : :native_device,
        cache_mode=cache.mode,
    )
end

function _print_rows(io, rows)
    fields = (
        :backend,
        :sample,
        :cold_prepare_ms,
        :warm_single_populate_ms,
        :warm_single_replay_ms,
        :multi_populate_ms,
        :multi_replay_ms,
        :rebuild_ms,
        :steps,
        :warm_single_replay_per_step_ms,
        :multi_populate_per_step_ms,
        :multi_replay_per_step_ms,
        :entries_after_prepare,
        :entries_after_single_populate,
        :entries_after_single_replay,
        :entries_after_multi_populate,
        :entries_after_multi_replay,
        :same_raycore_single,
        :same_raycore_multi,
        :raycore_tlas_source,
        :cache_mode,
    )
    println(io, join(string.(fields), '\t'))
    for row in rows
        println(io, join((string(getproperty(row, field)) for field in fields), '\t'))
    end
end

function main()
    workload = _cache_workload()
    backend = _cache_required_metal_backend()
    @info "Raycore static-cache lifecycle benchmark" config=workload.config samples=CACHE_BENCH_SAMPLES warmups=CACHE_BENCH_WARMUPS steps=CACHE_BENCH_STEPS
    for _ in 1:CACHE_BENCH_WARMUPS
        _measure_cache_lifecycle(workload, backend)
    end

    rows = NamedTuple[]
    for sample in 1:CACHE_BENCH_SAMPLES
        row = merge((sample=sample,), _measure_cache_lifecycle(workload, backend))
        push!(rows, row)
        @info "Raycore static-cache lifecycle sample" row
    end

    if isempty(CACHE_BENCH_OUTPUT)
        _print_rows(stdout, rows)
    else
        open(CACHE_BENCH_OUTPUT, "w") do io
            _print_rows(io, rows)
        end
        @info "Wrote Raycore static-cache lifecycle benchmark results" path=CACHE_BENCH_OUTPUT rows=length(rows)
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
