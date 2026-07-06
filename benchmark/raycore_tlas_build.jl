const _TLAS_BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_METAL", ""))
const _TLAS_BENCH_WANTS_METAL = _TLAS_BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _TLAS_BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using Adapt
using KernelAbstractions

const TLAS_BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_SAMPLES", "5"))
const TLAS_BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_WARMUPS", "2"))
const TLAS_BENCH_BACKENDS = [strip(lowercase(v)) for v in split(get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_BACKENDS", "raycore_ka_cpu,raycore_metal_gpu"), ',') if !isempty(strip(v))]
const TLAS_BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_OUTPUT", "")

function _flag_required(name::AbstractString)
    return lowercase(get(ENV, name, "")) in ("required", "force", "error")
end

function _optional_metal_backend()
    _TLAS_BENCH_WANTS_METAL || return nothing
    if !Sys.isapple() || Sys.ARCH != :aarch64
        msg = "Metal TLAS benchmark requires Apple Silicon macOS."
        _flag_required("ARCHIMEDLIGHT_TLAS_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end
    if Base.find_package("Metal") === nothing
        msg = "Metal.jl is not available. Run this benchmark with an environment that includes Metal."
        _flag_required("ARCHIMEDLIGHT_TLAS_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _benchmark_workload()
    config = get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_CONFIG", joinpath(@__DIR__, "..", "example_2", "config.yml"))
    paving = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_PAVING", "64"))
    pixel_size = parse(Float64, get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_PIXEL_SIZE", "0.04"))
    sectors = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_BENCH_SECTORS", "6"))
    options, scene, _meteo, models = ArchimedLight.read_config(config; plot_paving_override=paving)
    options = ArchimedLight.LightOptions(
        options;
        cache_radiation=true,
        scattering=true,
        meteo_range=nothing,
        pixel_size=pixel_size,
        turtle_sectors=sectors,
    )
    prepared = ArchimedLight._prepare_interception_data(
        scene,
        models,
        options;
        include_budget_maps=true,
        include_raycore_instancing=true,
    )
    return (; config, options, prepared)
end

function _backend_cases()
    cases = Dict{String,Any}(
        "raycore_ka_cpu" => KernelAbstractions.CPU(),
    )
    metal = _optional_metal_backend()
    metal !== nothing && (cases["raycore_metal_gpu"] = metal)
    selected = [(name=name, backend=cases[name]) for name in TLAS_BENCH_BACKENDS if haskey(cases, name)]
    isempty(selected) && error(
        "ARCHIMEDLIGHT_TLAS_BENCH_BACKENDS selected no available backends. " *
        "Use raycore_ka_cpu or raycore_metal_gpu; Metal also requires ARCHIMEDLIGHT_TLAS_BENCH_METAL=1.",
    )
    return selected
end

function _tlas_source(data)
    return data.tlas.backend isa KernelAbstractions.CPU ? :native_cpu : :native_device
end

function _measure_tlas_build(prepared, options, backend)
    config = ArchimedLight.RaycoreBackendConfig(backend=backend, validate=false)
    GC.gc()
    started = time_ns()
    data = ArchimedLight._raycore_scene_data(prepared, config; toricity=options.toricity)
    elapsed_ms = (time_ns() - started) / 1e6
    source = _tlas_source(data)
    if !(backend isa KernelAbstractions.CPU) && source != :native_device
        error("Raycore TLAS benchmark expected native device TLAS, got $source.")
    end
    return (
        build_ms=elapsed_ms,
        tlas_source=source,
        geometry_mode=data.geometry_mode,
        chunked_tlas=data.chunked_tlas,
        tlas_instances=length(data.tlas.instances),
        tlas_geometries=Raycore.n_geometries(data.tlas),
        static_tlas_ready=data.tlas.static_tlas !== nothing,
        kernel_tlas_same_object=data.kernel_tlas === data.tlas.static_tlas,
    )
end

function _measure_legacy_cpu_static_tlas_build(prepared, options, backend)
    config = ArchimedLight.RaycoreBackendConfig(backend=backend, validate=false)
    GC.gc()
    started = time_ns()
    cpu_config = ArchimedLight._raycore_config_with(config; backend=KernelAbstractions.CPU())
    cpu_data = ArchimedLight._raycore_scene_data(prepared, cpu_config; toricity=options.toricity)
    kernel_tlas = Adapt.adapt(backend, cpu_data.tlas.static_tlas)
    elapsed_ms = (time_ns() - started) / 1e6
    return (
        build_ms=elapsed_ms,
        tlas_source=:legacy_cpu_static,
        geometry_mode=cpu_data.geometry_mode,
        chunked_tlas=cpu_data.chunked_tlas,
        tlas_instances=length(cpu_data.tlas.instances),
        tlas_geometries=Raycore.n_geometries(cpu_data.tlas),
        static_tlas_ready=cpu_data.tlas.static_tlas !== nothing,
        kernel_tlas_same_object=kernel_tlas === cpu_data.tlas.static_tlas,
    )
end

function _print_rows(io, rows)
    fields = (
        :backend,
        :mode,
        :sample,
        :build_ms,
        :tlas_source,
        :geometry_mode,
        :chunked_tlas,
        :tlas_instances,
        :tlas_geometries,
        :static_tlas_ready,
        :kernel_tlas_same_object,
    )
    println(io, join(string.(fields), '\t'))
    for row in rows
        println(io, join((string(getproperty(row, field)) for field in fields), '\t'))
    end
end

function main()
    workload = _benchmark_workload()
    cases = _backend_cases()
    @info "Raycore TLAS build benchmark" config=workload.config samples=TLAS_BENCH_SAMPLES warmups=TLAS_BENCH_WARMUPS backends=[case.name for case in cases]
    rows = NamedTuple[]
    for case in cases
        modes =
            case.backend isa KernelAbstractions.CPU ?
            ((name=:native, measure=_measure_tlas_build),) :
            (
                (name=:legacy_cpu_static, measure=_measure_legacy_cpu_static_tlas_build),
                (name=:native_device, measure=_measure_tlas_build),
            )
        for mode in modes
            for _ in 1:TLAS_BENCH_WARMUPS
                mode.measure(workload.prepared, workload.options, case.backend)
            end
            for sample in 1:TLAS_BENCH_SAMPLES
                result = mode.measure(workload.prepared, workload.options, case.backend)
                row = merge((backend=case.name, mode=mode.name, sample=sample), result)
                push!(rows, row)
                @info "Raycore TLAS build sample" row
            end
        end
    end

    if isempty(TLAS_BENCH_OUTPUT)
        _print_rows(stdout, rows)
    else
        open(TLAS_BENCH_OUTPUT, "w") do io
            _print_rows(io, rows)
        end
        @info "Wrote Raycore TLAS build benchmark results" path=TLAS_BENCH_OUTPUT rows=length(rows)
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
