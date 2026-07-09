const _STAGE_BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_METAL", ""))
const _STAGE_BENCH_WANTS_METAL = _STAGE_BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _STAGE_BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using KernelAbstractions

const TLAS_STAGE_BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_SAMPLES", "5"))
const TLAS_STAGE_BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_WARMUPS", "2"))
const TLAS_STAGE_BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_OUTPUT", "")

function _stage_required_metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        error("Metal TLAS stage benchmark requires Apple Silicon macOS.")
    end
    if Base.find_package("Metal") === nothing
        error("Metal.jl is not available. Run this benchmark with an environment that includes Metal.")
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _stage_workload()
    config = get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_CONFIG", joinpath(@__DIR__, "..", "example_2", "config.yml"))
    paving = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_PAVING", "64"))
    pixel_size = parse(Float64, get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_PIXEL_SIZE", "0.04"))
    sectors = parse(Int, get(ENV, "ARCHIMEDLIGHT_TLAS_STAGE_BENCH_SECTORS", "6"))
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

function _measure_stages(workload, backend)
    prepared = workload.prepared
    config = ArchimedLight.RaycoreBackendConfig(backend=backend, validate=false)
    GC.gc()

    started = time_ns()
    tlas = Raycore.TLAS(config.backend)
    construct_ms = (time_ns() - started) / 1e6

    started = time_ns()
    mesh = ArchimedLight._raycore_mesh_from_geometry(prepared.geometry)
    mesh_ms = (time_ns() - started) / 1e6

    started = time_ns()
    if workload.options.toricity
        push!(tlas, mesh, ArchimedLight._raycore_toric_transforms(prepared.geometry; toricity=true))
    else
        push!(tlas, mesh)
    end
    push_mesh_ms = (time_ns() - started) / 1e6

    started = time_ns()
    Raycore.sync!(tlas)
    sync_ms = (time_ns() - started) / 1e6

    started = time_ns()
    data = ArchimedLight._raycore_scene_data(prepared, config; toricity=workload.options.toricity)
    scene_data_ms = (time_ns() - started) / 1e6

    return (
        backend=backend isa KernelAbstractions.CPU ? :raycore_ka_cpu : :raycore_metal_gpu,
        construct_ms=construct_ms,
        mesh_ms=mesh_ms,
        push_mesh_ms=push_mesh_ms,
        sync_ms=sync_ms,
        scene_data_ms=scene_data_ms,
        tlas_source=data.tlas.backend isa KernelAbstractions.CPU ? :native_cpu : :native_device,
        geometry_mode=data.geometry_mode,
        tlas_instances=length(data.tlas.instances),
        tlas_geometries=Raycore.n_geometries(data.tlas),
    )
end

function _print_rows(io, rows)
    fields = (
        :backend,
        :sample,
        :construct_ms,
        :mesh_ms,
        :push_mesh_ms,
        :sync_ms,
        :scene_data_ms,
        :tlas_source,
        :geometry_mode,
        :tlas_instances,
        :tlas_geometries,
    )
    println(io, join(string.(fields), '\t'))
    for row in rows
        println(io, join((string(getproperty(row, field)) for field in fields), '\t'))
    end
end

function main()
    workload = _stage_workload()
    backend = _stage_required_metal_backend()
    @info "Raycore TLAS stage benchmark" config=workload.config samples=TLAS_STAGE_BENCH_SAMPLES warmups=TLAS_STAGE_BENCH_WARMUPS
    for _ in 1:TLAS_STAGE_BENCH_WARMUPS
        _measure_stages(workload, backend)
    end

    rows = NamedTuple[]
    for sample in 1:TLAS_STAGE_BENCH_SAMPLES
        row = merge((sample=sample,), _measure_stages(workload, backend))
        push!(rows, row)
        @info "Raycore TLAS stage sample" row
    end

    if isempty(TLAS_STAGE_BENCH_OUTPUT)
        _print_rows(stdout, rows)
    else
        open(TLAS_STAGE_BENCH_OUTPUT, "w") do io
            _print_rows(io, rows)
        end
        @info "Wrote Raycore TLAS stage benchmark results" path=TLAS_STAGE_BENCH_OUTPUT rows=length(rows)
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
