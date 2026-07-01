using ArchimedLight
using KernelAbstractions
using Metal
using Raycore

const PREP_REUSE_SCENE_PATH = get(
    ENV,
    "ARCHIMEDLIGHT_PREP_REUSE_SCENE",
    "/Users/rvezy/Documents/cirad/Travail_AMAP/Projets/PalmStudio/simulations/Light_interception_vs_ages/archimed_simulations/sim.ops",
)
const PREP_REUSE_PROGRESS_OUTPUT = get(ENV, "ARCHIMEDLIGHT_PREP_REUSE_PROGRESS_OUTPUT", "")
const PREP_REUSE_ITERATIONS = parse(Int, get(ENV, "ARCHIMEDLIGHT_PREP_REUSE_ITERATIONS", "2"))

function _prep_reuse_progress!(event; kwargs...)
    isempty(PREP_REUSE_PROGRESS_OUTPUT) && return nothing
    open(PREP_REUSE_PROGRESS_OUTPUT, "a") do io
        print(io, time(), '\t', event)
        for (key, value) in kwargs
            print(io, '\t', key, '=', value)
        end
        println(io)
    end
    return nothing
end

function _prep_reuse_required_metal_backend()
    Sys.isapple() && Sys.ARCH == :aarch64 || error("Metal requires Apple Silicon macOS.")
    dev = Metal.MtlArray(zeros(Float32, 1))
    return KernelAbstractions.get_backend(dev)
end

function _prep_reuse_model_pairs(scene)
    pairs = Set{Tuple{String,String}}()
    node_ids = unique(scene.face2node)
    node_group, node_type = ArchimedLight._scene_group_type_maps(scene, node_ids)
    for nid in node_ids
        group = strip(get(node_group, nid, ""))
        type_name = strip(get(node_type, nid, ""))
        push!(pairs, (group, type_name))
    end
    return sort!(collect(pairs))
end

function _prep_reuse_default_models_for_scene(scene)
    by_group = Dict{String,Vector{Pair{String,Any}}}()
    for (group, type_name) in _prep_reuse_model_pairs(scene)
        material =
            occursin("soil", lowercase(group)) || occursin("ground", lowercase(group)) ||
            occursin("pavement", lowercase(group)) ?
            ArchimedLight.translucent(par=0.12, nir=0.60) :
            ArchimedLight.translucent(par=0.15, nir=0.90)
        push!(get!(by_group, group, Pair{String,Any}[]), type_name => material)
    end
    return ArchimedLight.models_for((group => Tuple(entries) for (group, entries) in by_group)...)
end

function _prep_reuse_elapsed!(f::F, event; kwargs...) where {F}
    started = time_ns()
    result = try
        f()
    finally
        _prep_reuse_progress!(event; ms=(time_ns() - started) / 1e6, kwargs...)
    end
    return result
end

function _prep_reuse_shape(data)
    data === nothing && return (
        tlas_source=:none,
        geometry_mode=:none,
        tlas_instances=0,
        tlas_geometries=0,
        max_blas_depth=0,
    )
    source = data.backend isa KernelAbstractions.CPU ? :native_cpu :
             data.tlas.backend isa KernelAbstractions.CPU ? :cpu_built_static : :native_device
    traversal = ArchimedLight._raycore_traversal_stack_depth_status(data)
    return (
        tlas_source=source,
        geometry_mode=data.geometry_mode,
        tlas_instances=length(data.tlas.instances),
        tlas_geometries=Raycore.n_geometries(data.tlas),
        max_blas_depth=traversal.blas_depth,
    )
end

function main()
    scene = _prep_reuse_elapsed!(:scene_read_done; path=PREP_REUSE_SCENE_PATH) do
        ArchimedLight.read_scene(PREP_REUSE_SCENE_PATH)
    end
    models = _prep_reuse_elapsed!(:models_done) do
        _prep_reuse_default_models_for_scene(scene)
    end
    options = ArchimedLight.LightOptions(
        turtle_sectors=16,
        pixel_size=0.05,
        scattering=true,
        cache_radiation=false,
        all_in_turtle=true,
        toricity=true,
        scattering_max_iter=12,
    )
    ib = ArchimedLight.RaycoreInterceptionBackend(
        backend=_prep_reuse_required_metal_backend(),
        max_hits_per_pixel=64,
        workgroupsize=256,
        validate=false,
    )
    sb = ArchimedLight.RaycoreScatteringBackend(ib)

    rows = NamedTuple[]
    for iteration in 1:PREP_REUSE_ITERATIONS
        GC.gc()
        sim = ArchimedLight.LightSimulation(
            scene,
            models;
            options=options,
            interception_backend=ib,
            scattering_backend=sb,
        )
        report = ArchimedLight.check_simulation(sim)
        isempty(report.errors) || error("Invalid benchmark simulation: $(join(report.errors, "; "))")
        sim.validation = report
        elapsed = @elapsed begin
            sim.cache = ArchimedLight.prepare_light_cache(
                sim.scene,
                sim.models,
                sim.options;
                interception_backend=sim.interception_backend,
                scattering_mode=sim.scattering_mode,
                scattering_backend=sim.scattering_backend,
                memory_limit_bytes=sim.memory_limit_bytes,
            )
        end
        shape = _prep_reuse_shape(sim.cache.raycore_data)
        row = (
            iteration=iteration,
            prepare_ms=1000elapsed,
            face_chunk_limit=ArchimedLight._raycore_face_chunk_limit(),
            shape...,
        )
        push!(rows, row)
        _prep_reuse_progress!(:prepare_done; row...)
    end
    return rows
end

main()
