using ArchimedLight
using KernelAbstractions
using Metal
using Raycore

const STEP_BREAKDOWN_SCENE_PATH = get(
    ENV,
    "ARCHIMEDLIGHT_STEP_BREAKDOWN_SCENE",
    "/Users/rvezy/Documents/cirad/Travail_AMAP/Projets/PalmStudio/simulations/Light_interception_vs_ages/archimed_simulations/sim.ops",
)
const STEP_BREAKDOWN_PROGRESS_OUTPUT = get(ENV, "ARCHIMEDLIGHT_STEP_BREAKDOWN_PROGRESS_OUTPUT", "")
const STEP_BREAKDOWN_ITERATIONS = parse(Int, get(ENV, "ARCHIMEDLIGHT_STEP_BREAKDOWN_ITERATIONS", "2"))

function _step_breakdown_progress!(event; kwargs...)
    isempty(STEP_BREAKDOWN_PROGRESS_OUTPUT) && return nothing
    open(STEP_BREAKDOWN_PROGRESS_OUTPUT, "a") do io
        print(io, time(), '\t', event)
        for (key, value) in kwargs
            print(io, '\t', key, '=', value)
        end
        println(io)
    end
    return nothing
end

function _step_breakdown_elapsed!(f::F, event; kwargs...) where {F}
    started = time_ns()
    result = try
        f()
    finally
        _step_breakdown_progress!(event; ms=(time_ns() - started) / 1e6, kwargs...)
    end
    return result
end

function _step_breakdown_required_metal_backend()
    Sys.isapple() && Sys.ARCH == :aarch64 || error("Metal requires Apple Silicon macOS.")
    dev = Metal.MtlArray(zeros(Float32, 1))
    return KernelAbstractions.get_backend(dev)
end

function _step_breakdown_model_pairs(scene)
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

function _step_breakdown_default_models_for_scene(scene)
    by_group = Dict{String,Vector{Pair{String,Any}}}()
    for (group, type_name) in _step_breakdown_model_pairs(scene)
        material =
            occursin("soil", lowercase(group)) || occursin("ground", lowercase(group)) ||
            occursin("pavement", lowercase(group)) ?
            ArchimedLight.translucent(par=0.12, nir=0.60) :
            ArchimedLight.translucent(par=0.15, nir=0.90)
        push!(get!(by_group, group, Pair{String,Any}[]), type_name => material)
    end
    return ArchimedLight.models_for((group => Tuple(entries) for (group, entries) in by_group)...)
end

function _step_breakdown_totals(first, scat)
    dense_first = first.dense
    dense_first === nothing && error("Raycore step breakdown expects dense first-order results.")
    node_ids = dense_first.node_ids
    first_par = dense_first.incident_power.par
    first_nir = dense_first.incident_power.nir
    scat_par = scat === nothing || scat.dense === nothing ? zeros(Float64, length(node_ids)) : scat.dense.added_power.par
    scat_nir = scat === nothing || scat.dense === nothing ? zeros(Float64, length(node_ids)) : scat.dense.added_power.nir
    return (
        first_par=sum(first_par),
        first_nir=sum(first_nir),
        scattering_par=sum(scat_par),
        scattering_nir=sum(scat_nir),
    )
end

function main()
    scene = _step_breakdown_elapsed!(:scene_read_done; path=STEP_BREAKDOWN_SCENE_PATH) do
        ArchimedLight.read_scene(STEP_BREAKDOWN_SCENE_PATH)
    end
    models = _step_breakdown_elapsed!(:models_done) do
        _step_breakdown_default_models_for_scene(scene)
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
    sky = ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30)
    ib = ArchimedLight.RaycoreInterceptionBackend(
        backend=_step_breakdown_required_metal_backend(),
        max_hits_per_pixel=64,
        workgroupsize=256,
        validate=false,
    )
    sb = ArchimedLight.RaycoreScatteringBackend(ib)
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
    sim.cache = _step_breakdown_elapsed!(:prepare_light_cache_done) do
        ArchimedLight.prepare_light_cache(
            sim.scene,
            sim.models,
            sim.options;
            interception_backend=sim.interception_backend,
            scattering_mode=sim.scattering_mode,
            scattering_backend=sim.scattering_backend,
            memory_limit_bytes=sim.memory_limit_bytes,
        )
    end
    data = sim.cache.raycore_data
    shape = ArchimedLight._raycore_scene_shape_summary(data)
    _step_breakdown_progress!(
        :raycore_shape;
        geometry_mode=shape.geometry_mode,
        tlas_instances=shape.tlas_instances,
        tlas_geometries=shape.tlas_geometries,
        max_blas_depth=shape.max_blas_depth,
        face_chunk_limit=ArchimedLight._raycore_face_chunk_limit(),
    )

    profile_turtle = ArchimedLight.build_turtle(options, sky)
    profile_fluxes = ArchimedLight.compute_directional_fluxes(sky, profile_turtle, options)
    profile_accumulate_edges = [sector.source != :sun for sector in profile_turtle.sectors]
    profile_needs_sector_area = [
        profile_fluxes.par[i] != 0.0 || profile_fluxes.nir[i] != 0.0
        for i in eachindex(profile_turtle.sectors)
    ]
    stack_profile = _step_breakdown_elapsed!(:raycore_stack_profile_done) do
        ArchimedLight._raycore_stack_profile(
            data,
            [sector.direction for sector in profile_turtle.sectors],
            options;
            accumulate_edges=profile_accumulate_edges,
            needs_sector_area=profile_needs_sector_area,
        )
    end
    _step_breakdown_progress!(:raycore_stack_profile; stack_profile...)

    for iteration in 1:STEP_BREAKDOWN_ITERATIONS
        GC.gc()
        turtle = _step_breakdown_elapsed!(:build_turtle_done; iteration=iteration) do
            ArchimedLight.build_turtle(options, sky)
        end
        fluxes = _step_breakdown_elapsed!(:compute_fluxes_done; iteration=iteration) do
            ArchimedLight.compute_directional_fluxes(sky, turtle, options)
        end
        first, topology = _step_breakdown_elapsed!(:raycore_first_order_topology_done; iteration=iteration) do
            ArchimedLight._stream_first_order_with_raycore_scattering_topology(data, scene, models, turtle, fluxes, options)
        end
        scat = _step_breakdown_elapsed!(:raycore_scattering_done; iteration=iteration) do
            ArchimedLight._compute_scattering_with_flags(
                scene,
                models,
                turtle,
                first,
                options;
                mode=sim.cache.scattering_mode,
                backend=sim.cache.scattering_backend,
                scattering_topology=topology,
                raycore_data=data,
            )
        end
        budget = _step_breakdown_elapsed!(:integrate_light_done; iteration=iteration) do
            ArchimedLight.integrate_light(
                scene,
                models,
                first,
                scat,
                options;
                step_duration_seconds=900.0,
                component_area_per_node=data.prepared.component_area_per_node,
                absorption_par_per_node=data.prepared.absorption_par_per_node,
                absorption_nir_per_node=data.prepared.absorption_nir_per_node,
            )
        end
        totals = _step_breakdown_totals(first, scat)
        _step_breakdown_progress!(
            :iteration_summary;
            iteration=iteration,
            first_par=totals.first_par,
            first_nir=totals.first_nir,
            scattering_par=totals.scattering_par,
            scattering_nir=totals.scattering_nir,
            budget_nodes=length(first.dense.node_ids),
        )
    end
    warm_profile = _step_breakdown_elapsed!(:raycore_warm_stack_profile_done) do
        ArchimedLight._raycore_stack_profile(
            data,
            [sector.direction for sector in profile_turtle.sectors],
            options;
            accumulate_edges=profile_accumulate_edges,
            needs_sector_area=profile_needs_sector_area,
        )
    end
    _step_breakdown_progress!(:raycore_warm_stack_profile; warm_profile...)
    return nothing
end

main()
