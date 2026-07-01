using ArchimedLight
using Adapt
using GeometryBasics
using KernelAbstractions
using Metal
using PlantGeom
using Raycore
using StaticArrays

const PREP_SCENE_PATH = get(
    ENV,
    "ARCHIMEDLIGHT_PREP_BENCH_SCENE",
    "/Users/rvezy/Documents/cirad/Travail_AMAP/Projets/PalmStudio/simulations/Light_interception_vs_ages/archimed_simulations/sim.ops",
)
const PREP_PROGRESS_OUTPUT = get(ENV, "ARCHIMEDLIGHT_PREP_BENCH_PROGRESS_OUTPUT", "")

function _prep_progress!(event; kwargs...)
    isempty(PREP_PROGRESS_OUTPUT) && return nothing
    open(PREP_PROGRESS_OUTPUT, "a") do io
        print(io, time(), '\t', event)
        for (key, value) in kwargs
            print(io, '\t', key, '=', value)
        end
        println(io)
    end
    return nothing
end

function _prep_required_metal_backend()
    Sys.isapple() && Sys.ARCH == :aarch64 || error("Metal requires Apple Silicon macOS.")
    dev = Metal.MtlArray(zeros(Float32, 1))
    return KernelAbstractions.get_backend(dev)
end

function _prep_model_pairs(scene)
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

function _prep_default_models_for_scene(scene)
    by_group = Dict{String,Vector{Pair{String,Any}}}()
    for (group, type_name) in _prep_model_pairs(scene)
        material =
            occursin("soil", lowercase(group)) || occursin("ground", lowercase(group)) ||
            occursin("pavement", lowercase(group)) ?
            ArchimedLight.translucent(par=0.12, nir=0.60) :
            ArchimedLight.translucent(par=0.15, nir=0.90)
        push!(get!(by_group, group, Pair{String,Any}[]), type_name => material)
    end
    return ArchimedLight.models_for((group => Tuple(entries) for (group, entries) in by_group)...)
end

function _prep_elapsed!(f::F, event; kwargs...) where {F}
    started = time_ns()
    result = try
        f()
    finally
        _prep_progress!(event; ms=(time_ns() - started) / 1e6, kwargs...)
    end
    return result
end

function _prepare_interception_data_with_breakdown(scene, models, options)
    geometry = _prep_elapsed!(:scene_geometry_for_interception_done) do
        _scene_geometry_for_interception_with_breakdown(scene, models, options)
    end
    _prep_progress!(
        :geometry_summary;
        faces=length(geometry.faces),
        nodes=length(geometry.node_ids),
        vertices=length(geometry.vertices),
        instancing=geometry.raycore_instanced_geometry !== nothing,
    )

    node_interception_by_index = _prep_elapsed!(:resolved_models_done) do
        ArchimedLight._resolved_node_interception_models(geometry, models)
    end
    virtual_nodes = _prep_elapsed!(:virtual_nodes_done) do
        ArchimedLight._virtual_sensor_node_ids(geometry, node_interception_by_index)
    end
    virtual_node_mask = _prep_elapsed!(:virtual_node_mask_done) do
        [nid in virtual_nodes for nid in geometry.node_ids]
    end
    upper_hit = ArchimedLight._use_upper_hit_pixel_table(models, options)
    cache_ctx = _prep_elapsed!(:projection_cache_context_done) do
        ArchimedLight._projection_cache_context(
            geometry.vertices,
            geometry.faces,
            geometry.face2node,
            geometry.plotbox,
            options,
        )
    end
    emit_par, emit_nir = _prep_elapsed!(:emitter_power_done) do
        ArchimedLight._emitter_power_per_node(scene, models)
    end
    emitter_par_power_by_index = zeros(Float64, length(geometry.node_ids))
    emitter_nir_power_by_index = zeros(Float64, length(geometry.node_ids))
    for (nid, power) in emit_par
        emitter_par_power_by_index[geometry.node_index[nid]] = power
    end
    for (nid, power) in emit_nir
        emitter_nir_power_by_index[geometry.node_index[nid]] = power
    end
    emitter_nodes = Set(union(keys(emit_par), keys(emit_nir)))
    emitter_node_mask = _prep_elapsed!(:emitter_node_mask_done) do
        [nid in emitter_nodes for nid in geometry.node_ids]
    end
    component_area_per_node = _prep_elapsed!(:component_area_done) do
        ArchimedLight._interception_area_per_node_from_geometry(geometry)
    end
    absorption_par_per_node, absorption_nir_per_node = _prep_elapsed!(:absorptance_maps_done) do
        ArchimedLight._node_absorptance_maps_from_geometry(
            options,
            geometry,
            virtual_node_mask,
            node_interception_by_index,
        )
    end

    return ArchimedLight.PreparedInterceptionData(
        geometry,
        node_interception_by_index,
        ArchimedLight._node_transparency_by_index(node_interception_by_index, virtual_node_mask),
        virtual_nodes,
        virtual_node_mask,
        upper_hit,
        cache_ctx,
        emit_par,
        emit_nir,
        emitter_par_power_by_index,
        emitter_nir_power_by_index,
        emitter_nodes,
        emitter_node_mask,
        component_area_per_node,
        absorption_par_per_node,
        absorption_nir_per_node,
    )
end

function _scene_geometry_for_interception_with_breakdown(scene, models, options)
    raw_vertices = _prep_elapsed!(:decompose_vertices_done) do
        GeometryBasics.decompose(GeometryBasics.Point3, scene.merged_mesh)
    end
    vertices = _prep_elapsed!(:convert_vertices_done) do
        [StaticArrays.SVector{3,Float64}(v[1], v[2], v[3]) for v in raw_vertices]
    end
    all_faces = _prep_elapsed!(:decompose_faces_done) do
        collect(GeometryBasics.decompose(PlantGeom.Face3, scene.merged_mesh))
    end
    all_face2node = _prep_elapsed!(:collect_face2node_done) do
        collect(scene.face2node)
    end

    ignored = _prep_elapsed!(:ignored_group_types_done) do
        ArchimedLight._ignored_group_types(models)
    end
    all_node_ids = _prep_elapsed!(:all_unique_node_ids_done) do
        unique(all_face2node)
    end
    all_node_group, all_node_type = _prep_elapsed!(:metadata_maps_done) do
        ArchimedLight._scene_group_type_maps(scene, all_node_ids)
    end
    ignored_nodes = _prep_elapsed!(:ignored_node_ids_done) do
        ArchimedLight._ignored_node_ids(all_node_ids, all_node_group, all_node_type, ignored)
    end
    _prep_elapsed!(:validate_scene_models_done) do
        missing = Set{Tuple{String,String}}()
        for nid in all_node_ids
            nid in ignored_nodes && continue
            group = strip(get(all_node_group, nid, ""))
            type_name = strip(get(all_node_type, nid, ""))
            ArchimedLight._type_model(models, group, type_name) === nothing && push!(missing, (group, type_name))
        end
        isempty(missing) || error("Missing models: $missing")
    end

    faces, face2node = _prep_elapsed!(:filter_ignored_faces_done) do
        if isempty(ignored)
            all_faces, all_face2node
        else
            kept_faces = PlantGeom.Face3[]
            kept_face2node = Int[]
            sizehint!(kept_faces, length(all_faces))
            sizehint!(kept_face2node, length(all_face2node))
            for i in eachindex(all_faces)
                node_id = all_face2node[i]
                node_id in ignored_nodes && continue
                push!(kept_faces, all_faces[i])
                push!(kept_face2node, node_id)
            end
            kept_faces, kept_face2node
        end
    end
    isempty(face2node) && error("No intercepting geometry left after applying ignore rules.")

    plotbox = _prep_elapsed!(:plotbox_done) do
        ArchimedLight._plotbox(scene, vertices, options.pixel_size)
    end
    node_ids = _prep_elapsed!(:unique_node_ids_done) do
        unique(face2node)
    end
    node_index = _prep_elapsed!(:node_index_done) do
        Dict{Int,Int}(nid => i for (i, nid) in enumerate(node_ids))
    end
    face2node_index = _prep_elapsed!(:face2node_index_done) do
        [node_index[nid] for nid in face2node]
    end
    node_group = _prep_elapsed!(:node_group_done) do
        Dict{Int,String}(nid => get(all_node_group, nid, "") for nid in node_ids)
    end
    node_group_by_index = _prep_elapsed!(:node_group_by_index_done) do
        [get(node_group, nid, "") for nid in node_ids]
    end
    pavement_node_mask = _prep_elapsed!(:pavement_mask_done) do
        [group == "pavement" for group in node_group_by_index]
    end
    node_type = _prep_elapsed!(:node_type_done) do
        Dict{Int,String}(nid => get(all_node_type, nid, "") for nid in node_ids)
    end
    node_type_by_index = _prep_elapsed!(:node_type_by_index_done) do
        [get(node_type, nid, "") for nid in node_ids]
    end
    raycore_instanced_geometry = _prep_elapsed!(:raycore_instancing_probe_done) do
        ArchimedLight._raycore_instanced_geometry_from_scene(
            scene,
            vertices,
            node_ids,
            node_index,
            faces,
            face2node,
            face2node_index;
            toricity=options.toricity,
        )
    end

    return ArchimedLight.InterceptionSceneData(
        vertices,
        faces,
        face2node,
        face2node_index,
        node_ids,
        node_index,
        plotbox,
        node_group,
        node_group_by_index,
        pavement_node_mask,
        node_type,
        node_type_by_index,
        raycore_instanced_geometry,
        Ref{Any}(nothing),
    )
end

function _raycore_scene_data_with_breakdown(
    prepared,
    config,
    options;
    face_chunk_limit::Int=ArchimedLight._raycore_face_chunk_limit(),
)
    _prep_progress!(:chunk_limit_start; face_chunk_limit=face_chunk_limit)
    tlas = Raycore.TLAS(config.backend)
    transforms = _prep_elapsed!(:toric_transforms_done) do
        ArchimedLight._raycore_toric_transforms(prepared.geometry; toricity=options.toricity)
    end
    chunk_count = Ref(0)
    push_ms = Ref(0.0)
    chunk_emit_ms = Ref(0.0)
    _prep_elapsed!(:chunked_push_geometry_done) do
        ArchimedLight._raycore_foreach_mesh_chunk_from_geometry(
            mesh -> begin
                chunk_count[] += 1
                started = time_ns()
                options.toricity ? push!(tlas, mesh, transforms) : push!(tlas, mesh)
                push_ms[] += (time_ns() - started) / 1e6
                return nothing
            end,
            prepared.geometry;
            toricity=options.toricity,
            face_chunk_limit=face_chunk_limit,
        )
    end
    _prep_progress!(
        :chunked_push_geometry_summary;
        face_chunk_limit=face_chunk_limit,
        chunks=chunk_count[],
        push_ms=push_ms[],
        emit_overhead_ms=chunk_emit_ms[],
        instances=length(tlas.instances),
    )
    instance_count = length(tlas.instances)
    metadata_stride = length(prepared.geometry.node_ids) + 1

    _prep_elapsed!(:raycore_sync_done) do
        Raycore.sync!(tlas)
    end
    kernel_tlas = _prep_elapsed!(:adapt_tlas_done) do
        Adapt.adapt(config.backend, tlas)
    end
    hit_decoder = _prep_elapsed!(:identity_hit_decoder_done) do
        ArchimedLight._raycore_identity_hit_decoder(prepared, config.backend, instance_count)
    end
    n_pixels = prepared.geometry.plotbox.nx * prepared.geometry.plotbox.ny
    stack_len = n_pixels * config.max_hits_per_pixel
    n_nodes = length(prepared.geometry.node_ids)
    dense_pairs_i128 = Int128(n_nodes) * Int128(n_nodes)
    dense_bytes_i128 = dense_pairs_i128 * Int128(sizeof(Int32))
    dense_matrix_available =
        dense_pairs_i128 <= typemax(Int) &&
        dense_bytes_i128 <= config.dense_edge_limit_bytes &&
        KernelAbstractions.supports_atomics(config.backend)

    buffers = _prep_elapsed!(:raycore_buffers_done) do
        top_nodes_dev = KernelAbstractions.allocate(config.backend, UInt32, n_pixels)
        stack_counts_dev = KernelAbstractions.allocate(config.backend, Int32, n_pixels)
        stack_nodes_dev = KernelAbstractions.allocate(config.backend, UInt32, stack_len)
        stack_instance_indices_dev = KernelAbstractions.allocate(config.backend, UInt32, stack_len)
        stack_heights_dev = KernelAbstractions.allocate(config.backend, Float32, stack_len)
        stack_overflow_dev = KernelAbstractions.allocate(config.backend, Bool, n_pixels)
        edge_key_capacity = n_pixels * max(0, 2 * (config.max_hits_per_pixel - 1))
        auto_dense_enabled = config.edge_accumulation == :auto && false && dense_matrix_available
        sparse_scratch_enabled =
            config.edge_accumulation == :sparse_host_reduce ||
            (config.edge_accumulation == :auto && !auto_dense_enabled)
        edge_keys_dev =
            sparse_scratch_enabled ? KernelAbstractions.allocate(config.backend, UInt64, edge_key_capacity) : nothing
        edge_key_counts_dev =
            sparse_scratch_enabled ? KernelAbstractions.allocate(config.backend, Int32, n_pixels) : nothing
        dense_scratch_enabled =
            (config.edge_accumulation == :dense_atomic || auto_dense_enabled) &&
            dense_matrix_available
        dense_edge_counts_dev =
            dense_scratch_enabled ?
            KernelAbstractions.allocate(config.backend, Int32, Int(dense_pairs_i128)) :
            nothing
        node_counts_dev =
            KernelAbstractions.supports_atomics(config.backend) ?
            KernelAbstractions.allocate(config.backend, Int32, n_nodes) :
            nothing
        node_transparency_dev = KernelAbstractions.allocate(config.backend, Float32, n_nodes)
        sector_area_dev =
            KernelAbstractions.supports_atomics(config.backend) ?
            KernelAbstractions.allocate(config.backend, Float32, n_nodes) :
            nothing
        virtual_node_mask_dev = KernelAbstractions.allocate(config.backend, Bool, length(prepared.virtual_node_mask))
        pavement_node_mask_dev =
            KernelAbstractions.allocate(config.backend, Bool, length(prepared.geometry.pavement_node_mask))
        node_ids_dev = KernelAbstractions.allocate(config.backend, Int, n_nodes)
        KernelAbstractions.copyto!(config.backend, node_transparency_dev, Float32.(prepared.node_transparency_by_index))
        KernelAbstractions.copyto!(config.backend, virtual_node_mask_dev, prepared.virtual_node_mask)
        KernelAbstractions.copyto!(config.backend, pavement_node_mask_dev, prepared.geometry.pavement_node_mask)
        KernelAbstractions.copyto!(config.backend, node_ids_dev, prepared.geometry.node_ids)
        (
            top_nodes_dev,
            stack_counts_dev,
            stack_nodes_dev,
            stack_instance_indices_dev,
            stack_heights_dev,
            stack_overflow_dev,
            edge_keys_dev,
            edge_key_counts_dev,
            dense_edge_counts_dev,
            node_counts_dev,
            node_transparency_dev,
            sector_area_dev,
            virtual_node_mask_dev,
            pavement_node_mask_dev,
            node_ids_dev,
        )
    end

    vertical_span = ArchimedLight._raycore_geometry_vertical_span(prepared.geometry)
    return ArchimedLight.RaycoreSceneData(
        prepared,
        tlas,
        kernel_tlas,
        config.backend,
        buffers[1],
        UInt32[],
        buffers[2],
        buffers[3],
        buffers[4],
        buffers[5],
        buffers[6],
        buffers[7],
        buffers[8],
        buffers[9],
        Int32[],
        buffers[10],
        buffers[11],
        buffers[12],
        Int32[],
        Float32[],
        buffers[13],
        buffers[14],
        buffers[15],
        Int32[],
        UInt32[],
        UInt32[],
        Float32[],
        Bool[],
        UInt64[],
        Int32[],
        UInt64[],
        Dict{Tuple{Float64,Float64,Float64,Bool},Float32}(),
        Dict{Tuple{Float64,Float64,Float64},Vector{Float64}}(),
        Dict{Tuple{Float64,Float64,Float64,Bool},Any}(),
        hit_decoder,
        Dict{Tuple{Symbol,Float64,Float64,Float64},Vector{Float64}}(),
        config.workgroupsize,
        config.max_hits_per_pixel,
        config.hit_epsilon,
        config.edge_accumulation,
        config.dense_edge_limit_bytes,
        config.validate,
        vertical_span,
        true,
        :chunked_merged_mesh,
    )
end

function _prep_face_chunk_limits()
    raw = get(ENV, "ARCHIMEDLIGHT_PREP_BENCH_FACE_CHUNK_LIMITS", "")
    isempty(strip(raw)) && return [ArchimedLight._raycore_face_chunk_limit()]
    limits = Int[]
    for item in split(raw, ',')
        value = parse(Int, strip(item))
        push!(limits, value <= 0 ? typemax(Int) : value)
    end
    return limits
end

function main()
    _prep_progress!(:scene_read_start; path=PREP_SCENE_PATH)
    scene = _prep_elapsed!(:scene_read_done) do
        ArchimedLight.read_scene(PREP_SCENE_PATH)
    end
    _prep_progress!(:scene_loaded; faces=length(scene.face2node))

    models = _prep_elapsed!(:models_done) do
        _prep_default_models_for_scene(scene)
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
        backend=_prep_required_metal_backend(),
        max_hits_per_pixel=64,
        workgroupsize=256,
        validate=false,
    )

    prepared = _prep_elapsed!(:prepare_interception_done) do
        _prepare_interception_data_with_breakdown(scene, models, options)
    end
    _prep_progress!(
        :prepared_summary;
        faces=length(prepared.geometry.faces),
        nodes=length(prepared.geometry.node_ids),
        instancing=prepared.geometry.raycore_instanced_geometry !== nothing,
    )

    status = _prep_elapsed!(:prechunk_status_done) do
        ArchimedLight._raycore_prechunk_instance_limit_status(prepared, ib.config; toricity=options.toricity)
    end
    _prep_progress!(:prechunk_status; exceeded=status.exceeded)

    data = nothing
    for face_chunk_limit in _prep_face_chunk_limits()
        candidate = _prep_elapsed!(:initial_scene_data_done; face_chunk_limit=face_chunk_limit) do
            _raycore_scene_data_with_breakdown(
                prepared,
                ib.config,
                options;
                face_chunk_limit=face_chunk_limit,
            )
        end
        data = candidate
        traversal = ArchimedLight._raycore_traversal_stack_depth_status(candidate)
        _prep_progress!(
            :raycore_data_summary;
            face_chunk_limit=face_chunk_limit,
            chunked=true,
            geometry_mode=candidate.geometry_mode,
            tlas_instances=length(candidate.tlas.instances),
            max_blas_depth=traversal.blas_depth,
            stack_capacity=traversal.capacity,
            exceeds_stack=traversal.exceeds,
        )

        validation = _prep_elapsed!(:vertical_validation_done; face_chunk_limit=face_chunk_limit) do
            ArchimedLight._raycore_vertical_trace_validation(candidate, options)
        end
        _prep_progress!(:vertical_validation; face_chunk_limit=face_chunk_limit, ok=validation.ok)

        stack_validation = _prep_elapsed!(:stack_validation_done; face_chunk_limit=face_chunk_limit) do
            ArchimedLight._raycore_stack_trace_validation(candidate, options)
        end
        _prep_progress!(:stack_validation; face_chunk_limit=face_chunk_limit, ok=stack_validation.ok)
    end
    return data
end

main()
