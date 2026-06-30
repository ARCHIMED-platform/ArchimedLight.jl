using KernelAbstractions
using Test
using ArchimedLight

const METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_METAL", ""))
const METAL_REQUESTED = METAL_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const METAL_REQUIRED = METAL_FLAG in ("required", "force", "error")
const RAYCORE_INVALID_NODE = ArchimedLight.RAYCORE_INVALID_NODE
const RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY = ArchimedLight._raycore_software_traversal_stack_capacity()
const THREE_WAY_AGGREGATE_RTOL = 5e-3
const THREE_WAY_AGGREGATE_ATOL = 1e-3

function _skip_or_fail(message::AbstractString)
    if METAL_REQUIRED
        error(message)
    else
        @test_skip message
        return nothing
    end
end

function _metal_backend()
    if !Sys.isapple()
        return _skip_or_fail("Metal validation requires macOS.")
    end
    if Sys.ARCH != :aarch64
        return _skip_or_fail("Metal validation requires Apple Silicon (Sys.ARCH == :aarch64).")
    end
    if Base.find_package("Metal") === nothing
        return _skip_or_fail("Metal.jl is not available in this optional test environment. Run `julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'`.")
    end

    metal_mod = try
        Base.require(Main, :Metal)
    catch err
        return _skip_or_fail("Metal.jl could not be imported: $(sprint(showerror, err))")
    end

    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = try
        Base.invokelatest(array_type, zeros(Float32, 1))
    catch err
        return _skip_or_fail("Metal.jl is installed, but creating a Metal array failed: $(sprint(showerror, err))")
    end
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _bvh_depth(nodes)
    isempty(nodes) && return 0
    function visit(node_idx::Int, depth::Int)
        1 <= node_idx <= length(nodes) || return depth
        node = nodes[node_idx]
        node.child0 == RAYCORE_INVALID_NODE && return depth
        return max(visit(Int(node.child0), depth + 1), visit(Int(node.child1), depth + 1))
    end
    return visit(1, 1)
end

function _max_blas_depth_and_nodes(data)
    static = data.tlas.static_tlas
    nodes = Array(static.all_blas_nodes)
    descriptors = Array(static.blas_descriptors)
    isempty(nodes) && return (depth=0, node_count=0)
    isempty(descriptors) && return (depth=_bvh_depth(nodes), node_count=length(nodes))

    max_depth = 0
    for idx in eachindex(descriptors)
        first_node = Int(descriptors[idx].nodes_offset) + 1
        last_node = idx == length(descriptors) ? length(nodes) : Int(descriptors[idx+1].nodes_offset)
        first_node <= last_node || continue
        max_depth = max(max_depth, _bvh_depth(@view nodes[first_node:last_node]))
    end
    return (depth=max_depth, node_count=length(nodes))
end

function _raycore_topology_summary(data)
    static = data.tlas.static_tlas
    blas = _max_blas_depth_and_nodes(data)
    return (
        tlas_depth=_bvh_depth(Array(static.nodes)),
        blas_depth=blas.depth,
        blas_node_count=blas.node_count,
    )
end

function _raycore_static_tlas_source(data)
    data.backend isa KernelAbstractions.CPU && return :native_cpu
    return data.tlas.backend isa KernelAbstractions.CPU ? :cpu_built_static : :native_device
end

function _raw_stack_sweep(
    scene,
    models,
    options,
    sky,
    backend;
    max_hits_values=(32, 64, 128, 256),
    hit_epsilon=1.0f-4,
    face_chunk_limit=nothing,
    candidate_static_tlas_from_reference::Bool=false,
    assert_ok::Bool,
)
    turtle = ArchimedLight.build_turtle(options, sky)
    sector_idx = findfirst(sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0, turtle.sectors)
    @test sector_idx !== nothing
    sector_idx === nothing && return NamedTuple[]

    prepared = ArchimedLight._prepare_interception_data(
        scene,
        models,
        options;
        include_budget_maps=true,
        include_raycore_instancing=true,
    )
    rows = NamedTuple[]
    direction = turtle.sectors[sector_idx].direction
    for max_hits in max_hits_values
        cpu_data = ArchimedLight._raycore_scene_data(
            prepared,
            ArchimedLight.RaycoreBackendConfig(
                backend=KernelAbstractions.CPU(),
                max_hits_per_pixel=max_hits,
                hit_epsilon=hit_epsilon,
                validate=true,
            );
            toricity=options.toricity,
            face_chunk_limit=face_chunk_limit,
        )
        metal_data = ArchimedLight._raycore_scene_data(
            prepared,
            ArchimedLight.RaycoreBackendConfig(
                backend=backend,
                max_hits_per_pixel=max_hits,
                hit_epsilon=hit_epsilon,
                validate=true,
            );
            toricity=options.toricity,
            face_chunk_limit=face_chunk_limit,
        )
        if candidate_static_tlas_from_reference
            metal_data = ArchimedLight._raycore_scene_data_with_static_tlas(metal_data, cpu_data)
        end
        cpu_top_hits_before = ArchimedLight._raycore_trace_direction_top_hits_direct(cpu_data, direction, options)
        metal_top_hits_before = ArchimedLight._raycore_trace_direction_top_hits_direct(metal_data, direction, options)
        raw_comparison = ArchimedLight._raycore_raw_stack_comparison(
            cpu_data,
            metal_data,
            direction,
            options,
        )
        cpu_topology = _raycore_topology_summary(cpu_data)
        metal_topology = _raycore_topology_summary(metal_data)
        mismatch_pixel = raw_comparison.mismatch !== nothing && hasproperty(raw_comparison.mismatch, :pixel) ?
                         raw_comparison.mismatch.pixel :
                         nothing
        mismatch_ray = mismatch_pixel === nothing ? nothing : ArchimedLight._raycore_sky_to_ground_ray(
            cpu_data.prepared.geometry.plotbox,
            mismatch_pixel,
            direction,
            ArchimedLight._raycore_projection_far(cpu_data, direction, options),
        )
        cpu_top_hit = mismatch_ray === nothing ? nothing : ArchimedLight._raycore_trace_top_hits(
            cpu_data,
            [mismatch_ray.o],
            [mismatch_ray.d];
            t_mins=[mismatch_ray.t_min],
            t_maxs=[mismatch_ray.t_max],
        )
        metal_top_hit = mismatch_ray === nothing ? nothing : ArchimedLight._raycore_trace_top_hits(
            metal_data,
            [mismatch_ray.o],
            [mismatch_ray.d];
            t_mins=[mismatch_ray.t_min],
            t_maxs=[mismatch_ray.t_max],
        )
        reference_local_ray = _raycore_reference_instance_local_ray(
            cpu_data,
            mismatch_ray,
            raw_comparison.mismatch,
        )
        candidate_local_ray = _raycore_reference_instance_local_ray(
            metal_data,
            mismatch_ray,
            raw_comparison.mismatch,
        )
        row = (
            max_hits=max_hits,
            hit_epsilon=hit_epsilon,
            ok=raw_comparison.ok,
            reference_hits=raw_comparison.reference_hits,
            candidate_hits=raw_comparison.candidate_hits,
            hit_ratio=raw_comparison.hit_ratio,
            reference_occupied=raw_comparison.reference_occupied,
            candidate_occupied=raw_comparison.candidate_occupied,
            occupied_ratio=raw_comparison.occupied_ratio,
            reference_overflow=raw_comparison.reference_overflow,
            candidate_overflow=raw_comparison.candidate_overflow,
            reference_max_seen=raw_comparison.reference_max_seen,
            candidate_max_seen=raw_comparison.candidate_max_seen,
            diagnosis=raw_comparison.diagnosis,
            reference_geometry_mode=cpu_data.geometry_mode,
            candidate_geometry_mode=metal_data.geometry_mode,
            candidate_static_tlas_source=
                candidate_static_tlas_from_reference ? :reference_cpu : _raycore_static_tlas_source(metal_data),
            face_chunk_limit=face_chunk_limit,
            reference_instances=length(cpu_data.tlas.instances),
            candidate_instances=length(metal_data.tlas.instances),
            reference_tlas_depth=cpu_topology.tlas_depth,
            candidate_tlas_depth=metal_topology.tlas_depth,
            reference_blas_depth=cpu_topology.blas_depth,
            candidate_blas_depth=metal_topology.blas_depth,
            raycore_traversal_stack_capacity=RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
            reference_blas_depth_over_stack_capacity=cpu_topology.blas_depth > RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
            candidate_blas_depth_over_stack_capacity=metal_topology.blas_depth > RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
            reference_blas_node_count=cpu_topology.blas_node_count,
            candidate_blas_node_count=metal_topology.blas_node_count,
            mismatch_kind=raw_comparison.mismatch === nothing ? :none : raw_comparison.mismatch.kind,
            ray_origin_x=mismatch_ray === nothing ? missing : Float32(mismatch_ray.o[1]),
            ray_origin_y=mismatch_ray === nothing ? missing : Float32(mismatch_ray.o[2]),
            ray_origin_z=mismatch_ray === nothing ? missing : Float32(mismatch_ray.o[3]),
            ray_direction_x=mismatch_ray === nothing ? missing : Float32(mismatch_ray.d[1]),
            ray_direction_y=mismatch_ray === nothing ? missing : Float32(mismatch_ray.d[2]),
            ray_direction_z=mismatch_ray === nothing ? missing : Float32(mismatch_ray.d[3]),
            ray_t_min=mismatch_ray === nothing ? missing : mismatch_ray.t_min,
            ray_t_max=mismatch_ray === nothing ? missing : mismatch_ray.t_max,
            reference_top_node_before=mismatch_pixel === nothing ? missing : cpu_top_hits_before.nodes[mismatch_pixel],
            candidate_top_node_before=mismatch_pixel === nothing ? missing : metal_top_hits_before.nodes[mismatch_pixel],
            reference_top_instance_before=mismatch_pixel === nothing ? missing : cpu_top_hits_before.instance_indices[mismatch_pixel],
            candidate_top_instance_before=mismatch_pixel === nothing ? missing : metal_top_hits_before.instance_indices[mismatch_pixel],
            reference_top_distance_before=mismatch_pixel === nothing ? missing : cpu_top_hits_before.distances[mismatch_pixel],
            candidate_top_distance_before=mismatch_pixel === nothing ? missing : metal_top_hits_before.distances[mismatch_pixel],
            reference_top_node=mismatch_ray === nothing ? missing : cpu_top_hit.nodes[1],
            candidate_top_node=mismatch_ray === nothing ? missing : metal_top_hit.nodes[1],
            reference_top_instance=mismatch_ray === nothing ? missing : cpu_top_hit.instance_indices[1],
            candidate_top_instance=mismatch_ray === nothing ? missing : metal_top_hit.instance_indices[1],
            reference_top_distance=mismatch_ray === nothing ? missing : cpu_top_hit.distances[1],
            candidate_top_distance=mismatch_ray === nothing ? missing : metal_top_hit.distances[1],
            reference_unmatched_instance=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_instance,
            reference_unmatched_blas_index=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_index,
            reference_unmatched_instance_id=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_instance_id,
            reference_unmatched_inv_tx=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_inv_tx,
            reference_unmatched_inv_ty=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_inv_ty,
            reference_unmatched_inv_tz=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_inv_tz,
            reference_unmatched_local_origin_x=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_origin_x,
            reference_unmatched_local_origin_y=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_origin_y,
            reference_unmatched_local_origin_z=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_origin_z,
            reference_unmatched_local_direction_x=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_direction_x,
            reference_unmatched_local_direction_y=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_direction_y,
            reference_unmatched_local_direction_z=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_local_direction_z,
            reference_unmatched_blas_min_x=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_min_x,
            reference_unmatched_blas_min_y=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_min_y,
            reference_unmatched_blas_min_z=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_min_z,
            reference_unmatched_blas_max_x=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_max_x,
            reference_unmatched_blas_max_y=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_max_y,
            reference_unmatched_blas_max_z=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_max_z,
            reference_unmatched_blas_entry_t=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_entry_t,
            reference_unmatched_blas_exit_t=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_exit_t,
            reference_unmatched_blas_intersects=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_blas_intersects,
            reference_unmatched_tlas_leaf_entry_t=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_tlas_leaf_entry_t,
            reference_unmatched_tlas_leaf_exit_t=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_tlas_leaf_exit_t,
            reference_unmatched_tlas_leaf_intersects=reference_local_ray === nothing ? missing : reference_local_ray.reference_unmatched_tlas_leaf_intersects,
            candidate_expected_instance=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_instance,
            candidate_expected_blas_index=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_blas_index,
            candidate_expected_inv_tx=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_inv_tx,
            candidate_expected_inv_ty=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_inv_ty,
            candidate_expected_inv_tz=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_inv_tz,
            candidate_expected_local_origin_x=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_origin_x,
            candidate_expected_local_origin_y=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_origin_y,
            candidate_expected_local_origin_z=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_origin_z,
            candidate_expected_local_direction_x=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_direction_x,
            candidate_expected_local_direction_y=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_direction_y,
            candidate_expected_local_direction_z=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_local_direction_z,
            candidate_expected_blas_entry_t=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_blas_entry_t,
            candidate_expected_blas_exit_t=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_blas_exit_t,
            candidate_expected_blas_intersects=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_blas_intersects,
            candidate_expected_tlas_leaf_entry_t=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_tlas_leaf_entry_t,
            candidate_expected_tlas_leaf_exit_t=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_tlas_leaf_exit_t,
            candidate_expected_tlas_leaf_intersects=candidate_local_ray === nothing ? missing : candidate_local_ray.reference_unmatched_tlas_leaf_intersects,
            mismatch=raw_comparison.mismatch,
        )
        push!(rows, row)
        if assert_ok
            @test row.ok
            @test row.reference_hits == row.candidate_hits
            @test row.hit_ratio == 1.0
            @test row.reference_occupied == row.candidate_occupied
            @test row.occupied_ratio == 1.0
            @test row.reference_overflow == row.candidate_overflow == 0
            @test row.diagnosis == :exact
            @test !row.candidate_blas_depth_over_stack_capacity
            @test row.mismatch === nothing
        end
    end
    @info "Raycore Metal raw stack max-hit sweep" assert_ok rows
    return rows
end

_raw_mismatch(row) =
    hasproperty(row, :mismatch) && getproperty(row, :mismatch) !== nothing ? getproperty(row, :mismatch) : nothing

_raw_mismatch_field(row, field::Symbol, default="missing") = begin
    mismatch = _raw_mismatch(row)
    mismatch !== nothing && hasproperty(mismatch, field) ? string(getproperty(mismatch, field)) : default
end

function _raw_mismatch_hit_field(row, hit_field::Symbol, value_field::Symbol, default="missing")
    mismatch = _raw_mismatch(row)
    mismatch === nothing && return default
    hasproperty(mismatch, hit_field) || return default
    hit = getproperty(mismatch, hit_field)
    hit === nothing && return default
    return hasproperty(hit, value_field) ? string(getproperty(hit, value_field)) : default
end

function _ray_bbox_interval(origin, direction, p_min, p_max, t_min::Float32, t_max::Float32)
    safe_inv(component) = 1.0f0 / (abs(component) > 1.0f-5 ? component : copysign(1.0f-5, component))
    inv_d = (safe_inv(direction[1]), safe_inv(direction[2]), safe_inv(direction[3]))
    tx1 = (p_min[1] - origin[1]) * inv_d[1]
    tx2 = (p_max[1] - origin[1]) * inv_d[1]
    ty1 = (p_min[2] - origin[2]) * inv_d[2]
    ty2 = (p_max[2] - origin[2]) * inv_d[2]
    tz1 = (p_min[3] - origin[3]) * inv_d[3]
    tz2 = (p_max[3] - origin[3]) * inv_d[3]
    entry = max(max(min(tx1, tx2), min(ty1, ty2)), max(min(tz1, tz2), t_min))
    exit = min(min(max(tx1, tx2), max(ty1, ty2)), min(max(tz1, tz2), t_max))
    return (entry=Float32(entry), exit=Float32(exit), intersects=entry <= exit)
end

function _tlas_leaf_bounds_for_instance(data, instance_idx::Int)
    nodes = Array(data.tlas.static_tlas.nodes)
    leaf = findfirst(node -> node.child0 == RAYCORE_INVALID_NODE && node.child1 == UInt32(instance_idx - 1), nodes)
    leaf === nothing && return nothing
    node = nodes[leaf]
    return (p_min=node.aabb0_min, p_max=node.aabb0_max)
end

function _raycore_reference_instance_local_ray(data, mismatch_ray, mismatch)
    mismatch_ray === nothing && return nothing
    mismatch === nothing && return nothing
    hasproperty(mismatch, :first_unmatched_reference) || return nothing
    hit = mismatch.first_unmatched_reference
    hit === nothing && return nothing
    instance_idx = Int(hit.instance_index)
    instance_idx > 0 || return nothing
    1 <= instance_idx <= length(data.tlas.instances) || return nothing

    raycore_mod = Base.require(Main, :Raycore)
    instance = Array(data.tlas.instances[instance_idx:instance_idx])[1]
    local_o = raycore_mod.transform_point(instance.inv_transform, mismatch_ray.o)
    local_d = raycore_mod.transform_direction(instance.inv_transform, mismatch_ray.d)
    descriptors = Array(data.tlas.static_tlas.blas_descriptors)
    blas_bounds = descriptors[Int(instance.blas_index)].root_aabb
    blas_interval = _ray_bbox_interval(
        local_o,
        local_d,
        blas_bounds.p_min,
        blas_bounds.p_max,
        mismatch_ray.t_min,
        mismatch_ray.t_max,
    )
    tlas_bounds = _tlas_leaf_bounds_for_instance(data, instance_idx)
    tlas_interval = tlas_bounds === nothing ? nothing : _ray_bbox_interval(
        mismatch_ray.o,
        mismatch_ray.d,
        tlas_bounds.p_min,
        tlas_bounds.p_max,
        mismatch_ray.t_min,
        mismatch_ray.t_max,
    )
    return (
        reference_unmatched_instance=instance_idx,
        reference_unmatched_blas_index=Int(instance.blas_index),
        reference_unmatched_instance_id=Int(instance.instance_id),
        reference_unmatched_inv_tx=Float32(instance.inv_transform[4, 1]),
        reference_unmatched_inv_ty=Float32(instance.inv_transform[4, 2]),
        reference_unmatched_inv_tz=Float32(instance.inv_transform[4, 3]),
        reference_unmatched_local_origin_x=Float32(local_o[1]),
        reference_unmatched_local_origin_y=Float32(local_o[2]),
        reference_unmatched_local_origin_z=Float32(local_o[3]),
        reference_unmatched_local_direction_x=Float32(local_d[1]),
        reference_unmatched_local_direction_y=Float32(local_d[2]),
        reference_unmatched_local_direction_z=Float32(local_d[3]),
        reference_unmatched_blas_min_x=Float32(blas_bounds.p_min[1]),
        reference_unmatched_blas_min_y=Float32(blas_bounds.p_min[2]),
        reference_unmatched_blas_min_z=Float32(blas_bounds.p_min[3]),
        reference_unmatched_blas_max_x=Float32(blas_bounds.p_max[1]),
        reference_unmatched_blas_max_y=Float32(blas_bounds.p_max[2]),
        reference_unmatched_blas_max_z=Float32(blas_bounds.p_max[3]),
        reference_unmatched_blas_entry_t=blas_interval.entry,
        reference_unmatched_blas_exit_t=blas_interval.exit,
        reference_unmatched_blas_intersects=blas_interval.intersects,
        reference_unmatched_tlas_leaf_entry_t=tlas_interval === nothing ? missing : tlas_interval.entry,
        reference_unmatched_tlas_leaf_exit_t=tlas_interval === nothing ? missing : tlas_interval.exit,
        reference_unmatched_tlas_leaf_intersects=tlas_interval === nothing ? missing : tlas_interval.intersects,
    )
end

function _append_raw_stack_smoke_rows(path::AbstractString, label::AbstractString, rows)
    fields = (
        :label,
        :max_hits,
        :hit_epsilon,
        :ok,
        :diagnosis,
        :reference_hits,
        :candidate_hits,
        :hit_ratio,
        :reference_occupied,
        :candidate_occupied,
        :occupied_ratio,
        :reference_overflow,
        :candidate_overflow,
        :reference_max_seen,
        :candidate_max_seen,
        :reference_geometry_mode,
        :candidate_geometry_mode,
        :candidate_static_tlas_source,
        :face_chunk_limit,
        :reference_instances,
        :candidate_instances,
        :reference_blas_depth,
        :candidate_blas_depth,
        :reference_blas_node_count,
        :candidate_blas_node_count,
        :raycore_traversal_stack_capacity,
        :candidate_blas_depth_over_stack_capacity,
        :mismatch_kind,
        :mismatch_pixel,
        :mismatch_reference_count,
        :mismatch_candidate_count,
        :unmatched_reference_count,
        :unmatched_candidate_count,
        :first_unmatched_reference_metadata,
        :first_unmatched_reference_instance,
        :first_unmatched_reference_distance,
        :first_unmatched_candidate_metadata,
        :first_unmatched_candidate_instance,
        :first_unmatched_candidate_distance,
        :ray_origin_x,
        :ray_origin_y,
        :ray_origin_z,
        :ray_direction_x,
        :ray_direction_y,
        :ray_direction_z,
        :ray_t_min,
        :ray_t_max,
        :reference_top_node_before,
        :candidate_top_node_before,
        :reference_top_instance_before,
        :candidate_top_instance_before,
        :reference_top_distance_before,
        :candidate_top_distance_before,
        :reference_top_node,
        :candidate_top_node,
        :reference_top_instance,
        :candidate_top_instance,
        :reference_top_distance,
        :candidate_top_distance,
        :reference_unmatched_instance,
        :reference_unmatched_blas_index,
        :reference_unmatched_instance_id,
        :reference_unmatched_inv_tx,
        :reference_unmatched_inv_ty,
        :reference_unmatched_inv_tz,
        :reference_unmatched_local_origin_x,
        :reference_unmatched_local_origin_y,
        :reference_unmatched_local_origin_z,
        :reference_unmatched_local_direction_x,
        :reference_unmatched_local_direction_y,
        :reference_unmatched_local_direction_z,
        :reference_unmatched_blas_min_x,
        :reference_unmatched_blas_min_y,
        :reference_unmatched_blas_min_z,
        :reference_unmatched_blas_max_x,
        :reference_unmatched_blas_max_y,
        :reference_unmatched_blas_max_z,
        :reference_unmatched_blas_entry_t,
        :reference_unmatched_blas_exit_t,
        :reference_unmatched_blas_intersects,
        :reference_unmatched_tlas_leaf_entry_t,
        :reference_unmatched_tlas_leaf_exit_t,
        :reference_unmatched_tlas_leaf_intersects,
        :candidate_expected_instance,
        :candidate_expected_blas_index,
        :candidate_expected_inv_tx,
        :candidate_expected_inv_ty,
        :candidate_expected_inv_tz,
        :candidate_expected_local_origin_x,
        :candidate_expected_local_origin_y,
        :candidate_expected_local_origin_z,
        :candidate_expected_local_direction_x,
        :candidate_expected_local_direction_y,
        :candidate_expected_local_direction_z,
        :candidate_expected_blas_entry_t,
        :candidate_expected_blas_exit_t,
        :candidate_expected_blas_intersects,
        :candidate_expected_tlas_leaf_entry_t,
        :candidate_expected_tlas_leaf_exit_t,
        :candidate_expected_tlas_leaf_intersects,
    )
    write_header = !isfile(path)
    open(path, "a") do io
        write_header && println(io, join(string.(fields), '\t'))
        for row in rows
            values = map(fields) do field
                if field === :label
                    label
                elseif field === :mismatch_pixel
                    _raw_mismatch_field(row, :pixel)
                elseif field === :mismatch_reference_count
                    _raw_mismatch_field(row, :reference_count)
                elseif field === :mismatch_candidate_count
                    _raw_mismatch_field(row, :candidate_count)
                elseif field === :unmatched_reference_count
                    _raw_mismatch_field(row, :unmatched_reference_count)
                elseif field === :unmatched_candidate_count
                    _raw_mismatch_field(row, :unmatched_candidate_count)
                elseif field === :first_unmatched_reference_metadata
                    _raw_mismatch_hit_field(row, :first_unmatched_reference, :metadata)
                elseif field === :first_unmatched_reference_instance
                    _raw_mismatch_hit_field(row, :first_unmatched_reference, :instance_index)
                elseif field === :first_unmatched_reference_distance
                    _raw_mismatch_hit_field(row, :first_unmatched_reference, :distance)
                elseif field === :first_unmatched_candidate_metadata
                    _raw_mismatch_hit_field(row, :first_unmatched_candidate, :metadata)
                elseif field === :first_unmatched_candidate_instance
                    _raw_mismatch_hit_field(row, :first_unmatched_candidate, :instance_index)
                elseif field === :first_unmatched_candidate_distance
                    _raw_mismatch_hit_field(row, :first_unmatched_candidate, :distance)
                elseif hasproperty(row, field)
                    string(getproperty(row, field))
                else
                    "missing"
                end
            end
            println(io, join(values, '\t'))
        end
    end
    return nothing
end

function _assert_raw_stack_sweep_ok(scene, models, options, sky, backend; max_hits_values=(32, 64, 128, 256))
    return _raw_stack_sweep(scene, models, options, sky, backend; max_hits_values=max_hits_values, assert_ok=true)
end

function _assert_unresolved_raw_failures_are_diagnosed(rows; require_missing_hit::Bool=false)
    failing = [row for row in rows if !row.ok]
    isempty(failing) && return nothing

    @test all(row -> row.diagnosis != :exact, failing)
    @test all(row -> row.mismatch !== nothing, failing)
    @test all(row -> hasproperty(row.mismatch, :diagnosis), failing)
    @test all(row -> row.reference_overflow == 0 && row.candidate_overflow == 0, failing)
    @test all(row -> row.diagnosis in (
        :candidate_missing_hits,
        :candidate_extra_hits,
        :raw_order_difference,
        :raw_stack_membership_difference,
        :instance_index_difference,
        :distance_difference,
    ), failing)

    missing_hit_rows = [row for row in failing if row.diagnosis == :candidate_missing_hits]
    if require_missing_hit
        @test !isempty(missing_hit_rows)
    end
    for row in missing_hit_rows
        @test row.mismatch.unmatched_reference_count > 0
        @test row.mismatch.unmatched_candidate_count == 0
        @test row.mismatch.first_unmatched_reference !== nothing
    end
    return nothing
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

function _run_single_sky(scene, models, options, sky; interception_backend, scattering_backend=nothing)
    sim = ArchimedLight.LightSimulation(
        scene,
        models;
        options=options,
        interception_backend=interception_backend,
        scattering_backend=scattering_backend,
    )
    step = ArchimedLight.run_light(sim, sky; step_duration_seconds=900.0)
    return (step=step, sim=sim, totals=_step_totals(step))
end

function _assert_totals_close(candidate, reference; label)
    for key in keys(reference)
        @test isapprox(
            candidate[key],
            reference[key];
            rtol=THREE_WAY_AGGREGATE_RTOL,
            atol=THREE_WAY_AGGREGATE_ATOL,
        )
    end
    @info "Raycore three-way aggregate parity" label rtol=THREE_WAY_AGGREGATE_RTOL atol=THREE_WAY_AGGREGATE_ATOL reference candidate
end

@testset "Raycore Metal validation" begin
    if !METAL_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to run optional Metal validation."
    else
        backend = _metal_backend()
        if backend !== nothing
            raw_stack_smoke_output = "/tmp/archimed_raycore_metal_raw_stack_smoke.tsv"
            isfile(raw_stack_smoke_output) && rm(raw_stack_smoke_output)

            @test ArchimedLight.RaycoreInterceptionBackend(backend=backend).config.workgroupsize == 256
            @test ArchimedLight.RaycoreInterceptionBackend(backend=backend, workgroupsize=64).config.workgroupsize == 64

            include(joinpath(@__DIR__, "..", "..", "scripts", "raycore_device_smoke.jl"))
            result = run_raycore_device_smoke(;
                backend=backend,
                first_order_area_atol=1e-5,
                first_order_area_rtol=1e-5,
                first_order_power_atol=1e-4,
                first_order_power_rtol=1e-5,
                scattering_atol=1e-4,
                scattering_rtol=1e-4,
            )
            @test result.scattering_eltype == Float32
            caps = result.reduction_capabilities
            if !caps.supports_atomics
                @test !caps.dense_edge_accumulation
                @test !caps.node_count_reduction
                @test !caps.sector_area_reduction
                @test !caps.fused_count_area_reduction
            end

            scene = _smoke_scene()
            models = _smoke_models()
            raw_options = ArchimedLight.LightOptions(
                turtle_sectors=6,
                all_in_turtle=false,
                scattering=true,
                pixel_size=0.01,
                toricity=false,
            )
            raw_sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
            raw_turtle = ArchimedLight.build_turtle(raw_options, raw_sky)
            raw_sector_idx = findfirst(sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0, raw_turtle.sectors)
            @test raw_sector_idx !== nothing
            if raw_sector_idx !== nothing
                prepared = ArchimedLight._prepare_interception_data(
                    scene,
                    models,
                    raw_options;
                    include_budget_maps=true,
                    include_raycore_instancing=true,
                )
                cpu_data = ArchimedLight._raycore_scene_data(
                    prepared,
                    ArchimedLight.RaycoreBackendConfig(
                        backend=KernelAbstractions.CPU(),
                        max_hits_per_pixel=8,
                        validate=true,
                    );
                    toricity=raw_options.toricity,
                )
                metal_data = ArchimedLight._raycore_scene_data(
                    prepared,
                    ArchimedLight.RaycoreBackendConfig(
                        backend=backend,
                        max_hits_per_pixel=8,
                        validate=true,
                    );
                    toricity=raw_options.toricity,
                )
                raw_comparison = ArchimedLight._raycore_raw_stack_comparison(
                    cpu_data,
                    metal_data,
                    raw_turtle.sectors[raw_sector_idx].direction,
                    raw_options,
                )
                @test raw_comparison.ok
                @test raw_comparison.reference_hits == raw_comparison.candidate_hits
                @test raw_comparison.hit_ratio == 1.0
                @test raw_comparison.reference_overflow == raw_comparison.candidate_overflow == 0
            end

            simple_config = joinpath(@__DIR__, "..", "fast_fixtures", "simpleplant_16_toric", "input", "config.yml")
            simple_options, simple_scene, _simple_meteo, simple_models = ArchimedLight.read_config(simple_config)
            simple_options = ArchimedLight.LightOptions(
                simple_options;
                cache_radiation=true,
                scattering=true,
                pixel_size=0.40,
                turtle_sectors=6,
            )
            simple_sky = ArchimedLight.SkyState(180.0, 55.0, 420.0, 220.0, 0.70, 0.30)
            normal_simple = _run_single_sky(
                simple_scene,
                simple_models,
                simple_options,
                simple_sky;
                interception_backend=:raster_cpu,
            )
            ka_interception = ArchimedLight.RaycoreInterceptionBackend(
                backend=KernelAbstractions.CPU(),
                max_hits_per_pixel=64,
                validate=true,
            )
            ka_simple = _run_single_sky(
                simple_scene,
                simple_models,
                simple_options,
                simple_sky;
                interception_backend=ka_interception,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ka_interception),
            )
            metal_interception = ArchimedLight.RaycoreInterceptionBackend(
                backend=backend,
                max_hits_per_pixel=64,
                validate=true,
            )
            metal_simple = _run_single_sky(
                simple_scene,
                simple_models,
                simple_options,
                simple_sky;
                interception_backend=metal_interception,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(metal_interception),
            )
            @test normal_simple.sim.cache.resolved_interception_backend isa ArchimedLight.RasterCPUBackend
            @test ka_simple.sim.cache.resolved_interception_backend isa ArchimedLight.RaycoreInterceptionBackend
            @test metal_simple.sim.cache.resolved_interception_backend isa ArchimedLight.RaycoreInterceptionBackend
            @test ka_simple.sim.cache.fallback_reason == :none
            @test metal_simple.sim.cache.fallback_reason == :none
            _assert_totals_close(ka_simple.totals, normal_simple.totals; label=:raycore_ka_cpu)
            _assert_totals_close(metal_simple.totals, normal_simple.totals; label=:raycore_metal_gpu)
            simple_rows = _assert_raw_stack_sweep_ok(
                simple_scene,
                simple_models,
                simple_options,
                simple_sky,
                backend,
            )
            @test [row.max_hits for row in simple_rows] == [32, 64, 128, 256]

            coffee_config = joinpath(@__DIR__, "..", "..", "example_2", "config.yml")
            coffee_paving = parse(Int, get(ENV, "ARCHIMEDLIGHT_TEST_GPU_COFFEE_RAW_PAVING", "64"))
            coffee_options, coffee_scene, coffee_meteo, coffee_models =
                ArchimedLight.read_config(coffee_config; plot_paving_override=coffee_paving)
            coffee_options = ArchimedLight.LightOptions(
                coffee_options;
                cache_radiation=true,
                scattering=true,
                meteo_range=nothing,
                pixel_size=parse(Float64, get(ENV, "ARCHIMEDLIGHT_TEST_GPU_COFFEE_RAW_PIXEL_SIZE", "0.04")),
                turtle_sectors=6,
            )
            coffee_rows = ArchimedLight.prepare_meteo(coffee_meteo, coffee_options).rows
            coffee_step = min(parse(Int, get(ENV, "ARCHIMEDLIGHT_TEST_GPU_COFFEE_RAW_STEP", "1")), length(coffee_rows))
            coffee_sky = ArchimedLight.compute_sky(coffee_rows[coffee_step], coffee_options)
            coffee_raw_rows = _raw_stack_sweep(
                coffee_scene,
                coffee_models,
                coffee_options,
                coffee_sky,
                backend;
                assert_ok=false,
            )
            @test [row.max_hits for row in coffee_raw_rows] == [32, 64, 128, 256]
            @test all(row -> row.reference_hits > 0, coffee_raw_rows)
            @test all(row -> row.candidate_hits > 0, coffee_raw_rows)
            @test all(row -> isfinite(row.hit_ratio), coffee_raw_rows)
            @test all(row -> isfinite(row.occupied_ratio), coffee_raw_rows)
            @test all(row -> row.reference_overflow >= 0 && row.candidate_overflow >= 0, coffee_raw_rows)
            _append_raw_stack_smoke_rows(raw_stack_smoke_output, "coffee_reduced", coffee_raw_rows)
            _assert_unresolved_raw_failures_are_diagnosed(
                coffee_raw_rows;
                require_missing_hit=any(row -> row.candidate_blas_depth_over_stack_capacity, coffee_raw_rows),
            )
            @test all(row -> row.raycore_traversal_stack_capacity == RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY, coffee_raw_rows)
            @test all(row -> row.ok, coffee_raw_rows)
            @test all(row -> row.diagnosis == :exact, coffee_raw_rows)
            @info "Raycore Metal coffee raw stack diagnostic" coffee_raw_rows

            coffee_cpu_static_rows = _raw_stack_sweep(
                coffee_scene,
                coffee_models,
                coffee_options,
                coffee_sky,
                backend;
                max_hits_values=(64,),
                candidate_static_tlas_from_reference=true,
                assert_ok=false,
            )
            @test length(coffee_cpu_static_rows) == 1
            @test only(coffee_cpu_static_rows).candidate_static_tlas_source == :reference_cpu
            @test only(coffee_cpu_static_rows).reference_geometry_mode == only(coffee_cpu_static_rows).candidate_geometry_mode
            @test only(coffee_cpu_static_rows).reference_blas_depth == only(coffee_cpu_static_rows).candidate_blas_depth
            @test only(coffee_cpu_static_rows).reference_blas_node_count == only(coffee_cpu_static_rows).candidate_blas_node_count
            _append_raw_stack_smoke_rows(raw_stack_smoke_output, "coffee_reduced_cpu_static_tlas", coffee_cpu_static_rows)
            @test only(coffee_cpu_static_rows).ok
            @test only(coffee_cpu_static_rows).diagnosis == :exact
            @info "Raycore Metal coffee raw stack diagnostic with CPU-built static TLAS" coffee_cpu_static_rows

            coffee_no_dedupe_rows = _raw_stack_sweep(
                coffee_scene,
                coffee_models,
                coffee_options,
                coffee_sky,
                backend;
                max_hits_values=(64,),
                hit_epsilon=-1.0f0,
                assert_ok=false,
            )
            @test length(coffee_no_dedupe_rows) == 1
            @test only(coffee_no_dedupe_rows).hit_epsilon < 0
            @test only(coffee_no_dedupe_rows).reference_hits > 0
            @test only(coffee_no_dedupe_rows).candidate_hits > 0
            _append_raw_stack_smoke_rows(raw_stack_smoke_output, "coffee_reduced_no_dedupe", coffee_no_dedupe_rows)
            @test only(coffee_no_dedupe_rows).ok
            @test only(coffee_no_dedupe_rows).diagnosis == :exact
            _assert_unresolved_raw_failures_are_diagnosed(coffee_no_dedupe_rows)
            @info "Raycore Metal coffee raw stack diagnostic without duplicate suppression" coffee_no_dedupe_rows

            coffee_chunked_rows = reduce(vcat, [
                _raw_stack_sweep(
                    coffee_scene,
                    coffee_models,
                    coffee_options,
                    coffee_sky,
                    backend;
                    max_hits_values=(64,),
                    face_chunk_limit=limit,
                    assert_ok=false,
                )
                for limit in (1536, 1024, 768, 512)
            ])
            @test [row.face_chunk_limit for row in coffee_chunked_rows] == [1536, 1024, 768, 512]
            @test all(row -> row.reference_geometry_mode == :chunked_merged_mesh, coffee_chunked_rows)
            @test all(row -> row.candidate_geometry_mode == :chunked_merged_mesh, coffee_chunked_rows)
            @test all(row -> row.reference_hits > 0, coffee_chunked_rows)
            @test all(row -> row.candidate_hits > 0, coffee_chunked_rows)
            _append_raw_stack_smoke_rows(raw_stack_smoke_output, "coffee_forced_chunk", coffee_chunked_rows)
            @test all(row -> row.ok, coffee_chunked_rows)
            @test all(row -> row.diagnosis == :exact, coffee_chunked_rows)
            @info "Raycore Metal coffee prechunked raw stack diagnostic" coffee_chunked_rows

            strict_interception = ArchimedLight.RaycoreInterceptionBackend(
                backend=backend,
                max_hits_per_pixel=64,
                validate=true,
            )
            strict_cache = ArchimedLight.prepare_light_cache(
                coffee_scene,
                coffee_models,
                coffee_options;
                interception_backend=strict_interception,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(strict_interception),
            )
            @test strict_cache.fallback_reason == :none
            @test strict_cache.resolved_interception_backend isa ArchimedLight.RaycoreInterceptionBackend
            @test strict_cache.raycore_data !== nothing
            @test strict_cache.raycore_data.geometry_mode in (:merged_mesh, :chunked_merged_mesh)
            strict_shape = ArchimedLight._raycore_scene_shape_summary(strict_cache.raycore_data)
            @test strict_shape.max_blas_depth <= strict_shape.raycore_traversal_stack_capacity
            @test !strict_shape.max_blas_depth_over_stack_capacity
            @test ArchimedLight._raycore_stack_trace_validation(strict_cache.raycore_data, coffee_options).ok
            @info "Raycore Metal strict coffee cache survived full-stack validation" geometry_mode=strict_cache.raycore_data.geometry_mode instances=length(strict_cache.raycore_data.tlas.instances) max_blas_depth=strict_shape.max_blas_depth traversal_stack_capacity=strict_shape.raycore_traversal_stack_capacity

            complexity_script = joinpath(@__DIR__, "..", "..", "benchmark", "scene_complexity_sweep.jl")
            complexity_env = Dict(
                "ARCHIMEDLIGHT_BENCH_METAL" => get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", nothing),
                "ARCHIMEDLIGHT_BENCH_VALIDATE" => get(ENV, "ARCHIMEDLIGHT_BENCH_VALIDATE", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_BACKENDS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_BACKENDS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_SCENES" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SCENES", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_SECTORS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SECTORS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_SAMPLES" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SAMPLES", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_SIMPLE_PIXELS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SIMPLE_PIXELS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PIXELS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PIXELS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PAVING" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PAVING", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_COMPLEX_PIXELS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_COMPLEX_PIXELS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_ROWS" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_ROWS", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_PLANTS_PER_ROW" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_PLANTS_PER_ROW", nothing),
                "ARCHIMEDLIGHT_COMPLEXITY_OUTPUT" => get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_OUTPUT", nothing),
                "ARCHIMEDLIGHT_BENCH_MAX_HITS" => get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", nothing),
            )
            try
                agripv_opf = joinpath(
                    get(
                        ENV,
                        "ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT",
                        "/Users/rvezy/Documents/cirad/Articles_Rapports_Communications/Conférences/2026_fspm/coutellier_agripv",
                    ),
                    "0_simulations",
                    "archicrop",
                    "wheat",
                    "plant_1995-06-24.opf",
                )
                complexity_scenes = isfile(agripv_opf) ? "simple,coffee,complex" : "simple,coffee"
                ENV["ARCHIMEDLIGHT_BENCH_METAL"] = "required"
                ENV["ARCHIMEDLIGHT_BENCH_VALIDATE"] = "1"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_BACKENDS"] = "raycore_metal_gpu"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_SCENES"] = complexity_scenes
                ENV["ARCHIMEDLIGHT_COMPLEXITY_SECTORS"] = "6"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_SAMPLES"] = "1"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_SIMPLE_PIXELS"] = "0.40"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PIXELS"] = "0.04"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PAVING"] = "2500"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_COMPLEX_PIXELS"] = "0.20"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_ROWS"] = "1"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_PLANTS_PER_ROW"] = "2"
                ENV["ARCHIMEDLIGHT_COMPLEXITY_OUTPUT"] = "/tmp/archimed_raycore_metal_strict_scene_complexity_smoke.tsv"
                ENV["ARCHIMEDLIGHT_BENCH_MAX_HITS"] = "64"

                complexity_rows = include(complexity_script)
                @test isfile(ENV["ARCHIMEDLIGHT_COMPLEXITY_OUTPUT"])
                @test length(complexity_rows) == (occursin("complex", complexity_scenes) ? 3 : 2)
                @test all(row -> row.backend == "raycore_metal_gpu", complexity_rows)
                @test all(row -> row.status == "ok", complexity_rows)
                @test all(row -> row.resolved_backend == "RaycoreInterceptionBackend", complexity_rows)
                @test all(row -> row.fallback_reason == :none, complexity_rows)
                @test all(row -> row.geometry_mode != :none, complexity_rows)
                @test all(row -> row.raycore_tlas_instances > 0, complexity_rows)
                coffee_complexity = only(row for row in complexity_rows if row.scene == "coffee_docs")
                @test coffee_complexity.pixel_size == 0.04
                @test coffee_complexity.faces > 0
                @test coffee_complexity.faces > 100_000
                @test coffee_complexity.raycore_chunked
                @test coffee_complexity.raycore_tlas_instances > coffee_complexity.raycore_tlas_geometries
                fake_validation = (
                    hit_ratio=0.5,
                    occupied_ratio=0.75,
                    reference_hits=10,
                    raycore_hits=5,
                    reference_occupied=4,
                    raycore_occupied=3,
                    reference_overflow=false,
                    raycore_overflow=false,
                    directions_tested=1,
                    direction_count=3,
                )
                fake_error = ArchimedLight.RaycoreValidationError(
                    :raycore_stack_trace_validation,
                    :complexity,
                    "fabricated strict validation failure",
                    fake_validation,
                )
                fake_row = _error_case_result(fake_error, (name="raycore_metal_gpu",))
                @test fake_row.status == "error"
                @test fake_row.resolved_backend == "RaycoreInterceptionBackend"
                @test fake_row.fallback_reason == :raycore_stack_trace_validation
                @test fake_row.error_reason == :raycore_stack_trace_validation
                @test fake_row.error_stage == :complexity
                @test fake_row.validation_hit_ratio == 0.5
                @test fake_row.validation_reference_hits == 10
                if occursin("complex", complexity_scenes)
                    @test any(row -> row.scene == "agripv_wheat_panel", complexity_rows)
                else
                    @test_skip "Agrivoltaics OPF not available; strict scene-complexity benchmark smoke covered simple and coffee scenes only."
                end
                @info "Raycore Metal strict scene-complexity benchmark smoke stayed on Raycore" scenes=[row.scene for row in complexity_rows] geometry_modes=[row.geometry_mode for row in complexity_rows] instances=[row.raycore_tlas_instances for row in complexity_rows]
            finally
                for (key, value) in complexity_env
                    if value === nothing
                        delete!(ENV, key)
                    else
                        ENV[key] = value
                    end
                end
            end

            raw_agripv_opf = joinpath(
                get(
                    ENV,
                    "ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT",
                    "/Users/rvezy/Documents/cirad/Articles_Rapports_Communications/Conférences/2026_fspm/coutellier_agripv",
                ),
                "0_simulations",
                "archicrop",
                "wheat",
                "plant_1995-06-24.opf",
            )
            agripv_raw_flag = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_GPU_AGRIPV_RAW", isfile(raw_agripv_opf) ? "1" : ""))
            agripv_raw_requested = agripv_raw_flag in ("1", "true", "yes", "on", "required", "force", "error")
            agripv_raw_required = agripv_raw_flag in ("required", "force", "error")
            if isfile(raw_agripv_opf) && agripv_raw_requested
                raw_script = joinpath(@__DIR__, "..", "..", "scripts", "raycore_raw_stack_compare.jl")
                raw_env = Dict(
                    "ARCHIMEDLIGHT_RAW_STACK_BACKEND" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_BACKEND", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_WORKLOAD" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_WORKLOAD", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_SECTOR" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_SECTOR", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_ROWS" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_ROWS", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_PLANTS_PER_ROW" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_PLANTS_PER_ROW", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY", nothing),
                    "ARCHIMEDLIGHT_RAW_STACK_OUTPUT" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_OUTPUT", nothing),
                )
                try
                    ENV["ARCHIMEDLIGHT_RAW_STACK_BACKEND"] = "metal"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_WORKLOAD"] = "agripv"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP"] = "64"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS"] = "auto"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_SECTOR"] = "first_scattering"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE"] = "0.20"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_AGRIPV_ROWS"] = "1"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_AGRIPV_PLANTS_PER_ROW"] = "2"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX"] = "6"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY"] = "8"
                    ENV["ARCHIMEDLIGHT_RAW_STACK_OUTPUT"] = ""

                    raw_mod = Module(:RaycoreRawAgripvMetalSmoke)
                    Base.include(raw_mod, raw_script)
                    agripv_raw_rows = try
                        Base.invokelatest(raw_mod.main)
                    catch err
                        @error "Raycore Metal agripv raw stack diagnostic failed" exception=(err, catch_backtrace())
                        rethrow()
                    end
                    @test length(agripv_raw_rows) == 1
                    agripv_raw = only(agripv_raw_rows)
                    @test agripv_raw.workload == "agripv"
                    @test agripv_raw.backend == "metal"
                    @test agripv_raw.max_hits == 64
                    @test isfinite(agripv_raw.hit_ratio)
                    @test isfinite(agripv_raw.occupied_ratio)
                    @test agripv_raw.reference_hits > 0
                    @test agripv_raw.candidate_hits > 0
                    @test agripv_raw.reference_overflow >= 0
                    @test agripv_raw.candidate_overflow >= 0
                    @test agripv_raw.reference_matches_candidate_topology
                    @test agripv_raw.reference_geometry_mode == agripv_raw.candidate_geometry_mode
                    @test agripv_raw.reference_blas_depth == agripv_raw.candidate_blas_depth
                    @test agripv_raw.reference_blas_node_count == agripv_raw.candidate_blas_node_count
                    @test agripv_raw.reference_tlas_depth == agripv_raw.candidate_tlas_depth
                    @test agripv_raw.reference_instances == agripv_raw.candidate_instances
                    @test agripv_raw.raycore_traversal_stack_capacity == RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY
                    _append_raw_stack_smoke_rows(raw_stack_smoke_output, "agripv_reduced", agripv_raw_rows)
                    if agripv_raw.ok
                        @test agripv_raw.mismatch === nothing
                        @test agripv_raw.diagnosis == :exact
                    else
                        @test agripv_raw.mismatch !== nothing
                        @test agripv_raw.mismatch_summary != "none"
                        _assert_unresolved_raw_failures_are_diagnosed(agripv_raw_rows)
                        @test_broken agripv_raw.ok
                    end
                    @info "Raycore Metal agripv raw stack diagnostic" agripv_raw_rows
                finally
                    for (key, value) in raw_env
                        if value === nothing
                            delete!(ENV, key)
                        else
                            ENV[key] = value
                        end
                    end
                end
            else
                if agripv_raw_required
                    error("ARCHIMEDLIGHT_TEST_GPU_AGRIPV_RAW was required, but the agrivoltaics OPF was not available.")
                end
                @test_skip "Set ARCHIMEDLIGHT_TEST_GPU_AGRIPV_RAW=1 to run the slower Metal agripv raw-stack diagnostic."
            end
        end
    end
end

include("local_realistic_tests.jl")
