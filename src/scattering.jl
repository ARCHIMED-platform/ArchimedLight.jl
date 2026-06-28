function _dict_zero(node_ids)
    Dict{Int,Float64}(nid => 0.0 for nid in node_ids)
end

struct ScatteringTopologyCache
    pair_counts::ScatteringPairCounts
    sun_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
end

function _copy_node_values(source::Dict{Int,Float64}, node_ids)
    Dict{Int,Float64}(nid => get(source, nid, 0.0) for nid in node_ids)
end

function _sum_dict_values(d::Dict{Int,Float64})
    s = 0.0
    for v in values(d)
        s += v
    end
    s
end

"""
    _pack_scattering_edge(to, from) -> UInt64

Pack one scattering edge `(to, from)` into a single 64-bit key.

We reserve 32 bits for each node id:

    [to ........ 32 bits][from ...... 32 bits]

The explicit `UInt32` conversion constrains each id to one half of the packed key, and the
`UInt64` widening makes the shift/OR operations safe. This is cheaper to hash in the topology
builder than a `Tuple{Int,Int}` key.
"""
@inline function _pack_scattering_edge(to::Int, from::Int)
    return (UInt64(UInt32(to)) << 32) | UInt64(UInt32(from))
end

# Recover the upper and lower 32-bit halves of a packed `(to, from)` edge key.
@inline _unpack_scattering_to(edge::UInt64) = Int(UInt32(edge >> 32))
@inline _unpack_scattering_from(edge::UInt64) = Int(UInt32(edge & 0xffffffff))

function _edge_counts_from_packed(edge_counts::Dict{UInt64,Int})
    isempty(edge_counts) && return ScatteringPairCounts(Int[], Int[], Int[])

    # Sort once when materializing the final graph so iteration order is deterministic.
    packed_edges = sort!(collect(keys(edge_counts)))
    to_nodes = Int[]
    from_nodes = Int[]
    counts = Int[]
    sizehint!(to_nodes, length(packed_edges))
    sizehint!(from_nodes, length(packed_edges))
    sizehint!(counts, length(packed_edges))

    for edge in packed_edges
        count = edge_counts[edge]
        if count > 0
            push!(to_nodes, _unpack_scattering_to(edge))
            push!(from_nodes, _unpack_scattering_from(edge))
            push!(counts, count)
        end
    end
    return ScatteringPairCounts(to_nodes, from_nodes, counts)
end

function _merge_packed_edge_keys!(edge_counts::Dict{UInt64,Int}, edge_keys::AbstractVector{UInt64})
    @inbounds for edge in edge_keys
        edge == 0 && continue
        edge_counts[edge] = get(edge_counts, edge, 0) + 1
    end
    return edge_counts
end

@inline function _add_packed_edge_count!(edge_counts::Dict{UInt64,Int}, to::Int, from::Int)
    edge = _pack_scattering_edge(to, from)
    edge_counts[edge] = get(edge_counts, edge, 0) + 1
    return nothing
end

function _accept_scattering_link(node_group::Dict{Int,String}, to::Int, from::Int)
    !(get(node_group, to, "") == "pavement" && get(node_group, from, "") == "pavement")
end

mutable struct ScatteringStackScratch
    nearest_nonvirtual_below::Vector{Int}
end

ScatteringStackScratch() = ScatteringStackScratch(Int[])

function _accumulate_sun_hits!(
    sun_hits::Dict{Int,Int},
    projection::DirectionProjectionResult,
    node_ids::Union{Nothing,Vector{Int}}=nothing,
)
    for (nid, h) in projection.node_hits
        sun_hits[nid] = get(sun_hits, nid, 0) + h
    end
    return nothing
end

function _accumulate_sun_hits!(
    sun_hits::Dict{Int,Int},
    projection::DenseDirectionProjectionResult,
    node_ids::Union{Nothing,Vector{Int}},
)
    node_ids === nothing && error("DenseDirectionProjectionResult sun-hit accumulation requires node_ids")
    @inbounds for i in eachindex(projection.node_hits)
        h = projection.node_hits[i]
        h == 0 && continue
        nid = node_ids[i]
        sun_hits[nid] = get(sun_hits, nid, 0) + h
    end
    return nothing
end

function _stack_transfer_pairs!(
    edge_counts::Dict{UInt64,Int},
    stack,
    virtual_nodes::Set{Int},
    scratch::ScatteringStackScratch,
    node_group::Dict{Int,String},
)
    n_hits = length(stack)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        nid = _stack_hit_node(stack, h)
        if !(nid in virtual_nodes)
            nearest_below = nid
        end
        below[h] = nearest_below
    end

    nearest_above = 0
    @inbounds for h in 1:(n_hits - 1)
        to_above = _stack_hit_node(stack, h)
        if !(to_above in virtual_nodes)
            nearest_above = to_above
        end

        from_below = below[h + 1]
        if from_below != 0 && _accept_scattering_link(node_group, to_above, from_below)
            _add_packed_edge_count!(edge_counts, to_above, from_below)
        end

        to_below = _stack_hit_node(stack, h + 1)
        from_above = nearest_above
        if from_above != 0 && _accept_scattering_link(node_group, to_below, from_above)
            _add_packed_edge_count!(edge_counts, to_below, from_above)
        end
    end
    return nothing
end

function _accumulate_scattering_counts!(
    edge_counts::Dict{UInt64,Int},
    sun_hits::Dict{Int,Int},
    sector::TurtleSector,
    projection::DirectionProjectionResult,
    virtual_nodes::Set{Int},
    node_group::Dict{Int,String},
    scratch::ScatteringStackScratch;
    node_ids::Union{Nothing,Vector{Int}}=nothing,
    stacks_sorted::Bool=false,
)
    if sector.source == :sun
        _accumulate_sun_hits!(sun_hits, projection, node_ids)
        return nothing
    end

    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _stack_transfer_pairs!(edge_counts, stack, virtual_nodes, scratch, node_group)
    end
    return nothing
end

@inline function _accumulate_scattering_counts_dense!(
    edge_counts::Dict{UInt64,Int},
    stack,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch,
    node_ids::Vector{Int},
)
    n_hits = length(stack)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        node_idx = _stack_hit_node(stack, h)
        if !virtual_node_mask[node_idx]
            nearest_below = node_idx
        end
        below[h] = nearest_below
    end

    nearest_above = 0
    @inbounds for h in 1:(n_hits - 1)
        to_above_idx = _stack_hit_node(stack, h)
        if !virtual_node_mask[to_above_idx]
            nearest_above = to_above_idx
        end

        from_below_idx = below[h + 1]
        if from_below_idx != 0
            to_above = node_ids[to_above_idx]
            from_below = node_ids[from_below_idx]
            if !(pavement_node_mask[to_above_idx] && pavement_node_mask[from_below_idx])
                _add_packed_edge_count!(edge_counts, to_above, from_below)
            end
        end

        to_below_idx = _stack_hit_node(stack, h + 1)
        from_above_idx = nearest_above
        if from_above_idx != 0
            to_below = node_ids[to_below_idx]
            from_above = node_ids[from_above_idx]
            if !(pavement_node_mask[to_below_idx] && pavement_node_mask[from_above_idx])
                _add_packed_edge_count!(edge_counts, to_below, from_above)
            end
        end
    end
    return nothing
end

@inline function _accumulate_scattering_counts_dense!(
    edge_counts::Dict{UInt64,Int},
    stack::FlatPixelHitStack,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch,
    node_ids::Vector{Int},
)
    n_hits = length(stack)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    start = stack.parent.starts[stack.pixel_idx]
    nodes = stack.parent.nodes

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        node_idx = Int(nodes[start + h - 1])
        if !virtual_node_mask[node_idx]
            nearest_below = node_idx
        end
        below[h] = nearest_below
    end

    @inbounds begin
        current_idx = Int(nodes[start])
        current_node = node_ids[current_idx]
        current_is_pavement = pavement_node_mask[current_idx]
        nearest_above = virtual_node_mask[current_idx] ? 0 : current_idx
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[start + h])
            next_node = node_ids[next_idx]
            next_is_pavement = pavement_node_mask[next_idx]

            from_above_idx = nearest_above
            from_below_idx = below[h + 1]
            if from_below_idx != 0
                from_below_is_pavement = pavement_node_mask[from_below_idx]
                if !(current_is_pavement && from_below_is_pavement)
                    from_below = node_ids[from_below_idx]
                    _add_packed_edge_count!(edge_counts, current_node, from_below)
                end
            end

            if from_above_idx != 0
                from_above_is_pavement = pavement_node_mask[from_above_idx]
                if !(next_is_pavement && from_above_is_pavement)
                    from_above = node_ids[from_above_idx]
                    _add_packed_edge_count!(edge_counts, next_node, from_above)
                end
            end

            if !virtual_node_mask[next_idx]
                nearest_above = next_idx
            end
            current_idx = next_idx
            current_node = next_node
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

KernelAbstractions.@kernel function _raycore_scattering_edge_keys_kernel!(
    edge_keys,
    counts,
    nodes,
    virtual_node_mask,
    pavement_node_mask,
    node_ids,
    max_hits::Int,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        max_edges = 2 * (max_hits - 1)
        out_base = (pixel_idx - 1) * max_edges
        for slot in 1:max_edges
            edge_keys[out_base+slot] = UInt64(0)
        end

        n_hits = Int(counts[pixel_idx])
        if n_hits > 1
            stack_base = (pixel_idx - 1) * max_hits
            nearest_above_idx = 0
            out_slot = 0
            for h in 1:(n_hits - 1)
                current_idx = Int(nodes[stack_base+h])
                if !virtual_node_mask[current_idx]
                    nearest_above_idx = current_idx
                end

                from_below_idx = 0
                for k in (h + 1):n_hits
                    candidate_idx = Int(nodes[stack_base+k])
                    if !virtual_node_mask[candidate_idx]
                        from_below_idx = candidate_idx
                        break
                    end
                end

                if from_below_idx != 0 && !(pavement_node_mask[current_idx] && pavement_node_mask[from_below_idx])
                    out_slot += 1
                    to_node = node_ids[current_idx]
                    from_node = node_ids[from_below_idx]
                    edge_keys[out_base+out_slot] = _pack_scattering_edge(to_node, from_node)
                end

                next_idx = Int(nodes[stack_base+h+1])
                if nearest_above_idx != 0 && !(pavement_node_mask[next_idx] && pavement_node_mask[nearest_above_idx])
                    out_slot += 1
                    to_node = node_ids[next_idx]
                    from_node = node_ids[nearest_above_idx]
                    edge_keys[out_base+out_slot] = _pack_scattering_edge(to_node, from_node)
                end
            end
        end
    end
end

function _raycore_scattering_edge_keys_from_traced_stacks(
    data::RaycoreSceneData,
    traced,
)
    n_pixels = length(traced.counts)
    max_hits = traced.max_hits
    max_edges = 2 * (max_hits - 1)
    (n_pixels == 0 || max_edges <= 0) && return UInt64[]

    backend = data.backend
    geometry = data.prepared.geometry
    counts_dev = KernelAbstractions.allocate(backend, Int32, n_pixels)
    nodes_dev = KernelAbstractions.allocate(backend, UInt32, length(traced.nodes))
    virtual_dev = KernelAbstractions.allocate(backend, Bool, length(data.prepared.virtual_node_mask))
    pavement_dev = KernelAbstractions.allocate(backend, Bool, length(geometry.pavement_node_mask))
    node_ids_dev = KernelAbstractions.allocate(backend, Int, length(geometry.node_ids))
    edge_keys_dev = KernelAbstractions.allocate(backend, UInt64, n_pixels * max_edges)

    KernelAbstractions.copyto!(backend, counts_dev, traced.counts)
    KernelAbstractions.copyto!(backend, nodes_dev, traced.nodes)
    KernelAbstractions.copyto!(backend, virtual_dev, data.prepared.virtual_node_mask)
    KernelAbstractions.copyto!(backend, pavement_dev, geometry.pavement_node_mask)
    KernelAbstractions.copyto!(backend, node_ids_dev, geometry.node_ids)

    kernel = _raycore_scattering_edge_keys_kernel!(backend, data.workgroupsize)
    kernel(
        edge_keys_dev,
        counts_dev,
        nodes_dev,
        virtual_dev,
        pavement_dev,
        node_ids_dev,
        max_hits;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)

    return Array(edge_keys_dev)
end

function _raycore_dense_edge_matrix_fits(data::RaycoreSceneData)
    n_nodes = length(data.prepared.geometry.node_ids)
    Int128(n_nodes) * Int128(n_nodes) <= typemax(Int) || return false
    bytes = Int128(n_nodes) * Int128(n_nodes) * Int128(sizeof(Int))
    return bytes <= data.dense_edge_limit_bytes
end

function _raycore_dense_edge_accumulation_supported(data::RaycoreSceneData)
    return _raycore_dense_edge_matrix_fits(data) && KernelAbstractions.supports_atomics(data.backend)
end

KernelAbstractions.@kernel function _raycore_scattering_dense_counts_kernel!(
    dense_counts,
    counts,
    nodes,
    virtual_node_mask,
    pavement_node_mask,
    max_hits::Int,
    n_nodes::Int,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        n_hits = Int(counts[pixel_idx])
        if n_hits > 1
            stack_base = (pixel_idx - 1) * max_hits
            nearest_above_idx = 0
            for h in 1:(n_hits - 1)
                current_idx = Int(nodes[stack_base+h])
                if !virtual_node_mask[current_idx]
                    nearest_above_idx = current_idx
                end

                from_below_idx = 0
                for k in (h + 1):n_hits
                    candidate_idx = Int(nodes[stack_base+k])
                    if !virtual_node_mask[candidate_idx]
                        from_below_idx = candidate_idx
                        break
                    end
                end

                if from_below_idx != 0 &&
                   !(pavement_node_mask[current_idx] && pavement_node_mask[from_below_idx])
                    dense_idx = (current_idx - 1) * n_nodes + from_below_idx
                    @atomic :monotonic dense_counts[dense_idx] += 1
                end

                next_idx = Int(nodes[stack_base+h+1])
                if nearest_above_idx != 0 &&
                   !(pavement_node_mask[next_idx] && pavement_node_mask[nearest_above_idx])
                    dense_idx = (next_idx - 1) * n_nodes + nearest_above_idx
                    @atomic :monotonic dense_counts[dense_idx] += 1
                end
            end
        end
    end
end

function _merge_dense_edge_counts!(
    edge_counts::Dict{UInt64,Int},
    dense_counts::AbstractVector{<:Integer},
    node_ids::AbstractVector{Int},
)
    n_nodes = length(node_ids)
    @inbounds for pair_idx in eachindex(dense_counts)
        count = Int(dense_counts[pair_idx])
        count == 0 && continue
        to_idx = ((pair_idx - 1) ÷ n_nodes) + 1
        from_idx = ((pair_idx - 1) % n_nodes) + 1
        edge = _pack_scattering_edge(node_ids[to_idx], node_ids[from_idx])
        edge_counts[edge] = get(edge_counts, edge, 0) + count
    end
    return edge_counts
end

function _raycore_scattering_dense_counts_from_traced_stacks(
    data::RaycoreSceneData,
    traced,
)
    n_nodes = length(data.prepared.geometry.node_ids)
    n_pixels = length(traced.counts)
    n_pairs = n_nodes * n_nodes
    (n_pixels == 0 || n_nodes == 0) && return Int[]

    backend = data.backend
    geometry = data.prepared.geometry
    counts_dev = KernelAbstractions.allocate(backend, Int32, n_pixels)
    nodes_dev = KernelAbstractions.allocate(backend, UInt32, length(traced.nodes))
    virtual_dev = KernelAbstractions.allocate(backend, Bool, length(data.prepared.virtual_node_mask))
    pavement_dev = KernelAbstractions.allocate(backend, Bool, length(geometry.pavement_node_mask))
    dense_counts_dev = KernelAbstractions.allocate(backend, Int, n_pairs)

    KernelAbstractions.copyto!(backend, counts_dev, traced.counts)
    KernelAbstractions.copyto!(backend, nodes_dev, traced.nodes)
    KernelAbstractions.copyto!(backend, virtual_dev, data.prepared.virtual_node_mask)
    KernelAbstractions.copyto!(backend, pavement_dev, geometry.pavement_node_mask)
    KernelAbstractions.copyto!(backend, dense_counts_dev, zeros(Int, n_pairs))

    kernel = _raycore_scattering_dense_counts_kernel!(backend, data.workgroupsize)
    kernel(
        dense_counts_dev,
        counts_dev,
        nodes_dev,
        virtual_dev,
        pavement_dev,
        traced.max_hits,
        n_nodes;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)

    return Array(dense_counts_dev)
end

function _accumulate_scattering_counts!(
    edge_counts::Dict{UInt64,Int},
    sun_hits::Dict{Int,Int},
    sector::TurtleSector,
    projection::DenseDirectionProjectionResult,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch;
    node_ids::Union{Nothing,Vector{Int}}=nothing,
    stacks_sorted::Bool=false,
)
    node_ids === nothing && error("DenseDirectionProjectionResult scattering accumulation requires node_ids")
    if sector.source == :sun
        _accumulate_sun_hits!(sun_hits, projection, node_ids)
        return nothing
    end

    for stack in values(projection.pixel_hits)
        n_hits = length(stack)
        n_hits <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _accumulate_scattering_counts_dense!(
            edge_counts,
            stack,
            virtual_node_mask,
            pavement_node_mask,
            scratch,
            node_ids,
        )
    end
    return nothing
end

function _pair_counts_for_scattering(scene::PlantGeom.SceneGeometry, models::LightModels, turtle::TurtleGrid, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    virtual_nodes = _virtual_sensor_node_ids(geometry.node_group, geometry.node_type, models)
    cache_ctx = _projection_cache_context(geometry.vertices, geometry.faces, geometry.face2node, geometry.plotbox, options)
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for sector in turtle.sectors
        projection =
            _direction_projection_cached(geometry.vertices, geometry.faces, geometry.face2node, sector.direction, options, geometry.plotbox, cache_ctx)
        _accumulate_scattering_counts!(
            edge_counts,
            sun_hits,
            sector,
            projection,
            virtual_nodes,
            geometry.node_group,
            scratch,
        )
    end

    return _edge_counts_from_packed(edge_counts), sun_hits, geometry.node_ids, geometry.node_group
end

function _pair_counts_from_projections(
    turtle::TurtleGrid,
    projections::AbstractVector,
    virtual_nodes::Set{Int},
    node_group::Dict{Int,String},
    node_ids::Union{Nothing,Vector{Int}}=nothing,
    pavement_node_mask::Union{Nothing,Vector{Bool}}=nothing,
    virtual_node_mask::Union{Nothing,Vector{Bool}}=nothing,
    stacks_sorted::Bool=false,
)
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for i in eachindex(turtle.sectors)
        projection = projections[i]
        if projection isa DenseDirectionProjectionResult
            _accumulate_scattering_counts!(
                edge_counts,
                sun_hits,
                turtle.sectors[i],
                projection,
                virtual_node_mask,
                pavement_node_mask,
                scratch;
                node_ids=node_ids,
                stacks_sorted=stacks_sorted,
            )
        else
            _accumulate_scattering_counts!(
                edge_counts,
                sun_hits,
                turtle.sectors[i],
                projection,
                virtual_nodes,
                node_group,
                scratch;
                node_ids=node_ids,
                stacks_sorted=stacks_sorted,
            )
        end
    end

    return _edge_counts_from_packed(edge_counts), sun_hits
end

function _pair_counts_from_streamed_projections(
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    options::LightOptions;
    stacks_sorted::Bool=false,
)
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for sector in turtle.sectors
        projection = _prepared_direction_projection(prepared, sector.direction, options)
        if projection isa DenseDirectionProjectionResult
            _accumulate_scattering_counts!(
                edge_counts,
                sun_hits,
                sector,
                projection,
                prepared.virtual_node_mask,
                prepared.geometry.pavement_node_mask,
                scratch;
                node_ids=prepared.geometry.node_ids,
                stacks_sorted=stacks_sorted,
            )
        else
            _accumulate_scattering_counts!(
                edge_counts,
                sun_hits,
                sector,
                projection,
                prepared.virtual_nodes,
                prepared.geometry.node_group,
                scratch;
                node_ids=prepared.geometry.node_ids,
                stacks_sorted=stacks_sorted,
            )
        end
    end

    return _edge_counts_from_packed(edge_counts), sun_hits
end

function _pair_counts_from_raycore_projections(
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    options::LightOptions;
    stacks_sorted::Bool=false,
)
    prepared = data.prepared
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for sector in turtle.sectors
        if sector.source != :sun &&
           !_raycore_use_raster_compat_projection(data) &&
           Float32(sector.direction[3]) > 0.0f0
            traced = _raycore_trace_direction_stacks(data, sector.direction, options)
            overflow_pixel = findfirst(traced.overflow)
            overflow_pixel === nothing || error(
                "Raycore max_hits_per_pixel=$(traced.max_hits) exceeded for pixel $overflow_pixel. " *
                "Increase RaycoreBackendConfig(max_hits_per_pixel=...).",
            )
            dense_supported = _raycore_dense_edge_accumulation_supported(data)
            if data.edge_accumulation == :dense_atomic && !KernelAbstractions.supports_atomics(data.backend)
                error(
                    "Raycore edge_accumulation=:dense_atomic requires a KernelAbstractions backend " *
                    "with atomic support. Use :auto or :sparse_host_reduce for this backend.",
                )
            end
            if data.edge_accumulation == :dense_atomic && !_raycore_dense_edge_matrix_fits(data)
                n_nodes = length(data.prepared.geometry.node_ids)
                error(
                    "Raycore edge_accumulation=:dense_atomic requires a dense $n_nodes x $n_nodes edge matrix " *
                    "larger than dense_edge_limit_bytes=$(data.dense_edge_limit_bytes). " *
                    "Increase RaycoreBackendConfig(dense_edge_limit_bytes=...) or use :sparse_host_reduce.",
                )
            end

            if data.edge_accumulation == :dense_atomic || (data.edge_accumulation == :auto && dense_supported)
                dense_counts = _raycore_scattering_dense_counts_from_traced_stacks(data, traced)
                _merge_dense_edge_counts!(edge_counts, dense_counts, data.prepared.geometry.node_ids)
            else
                edge_keys = _raycore_scattering_edge_keys_from_traced_stacks(data, traced)
                _merge_packed_edge_keys!(edge_counts, edge_keys)
            end
            continue
        end

        projection = _raycore_direction_projection(data, sector.direction, options)
        _accumulate_scattering_counts!(
            edge_counts,
            sun_hits,
            sector,
            projection,
            prepared.virtual_node_mask,
            prepared.geometry.pavement_node_mask,
            scratch;
            node_ids=prepared.geometry.node_ids,
            stacks_sorted=stacks_sorted,
        )
    end

    return _edge_counts_from_packed(edge_counts), sun_hits
end

function _all_dir_hits_for_scattering(first::FirstOrderResult, sun_hits::Dict{Int,Int}, options::LightOptions, node_ids)
    all_hits = Dict{Int,Int}(nid => get(first.hits_per_node, nid, 0) for nid in node_ids)
    if !options.all_in_turtle
        for (nid, hsun) in sun_hits
            all_hits[nid] = max(0, get(all_hits, nid, 0) - hsun)
        end
    end
    all_hits
end

function _group_optical_coeffs(models::LightModels)
    coeffs = Dict{Tuple{String,String},Dict{String,Float64}}()

    for group_model in values(models)
        group = group_model.group
        isempty(group) && continue

        isempty(group_model.types) && continue

        for (type_name, type_model) in group_model.types
            interception = type_model.interception
            interception === nothing && continue

            coeff =
                if _is_sensor_interception(interception)
                    Dict{String,Float64}("PAR" => 0.0, "NIR" => 0.0)
                else
                    props = interception.optical_properties
                    props === nothing && continue
                    band_coeffs = Dict{String,Float64}()
                    has_par = get(props.extras, "__has_par", true)
                    has_nir = get(props.extras, "__has_nir", true)
                    has_par && (band_coeffs["PAR"] = props.par)
                    has_nir && (band_coeffs["NIR"] = props.nir)
                    band_coeffs
                end

            coeffs[(group, type_name)] = coeff
            # Group-level fallback for nodes without explicit type metadata.
            if !haskey(coeffs, (group, "*"))
                coeffs[(group, "*")] = coeff
            end
        end
    end

    coeffs
end

struct ScatteringSceneContext
    pair_counts::ScatteringPairCounts
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
end

function _build_scattering_topology_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    pair_counts::ScatteringPairCounts,
    sun_hits::Dict{Int,Int},
)
    group_type_coeffs = _group_optical_coeffs(models)
    node_ids = copy(prepared.geometry.node_ids)
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(prepared.geometry.node_group, nid, _scene_group(scene, nid, ""))
        default_type = g == "pavement" ? "Cobblestone" : ""
        node_type[nid] = _scene_type(scene, nid, default_type)
    end
    return ScatteringTopologyCache(
        pair_counts,
        sun_hits,
        node_ids,
        prepared.geometry.node_group,
        node_type,
        group_type_coeffs,
    )
end

function _build_scattering_topology_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    projections::AbstractVector{DirectionProjectionResult},
    stacks_sorted::Bool=false,
)
    pair_counts, sun_hits = _pair_counts_from_projections(
        turtle,
        projections,
        prepared.virtual_nodes,
        prepared.geometry.node_group,
        prepared.geometry.node_ids,
        prepared.geometry.pavement_node_mask,
        prepared.virtual_node_mask,
        stacks_sorted,
    )
    return _build_scattering_topology_cache(scene, models, prepared, pair_counts, sun_hits)
end

function _build_scattering_topology_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    options::LightOptions;
    stacks_sorted::Bool=false,
)
    pair_counts, sun_hits = _pair_counts_from_streamed_projections(
        prepared,
        turtle,
        options;
        stacks_sorted=stacks_sorted,
    )
    return _build_scattering_topology_cache(scene, models, prepared, pair_counts, sun_hits)
end

function _node_ids_for_scattering(topology::ScatteringTopologyCache, first::FirstOrderResult)
    node_set = Set{Int}()
    for nid in topology.node_ids
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.par)
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.nir)
        push!(node_set, nid)
    end
    return collect(node_set)
end

function _transfer_graph_from_topology(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
)
    node_ids = _node_ids_for_scattering(topology, first)
    all_hits = _all_dir_hits_for_scattering(first, topology.sun_hits, options, node_ids)
    coeff_par, coeff_nir =
        _coeff_maps_by_node(node_ids, topology.node_group, topology.node_type, topology.group_type_coeffs, options)
    return ScatteringTransferGraph(
        topology.pair_counts,
        all_hits,
        node_ids,
        topology.node_group,
        topology.node_type,
        topology.group_type_coeffs,
        coeff_par,
        coeff_nir,
        options.scattering_coeff_par,
        options.scattering_coeff_nir,
    )
end

function _scattering_context(scene::PlantGeom.SceneGeometry, models::LightModels, turtle::TurtleGrid, first::FirstOrderResult, options::LightOptions)
    pair_counts, sun_hits, geom_node_ids, node_group = _pair_counts_for_scattering(scene, models, turtle, options)
    group_type_coeffs = _group_optical_coeffs(models)

    node_set = Set{Int}()
    for nid in geom_node_ids
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.par)
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.nir)
        push!(node_set, nid)
    end
    node_ids = collect(node_set)
    all_hits = _all_dir_hits_for_scattering(first, sun_hits, options, node_ids)
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(node_group, nid, _scene_group(scene, nid, ""))
        default_type = g == "pavement" ? "Cobblestone" : ""
        node_type[nid] = _scene_type(scene, nid, default_type)
    end
    return ScatteringSceneContext(pair_counts, all_hits, node_ids, node_group, node_type, group_type_coeffs)
end

function _scattering_backend_from_mode(mode::Symbol)
    mode == :raycast && return RaycastScatteringBackend()
    mode == :links && return LinksScatteringBackend()
    error("Unsupported scattering mode: $mode (supported: :raycast, :links)")
end

function _resolve_scattering_backend(mode::Symbol, backend::Nothing)
    return _scattering_backend_from_mode(mode)
end

function _resolve_scattering_backend(mode::Symbol, backend::ScatteringBackend)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode (supported: :raycast, :links)")
    return backend
end

function _resolve_scattering_backend(mode::Symbol, backend)
    error(
        "Unsupported scattering backend selector type: $(typeof(backend)). " *
        "Use `nothing`, `RaycastScatteringBackend()`, `LinksScatteringBackend()`, or `RaycoreScatteringBackend()`.",
    )
end

"""
    build_scattering_transfer_graph(scene, models, turtle, first, options; mode=:raycast, backend=nothing)::ScatteringTransferGraph

Build the scattering transfer graph (pair links and per-node hit normalization data)
independently from iterative scattering propagation.
"""
function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend,
)
    pair_counts, sun_hits = _pair_counts_from_raycore_projections(data, turtle, options; stacks_sorted=false)
    topology = _build_scattering_topology_cache(scene, models, data.prepared, pair_counts, sun_hits)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return build_scattering_transfer_graph(scene, models, turtle, first, options, b)
end

function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    prepared = _prepare_interception_data(scene, models, options)
    topology = _build_scattering_topology_cache(scene, models, prepared, turtle, options)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    # CPU reference currently uses the same transfer-graph construction for both modes.
    return build_scattering_transfer_graph(scene, models, turtle, first, options, RaycastScatteringBackend())
end

function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend,
)
    data = _prepare_raycore_interception_data(
        scene,
        models,
        options,
        RaycoreInterceptionBackend(backend.config),
    )
    pair_counts, sun_hits = _pair_counts_from_raycore_projections(data, turtle, options; stacks_sorted=false)
    topology = _build_scattering_topology_cache(scene, models, data.prepared, pair_counts, sun_hits)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return build_scattering_transfer_graph(topology, first, options, b)
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    return build_scattering_transfer_graph(topology, first, options, RaycastScatteringBackend())
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycoreScatteringBackend,
)
    return _transfer_graph_from_topology(topology, first, options)
end

function _default_band_coeff(options::LightOptions, band_key::String)
    bk = uppercase(band_key)
    bk == "NIR" && return options.scattering_coeff_nir
    return options.scattering_coeff_par
end

@inline function _group_type_band_coeff(
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}},
    group::String,
    type_name::String,
    band::String,
    default_coeff::Float64,
)
    coeffs = get(
        group_type_coeffs,
        (group, type_name),
        get(group_type_coeffs, (group, "*"), Dict{String,Float64}()),
    )
    return get(coeffs, band, default_coeff)
end

function _coeff_by_node(
    node_ids::Vector{Int},
    node_group::Dict{Int,String},
    node_type::Dict{Int,String},
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}},
    band_key::String,
    default_coeff::Float64,
)
    coeff_by_node = Dict{Int,Float64}()
    band = uppercase(band_key)
    for nid in node_ids
        g = get(node_group, nid, "")
        t = get(node_type, nid, "")
        coeff_by_node[nid] = _group_type_band_coeff(group_type_coeffs, g, t, band, default_coeff)
    end
    coeff_by_node
end

function _coeff_by_node(
    graph::ScatteringTransferGraph,
    band_key::String,
    default_coeff::Float64,
)
    band = uppercase(band_key)
    if band == "PAR"
        default_coeff == graph.default_coeff_par && return graph.coeff_par_by_node
    elseif band == "NIR"
        default_coeff == graph.default_coeff_nir && return graph.coeff_nir_by_node
    end
    return _coeff_by_node(
        graph.node_ids,
        graph.node_group,
        graph.node_type,
        graph.group_type_coeffs,
        band,
        default_coeff,
    )
end

function _coeff_maps_by_node(
    node_ids::Vector{Int},
    node_group::Dict{Int,String},
    node_type::Dict{Int,String},
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}},
    options::LightOptions,
)
    coeff_par = _coeff_by_node(
        node_ids,
        node_group,
        node_type,
        group_type_coeffs,
        "PAR",
        options.scattering_coeff_par,
    )
    coeff_nir = _coeff_by_node(
        node_ids,
        node_group,
        node_type,
        group_type_coeffs,
        "NIR",
        options.scattering_coeff_nir,
    )
    return coeff_par, coeff_nir
end

function _initial_scattering_power(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}},
    band::AbstractString,
)
    if initial_power_per_node !== nothing
        return _copy_node_values(initial_power_per_node, graph.node_ids)
    end
    source = uppercase(String(band)) == "NIR" ? first.incident_power.nir : first.incident_power.par
    return _copy_node_values(source, graph.node_ids)
end

function _propagate_scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    node_ids = graph.node_ids
    pair_counts = graph.pair_counts
    all_hits = graph.all_hits
    current = _copy_node_values(initial_power_per_node, node_ids)
    added = _dict_zero(node_ids)
    ref = _sum_dict_values(current)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    iterations = 0

    converged = false
    for it in 1:options.scattering_max_iter
        iterations = it

        hit_energy = _dict_zero(node_ids)
        for nid in node_ids
            nh = get(all_hits, nid, 0)
            if nh > 0
                hit_energy[nid] = get(current, nid, 0.0) * get(coeff_by_node, nid, default_coeff) / nh / 2.0
            end
        end

        next = _dict_zero(node_ids)
        for ((to, from), cnt) in pair_counts
            next[to] = get(next, to, 0.0) + cnt * get(hit_energy, from, 0.0)
        end

        total_next = _sum_dict_values(next)

        for nid in node_ids
            added[nid] = get(added, nid, 0.0) + get(next, nid, 0.0)
        end

        current = next

        if total_next < thr
            converged = true
            break
        end
    end
    return added, iterations, converged
end

function _scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    options::LightOptions,
    band_key::String,
    default_coeff::Float64,
)
    coeffs = _coeff_by_node(graph, band_key, default_coeff)
    return _propagate_scattering_one_band(initial_power_per_node, graph, coeffs, options, default_coeff)
end

KernelAbstractions.@kernel function _scattering_hit_energy_kernel!(
    hit_energy,
    current,
    coeff,
    all_hits,
)
    i = @index(Global, Linear)
    @inbounds begin
        nh = all_hits[i]
        if nh > 0
            T = eltype(hit_energy)
            hit_energy[i] = current[i] * coeff[i] / T(nh) / T(2)
        else
            hit_energy[i] = zero(eltype(hit_energy))
        end
    end
end

KernelAbstractions.@kernel function _scattering_transfer_by_node_kernel!(
    next,
    to_idx,
    from_idx,
    counts,
    hit_energy,
    n_edges::Int,
)
    node_idx = @index(Global, Linear)
    T = eltype(next)
    total = zero(T)
    @inbounds for edge_idx in 1:n_edges
        if to_idx[edge_idx] == node_idx
            total += T(counts[edge_idx]) * hit_energy[from_idx[edge_idx]]
        end
    end
    @inbounds next[node_idx] = total
end

KernelAbstractions.@kernel function _scattering_zero_kernel!(values)
    i = @index(Global, Linear)
    @inbounds values[i] = zero(eltype(values))
end

KernelAbstractions.@kernel function _scattering_transfer_by_edge_atomic_kernel!(
    next,
    to_idx,
    from_idx,
    counts,
    hit_energy,
)
    edge_idx = @index(Global, Linear)
    @inbounds begin
        to = to_idx[edge_idx]
        from = from_idx[edge_idx]
        T = eltype(next)
        @atomic :monotonic next[to] += T(counts[edge_idx]) * hit_energy[from]
    end
end

KernelAbstractions.@kernel function _scattering_accumulate_added_kernel!(added, next)
    i = @index(Global, Linear)
    @inbounds added[i] += next[i]
end

function _dense_scattering_arrays(
    graph::ScatteringTransferGraph,
    initial_power_per_node::Dict{Int,Float64},
    coeff_by_node::Dict{Int,Float64},
    default_coeff::Float64,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    node_ids = graph.node_ids
    node_index = Dict{Int,Int}(nid => i for (i, nid) in pairs(node_ids))
    n_nodes = length(node_ids)
    initial = zeros(T, n_nodes)
    coeff = zeros(T, n_nodes)
    all_hits = zeros(Int, n_nodes)
    @inbounds for (i, nid) in pairs(node_ids)
        initial[i] = get(initial_power_per_node, nid, 0.0)
        coeff[i] = get(coeff_by_node, nid, default_coeff)
        all_hits[i] = get(graph.all_hits, nid, 0)
    end

    n_edges = length(graph.pair_counts)
    to_idx = Vector{Int}(undef, n_edges)
    from_idx = Vector{Int}(undef, n_edges)
    counts = Vector{Int}(undef, n_edges)
    @inbounds for edge_idx in 1:n_edges
        to_idx[edge_idx] = node_index[graph.pair_counts.to_nodes[edge_idx]]
        from_idx[edge_idx] = node_index[graph.pair_counts.from_nodes[edge_idx]]
        counts[edge_idx] = graph.pair_counts.counts[edge_idx]
    end

    return (
        initial=initial,
        coeff=coeff,
        all_hits=all_hits,
        to_idx=to_idx,
        from_idx=from_idx,
        counts=counts,
    )
end

function _dense_vector_to_node_dict(node_ids::Vector{Int}, values::AbstractVector{<:Real})
    out = Dict{Int,Float64}()
    @inbounds for i in eachindex(node_ids)
        out[node_ids[i]] = Float64(values[i])
    end
    return out
end

function _propagate_scattering_one_band_device(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
)
    node_ids = graph.node_ids
    n_nodes = length(node_ids)
    n_nodes == 0 && return Dict{Int,Float64}(), 0, true

    arrays = _dense_scattering_arrays(graph, initial_power_per_node, coeff_by_node, default_coeff, scattering_eltype)
    ref = sum(x -> Float64(x), arrays.initial)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    n_edges = length(arrays.counts)

    current_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    next_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    added_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    hit_energy_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    coeff_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    all_hits_dev = KernelAbstractions.allocate(backend, Int, n_nodes)
    to_idx_dev = KernelAbstractions.allocate(backend, Int, n_edges)
    from_idx_dev = KernelAbstractions.allocate(backend, Int, n_edges)
    counts_dev = KernelAbstractions.allocate(backend, Int, n_edges)

    KernelAbstractions.copyto!(backend, current_dev, arrays.initial)
    KernelAbstractions.copyto!(backend, next_dev, zeros(scattering_eltype, n_nodes))
    KernelAbstractions.copyto!(backend, added_dev, zeros(scattering_eltype, n_nodes))
    KernelAbstractions.copyto!(backend, hit_energy_dev, zeros(scattering_eltype, n_nodes))
    KernelAbstractions.copyto!(backend, coeff_dev, arrays.coeff)
    KernelAbstractions.copyto!(backend, all_hits_dev, arrays.all_hits)
    KernelAbstractions.copyto!(backend, to_idx_dev, arrays.to_idx)
    KernelAbstractions.copyto!(backend, from_idx_dev, arrays.from_idx)
    KernelAbstractions.copyto!(backend, counts_dev, arrays.counts)

    hit_kernel = _scattering_hit_energy_kernel!(backend, workgroupsize)
    zero_kernel = _scattering_zero_kernel!(backend, workgroupsize)
    transfer_by_node_kernel = _scattering_transfer_by_node_kernel!(backend, workgroupsize)
    transfer_by_edge_kernel = _scattering_transfer_by_edge_atomic_kernel!(backend, workgroupsize)
    add_kernel = _scattering_accumulate_added_kernel!(backend, workgroupsize)
    use_atomic_transfer = n_edges > 0 && KernelAbstractions.supports_atomics(backend)

    iterations = 0
    converged = false
    for it in 1:options.scattering_max_iter
        iterations = it
        hit_kernel(hit_energy_dev, current_dev, coeff_dev, all_hits_dev; ndrange=n_nodes)
        if use_atomic_transfer
            zero_kernel(next_dev; ndrange=n_nodes)
            transfer_by_edge_kernel(next_dev, to_idx_dev, from_idx_dev, counts_dev, hit_energy_dev; ndrange=n_edges)
        else
            transfer_by_node_kernel(next_dev, to_idx_dev, from_idx_dev, counts_dev, hit_energy_dev, n_edges; ndrange=n_nodes)
        end
        add_kernel(added_dev, next_dev; ndrange=n_nodes)
        KernelAbstractions.synchronize(backend)

        total_next = sum(x -> Float64(x), Array(next_dev))
        current_dev, next_dev = next_dev, current_dev
        if total_next < thr
            converged = true
            break
        end
    end

    added = _dense_vector_to_node_dict(node_ids, Array(added_dev))
    return added, iterations, converged
end

"""
    compute_scattering_band(graph, first, options; mode=:raycast, backend=nothing, band="PAR", initial_power_per_node=nothing, default_coeff=nothing)

Run one-band iterative scattering from a pre-built transfer graph.
"""
function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return compute_scattering_band(
        graph,
        first,
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    initial = _initial_scattering_power(graph, first, initial_power_per_node, band)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    coeffs = _coeff_by_node(graph, String(band), dflt)
    added, iterations, converged = _propagate_scattering_one_band_device(
        initial,
        graph,
        coeffs,
        options,
        dflt,
        backend.config.backend,
        backend.config.workgroupsize,
        backend.config.scattering_eltype,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged)
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    initial = _initial_scattering_power(graph, first, initial_power_per_node, band)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    added, iterations, converged = _scattering_one_band(
        initial,
        graph,
        options,
        String(band),
        dflt,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged)
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering_band(
        graph,
        first,
        options,
        RaycastScatteringBackend();
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    graph = build_scattering_transfer_graph(scene, models, data, turtle, first, options, backend)
    return compute_scattering_band(
        graph,
        first,
        options,
        backend;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, models, turtle, first, options, b)
    return compute_scattering_band(
        graph,
        first,
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(topology, first, options, b)
    return compute_scattering_band(
        graph,
        first,
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

"""
    compute_scattering(scene, models, turtle, first, options; mode=:raycast, backend=nothing)::ScatteringResult
    compute_scattering(graph, first, options; mode=:raycast, backend=nothing)::ScatteringResult

Compute iterative multiple scattering for PAR and NIR from first-order incident power.
When a `ScatteringTransferGraph` is provided, transfer-link construction is skipped and only
iterative propagation is run.
"""
function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return compute_scattering(graph, first, options, b)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    initial_par = _initial_scattering_power(graph, first, nothing, "PAR")
    initial_nir = _initial_scattering_power(graph, first, nothing, "NIR")

    added_par, it_par, conv_par = _propagate_scattering_one_band(
        initial_par,
        graph,
        graph.coeff_par_by_node,
        options,
        graph.default_coeff_par,
    )
    added_nir, it_nir, conv_nir = _propagate_scattering_one_band(
        initial_nir,
        graph,
        graph.coeff_nir_by_node,
        options,
        graph.default_coeff_nir,
    )

    ScatteringResult(SpectralNodeValues(added_par, added_nir), max(it_par, it_nir), conv_par && conv_nir)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering(graph, first, options, RaycastScatteringBackend())
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend,
)
    par = compute_scattering_band(graph, first, options, backend; band="PAR")
    nir = compute_scattering_band(graph, first, options, backend; band="NIR")

    return ScatteringResult(
        SpectralNodeValues(par.added_power_per_node, nir.added_power_per_node),
        max(par.iterations, nir.iterations),
        par.converged && nir.converged,
    )
end

function compute_scattering(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RaycoreScatteringBackend,
)
    graph = build_scattering_transfer_graph(scene, models, data, turtle, first, options, backend)
    return compute_scattering(graph, first, options, backend)
end

function compute_scattering(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, models, turtle, first, options, b)
    return compute_scattering(graph, first, options, b)
end

function compute_scattering(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(topology, first, options, b)
    return compute_scattering(graph, first, options, b)
end
