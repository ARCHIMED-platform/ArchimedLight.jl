function _dict_zero(node_ids)
    out = Dict{Int,Float64}()
    sizehint!(out, length(node_ids))
    for nid in node_ids
        out[nid] = 0.0
    end
    return out
end

struct ScatteringTopologyCache
    pair_counts::ScatteringPairCounts
    sun_hits::Dict{Int,Int}
    sun_hits_by_node::Vector{Int}
    node_ids::Vector{Int}
    pair_to_idx::Vector{Int}
    pair_from_idx::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
    dense_static::Base.RefValue{Union{Nothing,DenseScatteringStaticGraph}}
end

function _copy_node_values(source::Dict{Int,Float64}, node_ids)
    out = Dict{Int,Float64}()
    sizehint!(out, length(node_ids))
    for nid in node_ids
        out[nid] = get(source, nid, 0.0)
    end
    return out
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

function _compact_nonzero_edge_keys!(edge_keys::Vector{UInt64})
    write_idx = 1
    @inbounds for edge in edge_keys
        edge == 0 && continue
        edge_keys[write_idx] = edge
        write_idx += 1
    end
    resize!(edge_keys, write_idx - 1)
    return edge_keys
end

function _merge_sorted_packed_edge_keys!(
    edge_counts::Dict{UInt64,Int},
    edge_keys::AbstractVector{UInt64},
    n_edges::Int=length(edge_keys),
)
    n_edges == 0 && return edge_counts
    sorted_edges = n_edges == length(edge_keys) ? edge_keys : view(edge_keys, 1:n_edges)
    sort!(sorted_edges)
    current = sorted_edges[1]
    count = 1
    @inbounds for i in 2:n_edges
        edge = sorted_edges[i]
        if edge == current
            count += 1
        else
            edge_counts[current] = get(edge_counts, current, 0) + count
            current = edge
            count = 1
        end
    end
    edge_counts[current] = get(edge_counts, current, 0) + count
    return edge_counts
end

function _merge_packed_edge_keys!(edge_counts::Dict{UInt64,Int}, edge_keys::Vector{UInt64})
    _compact_nonzero_edge_keys!(edge_keys)
    return _merge_sorted_packed_edge_keys!(edge_counts, edge_keys)
end

function _merge_packed_edge_keys!(edge_counts::Dict{UInt64,Int}, edge_keys::AbstractVector{UInt64})
    nonzero_edges = UInt64[]
    sizehint!(nonzero_edges, length(edge_keys))
    @inbounds for edge in edge_keys
        edge == 0 && continue
        push!(nonzero_edges, edge)
    end
    return _merge_sorted_packed_edge_keys!(edge_counts, nonzero_edges)
end

function _merge_counted_packed_edge_keys!(
    edge_counts::Dict{UInt64,Int},
    edge_keys::AbstractVector{UInt64},
    edge_key_counts::AbstractVector{<:Integer},
    max_edges::Int,
    compact_edges::Vector{UInt64}=Vector{UInt64}(undef, 0),
)
    n_edges = 0
    @inbounds for count in edge_key_counts
        n_edges += Int(count)
    end
    n_edges == 0 && return edge_counts

    length(compact_edges) < n_edges && resize!(compact_edges, n_edges)
    write_idx = 1
    @inbounds for pixel_idx in eachindex(edge_key_counts)
        out_count = Int(edge_key_counts[pixel_idx])
        out_count == 0 && continue
        out_base = (pixel_idx - 1) * max_edges
        for slot in 1:out_count
            compact_edges[write_idx] = edge_keys[out_base+slot]
            write_idx += 1
        end
    end
    return _merge_sorted_packed_edge_keys!(edge_counts, compact_edges, n_edges)
end

@inline function _add_packed_edge_count!(edge_counts::Dict{UInt64,Int}, to::Int, from::Int)
    edge = _pack_scattering_edge(to, from)
    edge_counts[edge] = get(edge_counts, edge, 0) + 1
    return nothing
end

@inline _dense_edge_index(to_idx::Int, from_idx::Int, n_nodes::Int) = (to_idx - 1) * n_nodes + from_idx
@inline _dense_edge_to_idx(edge::Int, n_nodes::Int) = ((edge - 1) ÷ n_nodes) + 1
@inline _dense_edge_from_idx(edge::Int, n_nodes::Int) = ((edge - 1) % n_nodes) + 1

@inline function _add_dense_index_edge_count!(
    edge_counts::Dict{Int,Int},
    to_idx::Int,
    from_idx::Int,
    n_nodes::Int,
)
    edge = _dense_edge_index(to_idx, from_idx, n_nodes)
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

function _accumulate_sun_hits!(
    sun_hits_by_node::Vector{Int},
    projection::DenseDirectionProjectionResult,
)
    @inbounds for i in eachindex(projection.node_hits, sun_hits_by_node)
        h = projection.node_hits[i]
        h == 0 && continue
        sun_hits_by_node[i] += h
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

function _stack_transfer_pairs_no_virtual!(
    edge_counts::Dict{UInt64,Int},
    stack,
    node_group::Dict{Int,String},
)
    n_hits = length(stack)
    n_hits <= 1 && return nothing

    @inbounds begin
        current_node = _stack_hit_node(stack, 1)
        for h in 1:(n_hits - 1)
            next_node = _stack_hit_node(stack, h + 1)
            if _accept_scattering_link(node_group, current_node, next_node)
                _add_packed_edge_count!(edge_counts, current_node, next_node)
                _add_packed_edge_count!(edge_counts, next_node, current_node)
            end
            current_node = next_node
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

    no_virtual_nodes = isempty(virtual_nodes)
    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        if no_virtual_nodes
            _stack_transfer_pairs_no_virtual!(edge_counts, stack, node_group)
        else
            _stack_transfer_pairs!(edge_counts, stack, virtual_nodes, scratch, node_group)
        end
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

@inline function _accumulate_scattering_counts_dense_no_virtual!(
    edge_counts::Dict{UInt64,Int},
    stack,
    pavement_node_mask::Vector{Bool},
    node_ids::Vector{Int},
)
    n_hits = length(stack)
    n_hits <= 1 && return nothing

    @inbounds begin
        current_idx = _stack_hit_node(stack, 1)
        current_node = node_ids[current_idx]
        current_is_pavement = pavement_node_mask[current_idx]
        for h in 1:(n_hits - 1)
            next_idx = _stack_hit_node(stack, h + 1)
            next_node = node_ids[next_idx]
            next_is_pavement = pavement_node_mask[next_idx]
            if !(current_is_pavement && next_is_pavement)
                _add_packed_edge_count!(edge_counts, current_node, next_node)
                _add_packed_edge_count!(edge_counts, next_node, current_node)
            end
            current_idx = next_idx
            current_node = next_node
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

@inline function _accumulate_scattering_counts_dense_no_virtual!(
    edge_counts::Dict{UInt64,Int},
    stack::FlatPixelHitStack,
    pavement_node_mask::Vector{Bool},
    node_ids::Vector{Int},
)
    n_hits = length(stack)
    n_hits <= 1 && return nothing

    start = _flat_stack_start(stack)
    nodes = stack.parent.nodes
    @inbounds begin
        current_idx = Int(nodes[start])
        current_node = node_ids[current_idx]
        current_is_pavement = pavement_node_mask[current_idx]
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[start + h])
            next_node = node_ids[next_idx]
            next_is_pavement = pavement_node_mask[next_idx]
            if !(current_is_pavement && next_is_pavement)
                _add_packed_edge_count!(edge_counts, current_node, next_node)
                _add_packed_edge_count!(edge_counts, next_node, current_node)
            end
            current_idx = next_idx
            current_node = next_node
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

@inline function _accumulate_scattering_counts_dense_from_nodes_no_virtual!(
    edge_counts::Dict{UInt64,Int},
    nodes,
    stack_base::Int,
    n_hits::Int,
    pavement_node_mask::Vector{Bool},
    node_ids::Vector{Int},
)
    n_hits <= 1 && return nothing

    @inbounds begin
        current_idx = Int(nodes[stack_base+1])
        current_node = node_ids[current_idx]
        current_is_pavement = pavement_node_mask[current_idx]
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[stack_base+h+1])
            next_node = node_ids[next_idx]
            next_is_pavement = pavement_node_mask[next_idx]
            if !(current_is_pavement && next_is_pavement)
                _add_packed_edge_count!(edge_counts, current_node, next_node)
                _add_packed_edge_count!(edge_counts, next_node, current_node)
            end
            current_idx = next_idx
            current_node = next_node
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

@inline function _accumulate_scattering_counts_dense_from_nodes!(
    edge_counts::Dict{UInt64,Int},
    nodes,
    stack_base::Int,
    n_hits::Int,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch,
    node_ids::Vector{Int},
)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        node_idx = Int(nodes[stack_base+h])
        if !virtual_node_mask[node_idx]
            nearest_below = node_idx
        end
        below[h] = nearest_below
    end

    @inbounds begin
        current_idx = Int(nodes[stack_base+1])
        current_node = node_ids[current_idx]
        current_is_pavement = pavement_node_mask[current_idx]
        nearest_above = virtual_node_mask[current_idx] ? 0 : current_idx
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[stack_base+h+1])
            next_node = node_ids[next_idx]
            next_is_pavement = pavement_node_mask[next_idx]

            from_below_idx = below[h + 1]
            if from_below_idx != 0
                from_below_is_pavement = pavement_node_mask[from_below_idx]
                if !(current_is_pavement && from_below_is_pavement)
                    _add_packed_edge_count!(edge_counts, current_node, node_ids[from_below_idx])
                end
            end

            if nearest_above != 0
                from_above_is_pavement = pavement_node_mask[nearest_above]
                if !(next_is_pavement && from_above_is_pavement)
                    _add_packed_edge_count!(edge_counts, next_node, node_ids[nearest_above])
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

@inline function _accumulate_scattering_index_counts_from_nodes_no_virtual!(
    edge_counts::Dict{Int,Int},
    nodes,
    stack_base::Int,
    n_hits::Int,
    pavement_node_mask::Vector{Bool},
    n_nodes::Int,
)
    n_hits <= 1 && return nothing

    @inbounds begin
        current_idx = Int(nodes[stack_base+1])
        current_is_pavement = pavement_node_mask[current_idx]
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[stack_base+h+1])
            next_is_pavement = pavement_node_mask[next_idx]
            if !(current_is_pavement && next_is_pavement)
                _add_dense_index_edge_count!(edge_counts, current_idx, next_idx, n_nodes)
                _add_dense_index_edge_count!(edge_counts, next_idx, current_idx, n_nodes)
            end
            current_idx = next_idx
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

@inline function _accumulate_scattering_index_counts_from_nodes!(
    edge_counts::Dict{Int,Int},
    nodes,
    stack_base::Int,
    n_hits::Int,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch,
    n_nodes::Int,
)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        node_idx = Int(nodes[stack_base+h])
        if !virtual_node_mask[node_idx]
            nearest_below = node_idx
        end
        below[h] = nearest_below
    end

    @inbounds begin
        current_idx = Int(nodes[stack_base+1])
        current_is_pavement = pavement_node_mask[current_idx]
        nearest_above = virtual_node_mask[current_idx] ? 0 : current_idx
        for h in 1:(n_hits - 1)
            next_idx = Int(nodes[stack_base+h+1])
            next_is_pavement = pavement_node_mask[next_idx]

            from_below_idx = below[h+1]
            if from_below_idx != 0
                from_below_is_pavement = pavement_node_mask[from_below_idx]
                if !(current_is_pavement && from_below_is_pavement)
                    _add_dense_index_edge_count!(edge_counts, current_idx, from_below_idx, n_nodes)
                end
            end

            if nearest_above != 0
                from_above_is_pavement = pavement_node_mask[nearest_above]
                if !(next_is_pavement && from_above_is_pavement)
                    _add_dense_index_edge_count!(edge_counts, next_idx, nearest_above, n_nodes)
                end
            end

            if !virtual_node_mask[next_idx]
                nearest_above = next_idx
            end
            current_idx = next_idx
            current_is_pavement = next_is_pavement
        end
    end
    return nothing
end

@inline _profile_flag_at(flag::Bool, ::Int) = flag
@inline _profile_flag_at(flag, i::Int) = Bool(flag[i])

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
    no_virtual_nodes::Union{Nothing,Bool}=nothing,
)
    node_ids === nothing && error("DenseDirectionProjectionResult scattering accumulation requires node_ids")
    if sector.source == :sun
        _accumulate_sun_hits!(sun_hits, projection, node_ids)
        return nothing
    end

    no_virtual = no_virtual_nodes === nothing ? !any(virtual_node_mask) : no_virtual_nodes
    for stack in values(projection.pixel_hits)
        n_hits = length(stack)
        n_hits <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        if no_virtual
            _accumulate_scattering_counts_dense_no_virtual!(
                edge_counts,
                stack,
                pavement_node_mask,
                node_ids,
            )
        else
            _accumulate_scattering_counts_dense!(
                edge_counts,
                stack,
                virtual_node_mask,
                pavement_node_mask,
                scratch,
                node_ids,
            )
        end
    end
    return nothing
end

function _accumulate_scattering_counts!(
    edge_counts::Dict{UInt64,Int},
    sun_hits_by_node::Vector{Int},
    sector::TurtleSector,
    projection::DenseDirectionProjectionResult,
    virtual_node_mask::Vector{Bool},
    pavement_node_mask::Vector{Bool},
    scratch::ScatteringStackScratch;
    node_ids::Union{Nothing,Vector{Int}}=nothing,
    stacks_sorted::Bool=false,
    no_virtual_nodes::Union{Nothing,Bool}=nothing,
)
    node_ids === nothing && error("DenseDirectionProjectionResult scattering accumulation requires node_ids")
    if sector.source == :sun
        _accumulate_sun_hits!(sun_hits_by_node, projection)
        return nothing
    end

    no_virtual = no_virtual_nodes === nothing ? !any(virtual_node_mask) : no_virtual_nodes
    for stack in values(projection.pixel_hits)
        n_hits = length(stack)
        n_hits <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        if no_virtual
            _accumulate_scattering_counts_dense_no_virtual!(
                edge_counts,
                stack,
                pavement_node_mask,
                node_ids,
            )
        else
            _accumulate_scattering_counts_dense!(
                edge_counts,
                stack,
                virtual_node_mask,
                pavement_node_mask,
                scratch,
                node_ids,
            )
        end
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
    no_virtual_nodes = isempty(virtual_nodes)

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
                no_virtual_nodes=no_virtual_nodes,
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
    no_virtual_nodes = isempty(prepared.virtual_nodes)

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
                no_virtual_nodes=no_virtual_nodes,
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

function _all_dir_hits_for_scattering(first::FirstOrderResult, sun_hits::Dict{Int,Int}, options::LightOptions, node_ids)
    dense = first.dense
    all_hits =
        if dense !== nothing && isempty(first.hits_per_node)
            if dense.node_ids == node_ids
                Dict{Int,Int}(nid => dense.hits_per_node[i] for (i, nid) in pairs(node_ids))
            else
                dense_index = Dict{Int,Int}(nid => i for (i, nid) in pairs(dense.node_ids))
                Dict{Int,Int}(nid => get(dense.hits_per_node, get(dense_index, nid, 0), 0) for nid in node_ids)
            end
        else
            Dict{Int,Int}(nid => get(first.hits_per_node, nid, 0) for nid in node_ids)
        end
    if !options.all_in_turtle
        for (nid, hsun) in sun_hits
            all_hits[nid] = max(0, get(all_hits, nid, 0) - hsun)
        end
    end
    return all_hits
end

function _dense_all_dir_hits_for_scattering(
    first::FirstOrderResult,
    sun_hits_by_node::Vector{Int},
    options::LightOptions,
    node_ids::Vector{Int},
    topology_node_ids::Vector{Int},
)
    dense = first.dense
    all_hits =
        if dense !== nothing && isempty(first.hits_per_node) && dense.node_ids == node_ids
            copy(dense.hits_per_node)
        else
            out = zeros(Int, length(node_ids))
            if dense !== nothing && isempty(first.hits_per_node)
                dense_index = Dict{Int,Int}(nid => i for (i, nid) in pairs(dense.node_ids))
                @inbounds for (i, nid) in pairs(node_ids)
                    out[i] = get(dense.hits_per_node, get(dense_index, nid, 0), 0)
                end
            else
                @inbounds for (i, nid) in pairs(node_ids)
                    out[i] = get(first.hits_per_node, nid, 0)
                end
            end
            out
        end
    if !options.all_in_turtle
        if node_ids == topology_node_ids
            @inbounds for i in eachindex(all_hits, sun_hits_by_node)
                all_hits[i] = max(0, all_hits[i] - sun_hits_by_node[i])
            end
        else
            topology_index = Dict{Int,Int}(nid => i for (i, nid) in pairs(topology_node_ids))
            @inbounds for (i, nid) in pairs(node_ids)
                sun_idx = get(topology_index, nid, 0)
                sun_idx == 0 && continue
                all_hits[i] = max(0, all_hits[i] - sun_hits_by_node[sun_idx])
            end
        end
    end
    return all_hits
end

function _dense_int_vector_to_node_dict(node_ids::Vector{Int}, values::AbstractVector{<:Integer})
    out = Dict{Int,Int}()
    sizehint!(out, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[node_ids[i]] = Int(values[i])
    end
    return out
end

function _dense_nonzero_int_vector_to_node_dict(node_ids::Vector{Int}, values::AbstractVector{<:Integer})
    out = Dict{Int,Int}()
    @inbounds for i in eachindex(node_ids, values)
        value = Int(values[i])
        value == 0 && continue
        out[node_ids[i]] = value
    end
    return out
end

function _merge_sun_hits_into_dense!(
    sun_hits_by_node::Vector{Int},
    sun_hits::Dict{Int,Int},
    node_index::Dict{Int,Int},
)
    isempty(sun_hits) && return sun_hits_by_node
    for (nid, h) in sun_hits
        idx = get(node_index, nid, 0)
        idx == 0 && continue
        sun_hits_by_node[idx] += h
    end
    return sun_hits_by_node
end

function _merge_sun_hits_into_dense!(
    sun_hits_by_node::Vector{Int},
    sun_hits::Dict{Int,Int},
    node_ids::Vector{Int},
)
    return _merge_sun_hits_into_dense!(
        sun_hits_by_node,
        sun_hits,
        Dict{Int,Int}(nid => i for (i, nid) in pairs(node_ids)),
    )
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

function _topology_dense_index_arrays(
    pair_counts::ScatteringPairCounts,
    sun_hits::Dict{Int,Int},
    node_ids::Vector{Int},
    node_index::Dict{Int,Int},
)
    n_edges = length(pair_counts)
    pair_to_idx = Vector{Int}(undef, n_edges)
    pair_from_idx = Vector{Int}(undef, n_edges)
    @inbounds for edge_idx in 1:n_edges
        pair_to_idx[edge_idx] = node_index[pair_counts.to_nodes[edge_idx]]
        pair_from_idx[edge_idx] = node_index[pair_counts.from_nodes[edge_idx]]
    end
    sun_hits_by_node = zeros(Int, length(node_ids))
    for (nid, h) in sun_hits
        idx = get(node_index, nid, 0)
        idx == 0 && continue
        sun_hits_by_node[idx] = h
    end
    return pair_to_idx, pair_from_idx, sun_hits_by_node
end

function _topology_dense_index_arrays(
    pair_counts::ScatteringPairCounts,
    sun_hits::Dict{Int,Int},
    node_ids::Vector{Int},
)
    return _topology_dense_index_arrays(
        pair_counts,
        sun_hits,
        node_ids,
        Dict{Int,Int}(nid => i for (i, nid) in pairs(node_ids)),
    )
end

function _merge_packed_edge_counts_as_indexed!(
    indexed_edge_counts::Dict{Int,Int},
    packed_edge_counts::Dict{UInt64,Int},
    node_ids::Vector{Int},
    node_index::Dict{Int,Int},
)
    isempty(packed_edge_counts) && return indexed_edge_counts
    n_nodes = length(node_ids)
    for (edge, count) in packed_edge_counts
        to_idx = get(node_index, _unpack_scattering_to(edge), 0)
        from_idx = get(node_index, _unpack_scattering_from(edge), 0)
        (to_idx == 0 || from_idx == 0 || count <= 0) && continue
        dense_edge = _dense_edge_index(to_idx, from_idx, n_nodes)
        indexed_edge_counts[dense_edge] = get(indexed_edge_counts, dense_edge, 0) + count
    end
    return indexed_edge_counts
end

function _merge_packed_edge_counts_as_indexed!(
    indexed_edge_counts::Dict{Int,Int},
    packed_edge_counts::Dict{UInt64,Int},
    node_ids::Vector{Int},
)
    return _merge_packed_edge_counts_as_indexed!(
        indexed_edge_counts,
        packed_edge_counts,
        node_ids,
        Dict{Int,Int}(nid => i for (i, nid) in pairs(node_ids)),
    )
end

function _pair_counts_from_indexed_edges(
    indexed_edge_counts::Dict{Int,Int},
    node_ids::Vector{Int},
)
    isempty(indexed_edge_counts) && return ScatteringPairCounts(Int[], Int[], Int[]), Int[], Int[]
    n_nodes = length(node_ids)
    dense_edges = sort!(collect(keys(indexed_edge_counts)))
    to_nodes = Int[]
    from_nodes = Int[]
    to_idx = Int[]
    from_idx = Int[]
    counts = Int[]
    sizehint!(to_nodes, length(dense_edges))
    sizehint!(from_nodes, length(dense_edges))
    sizehint!(to_idx, length(dense_edges))
    sizehint!(from_idx, length(dense_edges))
    sizehint!(counts, length(dense_edges))

    @inbounds for edge in dense_edges
        count = indexed_edge_counts[edge]
        count > 0 || continue
        t_idx = _dense_edge_to_idx(edge, n_nodes)
        f_idx = _dense_edge_from_idx(edge, n_nodes)
        push!(to_idx, t_idx)
        push!(from_idx, f_idx)
        push!(to_nodes, node_ids[t_idx])
        push!(from_nodes, node_ids[f_idx])
        push!(counts, count)
    end
    return ScatteringPairCounts(to_nodes, from_nodes, counts), to_idx, from_idx
end

function _build_scattering_topology_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    pair_counts::ScatteringPairCounts,
    sun_hits::Dict{Int,Int},
)
    return _build_scattering_topology_cache(
        scene,
        models,
        prepared,
        pair_counts,
        zeros(Int, length(prepared.geometry.node_ids)),
        sun_hits,
    )
end

function _build_scattering_topology_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    pair_counts::ScatteringPairCounts,
    sun_hits_by_node::Vector{Int},
    extra_sun_hits::Dict{Int,Int},
)
    group_type_coeffs = _group_optical_coeffs(models)
    node_ids = copy(prepared.geometry.node_ids)
    length(sun_hits_by_node) == length(node_ids) ||
        throw(ArgumentError("dense sun-hit vector length $(length(sun_hits_by_node)) does not match topology node count $(length(node_ids))"))
    node_index = prepared.geometry.node_index
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(prepared.geometry.node_group, nid, "")
        default_type = g == "pavement" ? "Cobblestone" : ""
        type_name = get(prepared.geometry.node_type, nid, default_type)
        node_type[nid] = isempty(type_name) ? default_type : type_name
    end
    pair_to_idx, pair_from_idx, _unused_sun_hits =
        _topology_dense_index_arrays(pair_counts, Dict{Int,Int}(), node_ids, node_index)
    dense_sun_hits = copy(sun_hits_by_node)
    _merge_sun_hits_into_dense!(dense_sun_hits, extra_sun_hits, node_index)
    sun_hits = _dense_nonzero_int_vector_to_node_dict(node_ids, dense_sun_hits)
    return ScatteringTopologyCache(
        pair_counts,
        sun_hits,
        dense_sun_hits,
        node_ids,
        pair_to_idx,
        pair_from_idx,
        prepared.geometry.node_group,
        node_type,
        group_type_coeffs,
        Ref{Union{Nothing,DenseScatteringStaticGraph}}(nothing),
    )
end

function _build_scattering_topology_cache_from_indexed_edges(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    indexed_edge_counts::Dict{Int,Int},
    packed_edge_counts::Dict{UInt64,Int},
    sun_hits::Dict{Int,Int},
)
    return _build_scattering_topology_cache_from_indexed_edges(
        scene,
        models,
        prepared,
        indexed_edge_counts,
        packed_edge_counts,
        zeros(Int, length(prepared.geometry.node_ids)),
        sun_hits,
    )
end

function _build_scattering_topology_cache_from_indexed_edges(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    indexed_edge_counts::Dict{Int,Int},
    packed_edge_counts::Dict{UInt64,Int},
    sun_hits_by_node::Vector{Int},
    extra_sun_hits::Dict{Int,Int},
)
    group_type_coeffs = _group_optical_coeffs(models)
    node_ids = copy(prepared.geometry.node_ids)
    length(sun_hits_by_node) == length(node_ids) ||
        throw(ArgumentError("dense sun-hit vector length $(length(sun_hits_by_node)) does not match topology node count $(length(node_ids))"))
    node_index = prepared.geometry.node_index
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(prepared.geometry.node_group, nid, "")
        default_type = g == "pavement" ? "Cobblestone" : ""
        type_name = get(prepared.geometry.node_type, nid, default_type)
        node_type[nid] = isempty(type_name) ? default_type : type_name
    end
    _merge_packed_edge_counts_as_indexed!(indexed_edge_counts, packed_edge_counts, node_ids, node_index)
    pair_counts, pair_to_idx, pair_from_idx = _pair_counts_from_indexed_edges(indexed_edge_counts, node_ids)
    dense_sun_hits = copy(sun_hits_by_node)
    _merge_sun_hits_into_dense!(dense_sun_hits, extra_sun_hits, node_index)
    sun_hits = _dense_nonzero_int_vector_to_node_dict(node_ids, dense_sun_hits)
    return ScatteringTopologyCache(
        pair_counts,
        sun_hits,
        dense_sun_hits,
        node_ids,
        pair_to_idx,
        pair_from_idx,
        prepared.geometry.node_group,
        node_type,
        group_type_coeffs,
        Ref{Union{Nothing,DenseScatteringStaticGraph}}(nothing),
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

function _dense_scattering_static_graph!(
    topology::ScatteringTopologyCache,
    coeff_par::Vector{Float64},
    coeff_nir::Vector{Float64},
)
    dense_static = topology.dense_static[]
    if dense_static === nothing
        dense_static = DenseScatteringStaticGraph(
            topology.pair_to_idx,
            topology.pair_from_idx,
            copy(topology.pair_counts.counts),
            coeff_par,
            coeff_nir,
            IdDict{Any,Any}(),
        )
        topology.dense_static[] = dense_static
    end
    return dense_static
end

function _dense_scattering_static_graph!(
    topology::ScatteringTopologyCache,
    coeff_par_by_node::Dict{Int,Float64},
    coeff_nir_by_node::Dict{Int,Float64},
)
    dense_static = topology.dense_static[]
    if dense_static === nothing
        dense_static = DenseScatteringStaticGraph(
            topology.pair_to_idx,
            topology.pair_from_idx,
            topology.pair_counts.counts,
            topology.node_ids,
            coeff_par_by_node,
            coeff_nir_by_node,
        )
        topology.dense_static[] = dense_static
    end
    return dense_static
end

function _node_ids_for_scattering(topology::ScatteringTopologyCache, first::FirstOrderResult)
    node_ids = copy(topology.node_ids)
    dense = first.dense
    if dense !== nothing &&
       isempty(first.incident_power.par) &&
       isempty(first.incident_power.nir) &&
       dense.node_ids == topology.node_ids
        return node_ids
    end
    node_set = Set{Int}(node_ids)
    dense_ids = dense === nothing ? () : dense.node_ids
    for nid in Iterators.flatten((keys(first.incident_power.par), keys(first.incident_power.nir), dense_ids))
        if !(nid in node_set)
            push!(node_ids, nid)
            push!(node_set, nid)
        end
    end
    return node_ids
end

function _transfer_graph_from_topology(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
)
    node_ids = _node_ids_for_scattering(topology, first)
    dense_all_hits = _dense_all_dir_hits_for_scattering(
        first,
        topology.sun_hits_by_node,
        options,
        node_ids,
        topology.node_ids,
    )
    all_hits = _dense_int_vector_to_node_dict(node_ids, dense_all_hits)
    coeff_par, coeff_nir, dense_static =
        if node_ids == topology.node_ids
            cp, cn, cp_dense, cn_dense =
                _coeff_maps_and_vectors_by_node(
                    node_ids,
                    topology.node_group,
                    topology.node_type,
                    topology.group_type_coeffs,
                    options,
                )
            cp, cn, _dense_scattering_static_graph!(topology, cp_dense, cn_dense)
        else
            cp, cn =
                _coeff_maps_by_node(
                    node_ids,
                    topology.node_group,
                    topology.node_type,
                    topology.group_type_coeffs,
                    options,
                )
            cp, cn, nothing
        end
    graph = ScatteringTransferGraph(
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
        dense_static,
    )
    if dense_static !== nothing
        graph.dense[] = DenseScatteringGraph(dense_all_hits, dense_static)
    end
    return graph
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
        "Use `nothing`, `RaycastScatteringBackend()`, `LinksScatteringBackend()`, or `RasterGPUScatteringBackend()`.",
    )
end

function _rastergpu_dense_edge_counts_fits(data::RasterGPUSceneData)
    return data.dense_edge_counts_dev !== nothing && data.dense_edge_counts_host !== nothing
end

function _rastergpu_auto_dense_edge_accumulation_supported(data::RasterGPUSceneData)
    return _rastergpu_dense_edge_counts_fits(data) && KernelAbstractions.supports_atomics(data.backend)
end

KernelAbstractions.@kernel function _rastergpu_scattering_edge_keys_kernel!(
    edge_keys,
    edge_key_counts,
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
        out_slot = 0
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
                    out_slot += 1
                    to_node = node_ids[current_idx]
                    from_node = node_ids[from_below_idx]
                    edge_keys[out_base+out_slot] = _pack_scattering_edge(to_node, from_node)
                end

                next_idx = Int(nodes[stack_base+h+1])
                if nearest_above_idx != 0 &&
                   !(pavement_node_mask[next_idx] && pavement_node_mask[nearest_above_idx])
                    out_slot += 1
                    to_node = node_ids[next_idx]
                    from_node = node_ids[nearest_above_idx]
                    edge_keys[out_base+out_slot] = _pack_scattering_edge(to_node, from_node)
                end
            end
        end
        edge_key_counts[pixel_idx] = Int32(out_slot)
    end
end

KernelAbstractions.@kernel function _rastergpu_scattering_dense_counts_kernel!(
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

KernelAbstractions.@kernel function _rastergpu_clear_dense_counts_kernel!(dense_counts)
    pair_idx = @index(Global, Linear)
    @inbounds dense_counts[pair_idx] = 0
end

function _rastergpu_scattering_edge_keys_from_device_stacks(data::RasterGPUSceneData)
    n_pixels = length(data.overflow_host)
    max_hits = data.max_hits_per_pixel
    max_edges = 2 * (max_hits - 1)
    (n_pixels == 0 || max_edges <= 0) && return (keys=UInt64[], counts=Int32[], max_edges=max_edges)

    edge_keys_dev = data.edge_keys_dev
    edge_key_counts_dev = data.edge_key_counts_dev
    if edge_keys_dev === nothing || edge_key_counts_dev === nothing
        error(
            "RasterGPU sparse edge-key scratch was not allocated. " *
            "Use edge_accumulation=:auto or :sparse_host_reduce.",
        )
    end

    kernel = _rastergpu_scattering_edge_keys_kernel!(data.backend, data.workgroupsize)
    kernel(
        edge_keys_dev,
        edge_key_counts_dev,
        data.counts_dev,
        data.nodes_dev,
        data.virtual_node_mask_dev,
        data.pavement_node_mask_dev,
        data.node_ids_dev,
        max_hits;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(data.backend)

    edge_key_len = n_pixels * max_edges
    length(data.edge_keys_host) < edge_key_len && resize!(data.edge_keys_host, edge_key_len)
    length(data.edge_key_counts_host) < n_pixels && resize!(data.edge_key_counts_host, n_pixels)
    copyto!(data.edge_keys_host, edge_keys_dev)
    copyto!(data.edge_key_counts_host, edge_key_counts_dev)

    return (
        keys=data.edge_keys_host,
        counts=data.edge_key_counts_host,
        max_edges=max_edges,
    )
end

function _rastergpu_accumulate_device_scattering_edges!(
    edge_counts::Dict{UInt64,Int},
    data::RasterGPUSceneData,
)
    dense_supported = _rastergpu_auto_dense_edge_accumulation_supported(data)
    if data.edge_accumulation == :dense_atomic && !KernelAbstractions.supports_atomics(data.backend)
        error(
            "RasterGPU edge_accumulation=:dense_atomic requires a KernelAbstractions backend " *
            "with atomic support. Use :auto or :sparse_host_reduce for this backend.",
        )
    end
    if data.edge_accumulation == :dense_atomic && !_rastergpu_dense_edge_counts_fits(data)
        n_nodes = length(data.prepared.geometry.node_ids)
        error(
            "RasterGPU edge_accumulation=:dense_atomic requires a dense $n_nodes x $n_nodes edge matrix " *
            "larger than dense_edge_limit_bytes=$(data.dense_edge_limit_bytes). " *
            "Increase RasterGPUBackendConfig(dense_edge_limit_bytes=...) or use :sparse_host_reduce.",
        )
    end

    if data.edge_accumulation != :dense_atomic && !(data.edge_accumulation == :auto && dense_supported)
        edge_keys = _rastergpu_scattering_edge_keys_from_device_stacks(data)
        _merge_counted_packed_edge_keys!(
            edge_counts,
            edge_keys.keys,
            edge_keys.counts,
            edge_keys.max_edges,
            data.edge_compact_host,
        )
        return edge_counts
    end

    n_nodes = length(data.prepared.geometry.node_ids)
    n_pairs = n_nodes * n_nodes
    n_pairs == 0 && return edge_counts
    if _rastergpu_use_fused_dense_edges(data)
        copyto!(data.dense_edge_counts_host, data.dense_edge_counts_dev)
        return _merge_dense_edge_counts!(
            edge_counts,
            data.dense_edge_counts_host,
            data.prepared.geometry.node_ids,
        )
    end

    clear_kernel = _rastergpu_clear_dense_counts_kernel!(data.backend, data.workgroupsize)
    clear_kernel(data.dense_edge_counts_dev; ndrange=n_pairs)
    KernelAbstractions.synchronize(data.backend)

    n_pixels = data.prepared.geometry.plotbox.nx * data.prepared.geometry.plotbox.ny
    if n_pixels > 0
        edge_kernel = _rastergpu_scattering_dense_counts_kernel!(data.backend, data.workgroupsize)
        edge_kernel(
            data.dense_edge_counts_dev,
            data.counts_dev,
            data.nodes_dev,
            data.virtual_node_mask_dev,
            data.pavement_node_mask_dev,
            data.max_hits_per_pixel,
            n_nodes;
            ndrange=n_pixels,
        )
        KernelAbstractions.synchronize(data.backend)
    end

    copyto!(data.dense_edge_counts_host, data.dense_edge_counts_dev)
    return _merge_dense_edge_counts!(
        edge_counts,
        data.dense_edge_counts_host,
        data.prepared.geometry.node_ids,
    )
end

function _pair_counts_from_rastergpu_projections(
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    options::LightOptions,
)
    prepared = data.prepared
    geometry = prepared.geometry
    edge_counts = Dict{UInt64,Int}()
    sun_hits_by_node = zeros(Int, length(geometry.node_ids))

    for sector in turtle.sectors
        if sector.source == :sun
            _rastergpu_project_direction!(data, sector.direction, options)
            copyto!(data.node_counts_host, data.node_counts_dev)
            @inbounds for idx in eachindex(sun_hits_by_node)
                sun_hits_by_node[idx] += Int(data.node_counts_host[idx])
            end
        else
            _rastergpu_project_direction!(data, sector.direction, options; accumulate_dense_edges=true)
            _rastergpu_accumulate_device_scattering_edges!(edge_counts, data)
        end
    end

    return _edge_counts_from_packed(edge_counts), sun_hits_by_node
end

"""
    build_scattering_transfer_graph(scene, models, turtle, first, options; mode=:raycast, backend=nothing)::ScatteringTransferGraph

Build the scattering transfer graph (pair links and per-node hit normalization data)
independently from iterative scattering propagation.
"""
function build_scattering_transfer_graph(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RasterGPUScatteringBackend,
)
    pair_counts, sun_hits_by_node = _pair_counts_from_rastergpu_projections(data, turtle, options)
    topology = _build_scattering_topology_cache(
        scene,
        models,
        data.prepared,
        pair_counts,
        sun_hits_by_node,
        Dict{Int,Int}(),
    )
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
    backend::RasterGPUScatteringBackend,
)
    prepared = _prepare_interception_data(scene, models, options)
    data = _rastergpu_scene_data(prepared, backend.config)
    return build_scattering_transfer_graph(scene, models, data, turtle, first, options, backend)
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
    ::RasterGPUScatteringBackend,
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
    key = (group, type_name)
    coeffs =
        if haskey(group_type_coeffs, key)
            group_type_coeffs[key]
        else
            fallback_key = (group, "*")
            haskey(group_type_coeffs, fallback_key) || return default_coeff
            group_type_coeffs[fallback_key]
        end
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

function _coeff_maps_and_vectors_by_node(
    node_ids::Vector{Int},
    node_group::Dict{Int,String},
    node_type::Dict{Int,String},
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}},
    options::LightOptions,
)
    coeff_par = Dict{Int,Float64}()
    coeff_nir = Dict{Int,Float64}()
    sizehint!(coeff_par, length(node_ids))
    sizehint!(coeff_nir, length(node_ids))
    coeff_par_values = Vector{Float64}(undef, length(node_ids))
    coeff_nir_values = Vector{Float64}(undef, length(node_ids))
    @inbounds for (i, nid) in pairs(node_ids)
        g = get(node_group, nid, "")
        t = get(node_type, nid, "")
        par = _group_type_band_coeff(group_type_coeffs, g, t, "PAR", options.scattering_coeff_par)
        nir = _group_type_band_coeff(group_type_coeffs, g, t, "NIR", options.scattering_coeff_nir)
        coeff_par[nid] = par
        coeff_nir[nid] = nir
        coeff_par_values[i] = par
        coeff_nir_values[i] = nir
    end
    return coeff_par, coeff_nir, coeff_par_values, coeff_nir_values
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
    dense = first.dense
    if dense !== nothing && dense.node_ids == graph.node_ids
        source = uppercase(String(band)) == "NIR" ? dense.incident_power.nir : dense.incident_power.par
        return _dense_vector_to_node_dict(graph.node_ids, source)
    end
    source = uppercase(String(band)) == "NIR" ? first.incident_power.nir : first.incident_power.par
    return _copy_node_values(source, graph.node_ids)
end

function _dense_initial_from_node_dict(
    node_ids::Vector{Int},
    values::Dict{Int,Float64},
    ::Type{T},
) where {T<:AbstractFloat}
    out = zeros(T, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[i] = T(get(values, node_ids[i], 0.0))
    end
    return out
end

function _dense_initial_scattering_power(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}},
    band::AbstractString,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    if initial_power_per_node !== nothing
        return _dense_initial_from_node_dict(graph.node_ids, initial_power_per_node, T)
    end
    dense = first.dense
    if dense !== nothing && dense.node_ids == graph.node_ids
        source = uppercase(String(band)) == "NIR" ? dense.incident_power.nir : dense.incident_power.par
        return T === Float64 ? source : T.(source)
    end
    source = uppercase(String(band)) == "NIR" ? first.incident_power.nir : first.incident_power.par
    return _dense_initial_from_node_dict(graph.node_ids, source, T)
end

function _propagate_scattering_one_band_dense_only(
    initial_power_per_node::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    n_nodes = length(graph.node_ids)
    n_nodes == 0 && return Float64[], 0, true

    dense = _dense_scattering_graph(graph)
    static = dense.static
    initial, coeff = _dense_scattering_band_arrays(graph, initial_power_per_node, coeff_by_node, default_coeff)
    current = copy(initial)
    next = zeros(Float64, n_nodes)
    added = zeros(Float64, n_nodes)
    hit_energy = zeros(Float64, n_nodes)
    ref = sum(initial)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    iterations = 0

    converged = false
    for it in 1:options.scattering_max_iter
        iterations = it

        @inbounds for i in 1:n_nodes
            nh = dense.all_hits[i]
            if nh > 0
                hit_energy[i] = current[i] * coeff[i] / nh / 2.0
            else
                hit_energy[i] = 0.0
            end
        end

        fill!(next, 0.0)
        @inbounds for edge_idx in eachindex(static.counts)
            next[static.to_idx[edge_idx]] += static.counts[edge_idx] * hit_energy[static.from_idx[edge_idx]]
        end

        total_next = sum(next)

        @inbounds for i in 1:n_nodes
            added[i] += next[i]
        end

        current, next = next, current

        if total_next < thr
            converged = true
            break
        end
    end
    return added, iterations, converged
end

function _propagate_scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    added, iterations, converged, _ = _propagate_scattering_one_band_dense(
        initial_power_per_node,
        graph,
        coeff_by_node,
        options,
        default_coeff,
    )
    return added, iterations, converged
end

function _propagate_scattering_one_band_dense(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    node_ids = graph.node_ids
    if isempty(node_ids)
        return Dict{Int,Float64}(), 0, true, Float64[]
    end
    added, iterations, converged = _propagate_scattering_one_band_dense_only(
        initial_power_per_node,
        graph,
        coeff_by_node,
        options,
        default_coeff,
    )
    return _dense_vector_to_node_dict(node_ids, added), iterations, converged, added
end

function _dense_coeff_vector(
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    default_coeff::Float64,
)
    dense = _dense_scattering_graph(graph)
    if coeff_by_node === graph.coeff_par_by_node && default_coeff == graph.default_coeff_par
        return dense.static.coeff_par
    elseif coeff_by_node === graph.coeff_nir_by_node && default_coeff == graph.default_coeff_nir
        return dense.static.coeff_nir
    end
    coeff = zeros(Float64, length(graph.node_ids))
    @inbounds for (i, nid) in pairs(graph.node_ids)
        coeff[i] = get(coeff_by_node, nid, default_coeff)
    end
    return coeff
end

function _propagate_scattering_band_matrix_dense(
    initial_power_by_rhs::AbstractMatrix{Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    n_rhs, n_nodes = size(initial_power_by_rhs)
    n_nodes == length(graph.node_ids) ||
        throw(ArgumentError("dense initial scattering matrix has $n_nodes nodes, expected $(length(graph.node_ids))"))
    n_rhs == 0 && return zeros(Float64, 0, n_nodes), Int[], Bool[]
    n_nodes == 0 && return zeros(Float64, n_rhs, 0), zeros(Int, n_rhs), trues(n_rhs)

    dense = _dense_scattering_graph(graph)
    static = dense.static
    coeff = _dense_coeff_vector(graph, coeff_by_node, default_coeff)
    current = copy(initial_power_by_rhs)
    next = zeros(Float64, n_rhs, n_nodes)
    added = zeros(Float64, n_rhs, n_nodes)
    hit_energy = zeros(Float64, n_rhs, n_nodes)
    thresholds = zeros(Float64, n_rhs)
    active = trues(n_rhs)
    iterations = zeros(Int, n_rhs)
    converged = falses(n_rhs)

    @inbounds for r in 1:n_rhs
        ref = 0.0
        for i in 1:n_nodes
            ref += current[r, i]
        end
        thresholds[r] = options.scattering_stop_ratio * max(ref, eps(Float64))
    end

    @inbounds for it in 1:options.scattering_max_iter
        any_active = false
        for r in 1:n_rhs
            if active[r]
                iterations[r] = it
                any_active = true
            end
        end
        any_active || break

        for i in 1:n_nodes
            nh = dense.all_hits[i]
            if nh > 0
                factor = coeff[i] / nh / 2.0
                for r in 1:n_rhs
                    active[r] && (hit_energy[r, i] = current[r, i] * factor)
                end
            else
                for r in 1:n_rhs
                    active[r] && (hit_energy[r, i] = 0.0)
                end
            end
        end

        fill!(next, 0.0)
        for edge_idx in eachindex(static.counts)
            to = static.to_idx[edge_idx]
            from = static.from_idx[edge_idx]
            count = static.counts[edge_idx]
            for r in 1:n_rhs
                active[r] && (next[r, to] += count * hit_energy[r, from])
            end
        end

        for r in 1:n_rhs
            active[r] || continue
            total_next = 0.0
            for i in 1:n_nodes
                v = next[r, i]
                total_next += v
                added[r, i] += v
            end
            if total_next < thresholds[r]
                converged[r] = true
                active[r] = false
            end
        end

        current, next = next, current
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

function _dense_scattering_graph(graph::ScatteringTransferGraph)
    dense = graph.dense[]
    if dense === nothing
        dense =
            graph.dense_static === nothing ?
            DenseScatteringGraph(
                graph.pair_counts,
                graph.all_hits,
                graph.node_ids,
                graph.coeff_par_by_node,
                graph.coeff_nir_by_node,
            ) :
            DenseScatteringGraph(graph.all_hits, graph.node_ids, graph.dense_static)
        graph.dense[] = dense
    end
    return dense
end

function _dense_scattering_band_arrays(
    graph::ScatteringTransferGraph,
    initial_power_per_node::Dict{Int,Float64},
    coeff_by_node::Dict{Int,Float64},
    default_coeff::Float64,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    node_ids = graph.node_ids
    n_nodes = length(node_ids)
    dense = _dense_scattering_graph(graph)
    initial = _dense_initial_from_node_dict(node_ids, initial_power_per_node, T)
    cached_coeffs =
        if coeff_by_node === graph.coeff_par_by_node && default_coeff == graph.default_coeff_par
            dense.static.coeff_par
        elseif coeff_by_node === graph.coeff_nir_by_node && default_coeff == graph.default_coeff_nir
            dense.static.coeff_nir
        else
            nothing
        end
    coeff = cached_coeffs === nothing || T !== Float64 ? zeros(T, n_nodes) : cached_coeffs
    @inbounds for (i, nid) in pairs(node_ids)
        if cached_coeffs === nothing || T !== Float64
            coeff[i] = get(coeff_by_node, nid, default_coeff)
        end
    end
    return initial, coeff
end

function _dense_scattering_band_arrays(
    graph::ScatteringTransferGraph,
    initial_power::AbstractVector{<:Real},
    coeff_by_node::Dict{Int,Float64},
    default_coeff::Float64,
    ::Type{T}=Float64,
) where {T<:AbstractFloat}
    node_ids = graph.node_ids
    n_nodes = length(node_ids)
    length(initial_power) == n_nodes ||
        throw(ArgumentError("dense initial scattering power length $(length(initial_power)) does not match graph node count $n_nodes"))
    dense = _dense_scattering_graph(graph)
    initial =
        if initial_power isa Vector{T}
            initial_power
        else
            T.(initial_power)
        end
    cached_coeffs =
        if coeff_by_node === graph.coeff_par_by_node && default_coeff == graph.default_coeff_par
            dense.static.coeff_par
        elseif coeff_by_node === graph.coeff_nir_by_node && default_coeff == graph.default_coeff_nir
            dense.static.coeff_nir
        else
            nothing
        end
    coeff = cached_coeffs === nothing || T !== Float64 ? zeros(T, n_nodes) : cached_coeffs
    if cached_coeffs === nothing || T !== Float64
        @inbounds for (i, nid) in pairs(node_ids)
            coeff[i] = get(coeff_by_node, nid, default_coeff)
        end
    end
    return initial, coeff
end

function _scattering_static_edge_device_arrays(static::DenseScatteringStaticGraph, backend)
    return get!(static.device_cache, backend) do
        n_edges = length(static.counts)
        to_idx_dev = KernelAbstractions.allocate(backend, Int, n_edges)
        from_idx_dev = KernelAbstractions.allocate(backend, Int, n_edges)
        counts_dev = KernelAbstractions.allocate(backend, Int, n_edges)
        KernelAbstractions.copyto!(backend, to_idx_dev, static.to_idx)
        KernelAbstractions.copyto!(backend, from_idx_dev, static.from_idx)
        KernelAbstractions.copyto!(backend, counts_dev, static.counts)
        (
            to_idx_dev=to_idx_dev,
            from_idx_dev=from_idx_dev,
            counts_dev=counts_dev,
            n_edges=n_edges,
        )
    end
end

function _scattering_graph_all_hits_device_array(dense::DenseScatteringGraph, backend)
    return get!(dense.device_cache, backend) do
        all_hits_dev = KernelAbstractions.allocate(backend, Int, length(dense.all_hits))
        KernelAbstractions.copyto!(backend, all_hits_dev, dense.all_hits)
        all_hits_dev
    end
end

function _cached_static_coeff_device_array(
    static::DenseScatteringStaticGraph,
    backend,
    coeff_values::Vector{Float64},
    coeff_key::Symbol,
    ::Type{T},
) where {T<:AbstractFloat}
    cache_key = (objectid(backend), coeff_key, T)
    return get!(static.device_cache, cache_key) do
        coeff_dev = KernelAbstractions.allocate(backend, T, length(coeff_values))
        host_values = T === Float64 ? coeff_values : T.(coeff_values)
        KernelAbstractions.copyto!(backend, coeff_dev, host_values)
        coeff_dev
    end
end

function _scattering_coeff_device_array(
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    default_coeff::Float64,
    backend,
    ::Type{T},
) where {T<:AbstractFloat}
    dense = _dense_scattering_graph(graph)
    static = dense.static
    if coeff_by_node === graph.coeff_par_by_node && default_coeff == graph.default_coeff_par
        return _cached_static_coeff_device_array(static, backend, static.coeff_par, :coeff_par, T)
    elseif coeff_by_node === graph.coeff_nir_by_node && default_coeff == graph.default_coeff_nir
        return _cached_static_coeff_device_array(static, backend, static.coeff_nir, :coeff_nir, T)
    end

    coeff = zeros(T, length(graph.node_ids))
    @inbounds for (i, nid) in pairs(graph.node_ids)
        coeff[i] = T(get(coeff_by_node, nid, default_coeff))
    end
    coeff_dev = KernelAbstractions.allocate(backend, T, length(coeff))
    KernelAbstractions.copyto!(backend, coeff_dev, coeff)
    return coeff_dev
end

function _copy_scattering_static_device_arrays(graph::ScatteringTransferGraph, backend)
    dense = _dense_scattering_graph(graph)
    all_hits_dev = _scattering_graph_all_hits_device_array(dense, backend)
    edge_arrays = _scattering_static_edge_device_arrays(dense.static, backend)
    return (
        all_hits_dev=all_hits_dev,
        to_idx_dev=edge_arrays.to_idx_dev,
        from_idx_dev=edge_arrays.from_idx_dev,
        counts_dev=edge_arrays.counts_dev,
        n_edges=edge_arrays.n_edges,
    )
end

function _dense_vector_to_node_dict(node_ids::Vector{Int}, values::AbstractVector{<:Real})
    out = Dict{Int,Float64}()
    sizehint!(out, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[node_ids[i]] = Float64(values[i])
    end
    return out
end

_dense_float_vector(values::Vector{Float64}) = values
_dense_float_vector(values::AbstractVector{<:Real}) = Float64.(values)

function _dense_vector_from_node_dict(node_ids::Vector{Int}, values::Dict{Int,Float64})
    out = zeros(Float64, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[i] = get(values, node_ids[i], 0.0)
    end
    return out
end

function _scattering_result_from_dense(
    node_ids::Vector{Int},
    added_par::AbstractVector{<:Real},
    added_nir::AbstractVector{<:Real},
    iterations::Int,
    converged::Bool,
)
    dense_par = _dense_float_vector(added_par)
    dense_nir = _dense_float_vector(added_nir)
    return ScatteringResult(
        SpectralNodeValues(
            _all_dense_float_node_map(node_ids, dense_par),
            _all_dense_float_node_map(node_ids, dense_nir),
        ),
        iterations,
        converged,
        DenseScatteringResult(node_ids, DenseSpectralNodeValues(dense_par, dense_nir)),
    )
end

function _propagate_scattering_one_band_device_dense_only_static(
    initial_power_per_node::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
    static_device_arrays,
)
    node_ids = graph.node_ids
    n_nodes = length(node_ids)
    n_nodes == 0 && return Float64[], 0, true

    initial =
        if initial_power_per_node isa AbstractVector
            length(initial_power_per_node) == n_nodes ||
                throw(ArgumentError("dense initial scattering power length $(length(initial_power_per_node)) does not match graph node count $n_nodes"))
            initial_power_per_node isa Vector{scattering_eltype} ?
            initial_power_per_node :
            scattering_eltype.(initial_power_per_node)
        else
            _dense_initial_from_node_dict(node_ids, initial_power_per_node, scattering_eltype)
        end
    ref = sum(x -> Float64(x), initial)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    n_edges = static_device_arrays.n_edges

    current_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    next_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    added_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    hit_energy_dev = KernelAbstractions.allocate(backend, scattering_eltype, n_nodes)
    coeff_dev = _scattering_coeff_device_array(graph, coeff_by_node, default_coeff, backend, scattering_eltype)
    zero_host = zeros(scattering_eltype, n_nodes)
    next_host = Vector{scattering_eltype}(undef, n_nodes)
    added_host = Vector{scattering_eltype}(undef, n_nodes)

    KernelAbstractions.copyto!(backend, current_dev, initial)
    KernelAbstractions.copyto!(backend, next_dev, zero_host)
    KernelAbstractions.copyto!(backend, added_dev, zero_host)
    KernelAbstractions.copyto!(backend, hit_energy_dev, zero_host)

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
        hit_kernel(hit_energy_dev, current_dev, coeff_dev, static_device_arrays.all_hits_dev; ndrange=n_nodes)
        if use_atomic_transfer
            zero_kernel(next_dev; ndrange=n_nodes)
            transfer_by_edge_kernel(
                next_dev,
                static_device_arrays.to_idx_dev,
                static_device_arrays.from_idx_dev,
                static_device_arrays.counts_dev,
                hit_energy_dev;
                ndrange=n_edges,
            )
        else
            transfer_by_node_kernel(
                next_dev,
                static_device_arrays.to_idx_dev,
                static_device_arrays.from_idx_dev,
                static_device_arrays.counts_dev,
                hit_energy_dev,
                n_edges;
                ndrange=n_nodes,
            )
        end
        add_kernel(added_dev, next_dev; ndrange=n_nodes)
        KernelAbstractions.synchronize(backend)

        copyto!(next_host, next_dev)
        total_next = sum(x -> Float64(x), next_host)
        current_dev, next_dev = next_dev, current_dev
        if total_next < thr
            converged = true
            break
        end
    end

    copyto!(added_host, added_dev)
    return _dense_float_vector(added_host), iterations, converged
end

function _propagate_scattering_one_band_device_static(
    initial_power_per_node::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
    static_device_arrays,
)
    dense, iterations, converged = _propagate_scattering_one_band_device_dense_only_static(
        initial_power_per_node,
        graph,
        coeff_by_node,
        options,
        default_coeff,
        backend,
        workgroupsize,
        scattering_eltype,
        static_device_arrays,
    )
    return _dense_vector_to_node_dict(graph.node_ids, dense), iterations, converged, dense
end

function _propagate_scattering_one_band_device_dense_only(
    initial_power_per_node::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
)
    static_device_arrays = _copy_scattering_static_device_arrays(graph, backend)
    return _propagate_scattering_one_band_device_dense_only_static(
        initial_power_per_node,
        graph,
        coeff_by_node,
        options,
        default_coeff,
        backend,
        workgroupsize,
        scattering_eltype,
        static_device_arrays,
    )
end

function _propagate_scattering_one_band_device(
    initial_power_per_node::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
)
    static_device_arrays = _copy_scattering_static_device_arrays(graph, backend)
    return _propagate_scattering_one_band_device_static(
        initial_power_per_node,
        graph,
        coeff_by_node,
        options,
        default_coeff,
        backend,
        workgroupsize,
        scattering_eltype,
        static_device_arrays,
    )
end

function _propagate_scattering_two_bands_device(
    initial_par::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    coeff_par::Dict{Int,Float64},
    default_par::Float64,
    initial_nir::Union{Dict{Int,Float64},AbstractVector{<:Real}},
    coeff_nir::Dict{Int,Float64},
    default_nir::Float64,
    graph::ScatteringTransferGraph,
    options::LightOptions,
    backend,
    workgroupsize::Int,
    scattering_eltype::Type{<:AbstractFloat},
)
    static_device_arrays = _copy_scattering_static_device_arrays(graph, backend)
    added_par, it_par, conv_par, dense_par = _propagate_scattering_one_band_device_static(
        initial_par,
        graph,
        coeff_par,
        options,
        default_par,
        backend,
        workgroupsize,
        scattering_eltype,
        static_device_arrays,
    )
    added_nir, it_nir, conv_nir, dense_nir = _propagate_scattering_one_band_device_static(
        initial_nir,
        graph,
        coeff_nir,
        options,
        default_nir,
        backend,
        workgroupsize,
        scattering_eltype,
        static_device_arrays,
    )
    return (
        added_par=added_par,
        added_nir=added_nir,
        dense_par=dense_par,
        dense_nir=dense_nir,
        iterations=max(it_par, it_nir),
        converged=conv_par && conv_nir,
    )
end

_rastergpu_scattering_eltype(backend::RasterGPUScatteringBackend) =
    backend.config.backend isa KernelAbstractions.CPU ? Float64 : Float32

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    initial_power::AbstractVector{<:Real},
    options::LightOptions,
    backend::RasterGPUScatteringBackend;
    band::AbstractString="PAR",
    coeff_by_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    coeffs = isnothing(coeff_by_node) ? _coeff_by_node(graph, String(band), dflt) : coeff_by_node
    return _propagate_scattering_one_band_device_dense_only(
        initial_power,
        graph,
        coeffs,
        options,
        dflt,
        backend.config.backend,
        backend.config.workgroupsize,
        _rastergpu_scattering_eltype(backend),
    )
end

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RasterGPUScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    T = _rastergpu_scattering_eltype(backend)
    initial = _dense_initial_scattering_power(graph, first, initial_power_per_node, band, T)
    coeffs = _coeff_by_node(graph, String(band), dflt)
    return _propagate_scattering_one_band_device_dense_only(
        initial,
        graph,
        coeffs,
        options,
        dflt,
        backend.config.backend,
        backend.config.workgroupsize,
        T,
    )
end

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    initial_power::AbstractVector{<:Real},
    options::LightOptions,
    ::RaycastScatteringBackend;
    band::AbstractString="PAR",
    coeff_by_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    coeffs = isnothing(coeff_by_node) ? _coeff_by_node(graph, String(band), dflt) : coeff_by_node
    return _propagate_scattering_one_band_dense_only(
        initial_power,
        graph,
        coeffs,
        options,
        dflt,
    )
end

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    initial = _dense_initial_scattering_power(graph, first, initial_power_per_node, band)
    return _propagate_scattering_one_band_dense_only(
        initial,
        graph,
        _coeff_by_node(graph, String(band), dflt),
        options,
        dflt,
    )
end

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    initial_power::AbstractVector{<:Real},
    options::LightOptions,
    ::LinksScatteringBackend;
    band::AbstractString="PAR",
    coeff_by_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    return _compute_scattering_band_dense(
        graph,
        initial_power,
        options,
        RaycastScatteringBackend();
        band=band,
        coeff_by_node=coeff_by_node,
        default_coeff=default_coeff,
    )
end

function _compute_scattering_band_dense(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    return _compute_scattering_band_dense(
        graph,
        first,
        options,
        RaycastScatteringBackend();
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
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
    backend::RasterGPUScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    T = _rastergpu_scattering_eltype(backend)
    initial = _dense_initial_scattering_power(graph, first, initial_power_per_node, band, T)
    coeffs = _coeff_by_node(graph, String(band), dflt)
    added, iterations, converged, dense = _propagate_scattering_one_band_device(
        initial,
        graph,
        coeffs,
        options,
        dflt,
        backend.config.backend,
        backend.config.workgroupsize,
        T,
    )
    return (
        added_power_per_node=added,
        iterations=iterations,
        converged=converged,
        node_ids=graph.node_ids,
        dense_added_power_per_node=dense,
    )
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
    added, iterations, converged, dense = _propagate_scattering_one_band_dense(
        initial,
        graph,
        _coeff_by_node(graph, String(band), dflt),
        options,
        dflt,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged, node_ids=graph.node_ids, dense_added_power_per_node=dense)
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
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RasterGPUScatteringBackend;
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

    added_par, it_par, conv_par, dense_par = _propagate_scattering_one_band_dense(
        initial_par,
        graph,
        graph.coeff_par_by_node,
        options,
        graph.default_coeff_par,
    )
    added_nir, it_nir, conv_nir, dense_nir = _propagate_scattering_one_band_dense(
        initial_nir,
        graph,
        graph.coeff_nir_by_node,
        options,
        graph.default_coeff_nir,
    )

    ScatteringResult(
        SpectralNodeValues(added_par, added_nir),
        max(it_par, it_nir),
        conv_par && conv_nir,
        DenseScatteringResult(graph.node_ids, DenseSpectralNodeValues(dense_par, dense_nir)),
    )
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
    backend::RasterGPUScatteringBackend,
)
    T = _rastergpu_scattering_eltype(backend)
    initial_par = _dense_initial_scattering_power(graph, first, nothing, "PAR", T)
    initial_nir = _dense_initial_scattering_power(graph, first, nothing, "NIR", T)
    propagated = _propagate_scattering_two_bands_device(
        initial_par,
        graph.coeff_par_by_node,
        graph.default_coeff_par,
        initial_nir,
        graph.coeff_nir_by_node,
        graph.default_coeff_nir,
        graph,
        options,
        backend.config.backend,
        backend.config.workgroupsize,
        T,
    )

    return ScatteringResult(
        SpectralNodeValues(propagated.added_par, propagated.added_nir),
        propagated.iterations,
        propagated.converged,
        DenseScatteringResult(
            graph.node_ids,
            DenseSpectralNodeValues(propagated.dense_par, propagated.dense_nir),
        ),
    )
end

function compute_scattering(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    backend::RasterGPUScatteringBackend,
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
