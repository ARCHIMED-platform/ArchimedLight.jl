struct ProjectionCacheContext
    cache_dir::String
    scene_key::UInt64
end

const HitRecord = Tuple{Float64,Int}
const FlatPoolInt = Int32

"""
    SmallHitStack

Compact per-pixel hit stack optimized for the common short-stack case.

The first two hits are stored directly in the struct, so pixels with 0-2 hits
do not allocate a separate heap vector. When a third hit arrives, the stack
"spills" to a regular `Vector{HitRecord}` stored in `spill`.

This is a tradeoff:
- small inline capacity keeps each occupied pixel stack compact
- larger inline capacity would reduce spills in dense canopies, but would make
  every stack heavier even when it only contains one or two hits

The current inline capacity of 2 is a conservative default, not a universal
optimum. Users can override the storage mode with
`LightOptions(pixel_hit_stack_mode=...)`.
"""
mutable struct SmallHitStack <: AbstractVector{HitRecord}
    len::Int32
    hit1::HitRecord
    hit2::HitRecord
    spill::Union{Nothing,Vector{HitRecord}}
end

SmallHitStack() = SmallHitStack(0, (0.0, 0), (0.0, 0), nothing)

"""
    DensePixelHits

Dense storage for per-pixel hit stacks. This avoids hashing pixel indices during
projection when the plotbox is small enough that a flat table is cheaper than a
`Dict`.
"""
struct DensePixelHits{S}
    stacks::Vector{Union{Nothing,S}}
end

struct DenseUpperPixelHits
    heights::Vector{Float64}
    nodes::Vector{Int}
    occupied::Vector{Int}
end

mutable struct FlatPixelHitBuilder
    counts::Vector{Int}
    occupied::Vector{Int}
    total_hits::Int
    pixels::Vector{FlatPoolInt}
    heights::Vector{Float32}
    nodes::Vector{FlatPoolInt}
end

FlatPixelHitBuilder(n_pixels::Int) =
    FlatPixelHitBuilder(zeros(Int, n_pixels), Int[], 0, FlatPoolInt[], Float32[], FlatPoolInt[])

mutable struct FlatPixelHits
    starts::Vector{Int}
    counts::Vector{Int}
    occupied::Vector{Int}
    sorted::BitVector
    heights::Vector{Float32}
    nodes::Vector{FlatPoolInt}
end

struct FlatPixelHitStack <: AbstractVector{HitRecord}
    parent::FlatPixelHits
    pixel_idx::Int
end

struct UpperHitStack <: AbstractVector{HitRecord}
    parent::DenseUpperPixelHits
    pixel_idx::Int
end

mutable struct RasterScanlineScratch
    minY::Vector{Int}
    maxY::Vector{Int}
end

RasterScanlineScratch() = RasterScanlineScratch(Int[], Int[])

struct DirectionProjectionResult{PH}
    pixel_hits::PH
    node_hits::Dict{Int,Int}
    projected_mesh_area::Dict{Int,Float64}
    projected_pixels_area::Dict{Int,Float64}
end

struct DenseDirectionProjectionResult{PH}
    pixel_hits::PH
    node_hits::Vector{Int}
    projected_mesh_area::Vector{Float64}
    projected_pixels_area::Vector{Float64}
end

struct InterceptionSceneData
    vertices
    faces::Vector{PlantGeom.Face3}
    face2node::Vector{Int}
    face2node_index::Vector{Int}
    node_ids::Vector{Int}
    node_index::Dict{Int,Int}
    plotbox
    node_group::Dict{Int,String}
    node_group_by_index::Vector{String}
    pavement_node_mask::Vector{Bool}
    node_type::Dict{Int,String}
    node_type_by_index::Vector{String}
    raycore_instanced_geometry::Any
end

struct RaycorePrototypeInstances
    mesh
    transforms::Vector{GeometryBasics.Mat4f}
    node_indices::Vector{UInt32}
end

struct RaycoreInstancedSceneGeometry
    prototypes::Vector{RaycorePrototypeInstances}
    fallback_mesh
    metadata_stride::Int
    prototype_face_count::Int
    prototype_node_count::Int
    fallback_face_count::Int
end

struct PreparedInterceptionData
    geometry::InterceptionSceneData
    node_interception_by_index::Vector{Union{Nothing,InterceptionModel}}
    node_transparency_by_index::Vector{Float64}
    virtual_nodes::Set{Int}
    virtual_node_mask::Vector{Bool}
    upper_hit::Bool
    cache_ctx::Union{Nothing,ProjectionCacheContext}
    emitter_par_power_per_node::Dict{Int,Float64}
    emitter_nir_power_per_node::Dict{Int,Float64}
    emitter_par_power_by_index::Vector{Float64}
    emitter_nir_power_by_index::Vector{Float64}
    emitter_nodes::Set{Int}
    emitter_node_mask::Vector{Bool}
    component_area_per_node::Union{Nothing,Dict{Int,Float64}}
    absorption_par_per_node::Union{Nothing,Dict{Int,Float64}}
    absorption_nir_per_node::Union{Nothing,Dict{Int,Float64}}
end

struct RaycoreHitDecoder{D}
    node_index_by_instance_metadata::Vector{UInt32}
    node_index_by_instance_metadata_dev::D
    metadata_stride::Int
    instance_count::Int
end

function _raycore_hit_decoder(
    node_index_by_instance_metadata::Vector{UInt32},
    metadata_stride::Int,
    instance_count::Int,
    backend,
)
    node_index_by_instance_metadata_dev =
        KernelAbstractions.allocate(backend, UInt32, length(node_index_by_instance_metadata))
    KernelAbstractions.copyto!(backend, node_index_by_instance_metadata_dev, node_index_by_instance_metadata)
    return RaycoreHitDecoder(
        node_index_by_instance_metadata,
        node_index_by_instance_metadata_dev,
        metadata_stride,
        instance_count,
    )
end

function _raycore_identity_hit_decoder(
    prepared::PreparedInterceptionData,
    backend,
    instance_count::Int,
)
    metadata_stride = length(prepared.geometry.node_ids) + 1
    table = zeros(UInt32, max(instance_count, 0) * metadata_stride)
    @inbounds for instance_idx in 1:instance_count
        base = (instance_idx - 1) * metadata_stride
        for node_idx in 1:(metadata_stride - 1)
            table[base+node_idx+1] = UInt32(node_idx)
        end
    end
    return _raycore_hit_decoder(table, metadata_stride, instance_count, backend)
end

@inline function _raycore_decode_node_index(
    node_index_by_instance_metadata,
    metadata_stride::Int,
    metadata::UInt32,
    instance_idx::UInt32,
)
    (metadata == 0x00000000 || instance_idx == 0x00000000) && return UInt32(0)
    table_idx = (Int(instance_idx) - 1) * metadata_stride + Int(metadata) + 1
    @inbounds return node_index_by_instance_metadata[table_idx]
end

@inline function _raycore_decode_node_index(
    decoder::RaycoreHitDecoder,
    metadata::UInt32,
    instance_idx::UInt32,
)
    return _raycore_decode_node_index(
        decoder.node_index_by_instance_metadata,
        decoder.metadata_stride,
        metadata,
        instance_idx,
    )
end

struct RaycoreSceneData{P,T,K,B,N,C,S,J,H,O,E,F,D,G,Q,A,V,M,I,U,R}
    prepared::P
    tlas::T
    kernel_tlas::K
    backend::B
    top_nodes_dev::N
    top_nodes_host::Vector{UInt32}
    stack_counts_dev::C
    stack_nodes_dev::S
    stack_instance_indices_dev::J
    stack_heights_dev::H
    stack_overflow_dev::O
    edge_keys_dev::E
    edge_key_counts_dev::F
    dense_edge_counts_dev::D
    dense_edge_counts_host::Vector{Int32}
    node_counts_dev::G
    node_transparency_dev::Q
    sector_area_dev::A
    node_counts_host::Vector{Int32}
    sector_area_host::Vector{Float32}
    virtual_node_mask_dev::V
    pavement_node_mask_dev::M
    node_ids_dev::I
    stack_counts_host::Vector{Int32}
    stack_nodes_host::Vector{UInt32}
    stack_instance_indices_host::Vector{UInt32}
    stack_heights_host::Vector{Float32}
    stack_overflow_host::Vector{Bool}
    edge_keys_host::Vector{UInt64}
    edge_key_counts_host::Vector{Int32}
    edge_compact_host::Vector{UInt64}
    projection_far_cache::Dict{Tuple{Float64,Float64,Float64,Bool},Float32}
    projected_mesh_area_cache::Dict{Tuple{Float64,Float64,Float64},Vector{Float64}}
    raster_compat_projection_cache::Dict{Tuple{Float64,Float64,Float64,Bool},Any}
    hit_decoder::U
    area_ratio_cache::R
    workgroupsize::Int
    max_hits_per_pixel::Int
    hit_epsilon::Float32
    edge_accumulation::Symbol
    dense_edge_limit_bytes::Int
    validate::Bool
    vertical_span::Float64
    chunked_tlas::Bool
    geometry_mode::Symbol
end

const RAYCORE_INVALID_NODE = UInt32(0xffffffff)
const RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY =
    isdefined(Raycore, :SOFTWARE_TRAVERSAL_STACK_CAPACITY) ?
    Int(getfield(Raycore, :SOFTWARE_TRAVERSAL_STACK_CAPACITY)) :
    32

_raycore_software_traversal_stack_capacity() = RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY

function _raycore_bvh_depth(nodes)
    isempty(nodes) && return 0
    function visit(node_idx::Int, depth::Int)
        1 <= node_idx <= length(nodes) || return depth
        node = nodes[node_idx]
        node.child0 == RAYCORE_INVALID_NODE && return depth
        return max(visit(Int(node.child0), depth + 1), visit(Int(node.child1), depth + 1))
    end
    return visit(1, 1)
end

function _raycore_max_blas_depth_and_nodes(data::RaycoreSceneData)
    static = data.tlas.static_tlas
    nodes = Array(static.all_blas_nodes)
    descriptors = Array(static.blas_descriptors)
    isempty(nodes) && return (depth=0, node_count=0)
    isempty(descriptors) && return (depth=_raycore_bvh_depth(nodes), node_count=length(nodes))

    max_depth = 0
    for idx in eachindex(descriptors)
        first_node = Int(descriptors[idx].nodes_offset) + 1
        last_node = idx == length(descriptors) ? length(nodes) : Int(descriptors[idx+1].nodes_offset)
        first_node <= last_node || continue
        max_depth = max(max_depth, _raycore_bvh_depth(@view nodes[first_node:last_node]))
    end
    return (depth=max_depth, node_count=length(nodes))
end

function _raycore_max_blas_depth(data::RaycoreSceneData)
    return _raycore_max_blas_depth_and_nodes(data).depth
end

function _raycore_traversal_stack_depth_status(data::RaycoreSceneData)
    blas = _raycore_max_blas_depth_and_nodes(data)
    capacity = _raycore_software_traversal_stack_capacity()
    return (
        capacity=capacity,
        blas_depth=blas.depth,
        blas_node_count=blas.node_count,
        exceeds=blas.depth > capacity,
    )
end

function _raycore_scene_shape_summary(data::Union{Nothing,RaycoreSceneData})
    data === nothing && return (
        geometry_mode=:none,
        chunked_tlas=false,
        tlas_instances=0,
        tlas_geometries=0,
        node_count=0,
        raycore_traversal_stack_capacity=_raycore_software_traversal_stack_capacity(),
        max_blas_depth=0,
        max_blas_depth_over_stack_capacity=false,
        blas_node_count=0,
        expanded_face_count=0,
        expanded_face_instance_upper_bound=0,
        reference_prototype_count=0,
        reference_prototype_node_count=0,
        reference_prototype_face_count=0,
        reference_fallback_face_count=0,
        reference_compact_face_count=0,
        dense_edge_pairs=0,
        edge_key_capacity=0,
    )

    geometry = data.prepared.geometry
    n_faces = length(geometry.faces)
    instanced = geometry.raycore_instanced_geometry
    reference_prototype_count = instanced === nothing ? 0 : length(instanced.prototypes)
    reference_prototype_node_count = instanced === nothing ? 0 : instanced.prototype_node_count
    reference_prototype_face_count = instanced === nothing ? 0 : instanced.prototype_face_count
    reference_fallback_face_count = instanced === nothing ? 0 : instanced.fallback_face_count
    reference_compact_face_count = reference_prototype_face_count + reference_fallback_face_count
    tlas_instances = length(data.tlas.instances)
    tlas_geometries = Raycore.n_geometries(data.tlas)
    traversal_stack = _raycore_traversal_stack_depth_status(data)

    return (
        geometry_mode=data.geometry_mode,
        chunked_tlas=data.chunked_tlas,
        tlas_instances=tlas_instances,
        tlas_geometries=tlas_geometries,
        node_count=length(geometry.node_ids),
        raycore_traversal_stack_capacity=traversal_stack.capacity,
        max_blas_depth=traversal_stack.blas_depth,
        max_blas_depth_over_stack_capacity=traversal_stack.exceeds,
        blas_node_count=traversal_stack.blas_node_count,
        expanded_face_count=n_faces,
        expanded_face_instance_upper_bound=tlas_instances * n_faces,
        reference_prototype_count=reference_prototype_count,
        reference_prototype_node_count=reference_prototype_node_count,
        reference_prototype_face_count=reference_prototype_face_count,
        reference_fallback_face_count=reference_fallback_face_count,
        reference_compact_face_count=reference_compact_face_count,
        dense_edge_pairs=data.dense_edge_counts_dev === nothing ? 0 : length(data.dense_edge_counts_dev),
        edge_key_capacity=data.edge_keys_dev === nothing ? 0 : length(data.edge_keys_dev),
    )
end

Base.IndexStyle(::Type{<:SmallHitStack}) = IndexLinear()
Base.size(stack::SmallHitStack) = (Int(stack.len),)
Base.length(stack::SmallHitStack) = Int(stack.len)

Base.IndexStyle(::Type{UpperHitStack}) = IndexLinear()
Base.size(stack::UpperHitStack) = (length(stack),)
Base.length(stack::UpperHitStack) = @inbounds(stack.parent.nodes[stack.pixel_idx] == 0 ? 0 : 1)

function Base.getindex(stack::UpperHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    @inbounds return (stack.parent.heights[stack.pixel_idx], stack.parent.nodes[stack.pixel_idx])
end

function Base.setindex!(stack::UpperHitStack, hit::HitRecord, i::Int)
    @boundscheck checkbounds(stack, i)
    @inbounds begin
        stack.parent.heights[stack.pixel_idx] = hit[1]
        stack.parent.nodes[stack.pixel_idx] = hit[2]
    end
    return stack
end

function Base.deleteat!(stack::UpperHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    @inbounds begin
        stack.parent.heights[stack.pixel_idx] = 0.0
        stack.parent.nodes[stack.pixel_idx] = 0
    end
    return stack
end

function Base.getindex(stack::SmallHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    if stack.spill === nothing
        return i == 1 ? stack.hit1 : stack.hit2
    end
    return stack.spill[i]
end

function Base.setindex!(stack::SmallHitStack, hit::HitRecord, i::Int)
    @boundscheck checkbounds(stack, i)
    if stack.spill === nothing
        if i == 1
            stack.hit1 = hit
        else
            stack.hit2 = hit
        end
    else
        stack.spill[i] = hit
    end
    return stack
end

function Base.push!(stack::SmallHitStack, hit::HitRecord)
    len = length(stack)
    if len == 0
        stack.hit1 = hit
    elseif len == 1
        stack.hit2 = hit
    elseif len == 2
        stack.spill = HitRecord[stack.hit1, stack.hit2, hit]
    else
        push!(stack.spill, hit)
    end
    stack.len = Int32(len + 1)
    return stack
end

function Base.deleteat!(stack::SmallHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    len = length(stack)
    if stack.spill === nothing
        if len == 1
            stack.hit1 = (0.0, 0)
        elseif i == 1
            stack.hit1 = stack.hit2
            stack.hit2 = (0.0, 0)
        else
            stack.hit2 = (0.0, 0)
        end
    else
        deleteat!(stack.spill, i)
        if len - 1 == 2
            stack.hit1 = stack.spill[1]
            stack.hit2 = stack.spill[2]
            stack.spill = nothing
        elseif len - 1 == 1
            stack.hit1 = stack.spill[1]
            stack.hit2 = (0.0, 0)
            stack.spill = nothing
        elseif len - 1 == 0
            stack.hit1 = (0.0, 0)
            stack.hit2 = (0.0, 0)
            stack.spill = nothing
        end
    end
    stack.len = Int32(len - 1)
    return stack
end

function Base.sort!(stack::SmallHitStack; by=identity, rev::Bool=false, alg=Base.DEFAULT_STABLE)
    len = length(stack)
    len <= 1 && return stack
    if stack.spill === nothing
        a = by(stack.hit1)
        b = by(stack.hit2)
        swap = rev ? (a < b) : (a > b)
        if swap
            stack.hit1, stack.hit2 = stack.hit2, stack.hit1
        end
        return stack
    end
    sort!(stack.spill; by=by, rev=rev, alg=alg)
    return stack
end

@inline _stack_hit(stack, i::Int) = stack[i]
@inline function _stack_hit(stack::SmallHitStack, i::Int)
    if stack.spill === nothing
        return i == 1 ? stack.hit1 : stack.hit2
    end
    return @inbounds stack.spill[i]
end

@inline _hit_sort_lt(a, b) = _hit_height(a) > _hit_height(b)

@inline function _sort_hit_stack!(stack::SmallHitStack)
    len = length(stack)
    len <= 1 && return stack
    if stack.spill === nothing
        if _hit_height(stack.hit1) < _hit_height(stack.hit2)
            stack.hit1, stack.hit2 = stack.hit2, stack.hit1
        end
        return stack
    end
    sort!(stack.spill; lt=_hit_sort_lt, alg=Base.Sort.MergeSort)
    return stack
end

@inline function _sort_small_hit_vector!(stack::Vector{HitRecord})
    @inbounds for i in 2:length(stack)
        x = stack[i]
        xh = _hit_height(x)
        j = i - 1
        while j >= 1 && _hit_height(stack[j]) < xh
            stack[j+1] = stack[j]
            j -= 1
        end
        stack[j+1] = x
    end
    return stack
end

@inline function _sort_hit_stack!(stack::Vector{HitRecord})
    len = length(stack)
    len <= 1 && return stack
    if len <= 8
        return _sort_small_hit_vector!(stack)
    end
    sort!(stack; lt=_hit_sort_lt, alg=Base.Sort.MergeSort)
    return stack
end

@inline function _sort_small_hit_stack!(stack::FlatPixelHitStack)
    len = length(stack)
    start = _flat_stack_start(stack)
    heights = stack.parent.heights
    nodes = stack.parent.nodes
    @inbounds for i in 2:len
        idx = start + i - 1
        xh = heights[idx]
        xn = nodes[idx]
        j = i - 1
        while j >= 1 && heights[start+j-1] < xh
            src = start + j - 1
            dst = src + 1
            heights[dst] = heights[src]
            nodes[dst] = nodes[src]
            j -= 1
        end
        dst = start + j
        heights[dst] = xh
        nodes[dst] = xn
    end
    return stack
end

@inline function _sort_hit_stack!(stack::FlatPixelHitStack)
    stack.parent.sorted[stack.pixel_idx] && return stack
    len = length(stack)
    if len <= 1
        stack.parent.sorted[stack.pixel_idx] = true
        return stack
    end
    if len <= 8
        _sort_small_hit_stack!(stack)
        stack.parent.sorted[stack.pixel_idx] = true
        return stack
    end
    _sort_small_hit_stack!(stack)
    stack.parent.sorted[stack.pixel_idx] = true
    return stack
end

@inline _sort_hit_stack!(stack) = sort!(stack, by=_hit_height, rev=true, alg=Base.Sort.MergeSort)
@inline _sort_hit_stack!(stack::UpperHitStack) = stack

function _sort_projection_pixel_stacks!(projection)
    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        _sort_hit_stack!(stack)
    end
    return projection
end

@inline _stack_hit_height(stack, i::Int) = _hit_height(_stack_hit(stack, i))
@inline _stack_hit_node(stack, i::Int) = _hit_node(_stack_hit(stack, i))
@inline function _stack_hit_height(stack::UpperHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    return @inbounds stack.parent.heights[stack.pixel_idx]
end
@inline function _stack_hit_node(stack::UpperHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    return @inbounds stack.parent.nodes[stack.pixel_idx]
end
@inline function _stack_hit_height(stack::FlatPixelHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    return @inbounds Float64(stack.parent.heights[_flat_stack_start(stack)+i-1])
end
@inline function _stack_hit_node(stack::FlatPixelHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    return @inbounds Int(stack.parent.nodes[_flat_stack_start(stack)+i-1])
end

@inline function _accumulate_emitter_transfer_counts_dense!(
    edge_counts::Dict{UInt64,Int},
    total_from::Dict{Int,Int},
    stack,
    emitter_node_mask::Vector{Bool},
    node_ids::Vector{Int},
)
    @inbounds for j in eachindex(stack)
        src_idx = _stack_hit_node(stack, j)
        emitter_node_mask[src_idx] || continue

        to_idx = 0
        for k in (j+1):length(stack)
            node_idx = _stack_hit_node(stack, k)
            emitter_node_mask[node_idx] && continue
            to_idx = node_idx
            break
        end
        to_idx == 0 && continue

        src = node_ids[src_idx]
        to = node_ids[to_idx]
        edge = _pack_emitter_edge(to, src)
        edge_counts[edge] = get(edge_counts, edge, 0) + 1
        total_from[src] = get(total_from, src, 0) + 1
    end
    return nothing
end

@inline function _accumulate_emitter_transfer_counts_dense!(
    edge_counts::Dict{UInt64,Int},
    total_from::Dict{Int,Int},
    stack::FlatPixelHitStack,
    emitter_node_mask::Vector{Bool},
    node_ids::Vector{Int},
)
    start = _flat_stack_start(stack)
    nodes = stack.parent.nodes
    n_hits = length(stack)
    @inbounds for j in 0:(n_hits-1)
        src_idx = Int(nodes[start+j])
        emitter_node_mask[src_idx] || continue

        to_idx = 0
        for k in (j+1):(n_hits-1)
            node_idx = Int(nodes[start+k])
            emitter_node_mask[node_idx] && continue
            to_idx = node_idx
            break
        end
        to_idx == 0 && continue

        src = node_ids[src_idx]
        to = node_ids[to_idx]
        edge = _pack_emitter_edge(to, src)
        edge_counts[edge] = get(edge_counts, edge, 0) + 1
        total_from[src] = get(total_from, src, 0) + 1
    end
    return nothing
end

@inline function _accumulate_visible_area_dense!(
    visible_area::Vector{Float64},
    projection::DenseDirectionProjectionResult,
    stack::UpperHitStack,
    options::LightOptions,
    pixel_area::Float64,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
)
    isempty(stack) && return nothing
    node_idx = _stack_hit_node(stack, 1)
    ratio = _projection_area_ratio(projection, options, node_idx)
    visible_area[node_idx] += pixel_area * ratio
    return nothing
end

@inline function _accumulate_visible_area_dense!(
    visible_area::Vector{Float64},
    projection::DenseDirectionProjectionResult,
    stack,
    options::LightOptions,
    pixel_area::Float64,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
)
    @inbounds for hit in stack
        node_idx = _hit_node(hit)
        if virtual_node_mask[node_idx]
            ratio = _projection_area_ratio(projection, options, node_idx)
            visible_area[node_idx] += pixel_area * ratio
            continue
        end

        transparency = node_transparency_by_index[node_idx]
        intercepted_fraction = 1.0 - transparency
        if intercepted_fraction > 0.0
            ratio = _projection_area_ratio(projection, options, node_idx)
            visible_area[node_idx] += pixel_area * intercepted_fraction * ratio
        end
        transparency > 0.0 || break
    end
    return nothing
end

@inline function _accumulate_visible_area_dense!(
    visible_area::Vector{Float64},
    projection::DenseDirectionProjectionResult,
    stack::FlatPixelHitStack,
    options::LightOptions,
    pixel_area::Float64,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
)
    start = _flat_stack_start(stack)
    nodes = stack.parent.nodes
    n_hits = length(stack)
    @inbounds for h in 0:(n_hits-1)
        node_idx = Int(nodes[start+h])
        if virtual_node_mask[node_idx]
            ratio = _projection_area_ratio(projection, options, node_idx)
            visible_area[node_idx] += pixel_area * ratio
            continue
        end

        transparency = node_transparency_by_index[node_idx]
        intercepted_fraction = 1.0 - transparency
        if intercepted_fraction > 0.0
            ratio = _projection_area_ratio(projection, options, node_idx)
            visible_area[node_idx] += pixel_area * intercepted_fraction * ratio
        end
        transparency > 0.0 || break
    end
    return nothing
end

# Dense tables avoid hash traffic on packed canopies, but above this size the empty-cell
# overhead starts to outweigh the benefit on sparse scenes.
const _DENSE_PIXEL_HITS_MAX_CELLS = 500_000
const _AUTO_VECTOR_PIXEL_HITS_MIN_CELLS = 300_000

DensePixelHits(::Type{S}, n::Int) where {S} = DensePixelHits{S}(fill(nothing, n))
function DenseUpperPixelHits(n::Int)
    occupied = Int[]
    # Dense upper-hit tables already commit to per-pixel arrays, so reserving
    # the occupied index list up front avoids repeated growth in the hot append path.
    sizehint!(occupied, n)
    return DenseUpperPixelHits(zeros(Float64, n), zeros(Int, n), occupied)
end

Base.get(pixel_hits::DensePixelHits, idx::Int, default) =
    (1 <= idx <= length(pixel_hits.stacks) && pixel_hits.stacks[idx] !== nothing) ? pixel_hits.stacks[idx] : default

Base.get(pixel_hits::DenseUpperPixelHits, idx::Int, default) =
    (1 <= idx <= length(pixel_hits.nodes) && pixel_hits.nodes[idx] != 0) ? UpperHitStack(pixel_hits, idx) : default

function Base.get!(f::F, pixel_hits::DensePixelHits, idx::Int) where {F}
    stack = pixel_hits.stacks[idx]
    if stack === nothing
        stack = f()
        pixel_hits.stacks[idx] = stack
    end
    return stack
end

function Base.delete!(pixel_hits::DensePixelHits, idx::Int)
    1 <= idx <= length(pixel_hits.stacks) && (pixel_hits.stacks[idx] = nothing)
    return pixel_hits
end

function Base.delete!(pixel_hits::DenseUpperPixelHits, idx::Int)
    if 1 <= idx <= length(pixel_hits.nodes)
        pixel_hits.heights[idx] = 0.0
        pixel_hits.nodes[idx] = 0
    end
    return pixel_hits
end

Base.values(pixel_hits::DensePixelHits) = (stack for stack in pixel_hits.stacks if stack !== nothing)
Base.values(pixel_hits::DenseUpperPixelHits) =
    (UpperHitStack(pixel_hits, idx) for idx in pixel_hits.occupied if pixel_hits.nodes[idx] != 0)

Base.IndexStyle(::Type{FlatPixelHitStack}) = IndexLinear()
Base.size(stack::FlatPixelHitStack) = (length(stack),)
Base.length(stack::FlatPixelHitStack) = stack.parent.counts[stack.pixel_idx]

@inline _flat_stack_start(stack::FlatPixelHitStack) = stack.parent.starts[stack.pixel_idx]

function Base.getindex(stack::FlatPixelHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    idx = _flat_stack_start(stack) + i - 1
    return (Float64(stack.parent.heights[idx]), Int(stack.parent.nodes[idx]))
end

function Base.setindex!(stack::FlatPixelHitStack, hit::HitRecord, i::Int)
    @boundscheck checkbounds(stack, i)
    idx = _flat_stack_start(stack) + i - 1
    stack.parent.heights[idx] = Float32(hit[1])
    stack.parent.nodes[idx] = FlatPoolInt(hit[2])
    stack.parent.sorted[stack.pixel_idx] = false
    return stack
end

function Base.deleteat!(stack::FlatPixelHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    count = length(stack)
    start = _flat_stack_start(stack)
    @inbounds for pos in i:(count-1)
        src = start + pos
        dst = start + pos - 1
        stack.parent.heights[dst] = stack.parent.heights[src]
        stack.parent.nodes[dst] = stack.parent.nodes[src]
    end
    stack.parent.counts[stack.pixel_idx] = count - 1
    return stack
end

Base.get(pixel_hits::FlatPixelHits, idx::Int, default) =
    (1 <= idx <= length(pixel_hits.counts) && pixel_hits.counts[idx] > 0) ? FlatPixelHitStack(pixel_hits, idx) : default

function Base.delete!(pixel_hits::FlatPixelHits, idx::Int)
    if 1 <= idx <= length(pixel_hits.counts)
        pixel_hits.counts[idx] = 0
        pixel_hits.sorted[idx] = false
    end
    return pixel_hits
end

Base.values(pixel_hits::FlatPixelHits) =
    (FlatPixelHitStack(pixel_hits, idx) for idx in pixel_hits.occupied if pixel_hits.counts[idx] > 0)

@inline _use_dense_pixel_hits(plotbox) = (plotbox.nx * plotbox.ny) <= _DENSE_PIXEL_HITS_MAX_CELLS

"""
    _pixel_hit_stack_mode(options)

Normalize the user-facing pixel-hit stack storage selector.

Accepted values are:
- `"auto"`: current validated optimized default
- `"small"`: force `SmallHitStack`
- `"vector"`: force the legacy `Vector{HitRecord}` representation

`"auto"` keeps `SmallHitStack` for the validated sparse and mid-density cases,
but switches to the legacy `Vector{HitRecord}` storage when the plot raster is
large enough that dense canopies are likely to spill `SmallHitStack` constantly.
This favors large toric canopies like the bundled coffee example without
changing the behavior for the smaller regression fixtures.
"""
function _pixel_hit_stack_mode(options::LightOptions)
    mode = lowercase(strip(options.pixel_hit_stack_mode))
    mode in ("auto", "small", "vector") || error(
        "Unsupported pixel_hit_stack_mode=$(repr(options.pixel_hit_stack_mode)); " *
        "supported values are \"auto\", \"small\", \"vector\".",
    )
    return mode
end

@inline function _pixel_hit_stack_type(options::LightOptions)
    mode = _pixel_hit_stack_mode(options)
    return mode == "vector" ? Vector{HitRecord} : SmallHitStack
end

@inline function _pixel_hit_stack_type(options::LightOptions, plotbox)
    mode = _pixel_hit_stack_mode(options)
    if mode == "auto"
        n_cells = plotbox.nx * plotbox.ny
        if n_cells >= _AUTO_VECTOR_PIXEL_HITS_MIN_CELLS && _use_dense_pixel_hits(plotbox)
            return Vector{HitRecord}
        end
        return SmallHitStack
    end
    return mode == "vector" ? Vector{HitRecord} : SmallHitStack
end

@inline _new_hit_stack(::Type{SmallHitStack}) = SmallHitStack()
@inline _new_hit_stack(::Type{Vector{HitRecord}}) = HitRecord[]

function _pixel_hits_table(::Type{S}, dense::Bool, plotbox) where {S}
    if dense
        return DensePixelHits(S, plotbox.nx * plotbox.ny)
    end
    return Dict{Int,S}()
end

function _pixel_hits_table(::Type{S}, dense::Bool, plotbox, upper_hit::Bool) where {S}
    if dense && upper_hit
        return DenseUpperPixelHits(plotbox.nx * plotbox.ny)
    end
    return _pixel_hits_table(S, dense, plotbox)
end

@inline function _append_hit!(pixel_hits, idx::Int, hit::HitRecord, stack_type::Type)
    h = get!(pixel_hits, idx) do
        _new_hit_stack(stack_type)
    end
    push!(h, hit)
    return nothing
end

@inline function _append_upper_hit!(pixel_hits, idx::Int, hit::HitRecord, stack_type::Type)
    h = get!(pixel_hits, idx) do
        _new_hit_stack(stack_type)
    end
    if isempty(h)
        push!(h, hit)
    elseif hit[1] > _stack_hit_height(h, 1)
        h[1] = hit
    end
    return nothing
end

@inline function _append_upper_hit!(pixel_hits::DenseUpperPixelHits, idx::Int, hit::HitRecord, ::Type)
    @inbounds begin
        if pixel_hits.nodes[idx] == 0
            pixel_hits.heights[idx] = hit[1]
            pixel_hits.nodes[idx] = hit[2]
            push!(pixel_hits.occupied, idx)
        elseif hit[1] > pixel_hits.heights[idx]
            pixel_hits.heights[idx] = hit[1]
            pixel_hits.nodes[idx] = hit[2]
        end
    end
    return nothing
end

@inline function _append_hit!(pixel_hits::FlatPixelHitBuilder, idx::Int, hit::HitRecord, ::Type)
    if pixel_hits.counts[idx] == 0
        push!(pixel_hits.occupied, idx)
    end
    push!(pixel_hits.pixels, FlatPoolInt(idx))
    push!(pixel_hits.heights, Float32(hit[1]))
    push!(pixel_hits.nodes, FlatPoolInt(hit[2]))
    pixel_hits.counts[idx] += 1
    pixel_hits.total_hits += 1
    return nothing
end

@inline function _append_upper_hit!(pixel_hits::FlatPixelHitBuilder, idx::Int, hit::HitRecord, ::Type)
    error("FlatPixelHitBuilder does not support upper-hit mode.")
end

function _finalize_flat_pixel_hits(builder::FlatPixelHitBuilder)
    n_pixels = length(builder.counts)
    total_hits = builder.total_hits
    starts = zeros(Int, n_pixels)
    heights = Vector{Float32}(undef, total_hits)
    nodes = Vector{FlatPoolInt}(undef, total_hits)

    cursor = 1
    counts = builder.counts
    occupied = builder.occupied
    @inbounds for idx in occupied
        count = counts[idx]
        starts[idx] = cursor
        counts[idx] = cursor
        cursor += count
    end

    pixels = builder.pixels
    @inbounds for pos in eachindex(pixels)
        idx = Int(pixels[pos])
        out = counts[idx]
        heights[out] = builder.heights[pos]
        nodes[out] = builder.nodes[pos]
        counts[idx] = out + 1
    end
    @inbounds for idx in occupied
        counts[idx] -= starts[idx]
    end
    return FlatPixelHits(starts, counts, occupied, falses(n_pixels), heights, nodes)
end

function _dense_float_node_map(node_ids::Vector{Int}, values::Vector{Float64})
    out = Dict{Int,Float64}()
    @inbounds for i in eachindex(node_ids)
        v = values[i]
        v == 0.0 && continue
        out[node_ids[i]] = v
    end
    return out
end

function _dense_int_node_map(node_ids::Vector{Int}, values::Vector{Int})
    out = Dict{Int,Int}()
    @inbounds for i in eachindex(node_ids)
        v = values[i]
        v == 0 && continue
        out[node_ids[i]] = v
    end
    return out
end

function _all_dense_float_node_map(node_ids::Vector{Int}, values::Vector{Float64})
    out = Dict{Int,Float64}()
    sizehint!(out, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[node_ids[i]] = values[i]
    end
    return out
end

function _all_dense_int_node_map(node_ids::Vector{Int}, values::Vector{Int})
    out = Dict{Int,Int}()
    sizehint!(out, length(node_ids))
    @inbounds for i in eachindex(node_ids)
        out[node_ids[i]] = values[i]
    end
    return out
end

function _first_order_result_from_dense(
    node_ids::Vector{Int},
    projected_area_per_node::Vector{Float64},
    incident_power_par::Vector{Float64},
    incident_power_nir::Vector{Float64},
    hits_per_node::Vector{Int},
    materialize_public::Bool=true,
)
    projected_area_map =
        materialize_public ? _all_dense_float_node_map(node_ids, projected_area_per_node) : Dict{Int,Float64}()
    incident_par_map =
        materialize_public ? _all_dense_float_node_map(node_ids, incident_power_par) : Dict{Int,Float64}()
    incident_nir_map =
        materialize_public ? _all_dense_float_node_map(node_ids, incident_power_nir) : Dict{Int,Float64}()
    hits_map =
        materialize_public ? _all_dense_int_node_map(node_ids, hits_per_node) : Dict{Int,Int}()
    return FirstOrderResult(
        projected_area_map,
        SpectralNodeValues(
            incident_par_map,
            incident_nir_map,
        ),
        hits_map,
        DenseFirstOrderResult(
            node_ids,
            projected_area_per_node,
            DenseSpectralNodeValues(incident_power_par, incident_power_nir),
            hits_per_node,
        ),
    )
end

function _materialize_first_order_result(first::FirstOrderResult)
    dense = first.dense
    dense === nothing && return first
    expected = length(dense.node_ids)
    if length(first.projected_area_per_node) == expected &&
       length(first.incident_power.par) == expected &&
       length(first.incident_power.nir) == expected &&
       length(first.hits_per_node) == expected
        return first
    end
    return _first_order_result_from_dense(
        dense.node_ids,
        dense.projected_area_per_node,
        dense.incident_power.par,
        dense.incident_power.nir,
        dense.hits_per_node,
        true,
    )
end

function _as_bool_local(x, default::Bool)
    x === nothing && return default
    x isa Bool && return x
    x isa Number && return x != 0
    x isa String && return lowercase(strip(x)) in ("1", "true", "yes", "y", "on")
    return default
end

function _normalize_group_name_local(x)
    s = strip(string(x))
    isempty(s) && return s
    s = replace(s, r"^#\[[^\]]+\]\s*" => "")
    strip(s)
end

function _cfg_toricity(options::LightOptions)
    options.toricity
end

function _cfg_debug_drop_leading_hit(options::LightOptions)
    options.debug_drop_leading_hit
end

function _apply_debug_drop_leading_hit!(pixel_hits, node_hits, projected_pixels_area, plotbox, options::LightOptions)
    spec = _cfg_debug_drop_leading_hit(options)
    spec === nothing && return
    idx = spec.x + 1 + spec.y * plotbox.nx
    stack = get(pixel_hits, idx, nothing)
    stack === nothing && return
    isempty(stack) && return

    _sort_hit_stack!(stack)
    top_nid = _stack_hit_node(stack, 1)
    top_nid == spec.node_id || return

    deleteat!(stack, 1)
    if isempty(stack)
        delete!(pixel_hits, idx)
    end

    node_hits[top_nid] = max(0, get(node_hits, top_nid, 0) - 1)
    projected_pixels_area[top_nid] = max(0.0, get(projected_pixels_area, top_nid, 0.0) - plotbox.pixel_area)
end

function _apply_debug_drop_leading_hit!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_pixels_area::Vector{Float64},
    plotbox,
    options::LightOptions,
    node_ids::Vector{Int},
)
    spec = _cfg_debug_drop_leading_hit(options)
    spec === nothing && return
    idx = spec.x + 1 + spec.y * plotbox.nx
    stack = get(pixel_hits, idx, nothing)
    stack === nothing && return
    isempty(stack) && return

    _sort_hit_stack!(stack)
    top_idx = _stack_hit_node(stack, 1)
    top_nid = node_ids[top_idx]
    top_nid == spec.node_id || return

    deleteat!(stack, 1)
    if isempty(stack)
        delete!(pixel_hits, idx)
    end

    node_hits[top_idx] = max(0, node_hits[top_idx] - 1)
    projected_pixels_area[top_idx] = max(0.0, projected_pixels_area[top_idx] - plotbox.pixel_area)
end

function _is_sensor_interception(interception::InterceptionModel)
    interception.sensor || return lowercase(strip(interception.model)) == "virtualsensor"
    return true
end

function _virtual_sensor_node_ids(node_group::Dict{Int,String}, node_type::Dict{Int,String}, models::LightModels)
    out = Set{Int}()
    for (nid, g) in node_group
        type_model = _type_model(models, _normalize_group_name_local(g), strip(get(node_type, nid, "")))
        type_model === nothing && continue
        interception = type_model.interception
        interception === nothing && continue
        _is_sensor_interception(interception) || continue
        push!(out, nid)
    end
    out
end

function _virtual_sensor_node_ids(geometry::InterceptionSceneData, models::LightModels)
    out = Set{Int}()
    @inbounds for i in eachindex(geometry.node_ids)
        g = geometry.node_group_by_index[i]
        type_model = _type_model(models, _normalize_group_name_local(g), strip(geometry.node_type_by_index[i]))
        type_model === nothing && continue
        interception = type_model.interception
        interception === nothing && continue
        _is_sensor_interception(interception) || continue
        push!(out, geometry.node_ids[i])
    end
    return out
end

function _virtual_sensor_node_ids(geometry::InterceptionSceneData, node_interception_by_index::Vector{Union{Nothing,InterceptionModel}})
    out = Set{Int}()
    @inbounds for i in eachindex(geometry.node_ids)
        interception = node_interception_by_index[i]
        interception === nothing && continue
        _is_sensor_interception(interception) || continue
        push!(out, geometry.node_ids[i])
    end
    return out
end

function _has_sensor_models(models::LightModels)
    for group_model in values(models)
        for type_model in values(group_model.types)
            interception = type_model.interception
            interception === nothing && continue
            _is_sensor_interception(interception) && return true
        end
    end
    return false
end

function _ignored_group_types(models::LightModels)
    out = Dict{String,Set{String}}()
    for group_model in values(models)
        group = _normalize_group_name_local(group_model.group)
        isempty(group) && continue
        for (type_name, type_model) in group_model.types
            interception = type_model.interception
            interception === nothing && continue
            lowercase(strip(interception.model)) == "ignore" || continue
            tname = strip(type_name)
            isempty(tname) && continue
            s = get!(out, group) do
                Set{String}()
            end
            push!(s, tname)
        end
    end
    out
end

function _is_ignored_node(node_id::Int, scene::PlantGeom.SceneGeometry, ignored::Dict{String,Set{String}})
    isempty(ignored) && return false
    g = _normalize_group_name_local(_scene_group(scene, node_id, ""))
    t = strip(_scene_type(scene, node_id, ""))
    return haskey(ignored, g) && (t in ignored[g])
end

function _group_light_emitters(models::LightModels)
    out = Dict{Tuple{String,String},NamedTuple{(:par, :nir),Tuple{Float64,Float64}}}()
    for group_model in values(models)
        group = strip(group_model.group)
        isempty(group) && continue
        for (type_name0, type_model) in group_model.types
            emitter = type_model.light_emitter
            emitter === nothing && continue
            type_name = strip(type_name0)
            lowercase(strip(emitter.model)) == "lambertianemitter" || continue

            radiance = emitter.radiance
            radiance > 0.0 || continue

            gpar = emitter.gamma.par
            gnir = emitter.gamma.nir
            gpar = max(gpar, 0.0)
            gnir = max(gnir, 0.0)
            gsum = gpar + gnir
            if gsum > 0.0
                gpar /= gsum
                gnir /= gsum
            else
                gpar = 0.48
                gnir = 0.52
            end

            key = (group, type_name)
            cur = get(out, key, (par=0.0, nir=0.0))
            out[key] = (par=cur.par + radiance * gpar, nir=cur.nir + radiance * gnir)
        end
    end
    out
end

function _resolved_type_key(group_model::GroupModel, type_name::AbstractString)
    key = String(type_name)
    haskey(group_model.types, key) && return key
    haskey(group_model.types, "*") && return "*"

    stripped = lowercase(strip(key))
    if length(group_model.types) == 1 && (isempty(stripped) || stripped == "mesh")
        return first(keys(group_model.types))
    end

    if stripped == "mesh"
        first_key = nothing
        first_sig = nothing
        equivalent = true
        for (type_key, type_model) in group_model.types
            sig = (
                interception=
                if type_model.interception === nothing
                    nothing
                else
                    interception = type_model.interception
                    (
                        interception.use,
                        interception.model,
                        interception.transparency,
                        interception.optical_properties === nothing ? nothing :
                        (
                            interception.optical_properties.par,
                            interception.optical_properties.nir,
                        ),
                    )
                end,
                emitter=
                if type_model.light_emitter === nothing
                    nothing
                else
                    emitter = type_model.light_emitter
                    (
                        emitter.model,
                        emitter.radiance,
                        emitter.gamma.par,
                        emitter.gamma.nir,
                    )
                end,
            )
            if first_sig === nothing
                first_key = type_key
                first_sig = sig
            elseif sig != first_sig
                equivalent = false
                break
            end
        end
        equivalent && first_key !== nothing && return first_key
    end

    return nothing
end

function _type_model(models::LightModels, group::AbstractString, type_name::AbstractString)
    for group_key in (String(group), "*")
        haskey(models, group_key) || continue
        group_model = models[group_key]
        resolved = _resolved_type_key(group_model, type_name)
        resolved === nothing || return group_model.types[resolved]
    end
    return nothing
end

function _resolved_node_interception_models(geometry::InterceptionSceneData, models::LightModels)
    out = Vector{Union{Nothing,InterceptionModel}}(undef, length(geometry.node_ids))
    @inbounds for i in eachindex(geometry.node_ids)
        group = geometry.node_group_by_index[i]
        type_name0 = geometry.node_type_by_index[i]
        type_name = isempty(type_name0) && group == "pavement" ? "Cobblestone" : type_name0
        type_model = _type_model(models, _normalize_group_name_local(group), strip(type_name))
        out[i] = type_model === nothing ? nothing : type_model.interception
    end
    return out
end

function _validate_scene_models(scene::PlantGeom.SceneGeometry, face2node::Vector{Int}, models::LightModels, ignored::Dict{String,Set{String}})
    missing = Set{Tuple{String,String}}()
    for nid in unique(face2node)
        _is_ignored_node(nid, scene, ignored) && continue
        group = strip(_scene_group(scene, nid, ""))
        type_name = strip(_scene_type(scene, nid, ""))
        _type_model(models, group, type_name) === nothing && push!(missing, (group, type_name))
    end
    isempty(missing) && return nothing
    details = join(["($(repr(g)), $(repr(t)))" for (g, t) in sort!(collect(missing))], ", ")
    error("Missing models for simulated geometry nodes: $details")
end

function _use_upper_hit_pixel_table(models::LightModels, options::LightOptions)
    # First-order interception uses the upper hit unless downstream transfer
    # logic requires the full per-pixel stack.
    if options.scattering
        return false
    end
    _has_sensor_models(models) && return false
    isempty(_group_light_emitters(models)) || return false
    return true
end

function _emitter_power_per_node(scene::PlantGeom.SceneGeometry, models::LightModels)
    by_group_type = _group_light_emitters(models)
    isempty(by_group_type) && return Dict{Int,Float64}(), Dict{Int,Float64}()

    par = Dict{Int,Float64}()
    nir = Dict{Int,Float64}()
    for ((group, type_name), pwr) in by_group_type
        nids = Int[
            nid for nid in keys(scene.nodes) if _scene_group(scene, nid, "") == group && _scene_type(scene, nid, "") == type_name
        ]
        if isempty(nids)
            # Fallback for scenes where type labels are unavailable.
            nids = Int[nid for nid in keys(scene.nodes) if _scene_group(scene, nid, "") == group]
        end
        isempty(nids) && continue

        atot = sum(_scene_area(scene, nid, 0.0) for nid in nids)
        if atot > 0.0
            for nid in nids
                w = _scene_area(scene, nid, 0.0) / atot
                par[nid] = get(par, nid, 0.0) + pwr.par * w
                nir[nid] = get(nir, nid, 0.0) + pwr.nir * w
            end
        else
            w = 1.0 / length(nids)
            for nid in nids
                par[nid] = get(par, nid, 0.0) + pwr.par * w
                nir[nid] = get(nir, nid, 0.0) + pwr.nir * w
            end
        end
    end
    return par, nir
end

@inline function _pack_emitter_edge(to::Int, from::Int)
    return (UInt64(UInt32(to)) << 32) | UInt64(UInt32(from))
end

@inline _unpack_emitter_to(edge::UInt64) = Int(UInt32(edge >> 32))
@inline _unpack_emitter_from(edge::UInt64) = Int(UInt32(edge & 0xffffffff))

function _emitter_weights_from_packed_counts(edge_counts::Dict{UInt64,Int}, total_from::Dict{Int,Int})
    weights = Dict{Tuple{Int,Int},Float64}()
    for (edge, count) in edge_counts
        to = _unpack_emitter_to(edge)
        src = _unpack_emitter_from(edge)
        n = get(total_from, src, 0)
        n > 0 || continue
        weights[(to, src)] = count / n
    end
    return weights
end

function _accumulate_emitter_transfer_counts!(
    edge_counts::Dict{UInt64,Int},
    total_from::Dict{Int,Int},
    projection::DirectionProjectionResult,
    emitter_nodes::Set{Int};
    stacks_sorted::Bool=false,
)
    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)

        for j in eachindex(stack)
            src = _stack_hit_node(stack, j)
            src in emitter_nodes || continue

            to = 0
            for k in (j+1):length(stack)
                nid = _stack_hit_node(stack, k)
                nid in emitter_nodes && continue
                to = nid
                break
            end
            to == 0 && continue

            edge = _pack_emitter_edge(to, src)
            edge_counts[edge] = get(edge_counts, edge, 0) + 1
            total_from[src] = get(total_from, src, 0) + 1
        end
    end
    return nothing
end

function _accumulate_emitter_transfer_counts!(
    edge_counts::Dict{UInt64,Int},
    total_from::Dict{Int,Int},
    projection::DenseDirectionProjectionResult,
    emitter_node_mask::Vector{Bool},
    node_ids::Vector{Int};
    stacks_sorted::Bool=false,
)
    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _accumulate_emitter_transfer_counts_dense!(edge_counts, total_from, stack, emitter_node_mask, node_ids)
    end
    return nothing
end

function _emitter_transfer_weights(
    vertices,
    faces,
    face2node,
    turtle::TurtleGrid,
    options::LightOptions,
    plotbox,
    emitter_nodes::Set{Int},
    cache_ctx,
)
    isempty(emitter_nodes) && return Dict{Tuple{Int,Int},Float64}()

    edge_counts = Dict{UInt64,Int}()
    total_from = Dict{Int,Int}()

    for sector in turtle.sectors
        sector.source == :sun && continue
        projection =
            _direction_projection_cached(vertices, faces, face2node, sector.direction, options, plotbox, cache_ctx, upper_hit=false)
        _accumulate_emitter_transfer_counts!(
            edge_counts,
            total_from,
            projection,
            emitter_nodes;
            stacks_sorted=false,
        )
    end

    return _emitter_weights_from_packed_counts(edge_counts, total_from)
end

function _emitter_transfer_weights_from_projections(
    projections::AbstractVector{DirectionProjectionResult},
    turtle::TurtleGrid,
    emitter_nodes::Set{Int},
    stacks_sorted::Bool=false,
)
    isempty(emitter_nodes) && return Dict{Tuple{Int,Int},Float64}()

    edge_counts = Dict{UInt64,Int}()
    total_from = Dict{Int,Int}()

    for i in eachindex(turtle.sectors)
        turtle.sectors[i].source == :sun && continue
        projection = projections[i]
        _accumulate_emitter_transfer_counts!(
            edge_counts,
            total_from,
            projection,
            emitter_nodes;
            stacks_sorted=stacks_sorted,
        )
    end

    return _emitter_weights_from_packed_counts(edge_counts, total_from)
end

function _emitter_transfer_weights_from_raycore(
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    options::LightOptions,
)
    prepared = data.prepared
    isempty(prepared.emitter_nodes) && return Dict{Tuple{Int,Int},Float64}()

    edge_counts = Dict{UInt64,Int}()
    total_from = Dict{Int,Int}()
    geometry = prepared.geometry

    for sector in turtle.sectors
        sector.source == :sun && continue
        projection = _raycore_direction_projection_full_stack(data, sector.direction, options)
        _accumulate_emitter_transfer_counts!(
            edge_counts,
            total_from,
            projection,
            prepared.emitter_node_mask,
            geometry.node_ids;
            stacks_sorted=false,
        )
    end

    return _emitter_weights_from_packed_counts(edge_counts, total_from)
end

function _cfg_cache_pixel_table(options::LightOptions)
    options.cache_pixel_table
end

function _projection_cache_dir(options::LightOptions)
    options.cache_pixel_table || return joinpath(tempdir(), "archimedlight-pixel_tables_cache")
    joinpath(tempdir(), "archimedlight-pixel_tables_cache")
end

const _FNV64_OFFSET_BASIS = UInt64(0xcbf29ce484222325)
const _FNV64_PRIME = UInt64(0x00000100000001b3)

@inline function _stable_mix_u64(h::UInt64, x::UInt64)
    return (h ⊻ x) * _FNV64_PRIME
end

@inline function _stable_mix_i64(h::UInt64, x::Int)
    return _stable_mix_u64(h, reinterpret(UInt64, Int64(x)))
end

@inline function _stable_mix_f64(h::UInt64, x::Float64)
    y = ifelse(x == 0.0, 0.0, x)
    return _stable_mix_u64(h, reinterpret(UInt64, y))
end

function _projection_scene_key(vertices, faces, face2node, plotbox, options::LightOptions)
    h = _FNV64_OFFSET_BASIS
    h = _stable_mix_i64(h, length(vertices))
    h = _stable_mix_i64(h, length(faces))
    h = _stable_mix_i64(h, length(face2node))
    h = _stable_mix_u64(h, _cfg_toricity(options) ? UInt64(1) : UInt64(0))
    h = _stable_mix_f64(h, plotbox.origin_x)
    h = _stable_mix_f64(h, plotbox.origin_y)
    h = _stable_mix_f64(h, plotbox.xdim)
    h = _stable_mix_f64(h, plotbox.ydim)
    h = _stable_mix_i64(h, plotbox.nx)
    h = _stable_mix_i64(h, plotbox.ny)
    h = _stable_mix_f64(h, plotbox.pix_x)
    h = _stable_mix_f64(h, plotbox.pix_y)
    h = _stable_mix_f64(h, plotbox.pixel_area)

    @inbounds for v in vertices
        h = _stable_mix_f64(h, Float64(v[1]))
        h = _stable_mix_f64(h, Float64(v[2]))
        h = _stable_mix_f64(h, Float64(v[3]))
    end

    @inbounds for i in eachindex(faces)
        f = faces[i]
        h = _stable_mix_i64(h, Int(f[1]))
        h = _stable_mix_i64(h, Int(f[2]))
        h = _stable_mix_i64(h, Int(f[3]))
        h = _stable_mix_i64(h, Int(face2node[i]))
    end

    h
end

function _projection_dir_key(direction)
    h = _FNV64_OFFSET_BASIS
    h = _stable_mix_f64(h, Float64(direction[1]))
    h = _stable_mix_f64(h, Float64(direction[2]))
    h = _stable_mix_f64(h, Float64(direction[3]))
    h
end

function _projection_cache_context(vertices, faces, face2node, plotbox, options::LightOptions)
    _cfg_cache_pixel_table(options) || return nothing
    return ProjectionCacheContext(
        _projection_cache_dir(options),
        _projection_scene_key(vertices, faces, face2node, plotbox, options),
    )
end

function _projection_cache_path(cache_ctx::ProjectionCacheContext, direction, upper_hit::Bool=false, strict_java_float::Bool=false)
    scene_hex = string(cache_ctx.scene_key, base=16, pad=16)
    dir_hex = string(_projection_dir_key(direction), base=16, pad=16)
    mode = (upper_hit ? "u1" : "u0") * (strict_java_float ? "_sf1" : "_sf0")
    joinpath(cache_ctx.cache_dir, "proj_" * scene_hex * "_" * dir_hex * "_" * mode * ".jls")
end

function _read_projection_cache(path::AbstractString)
    try
        open(path, "r") do io
            payload = Serialization.deserialize(io)
            payload isa NamedTuple || return nothing
            (:version in propertynames(payload) && getproperty(payload, :version) == 1) || return nothing
            req = (:pixel_hits, :node_hits, :projected_mesh_area, :projected_pixels_area)
            all(n -> n in propertynames(payload), req) || return nothing
            return DirectionProjectionResult(
                getproperty(payload, :pixel_hits),
                getproperty(payload, :node_hits),
                getproperty(payload, :projected_mesh_area),
                getproperty(payload, :projected_pixels_area),
            )
        end
    catch
        return nothing
    end
end

function _write_projection_cache(path::AbstractString, result::DirectionProjectionResult)
    mkpath(dirname(path))
    tmp = path * ".tmp-" * string(getpid()) * "-" * string(time_ns())
    payload = (
        version=1,
        pixel_hits=result.pixel_hits,
        node_hits=result.node_hits,
        projected_mesh_area=result.projected_mesh_area,
        projected_pixels_area=result.projected_pixels_area,
    )
    open(tmp, "w") do io
        Serialization.serialize(io, payload)
    end
    mv(tmp, path; force=true)
end

function _extract_scene_xy_bounds(scene::PlantGeom.SceneGeometry, vertices)
    if scene.scene_xy_bounds !== nothing
        return scene.scene_xy_bounds
    end

    attrs = MultiScaleTreeGraph.attributes(scene.mtg)
    if hasproperty(attrs, :scene_dimensions)
        dims = getproperty(attrs, :scene_dimensions)
        p0 = dims[1]
        p1 = dims[2]
        x0 = Float64(p0[1])
        y0 = Float64(p0[2])
        x1 = Float64(p1[1])
        y1 = Float64(p1[2])
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
    end

    xs = Float64[p[1] for p in vertices]
    ys = Float64[p[2] for p in vertices]
    (minimum(xs), minimum(ys), maximum(xs), maximum(ys))
end

function _plotbox(scene::PlantGeom.SceneGeometry, vertices, pixel_size::Float64)
    pixel_size > 0.0 || error("pixel_size must be > 0 m")
    # Java fixtures enforce a strict upper bound at 50 cm.
    pixel_size <= 0.5 || error("pixel_size must be <= 0.5 m (50 cm)")

    # Match Java PlotBox numeric path:
    # - scene corners are Point2f (Float32)
    # - table size uses double division from float-backed plot dimensions
    # - pixel dimensions are stored as float
    x0_raw, y0_raw, x1_raw, y1_raw = _extract_scene_xy_bounds(scene, vertices)
    x0 = Float32(min(x0_raw, x1_raw))
    y0 = Float32(min(y0_raw, y1_raw))
    x1 = Float32(max(x0_raw, x1_raw))
    y1 = Float32(max(y0_raw, y1_raw))

    xdim = max(Float64(x1 - x0), pixel_size)
    ydim = max(Float64(y1 - y0), pixel_size)

    nx = max(1, floor(Int, xdim / pixel_size))
    ny = max(1, floor(Int, ydim / pixel_size))

    pix_x = Float64(Float32(xdim / nx))
    pix_y = Float64(Float32(ydim / ny))
    plot_area = Float64(Float32(Float32(xdim) * Float32(ydim)))
    pixel_area = plot_area / (nx * ny)

    (
        origin_x=Float64(x0),
        origin_y=Float64(y0),
        xdim=xdim,
        ydim=ydim,
        nx=nx,
        ny=ny,
        pix_x=pix_x,
        pix_y=pix_y,
        pixel_area=pixel_area,
    )
end

function _project_point_ground(point, direction)
    dzden = Float32(direction[3])
    dzden == 0.0f0 && return nothing
    pz = Float32(point[3])
    dz = -pz / dzden
    x = Float32(point[1]) + Float32(direction[1]) * dz
    y = Float32(point[2]) + Float32(direction[2]) * dz
    # Keep source altitude for stack sorting as in Java.
    StaticArrays.SVector{3,Float64}(Float64(x), Float64(y), Float64(pz))
end

function _triangle_area_xy(p1, p2, p3)
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2])) * 0.5
end

function _compute_normal(points)
    n = length(points)
    n < 3 && return StaticArrays.SVector{3,Float64}(0.0, 0.0, 0.0)
    for k in 1:(n-2)
        p0 = points[k]
        p1 = points[k+1]
        p2 = points[k+2]

        v1x = Float32(p1[1] - p0[1])
        v1y = Float32(p1[2] - p0[2])
        v1z = Float32(p1[3] - p0[3])
        n1 = sqrt((v1x * v1x) + (v1y * v1y) + (v1z * v1z))
        n1 <= 0.0f0 && continue
        v1x /= n1
        v1y /= n1
        v1z /= n1

        v2x = Float32(p2[1] - p1[1])
        v2y = Float32(p2[2] - p1[2])
        v2z = Float32(p2[3] - p1[3])
        n2 = sqrt((v2x * v2x) + (v2y * v2y) + (v2z * v2z))
        n2 <= 0.0f0 && continue
        v2x /= n2
        v2y /= n2
        v2z /= n2

        nx = (v1y * v2z) - (v1z * v2y)
        ny = (v1z * v2x) - (v1x * v2z)
        nz = (v1x * v2y) - (v1y * v2x)
        nnorm = sqrt((nx * nx) + (ny * ny) + (nz * nz))
        nnorm <= 0.0f0 && continue
        return StaticArrays.SVector{3,Float64}(Float64(nx / nnorm), Float64(ny / nnorm), Float64(nz / nnorm))
    end
    StaticArrays.SVector{3,Float64}(0.0, 0.0, 0.0)
end

function _compute_normal_f32(points)
    length(points) < 2 && return StaticArrays.SVector{3,Float32}(0.0f0, 0.0f0, 0.0f0)
    @inbounds for n in 1:(length(points)-2)
        p0 = points[n]
        p1 = points[n+1]
        p2 = points[n+2]

        v1x = Float32(p1[1] - p0[1])
        v1y = Float32(p1[2] - p0[2])
        v1z = Float32(p1[3] - p0[3])
        n1 = sqrt((v1x * v1x) + (v1y * v1y) + (v1z * v1z))
        n1 <= 0.0f0 && continue
        v1x /= n1
        v1y /= n1
        v1z /= n1

        v2x = Float32(p2[1] - p1[1])
        v2y = Float32(p2[2] - p1[2])
        v2z = Float32(p2[3] - p1[3])
        n2 = sqrt((v2x * v2x) + (v2y * v2y) + (v2z * v2z))
        n2 <= 0.0f0 && continue
        v2x /= n2
        v2y /= n2
        v2z /= n2

        nx = (v1y * v2z) - (v1z * v2y)
        ny = (v1z * v2x) - (v1x * v2z)
        nz = (v1x * v2y) - (v1y * v2x)
        nnorm = sqrt((nx * nx) + (ny * ny) + (nz * nz))
        nnorm <= 0.0f0 && continue
        return StaticArrays.SVector{3,Float32}(nx / nnorm, ny / nnorm, nz / nnorm)
    end
    StaticArrays.SVector{3,Float32}(0.0f0, 0.0f0, 0.0f0)
end

@inline function _project_vertex_to_ground_pixel(v, dirx::Float32, diry::Float32, dirz::Float32, ox::Float32, oy::Float32, pxs::Float32, pys::Float32, u::Float32)
    pz = Float32(v[3]) * u
    dz = -pz / dirz
    xw = Float32(v[1]) * u + dirx * dz
    yw = Float32(v[2]) * u + diry * dz
    projected = StaticArrays.SVector{3,Float64}(Float64(xw / u), Float64(yw / u), 0.0)
    pix = StaticArrays.SVector{3,Float64}(Float64((xw - ox) / pxs), Float64((yw - oy) / pys), Float64(pz))
    return projected, pix
end

function _get_border_pixels(p1, p2, i_origin::Int, minY::Vector{Int}, maxY::Vector{Int})
    p_min, p_max = p1[1] < p2[1] ? (p1, p2) : (p2, p1)
    dx = Float32(p_max[1] - p_min[1])
    dx < 1e-6 && return

    slope = Float32((Float32(p_max[2]) - Float32(p_min[2])) / dx)
    i_min = ceil(Int, Float32(p_min[1]))
    i_max = floor(Int, Float32(p_max[1]))

    @inbounds for i in i_min:i_max
        i0 = i - i_origin
        yi = slope * Float32(i - Float32(p_min[1])) + Float32(p_min[2])
        j = trunc(Int, floor(Float64(yi)) + 0.5)
        if 0 <= i0 < length(minY)
            idx = i0 + 1
            minY[idx] = min(minY[idx], j)
            maxY[idx] = max(maxY[idx], j)
        end
    end
end

function _wrap_index(i::Int, n::Int)
    ii = i
    while ii < 0
        ii += n
    end
    while ii >= n
        ii -= n
    end
    ii
end

function _project_triangle!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64},
    stack_node::Int,
    node_idx::Int,
    p1,
    p2,
    p3,
    direction,
    origin_x::Float64,
    origin_y::Float64,
    x_pix::Float64,
    y_pix::Float64,
    pixel_area::Float64,
    nx::Int,
    ny::Int,
    toricity::Bool,
    upper_hit::Bool,
    strict_java_float::Bool,
    stack_type::Type,
)
    _project_triangle!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        stack_node,
        node_idx,
        p1,
        p2,
        p3,
        direction,
        origin_x,
        origin_y,
        x_pix,
        y_pix,
        pixel_area,
        nx,
        ny,
        toricity,
        upper_hit,
        strict_java_float,
        1.0f0,
        stack_type,
    )
end

function _project_triangle!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64},
    stack_node::Int,
    node_idx::Int,
    p1,
    p2,
    p3,
    direction,
    origin_x::Float64,
    origin_y::Float64,
    x_pix::Float64,
    y_pix::Float64,
    pixel_area::Float64,
    nx::Int,
    ny::Int,
    toricity::Bool,
    upper_hit::Bool,
    strict_java_float::Bool,
    unit_scale::Float32,
    stack_type::Type,
)
    u = unit_scale > 0.0f0 ? unit_scale : 1.0f0
    ox = Float32(origin_x) * u
    oy = Float32(origin_y) * u
    pxs = Float32(x_pix) * u
    pys = Float32(y_pix) * u
    dirx = Float32(direction[1])
    diry = Float32(direction[2])
    dirz = Float32(direction[3])
    dirz == 0.0f0 && return

    projected1, pix1 = _project_vertex_to_ground_pixel(p1, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected2, pix2 = _project_vertex_to_ground_pixel(p2, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected3, pix3 = _project_vertex_to_ground_pixel(p3, dirx, diry, dirz, ox, oy, pxs, pys, u)

    iMin = min(floor(Int, pix1[1]), floor(Int, pix2[1]), floor(Int, pix3[1]))
    iMax = max(ceil(Int, pix1[1]), ceil(Int, pix2[1]), ceil(Int, pix3[1]))
    jMin = min(floor(Int, pix1[2]), floor(Int, pix2[2]), floor(Int, pix3[2]))
    jMax = max(ceil(Int, pix1[2]), ceil(Int, pix2[2]), ceil(Int, pix3[2]))
    kMin = min(floor(Int, pix1[3]), floor(Int, pix2[3]), floor(Int, pix3[3]))
    kMax = max(ceil(Int, pix1[3]), ceil(Int, pix2[3]), ceil(Int, pix3[3]))

    iLength = (iMax - iMin) + 1
    iLength <= 0 && return

    minY = fill(jMax, iLength)
    maxY = fill(jMin, iLength)

    _get_border_pixels(pix1, pix2, iMin, minY, maxY)
    _get_border_pixels(pix2, pix3, iMin, minY, maxY)
    _get_border_pixels(pix3, pix1, iMin, minY, maxY)

    normal = strict_java_float ? _compute_normal_f32((pix1, pix2, pix3)) : _compute_normal((pix1, pix2, pix3))
    slopeX_f32, slopeY_f32 =
        if abs(normal[3]) > 1e-5
            (Float32(normal[1] / normal[3]), Float32(normal[2] / normal[3]))
        else
            # Java fallback uses signed direction.z when normal.z is almost zero.
            dz = Float32(direction[3])
            (dz * Float32(normal[1]), dz * Float32(normal[2]))
        end

    z0 = Float32(pix1[3])
    z0 += slopeX_f32 * (Float32(pix1[1]) - Float32(iMin))
    z0 += slopeY_f32 * (Float32(pix1[2]) - Float32(jMin))

    tri_proj_area = _triangle_area_xy(projected1, projected2, projected3)
    nb_hits = 0
    node_hit_count = node_hits[node_idx]
    @inbounds for i in iMin:(iMax-1)
        ni = i - iMin
        zi = z0 - slopeX_f32 * Float32(ni)
        ymin_i = minY[ni+1]
        ymax_i = maxY[ni+1]

        for j in ymin_i:(ymax_i-1)
            nj = j - jMin
            zpix = zi - slopeY_f32 * Float32(nj)
            zpix = clamp(zpix, Float32(kMin), Float32(kMax))
            nb_hits += 1

            ii, jj =
                if toricity
                    (_wrap_index(i, nx), _wrap_index(j, ny))
                else
                    (i, j)
                end

            if toricity || ((0 <= ii < nx) && (0 <= jj < ny))
                idx = ii + 1 + jj * nx
                zpix_f32 = Float64(zpix / u)
                hit = (zpix_f32, stack_node)
                if upper_hit
                    _append_upper_hit!(pixel_hits, idx, hit, stack_type)
                else
                    _append_hit!(pixel_hits, idx, hit, stack_type)
                end
                node_hit_count += 1
            end
        end
    end

    node_hits[node_idx] = node_hit_count
    projected_mesh_area[node_idx] += tri_proj_area
    projected_pixels_area[node_idx] += nb_hits * pixel_area
end

@inline function _scanline_bounds!(scratch::RasterScanlineScratch, iLength::Int, jMin::Int, jMax::Int)
    resize!(scratch.minY, iLength)
    resize!(scratch.maxY, iLength)
    fill!(scratch.minY, jMax)
    fill!(scratch.maxY, jMin)
    return scratch.minY, scratch.maxY
end

function _project_triangle!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64},
    stack_node::Int,
    node_idx::Int,
    p1,
    p2,
    p3,
    direction,
    origin_x::Float64,
    origin_y::Float64,
    x_pix::Float64,
    y_pix::Float64,
    pixel_area::Float64,
    nx::Int,
    ny::Int,
    toricity::Bool,
    upper_hit::Bool,
    strict_java_float::Bool,
    unit_scale::Float32,
    stack_type::Type,
    scratch::RasterScanlineScratch,
)
    u = unit_scale > 0.0f0 ? unit_scale : 1.0f0
    ox = Float32(origin_x) * u
    oy = Float32(origin_y) * u
    pxs = Float32(x_pix) * u
    pys = Float32(y_pix) * u
    dirx = Float32(direction[1])
    diry = Float32(direction[2])
    dirz = Float32(direction[3])
    dirz == 0.0f0 && return

    projected1, pix1 = _project_vertex_to_ground_pixel(p1, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected2, pix2 = _project_vertex_to_ground_pixel(p2, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected3, pix3 = _project_vertex_to_ground_pixel(p3, dirx, diry, dirz, ox, oy, pxs, pys, u)

    iMin = min(floor(Int, pix1[1]), floor(Int, pix2[1]), floor(Int, pix3[1]))
    iMax = max(ceil(Int, pix1[1]), ceil(Int, pix2[1]), ceil(Int, pix3[1]))
    jMin = min(floor(Int, pix1[2]), floor(Int, pix2[2]), floor(Int, pix3[2]))
    jMax = max(ceil(Int, pix1[2]), ceil(Int, pix2[2]), ceil(Int, pix3[2]))
    kMin = min(floor(Int, pix1[3]), floor(Int, pix2[3]), floor(Int, pix3[3]))
    kMax = max(ceil(Int, pix1[3]), ceil(Int, pix2[3]), ceil(Int, pix3[3]))

    iLength = (iMax - iMin) + 1
    iLength <= 0 && return

    minY, maxY = _scanline_bounds!(scratch, iLength, jMin, jMax)

    _get_border_pixels(pix1, pix2, iMin, minY, maxY)
    _get_border_pixels(pix2, pix3, iMin, minY, maxY)
    _get_border_pixels(pix3, pix1, iMin, minY, maxY)

    normal = strict_java_float ? _compute_normal_f32((pix1, pix2, pix3)) : _compute_normal((pix1, pix2, pix3))
    slopeX_f32, slopeY_f32 =
        if abs(normal[3]) > 1e-5
            (Float32(normal[1] / normal[3]), Float32(normal[2] / normal[3]))
        else
            dz = Float32(direction[3])
            (dz * Float32(normal[1]), dz * Float32(normal[2]))
        end

    z0 = Float32(pix1[3])
    z0 += slopeX_f32 * (Float32(pix1[1]) - Float32(iMin))
    z0 += slopeY_f32 * (Float32(pix1[2]) - Float32(jMin))

    tri_proj_area = _triangle_area_xy(projected1, projected2, projected3)
    nb_hits = 0
    node_hit_count = node_hits[node_idx]
    @inbounds for i in iMin:(iMax-1)
        ni = i - iMin
        zi = z0 - slopeX_f32 * Float32(ni)
        ymin_i = minY[ni+1]
        ymax_i = maxY[ni+1]

        for j in ymin_i:(ymax_i-1)
            nj = j - jMin
            zpix = zi - slopeY_f32 * Float32(nj)
            zpix = clamp(zpix, Float32(kMin), Float32(kMax))
            nb_hits += 1

            ii, jj =
                if toricity
                    (_wrap_index(i, nx), _wrap_index(j, ny))
                else
                    (i, j)
                end

            if toricity || ((0 <= ii < nx) && (0 <= jj < ny))
                idx = ii + 1 + jj * nx
                zpix_f32 = Float64(zpix / u)
                hit = (zpix_f32, stack_node)
                if upper_hit
                    _append_upper_hit!(pixel_hits, idx, hit, stack_type)
                else
                    _append_hit!(pixel_hits, idx, hit, stack_type)
                end
                node_hit_count += 1
            end
        end
    end

    node_hits[node_idx] = node_hit_count
    projected_mesh_area[node_idx] += tri_proj_area
    projected_pixels_area[node_idx] += nb_hits * pixel_area
end

@inline _hit_height(hit) = hit[1]
@inline _hit_node(hit) = hit[2]

function _node_transparency_by_index(
    node_interception_by_index::Vector{Union{Nothing,InterceptionModel}},
    virtual_node_mask::Vector{Bool},
)
    out = zeros(Float64, length(node_interception_by_index))
    @inbounds for i in eachindex(node_interception_by_index)
        if virtual_node_mask[i]
            out[i] = 1.0
            continue
        end
        interception = node_interception_by_index[i]
        out[i] = interception === nothing ? 0.0 : clamp(interception.transparency, 0.0, 1.0)
    end
    return out
end

@inline function _projection_area_ratio(
    projection::DirectionProjectionResult,
    options::LightOptions,
    node_id::Int,
)
    options.area_ratio || return 1.0
    projected_pixels = get(projection.projected_pixels_area, node_id, 0.0)
    projected_pixels > 0.0 || return 1.0
    return get(projection.projected_mesh_area, node_id, 0.0) / projected_pixels
end

@inline function _projection_area_ratio(
    projection::DenseDirectionProjectionResult,
    options::LightOptions,
    node_idx::Int,
)
    options.area_ratio || return 1.0
    projected_pixels = projection.projected_pixels_area[node_idx]
    projected_pixels > 0.0 || return 1.0
    return projection.projected_mesh_area[node_idx] / projected_pixels
end

function _dense_projection_hits(projection::DirectionProjectionResult, geometry::InterceptionSceneData)
    return _dense_sector_int(projection.node_hits, geometry).values
end

function _dense_projection_hits(projection::DenseDirectionProjectionResult, geometry::InterceptionSceneData)
    return projection.node_hits
end

function _copy_projection_hits!(
    hits_per_node::Vector{Int},
    projection::DirectionProjectionResult,
    geometry::InterceptionSceneData,
)
    fill!(hits_per_node, 0)
    for (nid, h) in projection.node_hits
        hits_per_node[geometry.node_index[nid]] = h
    end
    return hits_per_node
end

function _copy_projection_hits!(
    hits_per_node::Vector{Int},
    projection::DenseDirectionProjectionResult,
    geometry::InterceptionSceneData,
)
    copyto!(hits_per_node, projection.node_hits)
    return hits_per_node
end

function _accumulate_projection_hits!(hits_per_node::Vector{Int}, projection::DirectionProjectionResult, geometry::InterceptionSceneData)
    for (nid, h) in projection.node_hits
        hits_per_node[geometry.node_index[nid]] += h
    end
    return nothing
end

function _accumulate_projection_hits!(
    hits_per_node::Vector{Int},
    projection::DenseDirectionProjectionResult,
    geometry::InterceptionSceneData,
)
    @inbounds for i in eachindex(hits_per_node)
        hits_per_node[i] += projection.node_hits[i]
    end
    return nothing
end

function _visible_area_from_projection(
    projection::DirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_nodes::Set{Int},
    stacks_sorted::Bool=false,
)
    visible_area = Dict{Int,Float64}()
    sizehint!(visible_area, length(projection.node_hits))
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)

        non_virtual_seen = false
        first_non_virtual = 0
        for hit in stack
            nid = _hit_node(hit)
            if nid in virtual_nodes
                if !non_virtual_seen
                    ratio = _projection_area_ratio(projection, options, nid)
                    visible_area[nid] = get(visible_area, nid, 0.0) + pixel_area * ratio
                end
            else
                first_non_virtual = nid
                non_virtual_seen = true
                break
            end
        end
        if first_non_virtual != 0
            ratio = _projection_area_ratio(projection, options, first_non_virtual)
            visible_area[first_non_virtual] = get(visible_area, first_non_virtual, 0.0) + pixel_area * ratio
        end
    end

    return visible_area
end

function _visible_area_from_projection_dense(
    projection::DirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_nodes::Set{Int},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    visible_area = zeros(Float64, length(geometry.node_ids))
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)

        non_virtual_seen = false
        first_non_virtual = 0
        for hit in stack
            nid = _hit_node(hit)
            if nid in virtual_nodes
                if !non_virtual_seen
                    ratio = _projection_area_ratio(projection, options, nid)
                    visible_area[geometry.node_index[nid]] += pixel_area * ratio
                end
            else
                first_non_virtual = nid
                non_virtual_seen = true
                break
            end
        end
        if first_non_virtual != 0
            ratio = _projection_area_ratio(projection, options, first_non_virtual)
            visible_area[geometry.node_index[first_non_virtual]] += pixel_area * ratio
        end
    end

    return visible_area
end

function _visible_area_from_projection_dense(
    projection::DirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    visible_area = zeros(Float64, length(geometry.node_ids))
    return _visible_area_from_projection_dense!(
        visible_area,
        projection,
        options,
        plotbox,
        virtual_node_mask,
        node_transparency_by_index,
        geometry,
        stacks_sorted,
    )
end

function _visible_area_from_projection_dense!(
    visible_area::Vector{Float64},
    projection::DirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    fill!(visible_area, 0.0)
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)

        for hit in stack
            nid = _hit_node(hit)
            node_idx = geometry.node_index[nid]
            if virtual_node_mask[node_idx]
                ratio = _projection_area_ratio(projection, options, nid)
                visible_area[node_idx] += pixel_area * ratio
                continue
            end

            transparency = node_transparency_by_index[node_idx]
            intercepted_fraction = 1.0 - transparency
            if intercepted_fraction > 0.0
                ratio = _projection_area_ratio(projection, options, nid)
                visible_area[node_idx] += pixel_area * intercepted_fraction * ratio
            end
            transparency > 0.0 || break
        end
    end

    return visible_area
end

function _visible_area_from_projection_dense(
    projection::DenseDirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_nodes::Set{Int},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    visible_area = zeros(Float64, length(geometry.node_ids))
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)

        non_virtual_seen = false
        first_non_virtual_idx = 0
        for hit in stack
            node_idx = _hit_node(hit)
            if geometry.node_ids[node_idx] in virtual_nodes
                if !non_virtual_seen
                    ratio = _projection_area_ratio(projection, options, node_idx)
                    visible_area[node_idx] += pixel_area * ratio
                end
            else
                first_non_virtual_idx = node_idx
                non_virtual_seen = true
                break
            end
        end
        if first_non_virtual_idx != 0
            ratio = _projection_area_ratio(projection, options, first_non_virtual_idx)
            visible_area[first_non_virtual_idx] += pixel_area * ratio
        end
    end

    return visible_area
end

function _visible_area_from_projection_dense(
    projection::DenseDirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    visible_area = zeros(Float64, length(geometry.node_ids))
    return _visible_area_from_projection_dense!(
        visible_area,
        projection,
        options,
        plotbox,
        virtual_node_mask,
        node_transparency_by_index,
        geometry,
        stacks_sorted,
    )
end

function _visible_area_from_projection_dense!(
    visible_area::Vector{Float64},
    projection::DenseDirectionProjectionResult,
    options::LightOptions,
    plotbox,
    virtual_node_mask::Vector{Bool},
    node_transparency_by_index::Vector{Float64},
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    fill!(visible_area, 0.0)
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _accumulate_visible_area_dense!(
            visible_area,
            projection,
            stack,
            options,
            pixel_area,
            virtual_node_mask,
            node_transparency_by_index,
        )
    end

    return visible_area
end

function _prepared_direction_projection(
    prepared::PreparedInterceptionData,
    direction,
    options::LightOptions,
)
    geometry = prepared.geometry
    strict_java_float = _strict_java_float(options)
    return _direction_projection_cached(geometry, direction, options, prepared.cache_ctx; upper_hit=prepared.upper_hit, strict_java_float=strict_java_float)
end

function _build_direction_projections(
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    options::LightOptions,
)
    n = length(turtle.sectors)
    n == 0 && return Any[]
    first_projection = _prepared_direction_projection(prepared, turtle.sectors[1].direction, options)
    projections = Vector{typeof(first_projection)}(undef, n)
    projections[1] = first_projection
    for i in 2:n
        projections[i] = _prepared_direction_projection(prepared, turtle.sectors[i].direction, options)
    end
    return projections
end

function _rasterize_direction_java(
    vertices,
    faces,
    face2node,
    direction,
    options::LightOptions,
    plotbox;
    cache_ctx=nothing,
    virtual_nodes=Set{Int}(),
    upper_hit::Bool=false,
)
    strict_java_float = _strict_java_float(options)
    projection =
        _direction_projection_cached(vertices, faces, face2node, direction, options, plotbox, cache_ctx; upper_hit=upper_hit, strict_java_float=strict_java_float)
    return _visible_area_from_projection(projection, options, plotbox, virtual_nodes), projection.node_hits
end

function _accumulate_direction_projection!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64},
    vertices,
    faces,
    stack_nodes::Vector{Int},
    face2node_index::Vector{Int},
    direction,
    plotbox,
    toricity::Bool,
    use_upper_hit::Bool,
    strict_java_float::Bool,
    unit_scale::Float32,
    stack_type::Type,
)
    origin_x = plotbox.origin_x
    origin_y = plotbox.origin_y
    pix_x = plotbox.pix_x
    pix_y = plotbox.pix_y
    pixel_area = plotbox.pixel_area
    nx = plotbox.nx
    ny = plotbox.ny
    @inbounds for fi in eachindex(faces)
        f = faces[fi]
        _project_triangle!(
            pixel_hits,
            node_hits,
            projected_mesh_area,
            projected_pixels_area,
            stack_nodes[fi],
            face2node_index[fi],
            vertices[f[1]],
            vertices[f[2]],
            vertices[f[3]],
            direction,
            origin_x,
            origin_y,
            pix_x,
            pix_y,
            pixel_area,
            nx,
            ny,
            toricity,
            use_upper_hit,
            strict_java_float,
            unit_scale,
            stack_type,
        )
    end
    return nothing
end

function _accumulate_direction_projection_prepared!(
    pixel_hits,
    node_hits::Vector{Int},
    projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64},
    vertices,
    faces,
    face2node_index::Vector{Int},
    direction,
    plotbox,
    toricity::Bool,
    use_upper_hit::Bool,
    strict_java_float::Bool,
    unit_scale::Float32,
    stack_type::Type,
)
    origin_x = plotbox.origin_x
    origin_y = plotbox.origin_y
    pix_x = plotbox.pix_x
    pix_y = plotbox.pix_y
    pixel_area = plotbox.pixel_area
    nx = plotbox.nx
    ny = plotbox.ny
    scratch = RasterScanlineScratch()
    @inbounds for fi in eachindex(faces)
        f = faces[fi]
        node_idx = face2node_index[fi]
        _project_triangle!(
            pixel_hits,
            node_hits,
            projected_mesh_area,
            projected_pixels_area,
            node_idx,
            node_idx,
            vertices[f[1]],
            vertices[f[2]],
            vertices[f[3]],
            direction,
            origin_x,
            origin_y,
            pix_x,
            pix_y,
            pixel_area,
            nx,
            ny,
            toricity,
            use_upper_hit,
            strict_java_float,
            unit_scale,
            stack_type,
            scratch,
        )
    end
    return nothing
end

function _direction_projection_materialized(
    vertices,
    faces,
    face2node,
    face2node_index::Vector{Int},
    node_ids::Vector{Int},
    direction,
    options::LightOptions,
    plotbox;
    upper_hit::Union{Nothing,Bool}=nothing,
    dense_pixel_hits::Bool=true,
)
    toricity = _cfg_toricity(options)
    use_upper_hit = upper_hit === nothing ? false : Bool(upper_hit)
    strict_java_float = _strict_java_float(options)
    unit_scale = Float32(_projection_unit_scale(options))
    stack_type = _pixel_hit_stack_type(options, plotbox)
    use_dense_table = dense_pixel_hits && _use_dense_pixel_hits(plotbox)
    pixel_hits = _pixel_hits_table(stack_type, use_dense_table, plotbox, use_upper_hit)
    node_hits = zeros(Int, length(node_ids))
    projected_mesh_area = zeros(Float64, length(node_ids))
    projected_pixels_area = zeros(Float64, length(node_ids))
    _accumulate_direction_projection!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        vertices,
        faces,
        face2node,
        face2node_index,
        direction,
        plotbox,
        toricity,
        use_upper_hit,
        strict_java_float,
        unit_scale,
        stack_type,
    )
    node_hits_map = _dense_int_node_map(node_ids, node_hits)
    projected_mesh_area_map = _dense_float_node_map(node_ids, projected_mesh_area)
    projected_pixels_area_map = _dense_float_node_map(node_ids, projected_pixels_area)
    _apply_debug_drop_leading_hit!(pixel_hits, node_hits_map, projected_pixels_area_map, plotbox, options)
    return DirectionProjectionResult(pixel_hits, node_hits_map, projected_mesh_area_map, projected_pixels_area_map)
end

function _direction_projection_prepared(
    vertices,
    faces,
    face2node_index::Vector{Int},
    node_ids::Vector{Int},
    direction,
    options::LightOptions,
    plotbox;
    upper_hit::Union{Nothing,Bool}=nothing,
    dense_pixel_hits::Bool=true,
)
    toricity = _cfg_toricity(options)
    use_upper_hit = upper_hit === nothing ? false : Bool(upper_hit)
    strict_java_float = _strict_java_float(options)
    unit_scale = Float32(_projection_unit_scale(options))
    stack_type = _pixel_hit_stack_type(options, plotbox)
    use_dense_table = dense_pixel_hits && _use_dense_pixel_hits(plotbox)
    use_flat_hit_pool =
        use_dense_table &&
        !use_upper_hit &&
        (stack_type === Vector{HitRecord}) &&
        (plotbox.nx * plotbox.ny <= typemax(FlatPoolInt)) &&
        (length(node_ids) <= typemax(FlatPoolInt))
    pixel_hits =
        if use_flat_hit_pool
            FlatPixelHitBuilder(plotbox.nx * plotbox.ny)
        else
            _pixel_hits_table(stack_type, use_dense_table, plotbox, use_upper_hit)
        end
    node_hits = zeros(Int, length(node_ids))
    projected_mesh_area = zeros(Float64, length(node_ids))
    projected_pixels_area = zeros(Float64, length(node_ids))
    _accumulate_direction_projection_prepared!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        vertices,
        faces,
        face2node_index,
        direction,
        plotbox,
        toricity,
        use_upper_hit,
        strict_java_float,
        unit_scale,
        stack_type,
    )
    use_flat_hit_pool && (pixel_hits = _finalize_flat_pixel_hits(pixel_hits))
    _apply_debug_drop_leading_hit!(pixel_hits, node_hits, projected_pixels_area, plotbox, options, node_ids)
    return DenseDirectionProjectionResult(pixel_hits, node_hits, projected_mesh_area, projected_pixels_area)
end

function _direction_projection(vertices, faces, face2node, direction, options::LightOptions, plotbox; upper_hit::Union{Nothing,Bool}=nothing)
    node_ids = unique(face2node)
    node_index = Dict{Int,Int}(nid => i for (i, nid) in enumerate(node_ids))
    face2node_index = [node_index[nid] for nid in face2node]
    return _direction_projection_materialized(
        vertices,
        faces,
        face2node,
        face2node_index,
        node_ids,
        direction,
        options,
        plotbox;
        upper_hit=upper_hit,
        dense_pixel_hits=false,
    )
end

function _direction_projection(
    geometry::InterceptionSceneData,
    direction,
    options::LightOptions;
    upper_hit::Union{Nothing,Bool}=nothing,
)
    return _direction_projection_materialized(
        geometry.vertices,
        geometry.faces,
        geometry.face2node,
        geometry.face2node_index,
        geometry.node_ids,
        direction,
        options,
        geometry.plotbox;
        upper_hit=upper_hit,
        dense_pixel_hits=true,
    )
end

function _direction_projection_cached(vertices, faces, face2node, direction, options::LightOptions, plotbox, cache_ctx; upper_hit::Bool=false, strict_java_float::Bool=false)
    if cache_ctx === nothing
        return _direction_projection(vertices, faces, face2node, direction, options, plotbox; upper_hit=upper_hit)
    end

    path = _projection_cache_path(cache_ctx, direction, upper_hit, strict_java_float)
    if isfile(path)
        cached = _read_projection_cache(path)
        cached !== nothing && return cached
        rm(path; force=true)
    end

    result = _direction_projection(vertices, faces, face2node, direction, options, plotbox; upper_hit=upper_hit)
    _write_projection_cache(path, result)
    return result
end

function _direction_projection_cached(
    geometry::InterceptionSceneData,
    direction,
    options::LightOptions,
    cache_ctx;
    upper_hit::Bool=false,
    strict_java_float::Bool=false,
)
    if cache_ctx === nothing
        return _direction_projection_prepared(
            geometry.vertices,
            geometry.faces,
            geometry.face2node_index,
            geometry.node_ids,
            direction,
            options,
            geometry.plotbox;
            upper_hit=upper_hit,
            dense_pixel_hits=true,
        )
    end

    path = _projection_cache_path(cache_ctx, direction, upper_hit, strict_java_float)
    if isfile(path)
        cached = _read_projection_cache(path)
        cached !== nothing && return cached
        rm(path; force=true)
    end

    result = _direction_projection_materialized(
        geometry.vertices,
        geometry.faces,
        geometry.face2node,
        geometry.face2node_index,
        geometry.node_ids,
        direction,
        options,
        geometry.plotbox;
        upper_hit=upper_hit,
        dense_pixel_hits=false,
    )
    _write_projection_cache(path, result)
    return result
end

function _strict_java_float(options::LightOptions)
    false
end

function _projection_unit_scale(options::LightOptions)
    1.0
end

function _paving_mesh(plotbox, cobble_count::Int, first_node_id::Int)
    # Match Java BoxPaving behavior: compute paving in cm with float-based coordinates.
    x_min_m = Float32(plotbox.origin_x)
    y_min_m = Float32(plotbox.origin_y)
    x_max_m = Float32(plotbox.origin_x + plotbox.xdim)
    y_max_m = Float32(plotbox.origin_y + plotbox.ydim)

    plot_x_m = Float64(Float32(x_max_m - x_min_m))
    plot_y_m = Float64(Float32(y_max_m - y_min_m))
    plot_area_m2 = plot_x_m * plot_y_m
    cobble_area_m2 = plot_area_m2 / max(cobble_count, 1)
    cobble_edge_m = sqrt(cobble_area_m2)

    nx = max(1, floor(Int, plot_x_m / cobble_edge_m))
    ny = max(1, floor(Int, plot_y_m / cobble_edge_m))
    cobble_x_cm = (plot_x_m / nx) * 100.0
    cobble_y_cm = (plot_y_m / ny) * 100.0

    x_min_cm = Float64(Float32(x_min_m * 100.0f0))
    y_min_cm = Float64(Float32(y_min_m * 100.0f0))
    x_max_cm = Float64(Float32(x_max_m * 100.0f0))
    y_max_cm = Float64(Float32(y_max_m * 100.0f0))

    min_size_cm = 1e-4
    z_cm = Float32(0.5)

    vertices = StaticArrays.SVector{3,Float64}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    node_area = Dict{Int,Float64}()

    node_id = first_node_id
    x_cm = x_min_cm
    while x_cm < x_max_cm
        x_size_cm = (x_cm + cobble_x_cm > x_max_cm) ? (x_max_cm - x_cm) : cobble_x_cm
        if x_size_cm > min_size_cm
            y_cm = y_min_cm
            while y_cm < y_max_cm
                y_size_cm = (y_cm + cobble_y_cm > y_max_cm) ? (y_max_cm - y_cm) : cobble_y_cm
                if y_size_cm > min_size_cm
                    x_center_cm = x_cm + (x_size_cm / 2.0)
                    y_center_cm = y_cm + (y_size_cm / 2.0)

                    x0 = Float32(-x_size_cm / 2.0)
                    x1 = Float32(x_size_cm / 2.0)
                    y0 = Float32(-y_size_cm / 2.0)
                    y1 = Float32(y_size_cm / 2.0)
                    xc = Float32(x_center_cm)
                    yc = Float32(y_center_cm)
                    # Match Java cm->m conversion path (float scaling by 0.01f).
                    p1x = Float32(x0 + xc) * 0.01f0
                    p2x = Float32(x1 + xc) * 0.01f0
                    p3x = Float32(x1 + xc) * 0.01f0
                    p4x = Float32(x0 + xc) * 0.01f0
                    p1y = Float32(y0 + yc) * 0.01f0
                    p2y = Float32(y0 + yc) * 0.01f0
                    p3y = Float32(y1 + yc) * 0.01f0
                    p4y = Float32(y1 + yc) * 0.01f0
                    z_m = z_cm * 0.01f0

                    p1 = StaticArrays.SVector{3,Float64}(Float64(p1x), Float64(p1y), Float64(z_m))
                    p2 = StaticArrays.SVector{3,Float64}(Float64(p2x), Float64(p2y), Float64(z_m))
                    p3 = StaticArrays.SVector{3,Float64}(Float64(p3x), Float64(p3y), Float64(z_m))
                    p4 = StaticArrays.SVector{3,Float64}(Float64(p4x), Float64(p4y), Float64(z_m))

                    base = length(vertices)
                    push!(vertices, p1, p2, p3, p4)
                    push!(faces, PlantGeom.Face3(base + 1, base + 2, base + 3))
                    push!(faces, PlantGeom.Face3(base + 1, base + 3, base + 4))
                    push!(face2node, node_id, node_id)
                    node_area[node_id] = (x_size_cm * y_size_cm) / 10000.0
                    node_id += 1
                end
                y_cm += cobble_y_cm
            end
        end
        x_cm += cobble_x_cm
    end

    return vertices, faces, face2node, node_area
end

function _scene_geometry_for_interception(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions;
    include_raycore_instancing::Bool=false,
)
    raw_vertices = GeometryBasics.decompose(GeometryBasics.Point3, scene.merged_mesh)
    vertices = [StaticArrays.SVector{3,Float64}(v[1], v[2], v[3]) for v in raw_vertices]
    all_faces = collect(GeometryBasics.decompose(PlantGeom.Face3, scene.merged_mesh))
    all_face2node = collect(scene.face2node)

    ignored = _ignored_group_types(models)
    _validate_scene_models(scene, all_face2node, models, ignored)
    faces = PlantGeom.Face3[]
    face2node = Int[]
    for i in eachindex(all_faces)
        node_id = all_face2node[i]
        _is_ignored_node(node_id, scene, ignored) && continue
        push!(faces, all_faces[i])
        push!(face2node, node_id)
    end
    isempty(face2node) && error("No intercepting geometry left after applying ignore rules.")

    plotbox = _plotbox(scene, vertices, options.pixel_size)

    node_ids = unique(face2node)
    node_index = Dict{Int,Int}(nid => i for (i, nid) in enumerate(node_ids))
    face2node_index = [node_index[nid] for nid in face2node]
    node_group = Dict{Int,String}(nid => _scene_group(scene, nid, "") for nid in node_ids)
    node_group_by_index = [get(node_group, nid, "") for nid in node_ids]
    pavement_node_mask = [group == "pavement" for group in node_group_by_index]
    node_type = Dict{Int,String}(nid => _scene_type(scene, nid, "") for nid in node_ids)
    node_type_by_index = [get(node_type, nid, "") for nid in node_ids]
    raycore_instanced_geometry =
        include_raycore_instancing ?
        _raycore_instanced_geometry_from_scene(
            scene,
            vertices,
            node_ids,
            node_index,
            faces,
            face2node,
            face2node_index;
            toricity=options.toricity,
        ) :
        nothing

    return InterceptionSceneData(
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
    )
end

function _raycore_reference_instancing_enabled()
    raw = lowercase(strip(get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING", "1")))
    return !(raw in ("0", "false", "no", "off"))
end

function _raycore_reference_instance_limit()
    raw = get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT", "")
    isempty(strip(raw)) && return 4096
    value = parse(Int, raw)
    return value <= 0 ? typemax(Int) : value
end

function _raycore_reference_min_face_savings_ratio()
    raw = get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_MIN_FACE_SAVINGS_RATIO", "")
    isempty(strip(raw)) && return 0.20
    return clamp(parse(Float64, raw), 0.0, 1.0)
end

function _raycore_mat4f_from_matrix(mat::AbstractMatrix)
    size(mat) == (4, 4) || return nothing
    all(isfinite, mat) || return nothing
    return GeometryBasics.Mat4f((Float32(mat[i]) for i in eachindex(mat))...)
end

function _raycore_transform_matrix4(transformation)
    mat = try
        PlantGeom.transformation_matrix4(transformation)
    catch
        return nothing
    end
    return _raycore_mat4f_from_matrix(mat)
end

function _raycore_reference_mesh_identity_key(ref_mesh)
    return (String(ref_mesh.name), hash(ref_mesh.mesh), ref_mesh.taper)
end

function _raycore_geometry_node_candidate_for_instancing(node)
    PlantGeom.has_geometry(node) || return false
    geom = node[:geometry]
    geom isa PlantGeom.Geometry || return false
    geom.ref_mesh.taper && return false
    return true
end

function _raycore_geometry_node_supported_for_instancing(node)
    _raycore_geometry_node_candidate_for_instancing(node) || return false
    geom = node[:geometry]
    _raycore_transform_matrix4(geom.transformation) === nothing && return false
    return true
end

function _raycore_mesh_from_reference_mesh(ref_mesh)
    raw_vertices = GeometryBasics.decompose(GeometryBasics.Point3, ref_mesh.mesh)
    raw_faces = collect(GeometryBasics.decompose(PlantGeom.Face3, ref_mesh.mesh))
    (isempty(raw_vertices) || isempty(raw_faces)) && return nothing
    points = Vector{GeometryBasics.Point3f}(undef, length(raw_vertices))
    faces = Vector{GeometryBasics.TriangleFace{Int}}(undef, length(raw_faces))
    face_meta = Vector{UInt32}(undef, length(raw_faces))
    @inbounds for (i, p) in pairs(raw_vertices)
        points[i] = GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3]))
    end
    @inbounds for (i, f) in pairs(raw_faces)
        faces[i] = GeometryBasics.TriangleFace{Int}(f[1], f[2], f[3])
        # The prototype-local metadata is intentionally independent of the
        # final Archimed node. All instances decode this local surface id
        # through the per-instance table.
        face_meta[i] = UInt32(1)
    end
    return GeometryBasics.Mesh(points, faces; face_meta=GeometryBasics.per_face(face_meta, faces))
end

function _raycore_fallback_mesh_from_geometry(
    vertices,
    faces::Vector{PlantGeom.Face3},
    face2node_index::Vector{Int},
    fallback_node_ids::Set{Int},
    face2node::Vector{Int},
)
    fallback_faces = Int[]
    @inbounds for i in eachindex(faces)
        face2node[i] in fallback_node_ids || continue
        push!(fallback_faces, i)
    end
    isempty(fallback_faces) && return nothing

    points = Vector{GeometryBasics.Point3f}(undef, length(vertices))
    out_faces = Vector{GeometryBasics.TriangleFace{Int}}(undef, length(fallback_faces))
    face_meta = Vector{UInt32}(undef, length(fallback_faces))
    @inbounds for (i, p) in pairs(vertices)
        points[i] = GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3]))
    end
    @inbounds for (out_i, src_i) in pairs(fallback_faces)
        f = faces[src_i]
        out_faces[out_i] = GeometryBasics.TriangleFace{Int}(f[1], f[2], f[3])
        face_meta[out_i] = UInt32(face2node_index[src_i])
    end
    return GeometryBasics.Mesh(points, out_faces; face_meta=GeometryBasics.per_face(face_meta, out_faces))
end

function _raycore_reference_instancing_stats_from_counts(;
    expanded_face_count::Int,
    prototype_face_count::Int,
    fallback_face_count::Int,
    prototype_node_count::Int,
    fallback_present::Bool,
    metadata_stride::Int,
    toricity::Bool=false,
    prototype_count::Int=0,
)
    toric_copies = _raycore_toric_copy_count(toricity)
    fallback_instance_count = fallback_present ? 1 : 0
    instance_count = (prototype_node_count + fallback_instance_count) * toric_copies
    compact_face_count = prototype_face_count + fallback_face_count
    saved_faces = max(expanded_face_count - compact_face_count, 0)
    savings_ratio = expanded_face_count == 0 ? 0.0 : saved_faces / expanded_face_count
    decoder_entries = instance_count * metadata_stride
    return (
        instance_count=instance_count,
        prototype_count=prototype_count,
        prototype_face_count=prototype_face_count,
        prototype_node_count=prototype_node_count,
        fallback_face_count=fallback_face_count,
        compact_face_count=compact_face_count,
        expanded_face_count=expanded_face_count,
        saved_faces=saved_faces,
        savings_ratio=savings_ratio,
        decoder_entries=decoder_entries,
    )
end

function _raycore_reference_instancing_stats_pass(stats)
    stats.saved_faces > 0 || return false
    stats.savings_ratio >= _raycore_reference_min_face_savings_ratio() || return false
    stats.instance_count <= _raycore_reference_instance_limit() || return false
    return true
end

function _raycore_reference_instancing_diagnostics(
    scene::PlantGeom.SceneGeometry,
    geometry::InterceptionSceneData;
    toricity::Bool=false,
)
    enabled = _raycore_reference_instancing_enabled()
    base = (
        enabled=enabled,
        status=enabled ? :unknown : :disabled,
        geometry_nodes=0,
        candidate_nodes=0,
        supported_nodes=0,
        tapered_nodes=0,
        non_geometry_nodes=0,
        no_transform_nodes=0,
        unique_refs=0,
        reusable_refs=0,
        reusable_nodes=0,
        fallback_nodes=length(geometry.node_ids),
        instance_count=0,
        prototype_count=0,
        prototype_face_count=0,
        fallback_face_count=length(geometry.faces),
        compact_face_count=length(geometry.faces),
        expanded_face_count=length(geometry.faces),
        saved_faces=0,
        savings_ratio=0.0,
        decoder_entries=0,
    )
    enabled || return base
    scene.mtg === nothing && return merge(base, (; status=:no_mtg))

    nodes = Any[]
    seen_node_ids = Set{Int}()
    MultiScaleTreeGraph.traverse!(scene.mtg) do node
        PlantGeom.has_geometry(node) || return nothing
        haskey(node, :scene_dimensions) && return nothing
        nid = MultiScaleTreeGraph.node_id(node)
        haskey(geometry.node_index, nid) || return nothing
        nid in seen_node_ids && return nothing
        push!(seen_node_ids, nid)
        push!(nodes, node)
        return nothing
    end
    isempty(nodes) && return merge(base, (; status=:no_geometry_nodes))

    ref_counts = Dict{Any,Int}()
    ref_face_counts = Dict{Any,Int}()
    candidate_nodes = 0
    supported_nodes = 0
    tapered_nodes = 0
    non_geometry_nodes = 0
    no_transform_nodes = 0
    for node in nodes
        geom = node[:geometry]
        if !(geom isa PlantGeom.Geometry)
            non_geometry_nodes += 1
            continue
        end
        if geom.ref_mesh.taper
            tapered_nodes += 1
            continue
        end
        candidate_nodes += 1
        key = _raycore_reference_mesh_identity_key(geom.ref_mesh)
        ref_counts[key] = get(ref_counts, key, 0) + 1
        if !haskey(ref_face_counts, key)
            ref_face_counts[key] = length(GeometryBasics.decompose(PlantGeom.Face3, geom.ref_mesh.mesh))
        end
        if _raycore_transform_matrix4(geom.transformation) === nothing
            no_transform_nodes += 1
        else
            supported_nodes += 1
        end
    end

    reusable_keys = Set{Any}(key for (key, count) in ref_counts if count >= 2)
    reusable_node_ids = Set{Int}()
    reusable_nodes = 0
    for node in nodes
        geom = node[:geometry]
        geom isa PlantGeom.Geometry || continue
        geom.ref_mesh.taper && continue
        key = _raycore_reference_mesh_identity_key(geom.ref_mesh)
        key in reusable_keys || continue
        reusable_nodes += 1
        push!(reusable_node_ids, MultiScaleTreeGraph.node_id(node))
    end

    fallback_node_ids = Set{Int}(nid for nid in geometry.node_ids if !(nid in reusable_node_ids))
    fallback_face_count = count(nid -> nid in fallback_node_ids, geometry.face2node)
    prototype_face_count = sum(ref_face_counts[key] for key in reusable_keys; init=0)
    metadata_stride = max(length(geometry.node_ids) + 1, 2)
    stats = _raycore_reference_instancing_stats_from_counts(;
        expanded_face_count=length(geometry.faces),
        prototype_face_count=prototype_face_count,
        fallback_face_count=fallback_face_count,
        prototype_node_count=reusable_nodes,
        fallback_present=!isempty(fallback_node_ids),
        metadata_stride=metadata_stride,
        toricity=toricity,
        prototype_count=length(reusable_keys),
    )
    status =
        isempty(reusable_keys) ? :no_reusable_references :
        no_transform_nodes > 0 ? :unsupported_transform :
        _raycore_reference_instancing_stats_pass(stats) ? :eligible :
        stats.savings_ratio < _raycore_reference_min_face_savings_ratio() ? :insufficient_savings :
        stats.instance_count > _raycore_reference_instance_limit() ? :instance_limit :
        :not_beneficial

    return (
        enabled=enabled,
        status=status,
        geometry_nodes=length(nodes),
        candidate_nodes=candidate_nodes,
        supported_nodes=supported_nodes,
        tapered_nodes=tapered_nodes,
        non_geometry_nodes=non_geometry_nodes,
        no_transform_nodes=no_transform_nodes,
        unique_refs=length(ref_counts),
        reusable_refs=length(reusable_keys),
        reusable_nodes=reusable_nodes,
        fallback_nodes=length(fallback_node_ids),
        instance_count=stats.instance_count,
        prototype_count=stats.prototype_count,
        prototype_face_count=stats.prototype_face_count,
        fallback_face_count=stats.fallback_face_count,
        compact_face_count=stats.compact_face_count,
        expanded_face_count=stats.expanded_face_count,
        saved_faces=stats.saved_faces,
        savings_ratio=stats.savings_ratio,
        decoder_entries=stats.decoder_entries,
    )
end

function _raycore_reference_instancing_diagnostics(
    scene::PlantGeom.SceneGeometry,
    prepared::PreparedInterceptionData,
    options::LightOptions,
)
    return _raycore_reference_instancing_diagnostics(scene, prepared.geometry; toricity=options.toricity)
end

function _raycore_instanced_geometry_from_scene(
    scene::PlantGeom.SceneGeometry,
    vertices,
    node_ids::Vector{Int},
    node_index::Dict{Int,Int},
    faces::Vector{PlantGeom.Face3},
    face2node::Vector{Int},
    face2node_index::Vector{Int},
    ;
    toricity::Bool=false,
)
    _raycore_reference_instancing_enabled() || return nothing
    scene.mtg === nothing && return nothing

    nodes = Any[]
    seen_node_ids = Set{Int}()
    MultiScaleTreeGraph.traverse!(scene.mtg) do node
        PlantGeom.has_geometry(node) || return nothing
        haskey(node, :scene_dimensions) && return nothing
        nid = MultiScaleTreeGraph.node_id(node)
        haskey(node_index, nid) || return nothing
        nid in seen_node_ids && return nothing
        push!(seen_node_ids, nid)
        push!(nodes, node)
        return nothing
    end
    isempty(nodes) && return nothing

    ref_counts = Dict{Any,Int}()
    supported = Dict{Int,Bool}()
    ref_face_counts = Dict{Any,Int}()
    for node in nodes
        nid = MultiScaleTreeGraph.node_id(node)
        ok = _raycore_geometry_node_candidate_for_instancing(node)
        supported[nid] = ok
        ok || continue
        ref_mesh = node[:geometry].ref_mesh
        key = _raycore_reference_mesh_identity_key(ref_mesh)
        ref_counts[key] = get(ref_counts, key, 0) + 1
        if !haskey(ref_face_counts, key)
            ref_face_counts[key] = length(GeometryBasics.decompose(PlantGeom.Face3, ref_mesh.mesh))
        end
    end

    reusable_nodes = Set{Int}()
    reusable_keys = Set{Any}()
    for node in nodes
        nid = MultiScaleTreeGraph.node_id(node)
        get(supported, nid, false) || continue
        key = _raycore_reference_mesh_identity_key(node[:geometry].ref_mesh)
        get(ref_counts, key, 0) >= 2 || continue
        push!(reusable_nodes, nid)
        push!(reusable_keys, key)
    end
    isempty(reusable_nodes) && return nothing

    prototype_face_count = sum(ref_face_counts[key] for key in reusable_keys; init=0)
    fallback_node_ids = Set{Int}(nid for nid in node_ids if !(nid in reusable_nodes))
    fallback_face_count = count(nid -> nid in fallback_node_ids, face2node)
    metadata_stride = max(length(node_ids) + 1, 2)
    early_stats = _raycore_reference_instancing_stats_from_counts(;
        expanded_face_count=length(faces),
        prototype_face_count=prototype_face_count,
        fallback_face_count=fallback_face_count,
        prototype_node_count=length(reusable_nodes),
        fallback_present=!isempty(fallback_node_ids),
        metadata_stride=metadata_stride,
        toricity=toricity,
        prototype_count=length(reusable_keys),
    )
    _raycore_reference_instancing_stats_pass(early_stats) || return nothing

    prototype_index = Dict{Any,Int}()
    prototype_ref_meshes = Any[]
    prototype_transforms = Vector{Vector{GeometryBasics.Mat4f}}()
    prototype_node_indices = Vector{Vector{UInt32}}()
    for node in nodes
        nid = MultiScaleTreeGraph.node_id(node)
        nid in reusable_nodes || continue
        geom = node[:geometry]
        key = _raycore_reference_mesh_identity_key(geom.ref_mesh)
        proto_idx = get!(prototype_index, key) do
            push!(prototype_ref_meshes, geom.ref_mesh)
            push!(prototype_transforms, GeometryBasics.Mat4f[])
            push!(prototype_node_indices, UInt32[])
            length(prototype_ref_meshes)
        end
        transform = _raycore_transform_matrix4(geom.transformation)
        transform === nothing && return nothing
        push!(prototype_transforms[proto_idx], transform)
        push!(prototype_node_indices[proto_idx], UInt32(node_index[nid]))
    end

    prototypes = RaycorePrototypeInstances[]
    for i in eachindex(prototype_ref_meshes)
        mesh = _raycore_mesh_from_reference_mesh(prototype_ref_meshes[i])
        mesh === nothing && return nothing
        push!(prototypes, RaycorePrototypeInstances(mesh, prototype_transforms[i], prototype_node_indices[i]))
    end

    fallback_mesh = _raycore_fallback_mesh_from_geometry(vertices, faces, face2node_index, fallback_node_ids, face2node)
    prototype_node_count = sum(length(proto.node_indices) for proto in prototypes)
    return RaycoreInstancedSceneGeometry(
        prototypes,
        fallback_mesh,
        metadata_stride,
        prototype_face_count,
        prototype_node_count,
        fallback_face_count,
    )
end

_light_render_geometry(geometry::InterceptionSceneData) =
    LightRenderGeometry(geometry.vertices, geometry.faces, geometry.face2node)

function _interception_area_per_node_from_geometry(geometry::InterceptionSceneData)
    area = zeros(Float64, length(geometry.node_ids))
    @inbounds for i in eachindex(geometry.faces)
        f = geometry.faces[i]
        a = _triangle_area3d(geometry.vertices[f[1]], geometry.vertices[f[2]], geometry.vertices[f[3]])
        area[geometry.face2node_index[i]] += a
    end
    return _all_dense_float_node_map(geometry.node_ids, area)
end

function _node_absorptance_maps_from_geometry(
    options::LightOptions,
    geometry::InterceptionSceneData,
    virtual_node_mask::Vector{Bool},
    node_interception_by_index::Vector{Union{Nothing,InterceptionModel}},
)
    sf_default_par = _default_scattering_factor_local(options, "PAR")
    sf_default_nir = _default_scattering_factor_local(options, "NIR")
    par = zeros(Float64, length(geometry.node_ids))
    nir = zeros(Float64, length(geometry.node_ids))
    @inbounds for i in eachindex(geometry.node_ids)
        if virtual_node_mask[i]
            par[i] = 0.0
            nir[i] = 0.0
            continue
        end
        interception = node_interception_by_index[i]
        props = interception === nothing ? nothing : interception.optical_properties
        if props === nothing
            par[i] = clamp(1.0 - sf_default_par, 0.0, 1.0)
            nir[i] = clamp(1.0 - sf_default_nir, 0.0, 1.0)
            continue
        end
        has_par = get(props.extras, "__has_par", true)
        has_nir = get(props.extras, "__has_nir", true)
        par_sf = has_par ? props.par : sf_default_par
        nir_sf = has_nir ? props.nir : sf_default_nir
        par[i] = clamp(1.0 - par_sf, 0.0, 1.0)
        nir[i] = clamp(1.0 - nir_sf, 0.0, 1.0)
    end
    return _all_dense_float_node_map(geometry.node_ids, par), _all_dense_float_node_map(geometry.node_ids, nir)
end

function _node_absorptance_per_band_from_geometry(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions,
    geometry::InterceptionSceneData,
    virtual_nodes::Set{Int},
    band::String,
)
    virtual_node_mask = [nid in virtual_nodes for nid in geometry.node_ids]
    node_interception_by_index = _resolved_node_interception_models(geometry, models)
    par, nir = _node_absorptance_maps_from_geometry(options, geometry, virtual_node_mask, node_interception_by_index)
    return uppercase(band) == "NIR" ? nir : par
end

function _prepare_interception_data(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions;
    include_budget_maps::Bool=false,
    include_raycore_instancing::Bool=false,
)
    geometry = _scene_geometry_for_interception(
        scene,
        models,
        options;
        include_raycore_instancing=include_raycore_instancing,
    )
    node_interception_by_index = _resolved_node_interception_models(geometry, models)
    virtual_nodes = _virtual_sensor_node_ids(geometry, node_interception_by_index)
    virtual_node_mask = [nid in virtual_nodes for nid in geometry.node_ids]
    upper_hit = _use_upper_hit_pixel_table(models, options)
    cache_ctx = _projection_cache_context(geometry.vertices, geometry.faces, geometry.face2node, geometry.plotbox, options)
    emit_par, emit_nir = _emitter_power_per_node(scene, models)
    emitter_par_power_by_index = zeros(Float64, length(geometry.node_ids))
    emitter_nir_power_by_index = zeros(Float64, length(geometry.node_ids))
    for (nid, power) in emit_par
        emitter_par_power_by_index[geometry.node_index[nid]] = power
    end
    for (nid, power) in emit_nir
        emitter_nir_power_by_index[geometry.node_index[nid]] = power
    end
    emitter_nodes = Set(union(keys(emit_par), keys(emit_nir)))
    emitter_node_mask = [nid in emitter_nodes for nid in geometry.node_ids]

    component_area_per_node =
        include_budget_maps ? _interception_area_per_node_from_geometry(geometry) : nothing
    absorption_par_per_node = nothing
    absorption_nir_per_node = nothing
    if include_budget_maps
        absorption_par_per_node, absorption_nir_per_node =
            _node_absorptance_maps_from_geometry(options, geometry, virtual_node_mask, node_interception_by_index)
    end

    return PreparedInterceptionData(
        geometry,
        node_interception_by_index,
        _node_transparency_by_index(node_interception_by_index, virtual_node_mask),
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

function _raycore_point3f(p)
    return GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3]))
end

function _raycore_vec3f(v)
    return GeometryBasics.Vec3f(Float32(v[1]), Float32(v[2]), Float32(v[3]))
end

function _raycore_toric_shifts(geometry::InterceptionSceneData; toricity::Bool=false)
    if toricity
        pb = geometry.plotbox
        return [(Float64(ix) * pb.xdim, Float64(iy) * pb.ydim, 0.0) for ix in -1:1 for iy in -1:1]
    end
    return [(0.0, 0.0, 0.0)]
end

function _raycore_translation_matrix(shift)
    sx, sy, sz = shift
    return GeometryBasics.Mat4f(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        Float32(sx), Float32(sy), Float32(sz), 1,
    )
end

function _raycore_toric_transforms(geometry::InterceptionSceneData; toricity::Bool=false)
    return [_raycore_translation_matrix(shift) for shift in _raycore_toric_shifts(geometry; toricity=toricity)]
end

function _raycore_matrix_from_mat4f(mat::GeometryBasics.Mat4f)
    out = Matrix{Float64}(undef, 4, 4)
    @inbounds for i in 1:4, j in 1:4
        out[i, j] = Float64(mat[i, j])
    end
    return out
end

function _raycore_shifted_transform(transform::GeometryBasics.Mat4f, shift)
    sx, sy, sz = shift
    shift_mat = zeros(Float64, 4, 4)
    @inbounds for i in 1:4
        shift_mat[i, i] = 1.0
    end
    shift_mat[1, 4] = Float64(sx)
    shift_mat[2, 4] = Float64(sy)
    shift_mat[3, 4] = Float64(sz)
    return _raycore_mat4f_from_matrix(shift_mat * _raycore_matrix_from_mat4f(transform))
end

function _raycore_mesh_from_geometry(geometry::InterceptionSceneData)
    n_vertices = length(geometry.vertices)
    n_faces = length(geometry.faces)
    points = Vector{GeometryBasics.Point3f}(undef, n_vertices)
    faces = Vector{GeometryBasics.TriangleFace{Int}}(undef, n_faces)
    face_node_indices = Vector{UInt32}(undef, n_faces)

    @inbounds for (i, p) in pairs(geometry.vertices)
        points[i] = GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3]))
    end
    @inbounds for (i, f) in pairs(geometry.faces)
        faces[i] = GeometryBasics.TriangleFace{Int}(f[1], f[2], f[3])
        face_node_indices[i] = UInt32(geometry.face2node_index[i])
    end

    face_meta = GeometryBasics.per_face(face_node_indices, faces)
    return GeometryBasics.Mesh(points, faces; face_meta=face_meta)
end

function _raycore_push_geometry!(tlas, geometry::InterceptionSceneData; toricity::Bool=false)
    mesh = _raycore_mesh_from_geometry(geometry)
    if toricity
        push!(tlas, mesh, _raycore_toric_transforms(geometry; toricity=true))
    else
        push!(tlas, mesh)
    end
    return toricity ? 9 : 1
end

function _raycore_decoder_row!(table::Vector{UInt32}, instance_idx::Int, metadata_stride::Int)
    offset = (instance_idx - 1) * metadata_stride
    length(table) >= offset + metadata_stride || resize!(table, offset + metadata_stride)
    @inbounds for i in 1:metadata_stride
        table[offset+i] = UInt32(0)
    end
    return offset
end

function _raycore_push_instanced_geometry!(
    tlas,
    instanced::RaycoreInstancedSceneGeometry,
    geometry::InterceptionSceneData;
    toricity::Bool=false,
)
    shifts = _raycore_toric_shifts(geometry; toricity=toricity)
    decoder_table = UInt32[]
    instance_idx = 0
    metadata_stride = instanced.metadata_stride

    for prototype in instanced.prototypes
        transforms = GeometryBasics.Mat4f[]
        node_indices = UInt32[]
        for (base_transform, node_idx) in zip(prototype.transforms, prototype.node_indices)
            for shift in shifts
                transform = _raycore_shifted_transform(base_transform, shift)
                transform === nothing && return nothing
                push!(transforms, transform)
                push!(node_indices, node_idx)
            end
        end
        isempty(transforms) && continue
        push!(tlas, prototype.mesh, transforms)
        for node_idx in node_indices
            instance_idx += 1
            offset = _raycore_decoder_row!(decoder_table, instance_idx, metadata_stride)
            decoder_table[offset+2] = node_idx
        end
    end

    if instanced.fallback_mesh !== nothing
        fallback_transforms = _raycore_toric_transforms(geometry; toricity=toricity)
        push!(tlas, instanced.fallback_mesh, fallback_transforms)
        for _ in fallback_transforms
            instance_idx += 1
            offset = _raycore_decoder_row!(decoder_table, instance_idx, metadata_stride)
            @inbounds for node_idx in eachindex(geometry.node_ids)
                decoder_table[offset+node_idx+1] = UInt32(node_idx)
            end
        end
    end

    instance_idx == 0 && return nothing
    return (instance_count=instance_idx, metadata_stride=metadata_stride, decoder_table=decoder_table)
end

function _raycore_reference_instancing_stats(
    instanced::RaycoreInstancedSceneGeometry,
    geometry::InterceptionSceneData;
    toricity::Bool=false,
)
    return _raycore_reference_instancing_stats_from_counts(;
        expanded_face_count=length(geometry.faces),
        prototype_face_count=instanced.prototype_face_count,
        fallback_face_count=instanced.fallback_face_count,
        prototype_node_count=instanced.prototype_node_count,
        fallback_present=instanced.fallback_mesh !== nothing,
        metadata_stride=instanced.metadata_stride,
        toricity=toricity,
        prototype_count=length(instanced.prototypes),
    )
end

function _raycore_should_use_reference_instancing(
    instanced::RaycoreInstancedSceneGeometry,
    geometry::InterceptionSceneData;
    toricity::Bool=false,
)
    stats = _raycore_reference_instancing_stats(instanced, geometry; toricity=toricity)
    return _raycore_reference_instancing_stats_pass(stats)
end

function _raycore_face_chunk_limit()
    raw = get(ENV, "ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT", "")
    isempty(strip(raw)) && return 1536
    value = parse(Int, raw)
    return value <= 0 ? typemax(Int) : value
end

function _raycore_prechunk_face_threshold()
    raw = get(ENV, "ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD", "")
    isempty(strip(raw)) && return 500_000
    value = parse(Int, raw)
    return value <= 0 ? typemax(Int) : value
end

_raycore_toric_copy_count(toricity::Bool) = toricity ? 9 : 1

function _raycore_effective_face_count(prepared::PreparedInterceptionData; toricity::Bool=false)
    return length(prepared.geometry.faces) * _raycore_toric_copy_count(toricity)
end

function _raycore_estimated_prechunk_instances(
    prepared::PreparedInterceptionData;
    toricity::Bool=false,
    face_chunk_limit::Int=_raycore_face_chunk_limit(),
)
    effective_faces = _raycore_effective_face_count(prepared; toricity=toricity)
    face_chunk_limit == typemax(Int) && return effective_faces == 0 ? 0 : 1
    return cld(effective_faces, max(face_chunk_limit, 1))
end

function _raycore_should_prechunk_scene(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig;
    toricity::Bool=false,
)
    config.backend isa KernelAbstractions.CPU && return false
    if prepared.geometry.raycore_instanced_geometry !== nothing &&
       _raycore_should_use_reference_instancing(
           prepared.geometry.raycore_instanced_geometry,
           prepared.geometry;
           toricity=toricity,
       )
        return false
    end
    _raycore_face_chunk_limit() == typemax(Int) && return false
    threshold = _raycore_prechunk_face_threshold()
    threshold == typemax(Int) && return false
    return _raycore_effective_face_count(prepared; toricity=toricity) > threshold
end

function _raycore_prechunk_instance_limit_status(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig;
    toricity::Bool=false,
)
    face_chunk_limit = _raycore_face_chunk_limit()
    max_instances = config.max_prechunk_instances
    estimated_instances =
        _raycore_estimated_prechunk_instances(prepared; toricity=toricity, face_chunk_limit=face_chunk_limit)
    should_prechunk = _raycore_should_prechunk_scene(prepared, config; toricity=toricity)
    exceeded =
        should_prechunk &&
        max_instances != typemax(Int) &&
        estimated_instances > max_instances
    return (
        exceeded=exceeded,
        should_prechunk=should_prechunk,
        estimated_instances=estimated_instances,
        max_instances=max_instances,
        face_chunk_limit=face_chunk_limit,
        effective_faces=_raycore_effective_face_count(prepared; toricity=toricity),
    )
end

function _raycore_prechunk_instance_limit_message(status)
    return "Raycore backend would require about $(status.estimated_instances) prechunked BLAS instances " *
           "for $(status.effective_faces) effective faces with face_chunk_limit=$(status.face_chunk_limit), " *
           "exceeding max_prechunk_instances=$(status.max_instances)."
end

function _raycore_is_valid_face(geometry::InterceptionSceneData, face)
    p1 = geometry.vertices[face[1]]
    p2 = geometry.vertices[face[2]]
    p3 = geometry.vertices[face[3]]
    return _triangle_area3d(p1, p2, p3) > 0.0
end

function _raycore_foreach_mesh_chunk_from_geometry(
    emit::F,
    geometry::InterceptionSceneData;
    toricity::Bool=false,
    face_chunk_limit::Int=_raycore_face_chunk_limit(),
) where {F}
    if face_chunk_limit <= 0
        emit(_raycore_mesh_from_geometry(geometry))
        return 1
    end

    source_faces = geometry.faces
    vertex_map = Dict{Int,Int}()
    points = GeometryBasics.Point3f[]
    faces = GeometryBasics.TriangleFace{Int}[]
    face_node_indices = UInt32[]
    chunk_count = 0

    function flush_chunk!()
        isempty(faces) && return nothing
        face_meta = GeometryBasics.per_face(copy(face_node_indices), copy(faces))
        emit(GeometryBasics.Mesh(copy(points), copy(faces); face_meta=face_meta))
        chunk_count += 1
        empty!(vertex_map)
        empty!(points)
        empty!(faces)
        empty!(face_node_indices)
        return nothing
    end

    @inbounds for (face_idx, face) in pairs(source_faces)
        _raycore_is_valid_face(geometry, face) || continue
        if length(faces) >= face_chunk_limit
            flush_chunk!()
        end
        local_vertices = ntuple(j -> begin
            vertex_idx = face[j]
            get!(vertex_map, vertex_idx) do
                p = geometry.vertices[vertex_idx]
                push!(points, GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])))
                length(points)
            end
        end, Val(3))
        push!(faces, GeometryBasics.TriangleFace{Int}(local_vertices...))
        push!(face_node_indices, UInt32(geometry.face2node_index[face_idx]))
    end
    flush_chunk!()
    return toricity ? 9 * chunk_count : chunk_count
end

function _raycore_mesh_chunks_from_geometry(
    geometry::InterceptionSceneData;
    toricity::Bool=false,
    face_chunk_limit::Int=_raycore_face_chunk_limit(),
)
    chunks = GeometryBasics.Mesh[]
    _raycore_foreach_mesh_chunk_from_geometry(
        mesh -> push!(chunks, mesh),
        geometry;
        toricity=toricity,
        face_chunk_limit=face_chunk_limit,
    )
    return chunks
end

function _raycore_stack_safe_scene_data(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig,
    options::LightOptions,
    data::RaycoreSceneData;
    toricity::Bool=options.toricity,
    skip_face_chunk_limit::Union{Nothing,Int}=data.chunked_tlas ? _raycore_face_chunk_limit() : nothing,
)
    config.backend isa KernelAbstractions.CPU && return data, false

    status = _raycore_traversal_stack_depth_status(data)
    status.exceeds || return data, false

    for limit in _raycore_retry_chunk_limits(skip_face_chunk_limit)
        estimated_instances = _raycore_estimated_prechunk_instances(prepared; toricity=toricity, face_chunk_limit=limit)
        if config.max_prechunk_instances != typemax(Int) && estimated_instances > config.max_prechunk_instances
            @info "Skipping Raycore proactive BLAS chunking because the chunk limit would exceed max_prechunk_instances." face_chunk_limit=limit estimated_instances max_prechunk_instances=config.max_prechunk_instances max_blas_depth=status.blas_depth traversal_stack_capacity=status.capacity
            continue
        end
        chunked = _raycore_scene_data(prepared, config; toricity=toricity, face_chunk_limit=limit)
        chunked_status = _raycore_traversal_stack_depth_status(chunked)
        if !chunked_status.exceeds
            @info "Raycore backend proactively chunked BLAS geometry to avoid exceeding the software traversal stack." face_chunk_limit=limit previous_blas_depth=status.blas_depth max_blas_depth=chunked_status.blas_depth traversal_stack_capacity=chunked_status.capacity
            return chunked, true
        end
    end

    @info "Raycore backend kept BLAS geometry whose depth exceeds the software traversal stack because no allowed chunk limit reduced it below capacity." max_blas_depth=status.blas_depth traversal_stack_capacity=status.capacity
    return data, false
end

function _raycore_initial_scene_data(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig,
    options::LightOptions;
    toricity::Bool=options.toricity,
)
    if _raycore_should_prechunk_scene(prepared, config; toricity=toricity)
        face_chunk_limit = _raycore_face_chunk_limit()
        data = _raycore_scene_data(prepared, config; toricity=toricity, face_chunk_limit=face_chunk_limit)
        @info "Raycore backend prechunked BLAS geometry before validation." face_chunk_limit effective_faces=_raycore_effective_face_count(prepared; toricity=toricity)
        safe_data, safe_chunked = _raycore_stack_safe_scene_data(
            prepared,
            config,
            options,
            data;
            toricity=toricity,
            skip_face_chunk_limit=face_chunk_limit,
        )
        return safe_data, true
    end
    data = _raycore_scene_data(prepared, config; toricity=toricity)
    return _raycore_stack_safe_scene_data(prepared, config, options, data; toricity=toricity)
end

function _raycore_scene_data(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig;
    toricity::Bool=false,
    face_chunk_limit::Union{Nothing,Int}=nothing,
)
    tlas = Raycore.TLAS(config.backend)
    chunked_tlas = face_chunk_limit !== nothing
    geometry_mode = :merged_mesh
    instanced_build =
        if face_chunk_limit === nothing &&
           prepared.geometry.raycore_instanced_geometry !== nothing &&
           _raycore_should_use_reference_instancing(
               prepared.geometry.raycore_instanced_geometry,
               prepared.geometry;
               toricity=toricity,
           )
            _raycore_push_instanced_geometry!(
                tlas,
                prepared.geometry.raycore_instanced_geometry,
                prepared.geometry;
                toricity=toricity,
            )
        else
            nothing
        end
    instance_count, metadata_stride, decoder_table =
        if instanced_build !== nothing
            geometry_mode = :reference_instances
            instanced_build.instance_count, instanced_build.metadata_stride, instanced_build.decoder_table
        elseif face_chunk_limit === nothing
            _raycore_push_geometry!(tlas, prepared.geometry; toricity=toricity),
            length(prepared.geometry.node_ids) + 1,
            UInt32[]
        else
            geometry_mode = :chunked_merged_mesh
            transforms = _raycore_toric_transforms(prepared.geometry; toricity=toricity)
            count = _raycore_foreach_mesh_chunk_from_geometry(
                mesh -> (toricity ? push!(tlas, mesh, transforms) : push!(tlas, mesh)),
                prepared.geometry;
                toricity=toricity,
                face_chunk_limit=face_chunk_limit,
            )
            count, length(prepared.geometry.node_ids) + 1, UInt32[]
        end
    actual_instance_count = length(tlas.instances)
    if isempty(decoder_table)
        instance_count = actual_instance_count
    elseif actual_instance_count != instance_count
        error(
            "Raycore TLAS instance count mismatch: decoder has $instance_count rows but TLAS has $actual_instance_count instances.",
        )
    end
    Raycore.sync!(tlas)
    kernel_tlas = Adapt.adapt(config.backend, tlas)
    hit_decoder =
        isempty(decoder_table) ?
        _raycore_identity_hit_decoder(prepared, config.backend, instance_count) :
        _raycore_hit_decoder(decoder_table, metadata_stride, instance_count, config.backend)
    n_pixels = prepared.geometry.plotbox.nx * prepared.geometry.plotbox.ny
    stack_len = n_pixels * config.max_hits_per_pixel
    top_nodes_dev = KernelAbstractions.allocate(config.backend, UInt32, n_pixels)
    stack_counts_dev = KernelAbstractions.allocate(config.backend, Int32, n_pixels)
    stack_nodes_dev = KernelAbstractions.allocate(config.backend, UInt32, stack_len)
    stack_instance_indices_dev = KernelAbstractions.allocate(config.backend, UInt32, stack_len)
    stack_heights_dev = KernelAbstractions.allocate(config.backend, Float32, stack_len)
    stack_overflow_dev = KernelAbstractions.allocate(config.backend, Bool, n_pixels)
    edge_key_capacity = n_pixels * max(0, 2 * (config.max_hits_per_pixel - 1))
    n_nodes = length(prepared.geometry.node_ids)
    dense_pairs_i128 = Int128(n_nodes) * Int128(n_nodes)
    dense_bytes_i128 = dense_pairs_i128 * Int128(sizeof(Int32))
    dense_matrix_available =
        dense_pairs_i128 <= typemax(Int) &&
        dense_bytes_i128 <= config.dense_edge_limit_bytes &&
        KernelAbstractions.supports_atomics(config.backend)
    auto_dense_enabled = config.edge_accumulation == :auto && !chunked_tlas && dense_matrix_available
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
        dense_scratch_enabled ? KernelAbstractions.allocate(config.backend, Int32, Int(dense_pairs_i128)) : nothing
    dense_edge_counts_host = Int32[]
    node_counts_dev =
        KernelAbstractions.supports_atomics(config.backend) ?
        KernelAbstractions.allocate(config.backend, Int32, n_nodes) :
        nothing
    node_transparency_dev = KernelAbstractions.allocate(config.backend, Float32, n_nodes)
    sector_area_dev =
        KernelAbstractions.supports_atomics(config.backend) ?
        KernelAbstractions.allocate(config.backend, Float32, n_nodes) :
        nothing
    node_counts_host = Int32[]
    sector_area_host = Float32[]
    virtual_node_mask_dev = KernelAbstractions.allocate(config.backend, Bool, length(prepared.virtual_node_mask))
    pavement_node_mask_dev = KernelAbstractions.allocate(config.backend, Bool, length(prepared.geometry.pavement_node_mask))
    node_ids_dev = KernelAbstractions.allocate(config.backend, Int, length(prepared.geometry.node_ids))
    KernelAbstractions.copyto!(config.backend, node_transparency_dev, Float32.(prepared.node_transparency_by_index))
    KernelAbstractions.copyto!(config.backend, virtual_node_mask_dev, prepared.virtual_node_mask)
    KernelAbstractions.copyto!(config.backend, pavement_node_mask_dev, prepared.geometry.pavement_node_mask)
    KernelAbstractions.copyto!(config.backend, node_ids_dev, prepared.geometry.node_ids)
    top_nodes_host = UInt32[]
    stack_counts_host = Int32[]
    stack_nodes_host = UInt32[]
    stack_instance_indices_host = UInt32[]
    stack_heights_host = Float32[]
    stack_overflow_host = Bool[]
    edge_keys_host = UInt64[]
    edge_key_counts_host = Int32[]
    edge_compact_host = UInt64[]
    vertical_span = _raycore_geometry_vertical_span(prepared.geometry)
    return RaycoreSceneData(
        prepared,
        tlas,
        kernel_tlas,
        config.backend,
        top_nodes_dev,
        top_nodes_host,
        stack_counts_dev,
        stack_nodes_dev,
        stack_instance_indices_dev,
        stack_heights_dev,
        stack_overflow_dev,
        edge_keys_dev,
        edge_key_counts_dev,
        dense_edge_counts_dev,
        dense_edge_counts_host,
        node_counts_dev,
        node_transparency_dev,
        sector_area_dev,
        node_counts_host,
        sector_area_host,
        virtual_node_mask_dev,
        pavement_node_mask_dev,
        node_ids_dev,
        stack_counts_host,
        stack_nodes_host,
        stack_instance_indices_host,
        stack_heights_host,
        stack_overflow_host,
        edge_keys_host,
        edge_key_counts_host,
        edge_compact_host,
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
        chunked_tlas,
        geometry_mode,
    )
end

function _prepare_raycore_interception_data(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions,
    backend::RaycoreInterceptionBackend;
    include_budget_maps::Bool=false,
)
    prepared = _prepare_interception_data(
        scene,
        models,
        options;
        include_budget_maps=include_budget_maps,
        include_raycore_instancing=true,
    )
    data, _chunked = _raycore_initial_scene_data(prepared, backend.config, options; toricity=options.toricity)
    return data
end

function _raycore_static_tlas(data::RaycoreSceneData)
    return data.tlas.static_tlas
end

@inline _raycore_kernel_tlas(data::RaycoreSceneData) = data.kernel_tlas

function _raycore_all_hits_available()
    return isdefined(Raycore, Symbol("all_hits!"))
end

function _raycore_require_all_hits!(context::AbstractString)
    _raycore_all_hits_available() && return nothing
    error(
        "$context requires Raycore.all_hits!, which is not available in the loaded Raycore version. " *
        "Use the Raycore branch or release that includes the all_hits! TLAS traversal API.",
    )
end

function _raycore_closest_hit_node_index(data::RaycoreSceneData, origin, direction)
    ray = Raycore.Ray(; o=_raycore_point3f(origin), d=_raycore_vec3f(direction))
    hit, tri, _distance, _bary, instance_idx = Raycore.closest_hit(_raycore_static_tlas(data), ray)
    return hit ? Int(_raycore_decode_node_index(data.hit_decoder, UInt32(tri.metadata), UInt32(instance_idx))) : 0
end

KernelAbstractions.@kernel function _raycore_trace_top_hits_kernel!(
    nodes,
    heights,
    distances,
    instance_indices,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origins,
    directions,
    t_mins,
    t_maxs,
)
    i = @index(Global, Linear)
    @inbounds begin
        ray = Raycore.Ray(;
            o=origins[i],
            d=directions[i],
            t_min=t_mins[i],
            t_max=t_maxs[i],
        )
        hit, tri, distance, _bary, instance_idx = Raycore.closest_hit(static_tlas, ray)
        if hit
            instance_indices[i] = UInt32(instance_idx)
            nodes[i] = _raycore_decode_node_index(
                node_index_by_instance_metadata,
                metadata_stride,
                UInt32(tri.metadata),
                UInt32(instance_idx),
            )
            hit_point = ray.o + ray.d * distance
            heights[i] = hit_point[3]
            distances[i] = distance
        else
            nodes[i] = UInt32(0)
            instance_indices[i] = UInt32(0)
            heights[i] = -Inf32
            distances[i] = Inf32
        end
    end
end

KernelAbstractions.@kernel function _raycore_trace_direction_top_hits_kernel!(
    nodes,
    heights,
    distances,
    instance_indices,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    dir_x::Float32,
    dir_y::Float32,
    dir_z::Float32,
    far::Float32,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        zero_idx = pixel_idx - 1
        ix = zero_idx % nx
        iy = zero_idx ÷ nx
        ground_x = origin_x + (Float32(ix) + 0.5f0) * pix_x
        ground_y = origin_y + (Float32(iy) + 0.5f0) * pix_y
        ray = Raycore.Ray(;
            o=GeometryBasics.Point3f(
                ground_x + dir_x * far,
                ground_y + dir_y * far,
                dir_z * far,
            ),
            d=GeometryBasics.Vec3f(-dir_x, -dir_y, -dir_z),
            t_min=0.0f0,
            t_max=2.0f0 * far,
        )
        hit, tri, distance, _bary, instance_idx = Raycore.closest_hit(static_tlas, ray)
        if hit
            instance_indices[pixel_idx] = UInt32(instance_idx)
            nodes[pixel_idx] = _raycore_decode_node_index(
                node_index_by_instance_metadata,
                metadata_stride,
                UInt32(tri.metadata),
                UInt32(instance_idx),
            )
            hit_point = ray.o + ray.d * distance
            heights[pixel_idx] = hit_point[3]
            distances[pixel_idx] = distance
        else
            nodes[pixel_idx] = UInt32(0)
            instance_indices[pixel_idx] = UInt32(0)
            heights[pixel_idx] = -Inf32
            distances[pixel_idx] = Inf32
        end
    end
end

KernelAbstractions.@kernel function _raycore_trace_direction_top_nodes_kernel!(
    nodes,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    dir_x::Float32,
    dir_y::Float32,
    dir_z::Float32,
    far::Float32,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        zero_idx = pixel_idx - 1
        ix = zero_idx % nx
        iy = zero_idx ÷ nx
        ground_x = origin_x + (Float32(ix) + 0.5f0) * pix_x
        ground_y = origin_y + (Float32(iy) + 0.5f0) * pix_y
        ray = Raycore.Ray(;
            o=GeometryBasics.Point3f(
                ground_x + dir_x * far,
                ground_y + dir_y * far,
                dir_z * far,
            ),
            d=GeometryBasics.Vec3f(-dir_x, -dir_y, -dir_z),
            t_min=0.0f0,
            t_max=2.0f0 * far,
        )
        hit, tri, _distance, _bary, instance_idx = Raycore.closest_hit(static_tlas, ray)
        nodes[pixel_idx] = hit ? _raycore_decode_node_index(
            node_index_by_instance_metadata,
            metadata_stride,
            UInt32(tri.metadata),
            UInt32(instance_idx),
        ) : UInt32(0)
    end
end

function _raycore_trace_top_hits(
    data::RaycoreSceneData,
    origins::AbstractVector,
    directions::AbstractVector;
    t_mins=nothing,
    t_maxs=nothing,
)
    n = length(origins)
    length(directions) == n || throw(ArgumentError("Raycore trace origins and directions must have the same length."))
    t_mins === nothing || length(t_mins) == n ||
        throw(ArgumentError("Raycore trace t_mins must have length $n."))
    t_maxs === nothing || length(t_maxs) == n ||
        throw(ArgumentError("Raycore trace t_maxs must have length $n."))

    backend = data.backend
    origins_host = GeometryBasics.Point3f[_raycore_point3f(origin) for origin in origins]
    directions_host = GeometryBasics.Vec3f[_raycore_vec3f(direction) for direction in directions]
    t_mins_host = t_mins === nothing ? fill(0.0f0, n) : Float32.(t_mins)
    t_maxs_host = t_maxs === nothing ? fill(Inf32, n) : Float32.(t_maxs)

    static_tlas = _raycore_kernel_tlas(data)

    origins_dev = KernelAbstractions.allocate(backend, GeometryBasics.Point3f, n)
    directions_dev = KernelAbstractions.allocate(backend, GeometryBasics.Vec3f, n)
    t_mins_dev = KernelAbstractions.allocate(backend, Float32, n)
    t_maxs_dev = KernelAbstractions.allocate(backend, Float32, n)
    nodes_dev = KernelAbstractions.allocate(backend, UInt32, n)
    heights_dev = KernelAbstractions.allocate(backend, Float32, n)
    distances_dev = KernelAbstractions.allocate(backend, Float32, n)
    instance_indices_dev = KernelAbstractions.allocate(backend, UInt32, n)

    KernelAbstractions.copyto!(backend, origins_dev, origins_host)
    KernelAbstractions.copyto!(backend, directions_dev, directions_host)
    KernelAbstractions.copyto!(backend, t_mins_dev, t_mins_host)
    KernelAbstractions.copyto!(backend, t_maxs_dev, t_maxs_host)

    kernel = _raycore_trace_top_hits_kernel!(backend, data.workgroupsize)
    kernel(
        nodes_dev,
        heights_dev,
        distances_dev,
        instance_indices_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        origins_dev,
        directions_dev,
        t_mins_dev,
        t_maxs_dev;
        ndrange=n,
    )
    KernelAbstractions.synchronize(backend)

    return (
        nodes=Array(nodes_dev),
        heights=Array(heights_dev),
        distances=Array(distances_dev),
        instance_indices=Array(instance_indices_dev),
    )
end

function _raycore_trace_direction_top_hits_direct(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)

    nodes_dev = KernelAbstractions.allocate(backend, UInt32, n_pixels)
    heights_dev = KernelAbstractions.allocate(backend, Float32, n_pixels)
    distances_dev = KernelAbstractions.allocate(backend, Float32, n_pixels)
    instance_indices_dev = KernelAbstractions.allocate(backend, UInt32, n_pixels)

    kernel = _raycore_trace_direction_top_hits_kernel!(backend, data.workgroupsize)
    kernel(
        nodes_dev,
        heights_dev,
        distances_dev,
        instance_indices_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)

    return (
        nodes=Array(nodes_dev),
        heights=Array(heights_dev),
        distances=Array(distances_dev),
        instance_indices=Array(instance_indices_dev),
    )
end

function _raycore_trace_direction_top_nodes_direct(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)

    nodes_dev = data.top_nodes_dev
    kernel = _raycore_trace_direction_top_nodes_kernel!(backend, data.workgroupsize)
    kernel(
        nodes_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)
    length(data.top_nodes_host) == n_pixels || resize!(data.top_nodes_host, n_pixels)
    copyto!(data.top_nodes_host, nodes_dev)
    return data.top_nodes_host
end

function _raycore_backend_requires_trace_validation(data::RaycoreSceneData)
    return data.validate && !(data.backend isa KernelAbstractions.CPU)
end

function _raycore_validation_message(reason::Symbol, validation)
    reason == :raycore_prechunk_instance_cap && return _raycore_prechunk_instance_limit_message(validation)
    reason == :raycore_trace_validation && return _raycore_trace_validation_message(validation)
    reason == :raycore_stack_trace_validation && return _raycore_stack_trace_validation_message(validation)
    return "Raycore backend failed validation for reason $reason."
end

function _raycore_throw_validation_error(
    config::RaycoreBackendConfig,
    reason::Symbol,
    stage::Symbol,
    validation,
)
    throw(RaycoreValidationError(reason, stage, _raycore_validation_message(reason, validation), validation))
end

function _raycore_reference_top_pixel_count(
    prepared::PreparedInterceptionData,
    direction,
    options::LightOptions,
)
    geometry = prepared.geometry
    projection = _direction_projection_cached(
        geometry,
        direction,
        options,
        prepared.cache_ctx;
        upper_hit=true,
        strict_java_float=_strict_java_float(options),
    )
    return count(_ -> true, values(projection.pixel_hits))
end

function _raycore_vertical_trace_validation(
    data::RaycoreSceneData,
    options::LightOptions;
    min_reference_pixels::Int=1024,
    min_reference_fraction::Float64=0.05,
    min_ratio::Float64=0.95,
)
    _raycore_backend_requires_trace_validation(data) || return (
        ok=true,
        required=false,
        reference_pixels=0,
        raycore_pixels=0,
        ratio=1.0,
        n_pixels=data.prepared.geometry.plotbox.nx * data.prepared.geometry.plotbox.ny,
    )

    direction = (0.0, 0.0, 1.0)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    pixel_area = plotbox.pixel_area
    pixel_area > 0.0 || return (
        ok=true,
        required=true,
        reference_pixels=0,
        raycore_pixels=0,
        ratio=1.0,
        n_pixels=n_pixels,
    )

    reference_pixels = _raycore_reference_top_pixel_count(data.prepared, direction, options)
    required_reference = max(min_reference_pixels, ceil(Int, min_reference_fraction * n_pixels))
    if reference_pixels < required_reference
        return (
            ok=true,
            required=true,
            reference_pixels=reference_pixels,
            raycore_pixels=0,
            ratio=1.0,
            n_pixels=n_pixels,
        )
    end

    raycore_nodes = _raycore_trace_direction_top_nodes_direct(data, direction, options)
    raycore_pixels = count(!=(0), raycore_nodes)
    ratio = raycore_pixels / max(reference_pixels, 1)
    return (
        ok=ratio >= min_ratio,
        required=true,
        reference_pixels=reference_pixels,
        raycore_pixels=raycore_pixels,
        ratio=ratio,
        n_pixels=n_pixels,
    )
end

function _raycore_trace_validation_message(validation)
    return "Raycore backend failed vertical trace validation " *
           "(raycore_pixels=$(validation.raycore_pixels), " *
           "reference_pixels=$(validation.reference_pixels), " *
           "ratio=$(round(validation.ratio; digits=4)), " *
           "n_pixels=$(validation.n_pixels))."
end

function _raycore_env_int(name::AbstractString, default::Int)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    try
        return parse(Int, raw)
    catch
        error("$name must be an integer, got $(repr(raw)).")
    end
end

function _raycore_env_float(name::AbstractString, default::Float64)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    try
        return parse(Float64, raw)
    catch
        error("$name must be a number, got $(repr(raw)).")
    end
end

function _raycore_stack_validation_min_reference_hits()
    return max(1, _raycore_env_int("ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_HITS", 512))
end

function _raycore_stack_validation_min_reference_occupied()
    return max(1, _raycore_env_int("ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_OCCUPIED", 128))
end

function _raycore_stack_validation_min_hit_ratio()
    return clamp(_raycore_env_float("ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_HIT_RATIO", 0.95), 0.0, 1.0)
end

function _raycore_stack_validation_min_occupied_ratio()
    return clamp(_raycore_env_float("ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_MIN_OCCUPIED_RATIO", 0.95), 0.0, 1.0)
end

function _raycore_stack_validation_direction_count()
    return clamp(_raycore_env_int("ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS", 3), 1, 3)
end

function _raycore_stack_validation_candidate_directions()
    return (
        (0.0, 0.0, 1.0),
        (0.5, 0.0, 0.8660254037844386),
        (-0.35355339059327373, 0.35355339059327373, 0.8660254037844386),
    )
end

function _raycore_stack_validation_directions()
    candidates = _raycore_stack_validation_candidate_directions()
    return candidates[1:_raycore_stack_validation_direction_count()]
end

function _raycore_stack_validation_cpu_data(data::RaycoreSceneData, options::LightOptions)
    config = RaycoreBackendConfig(
        ;
        backend=KernelAbstractions.CPU(),
        workgroupsize=data.workgroupsize,
        max_hits_per_pixel=data.max_hits_per_pixel,
        hit_epsilon=data.hit_epsilon,
        edge_accumulation=data.edge_accumulation,
        dense_edge_limit_bytes=data.dense_edge_limit_bytes,
        max_prechunk_instances=0,
        validate=false,
    )
    return _raycore_scene_data(data.prepared, config; toricity=options.toricity)
end

function _raycore_stack_trace_summary(data::RaycoreSceneData, direction, options::LightOptions)
    traced = _raycore_trace_direction_stack_nodes_device(data, direction, options)
    host = _raycore_copy_direction_stack_nodes_host!(data, traced)
    total_hits = 0
    occupied_pixels = 0
    max_seen = 0
    @inbounds for count in host.counts
        c = Int(count)
        total_hits += c
        occupied_pixels += c > 0
        max_seen = max(max_seen, c)
    end
    return (
        total_hits=total_hits,
        occupied_pixels=occupied_pixels,
        max_seen=max_seen,
        overflow=any(host.overflow),
    )
end

function _raycore_stack_trace_validation(
    data::RaycoreSceneData,
    options::LightOptions;
    min_reference_hits::Int=_raycore_stack_validation_min_reference_hits(),
    min_reference_occupied::Int=_raycore_stack_validation_min_reference_occupied(),
    min_hit_ratio::Float64=_raycore_stack_validation_min_hit_ratio(),
    min_occupied_ratio::Float64=_raycore_stack_validation_min_occupied_ratio(),
)
    _raycore_backend_requires_trace_validation(data) || return (
        ok=true,
        required=false,
        directions_tested=0,
        direction_count=0,
        min_reference_hits=min_reference_hits,
        min_reference_occupied=min_reference_occupied,
        min_hit_ratio=min_hit_ratio,
        min_occupied_ratio=min_occupied_ratio,
        reference_hits=0,
        raycore_hits=0,
        hit_ratio=1.0,
        reference_occupied=0,
        raycore_occupied=0,
        occupied_ratio=1.0,
        reference_overflow=false,
        raycore_overflow=false,
    )

    directions = _raycore_stack_validation_directions()
    reference = _raycore_stack_validation_cpu_data(data, options)
    directions_tested = 0
    reference_hits = 0
    raycore_hits = 0
    reference_occupied = 0
    raycore_occupied = 0
    reference_overflow = false
    raycore_overflow = false

    for direction in directions
        ref = _raycore_stack_trace_summary(reference, direction, options)
        ref.total_hits >= min_reference_hits || ref.occupied_pixels >= min_reference_occupied || continue

        got = _raycore_stack_trace_summary(data, direction, options)
        directions_tested += 1
        reference_hits += ref.total_hits
        raycore_hits += got.total_hits
        reference_occupied += ref.occupied_pixels
        raycore_occupied += got.occupied_pixels
        reference_overflow |= ref.overflow
        raycore_overflow |= got.overflow
    end

    directions_tested == 0 && return (
        ok=true,
        required=true,
        directions_tested=0,
        direction_count=length(directions),
        min_reference_hits=min_reference_hits,
        min_reference_occupied=min_reference_occupied,
        min_hit_ratio=min_hit_ratio,
        min_occupied_ratio=min_occupied_ratio,
        reference_hits=reference_hits,
        raycore_hits=raycore_hits,
        hit_ratio=1.0,
        reference_occupied=reference_occupied,
        raycore_occupied=raycore_occupied,
        occupied_ratio=1.0,
        reference_overflow=reference_overflow,
        raycore_overflow=raycore_overflow,
    )

    hit_ratio = raycore_hits / max(reference_hits, 1)
    occupied_ratio = raycore_occupied / max(reference_occupied, 1)
    return (
        ok=hit_ratio >= min_hit_ratio && occupied_ratio >= min_occupied_ratio,
        required=true,
        directions_tested=directions_tested,
        direction_count=length(directions),
        min_reference_hits=min_reference_hits,
        min_reference_occupied=min_reference_occupied,
        min_hit_ratio=min_hit_ratio,
        min_occupied_ratio=min_occupied_ratio,
        reference_hits=reference_hits,
        raycore_hits=raycore_hits,
        hit_ratio=hit_ratio,
        reference_occupied=reference_occupied,
        raycore_occupied=raycore_occupied,
        occupied_ratio=occupied_ratio,
        reference_overflow=reference_overflow,
        raycore_overflow=raycore_overflow,
    )
end

function _raycore_stack_trace_validation_message(validation)
    return "Raycore backend failed full-stack trace validation " *
           "(raycore_hits=$(validation.raycore_hits), " *
           "reference_hits=$(validation.reference_hits), " *
           "hit_ratio=$(round(validation.hit_ratio; digits=4)), " *
           "raycore_occupied=$(validation.raycore_occupied), " *
           "reference_occupied=$(validation.reference_occupied), " *
           "occupied_ratio=$(round(validation.occupied_ratio; digits=4)), " *
           "min_hit_ratio=$(validation.min_hit_ratio), " *
           "min_occupied_ratio=$(validation.min_occupied_ratio), " *
           "directions_tested=$(validation.directions_tested), " *
           "direction_count=$(validation.direction_count))."
end

function _raycore_retry_chunk_limits(skip_face_chunk_limit::Union{Nothing,Int}=nothing)
    face_chunk_limit = _raycore_face_chunk_limit()
    face_chunk_limit == typemax(Int) && return Int[]

    candidate_limits = Int[]
    for limit in (face_chunk_limit, min(face_chunk_limit, 1024), min(face_chunk_limit, 768), min(face_chunk_limit, 512))
        limit > 0 || continue
        limit == skip_face_chunk_limit && continue
        limit in candidate_limits && continue
        push!(candidate_limits, limit)
    end
    return candidate_limits
end

function _raycore_retry_chunked_scene_data(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig,
    options::LightOptions,
    validation;
    toricity::Bool=options.toricity,
    skip_face_chunk_limit::Union{Nothing,Int}=nothing,
)
    config.backend isa KernelAbstractions.CPU && return nothing, validation
    candidate_limits = _raycore_retry_chunk_limits(skip_face_chunk_limit)

    last_validation = validation
    for limit in candidate_limits
        estimated_instances = _raycore_estimated_prechunk_instances(prepared; toricity=toricity, face_chunk_limit=limit)
        if config.max_prechunk_instances != typemax(Int) && estimated_instances > config.max_prechunk_instances
            @info "Skipping Raycore validation retry because the chunk limit would exceed max_prechunk_instances." face_chunk_limit=limit estimated_instances max_prechunk_instances=config.max_prechunk_instances
            continue
        end
        chunked = _raycore_scene_data(prepared, config; toricity=toricity, face_chunk_limit=limit)
        chunked_validation = _raycore_vertical_trace_validation(chunked, options)
        last_validation = chunked_validation
        if chunked_validation.ok
            @info "Raycore backend passed vertical trace validation after chunking BLAS geometry." face_chunk_limit=limit
            return chunked, chunked_validation
        end
    end
    return nothing, last_validation
end

function _raycore_retry_stack_chunked_scene_data(
    prepared::PreparedInterceptionData,
    config::RaycoreBackendConfig,
    options::LightOptions,
    validation;
    toricity::Bool=options.toricity,
    skip_face_chunk_limit::Union{Nothing,Int}=nothing,
)
    config.backend isa KernelAbstractions.CPU && return nothing, validation
    candidate_limits = _raycore_retry_chunk_limits(skip_face_chunk_limit)

    last_validation = validation
    for limit in candidate_limits
        estimated_instances = _raycore_estimated_prechunk_instances(prepared; toricity=toricity, face_chunk_limit=limit)
        if config.max_prechunk_instances != typemax(Int) && estimated_instances > config.max_prechunk_instances
            @info "Skipping Raycore full-stack validation retry because the chunk limit would exceed max_prechunk_instances." face_chunk_limit=limit estimated_instances max_prechunk_instances=config.max_prechunk_instances
            continue
        end
        chunked = _raycore_scene_data(prepared, config; toricity=toricity, face_chunk_limit=limit)
        chunked_validation = _raycore_stack_trace_validation(chunked, options)
        last_validation = chunked_validation
        if chunked_validation.ok
            @info "Raycore backend passed full-stack trace validation after chunking BLAS geometry." face_chunk_limit=limit hit_ratio=chunked_validation.hit_ratio occupied_ratio=chunked_validation.occupied_ratio
            return chunked, chunked_validation
        end
    end
    return nothing, last_validation
end

KernelAbstractions.@kernel function _raycore_trace_stacks_kernel!(
    counts,
    nodes,
    heights,
    instance_indices,
    overflow,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origins,
    directions,
    t_mins,
    t_maxs,
    max_hits::Int,
    hit_epsilon::Float32,
)
    i = @index(Global, Linear)
    @inbounds begin
        base_origin = origins[i]
        base_direction = directions[i]
        t_min = t_mins[i]
        t_max = t_maxs[i]
        ray = Raycore.Ray(;
            o=base_origin,
            d=base_direction,
            t_min=t_min,
            t_max=t_max,
        )
        offset = (i - 1) * max_hits
        n_unique, overflow_hit = Raycore.all_hits!(
            nodes,
            heights,
            instance_indices,
            static_tlas,
            ray,
            offset,
            max_hits,
            hit_epsilon,
        )
        for slot in 1:Int(n_unique)
            out_idx = offset + slot
            distance = heights[out_idx]
            nodes[out_idx] = _raycore_decode_node_index(
                node_index_by_instance_metadata,
                metadata_stride,
                nodes[out_idx],
                instance_indices[out_idx],
            )
            hit_point = ray.o + ray.d * distance
            heights[out_idx] = hit_point[3]
        end
        counts[i] = n_unique
        overflow[i] = overflow_hit
    end
end

KernelAbstractions.@kernel function _raycore_trace_direction_stacks_kernel!(
    counts,
    nodes,
    heights,
    instance_indices,
    overflow,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    dir_x::Float32,
    dir_y::Float32,
    dir_z::Float32,
    far::Float32,
    max_hits::Int,
    hit_epsilon::Float32,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        zero_idx = pixel_idx - 1
        ix = zero_idx % nx
        iy = zero_idx ÷ nx
        ground_x = origin_x + (Float32(ix) + 0.5f0) * pix_x
        ground_y = origin_y + (Float32(iy) + 0.5f0) * pix_y
        base_origin = GeometryBasics.Point3f(
            ground_x + dir_x * far,
            ground_y + dir_y * far,
            dir_z * far,
        )
        base_direction = GeometryBasics.Vec3f(-dir_x, -dir_y, -dir_z)
        t_max = 2.0f0 * far
        ray = Raycore.Ray(;
            o=base_origin,
            d=base_direction,
            t_min=0.0f0,
            t_max=t_max,
        )
        offset = (pixel_idx - 1) * max_hits
        n_unique, overflow_hit = Raycore.all_hits!(
            nodes,
            heights,
            instance_indices,
            static_tlas,
            ray,
            offset,
            max_hits,
            hit_epsilon,
        )
        for slot in 1:Int(n_unique)
            out_idx = offset + slot
            distance = heights[out_idx]
            nodes[out_idx] = _raycore_decode_node_index(
                node_index_by_instance_metadata,
                metadata_stride,
                nodes[out_idx],
                instance_indices[out_idx],
            )
            hit_point = ray.o + ray.d * distance
            heights[out_idx] = hit_point[3]
        end
        counts[pixel_idx] = n_unique
        overflow[pixel_idx] = overflow_hit
    end
end

KernelAbstractions.@kernel function _raycore_trace_direction_stack_nodes_kernel!(
    counts,
    nodes,
    distances,
    instance_indices,
    overflow,
    node_index_by_instance_metadata,
    metadata_stride::Int,
    static_tlas,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    dir_x::Float32,
    dir_y::Float32,
    dir_z::Float32,
    far::Float32,
    max_hits::Int,
    hit_epsilon::Float32,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        zero_idx = pixel_idx - 1
        ix = zero_idx % nx
        iy = zero_idx ÷ nx
        ground_x = origin_x + (Float32(ix) + 0.5f0) * pix_x
        ground_y = origin_y + (Float32(iy) + 0.5f0) * pix_y
        base_origin = GeometryBasics.Point3f(
            ground_x + dir_x * far,
            ground_y + dir_y * far,
            dir_z * far,
        )
        base_direction = GeometryBasics.Vec3f(-dir_x, -dir_y, -dir_z)
        t_max = 2.0f0 * far
        ray = Raycore.Ray(;
            o=base_origin,
            d=base_direction,
            t_min=0.0f0,
            t_max=t_max,
        )
        offset = (pixel_idx - 1) * max_hits
        n_unique, overflow_hit = Raycore.all_hits!(
            nodes,
            distances,
            instance_indices,
            static_tlas,
            ray,
            offset,
            max_hits,
            hit_epsilon,
        )
        for slot in 1:Int(n_unique)
            out_idx = offset + slot
            nodes[out_idx] = _raycore_decode_node_index(
                node_index_by_instance_metadata,
                metadata_stride,
                nodes[out_idx],
                instance_indices[out_idx],
            )
        end
        counts[pixel_idx] = n_unique
        overflow[pixel_idx] = overflow_hit
    end
end

KernelAbstractions.@kernel function _raycore_trace_direction_raw_stacks_kernel!(
    counts,
    metadata,
    distances,
    instance_indices,
    overflow,
    static_tlas,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    dir_x::Float32,
    dir_y::Float32,
    dir_z::Float32,
    far::Float32,
    max_hits::Int,
    hit_epsilon::Float32,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        zero_idx = pixel_idx - 1
        ix = zero_idx % nx
        iy = zero_idx ÷ nx
        ground_x = origin_x + (Float32(ix) + 0.5f0) * pix_x
        ground_y = origin_y + (Float32(iy) + 0.5f0) * pix_y
        base_origin = GeometryBasics.Point3f(
            ground_x + dir_x * far,
            ground_y + dir_y * far,
            dir_z * far,
        )
        base_direction = GeometryBasics.Vec3f(-dir_x, -dir_y, -dir_z)
        ray = Raycore.Ray(;
            o=base_origin,
            d=base_direction,
            t_min=0.0f0,
            t_max=2.0f0 * far,
        )
        offset = (pixel_idx - 1) * max_hits
        n_unique, overflow_hit = Raycore.all_hits!(
            metadata,
            distances,
            instance_indices,
            static_tlas,
            ray,
            offset,
            max_hits,
            hit_epsilon,
        )
        counts[pixel_idx] = n_unique
        overflow[pixel_idx] = overflow_hit
    end
end

function _raycore_trace_stacks(
    data::RaycoreSceneData,
    origins::AbstractVector,
    directions::AbstractVector;
    t_mins=nothing,
    t_maxs=nothing,
)
    _raycore_require_all_hits!("Raycore stack tracing")
    n = length(origins)
    length(directions) == n || throw(ArgumentError("Raycore stack trace origins and directions must have the same length."))
    t_mins === nothing || length(t_mins) == n ||
        throw(ArgumentError("Raycore stack trace t_mins must have length $n."))
    t_maxs === nothing || length(t_maxs) == n ||
        throw(ArgumentError("Raycore stack trace t_maxs must have length $n."))

    backend = data.backend
    max_hits = data.max_hits_per_pixel
    origins_host = GeometryBasics.Point3f[_raycore_point3f(origin) for origin in origins]
    directions_host = GeometryBasics.Vec3f[_raycore_vec3f(direction) for direction in directions]
    t_mins_host = t_mins === nothing ? fill(0.0f0, n) : Float32.(t_mins)
    t_maxs_host = t_maxs === nothing ? fill(Inf32, n) : Float32.(t_maxs)

    static_tlas = _raycore_kernel_tlas(data)

    origins_dev = KernelAbstractions.allocate(backend, GeometryBasics.Point3f, n)
    directions_dev = KernelAbstractions.allocate(backend, GeometryBasics.Vec3f, n)
    t_mins_dev = KernelAbstractions.allocate(backend, Float32, n)
    t_maxs_dev = KernelAbstractions.allocate(backend, Float32, n)
    counts_dev = KernelAbstractions.allocate(backend, Int32, n)
    nodes_dev = KernelAbstractions.allocate(backend, UInt32, n * max_hits)
    heights_dev = KernelAbstractions.allocate(backend, Float32, n * max_hits)
    instance_indices_dev = KernelAbstractions.allocate(backend, UInt32, n * max_hits)
    overflow_dev = KernelAbstractions.allocate(backend, Bool, n)

    KernelAbstractions.copyto!(backend, origins_dev, origins_host)
    KernelAbstractions.copyto!(backend, directions_dev, directions_host)
    KernelAbstractions.copyto!(backend, t_mins_dev, t_mins_host)
    KernelAbstractions.copyto!(backend, t_maxs_dev, t_maxs_host)

    kernel = _raycore_trace_stacks_kernel!(backend, data.workgroupsize)
    kernel(
        counts_dev,
        nodes_dev,
        heights_dev,
        instance_indices_dev,
        overflow_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        origins_dev,
        directions_dev,
        t_mins_dev,
        t_maxs_dev,
        max_hits,
        data.hit_epsilon;
        ndrange=n,
    )
    KernelAbstractions.synchronize(backend)

    return (
        counts=Array(counts_dev),
        nodes=Array(nodes_dev),
        instance_indices=Array(instance_indices_dev),
        heights=Array(heights_dev),
        overflow=Array(overflow_dev),
        max_hits=max_hits,
    )
end

function _raycore_trace_direction_stack_nodes(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    _raycore_require_all_hits!("Raycore direction stack tracing")
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    max_hits = data.max_hits_per_pixel
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)
    counts_dev = data.stack_counts_dev
    nodes_dev = data.stack_nodes_dev
    distances_dev = data.stack_heights_dev
    instance_indices_dev = data.stack_instance_indices_dev
    overflow_dev = data.stack_overflow_dev

    kernel = _raycore_trace_direction_stack_nodes_kernel!(backend, data.workgroupsize)
    kernel(
        counts_dev,
        nodes_dev,
        distances_dev,
        instance_indices_dev,
        overflow_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far,
        max_hits,
        data.hit_epsilon;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)
    length(data.stack_counts_host) == n_pixels || resize!(data.stack_counts_host, n_pixels)
    length(data.stack_nodes_host) == n_pixels * max_hits || resize!(data.stack_nodes_host, n_pixels * max_hits)
    length(data.stack_instance_indices_host) == n_pixels * max_hits ||
        resize!(data.stack_instance_indices_host, n_pixels * max_hits)
    length(data.stack_overflow_host) == n_pixels || resize!(data.stack_overflow_host, n_pixels)
    copyto!(data.stack_counts_host, counts_dev)
    copyto!(data.stack_nodes_host, nodes_dev)
    copyto!(data.stack_instance_indices_host, instance_indices_dev)
    copyto!(data.stack_overflow_host, overflow_dev)

    return (
        counts=data.stack_counts_host,
        nodes=data.stack_nodes_host,
        instance_indices=data.stack_instance_indices_host,
        overflow=data.stack_overflow_host,
        max_hits=max_hits,
    )
end

function _raycore_trace_direction_stack_nodes_device(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    _raycore_require_all_hits!("Raycore device direction stack tracing")
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    max_hits = data.max_hits_per_pixel
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)
    counts_dev = data.stack_counts_dev
    nodes_dev = data.stack_nodes_dev
    distances_dev = data.stack_heights_dev
    instance_indices_dev = data.stack_instance_indices_dev
    overflow_dev = data.stack_overflow_dev

    kernel = _raycore_trace_direction_stack_nodes_kernel!(backend, data.workgroupsize)
    kernel(
        counts_dev,
        nodes_dev,
        distances_dev,
        instance_indices_dev,
        overflow_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far,
        max_hits,
        data.hit_epsilon;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)
    length(data.stack_overflow_host) == n_pixels || resize!(data.stack_overflow_host, n_pixels)
    copyto!(data.stack_overflow_host, overflow_dev)

    return (
        counts_dev=counts_dev,
        nodes_dev=nodes_dev,
        instance_indices_dev=instance_indices_dev,
        overflow=data.stack_overflow_host,
        max_hits=max_hits,
    )
end

function _raycore_trace_direction_raw_stacks(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    _raycore_require_all_hits!("Raycore raw direction stack tracing")
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    max_hits = data.max_hits_per_pixel
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)
    counts_dev = data.stack_counts_dev
    metadata_dev = data.stack_nodes_dev
    distances_dev = data.stack_heights_dev
    instance_indices_dev = data.stack_instance_indices_dev
    overflow_dev = data.stack_overflow_dev

    kernel = _raycore_trace_direction_raw_stacks_kernel!(backend, data.workgroupsize)
    kernel(
        counts_dev,
        metadata_dev,
        distances_dev,
        instance_indices_dev,
        overflow_dev,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far,
        max_hits,
        data.hit_epsilon;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)

    length(data.stack_counts_host) == n_pixels || resize!(data.stack_counts_host, n_pixels)
    length(data.stack_nodes_host) == n_pixels * max_hits || resize!(data.stack_nodes_host, n_pixels * max_hits)
    length(data.stack_instance_indices_host) == n_pixels * max_hits ||
        resize!(data.stack_instance_indices_host, n_pixels * max_hits)
    length(data.stack_heights_host) == n_pixels * max_hits || resize!(data.stack_heights_host, n_pixels * max_hits)
    length(data.stack_overflow_host) == n_pixels || resize!(data.stack_overflow_host, n_pixels)
    copyto!(data.stack_counts_host, counts_dev)
    copyto!(data.stack_nodes_host, metadata_dev)
    copyto!(data.stack_instance_indices_host, instance_indices_dev)
    copyto!(data.stack_heights_host, distances_dev)
    copyto!(data.stack_overflow_host, overflow_dev)

    return (
        counts=data.stack_counts_host,
        metadata=data.stack_nodes_host,
        instance_indices=data.stack_instance_indices_host,
        distances=data.stack_heights_host,
        overflow=data.stack_overflow_host,
        max_hits=max_hits,
    )
end

@inline function _raycore_stack_slot(trace, pixel_idx::Int, slot::Int)
    return (pixel_idx - 1) * trace.max_hits + slot
end

@inline function _raycore_raw_hit(trace, pixel_idx::Int, slot::Int)
    idx = _raycore_stack_slot(trace, pixel_idx, slot)
    return (
        metadata=trace.metadata[idx],
        instance_index=trace.instance_indices[idx],
        distance=trace.distances[idx],
    )
end

@inline function _raycore_same_raw_hit(a, b, distance_atol::Float32)
    return a.metadata == b.metadata &&
           a.instance_index == b.instance_index &&
           isapprox(a.distance, b.distance; atol=distance_atol, rtol=0.0f0)
end

function _raycore_unmatched_raw_hits(reference, candidate, pixel_idx::Int, distance_atol::Float32)
    ref_count = Int(reference.counts[pixel_idx])
    got_count = Int(candidate.counts[pixel_idx])
    ref_hits = [_raycore_raw_hit(reference, pixel_idx, slot) for slot in 1:ref_count]
    got_hits = [_raycore_raw_hit(candidate, pixel_idx, slot) for slot in 1:got_count]
    matched_candidate = falses(length(got_hits))
    unmatched_reference = typeof(ref_hits)()
    for ref_hit in ref_hits
        matched = findfirst(eachindex(got_hits)) do got_idx
            !matched_candidate[got_idx] && _raycore_same_raw_hit(ref_hit, got_hits[got_idx], distance_atol)
        end
        if matched === nothing
            push!(unmatched_reference, ref_hit)
        else
            matched_candidate[matched] = true
        end
    end
    unmatched_candidate = got_hits[.!matched_candidate]
    return (
        unmatched_reference_count=length(unmatched_reference),
        unmatched_candidate_count=length(unmatched_candidate),
        first_unmatched_reference=isempty(unmatched_reference) ? nothing : first(unmatched_reference),
        first_unmatched_candidate=isempty(unmatched_candidate) ? nothing : first(unmatched_candidate),
    )
end

function _raycore_same_identity_membership(reference, candidate, pixel_idx::Int)
    ref_count = Int(reference.counts[pixel_idx])
    got_count = Int(candidate.counts[pixel_idx])
    ref_count == got_count || return false
    got_hits = [_raycore_raw_hit(candidate, pixel_idx, slot) for slot in 1:got_count]
    matched_candidate = falses(length(got_hits))
    for slot in 1:ref_count
        ref_hit = _raycore_raw_hit(reference, pixel_idx, slot)
        matched = findfirst(eachindex(got_hits)) do got_idx
            !matched_candidate[got_idx] &&
                ref_hit.metadata == got_hits[got_idx].metadata &&
                ref_hit.instance_index == got_hits[got_idx].instance_index
        end
        matched === nothing && return false
        matched_candidate[matched] = true
    end
    return true
end

function _raycore_same_metadata_membership(reference, candidate, pixel_idx::Int)
    ref_count = Int(reference.counts[pixel_idx])
    got_count = Int(candidate.counts[pixel_idx])
    ref_count == got_count || return false
    ref_metadata = [_raycore_raw_hit(reference, pixel_idx, slot).metadata for slot in 1:ref_count]
    got_metadata = [_raycore_raw_hit(candidate, pixel_idx, slot).metadata for slot in 1:got_count]
    sort!(ref_metadata)
    sort!(got_metadata)
    return ref_metadata == got_metadata
end

function _raycore_raw_mismatch_details(reference, candidate, pixel_idx::Int, distance_atol::Float32)
    unmatched = _raycore_unmatched_raw_hits(reference, candidate, pixel_idx, distance_atol)
    ref_count = Int(reference.counts[pixel_idx])
    got_count = Int(candidate.counts[pixel_idx])
    diagnosis =
        if unmatched.unmatched_reference_count == 0 && unmatched.unmatched_candidate_count == 0
            :raw_order_difference
        elseif ref_count > got_count && unmatched.unmatched_candidate_count == 0
            :candidate_missing_hits
        elseif got_count > ref_count && unmatched.unmatched_reference_count == 0
            :candidate_extra_hits
        elseif ref_count == got_count && _raycore_same_identity_membership(reference, candidate, pixel_idx)
            :distance_difference
        elseif ref_count == got_count && _raycore_same_metadata_membership(reference, candidate, pixel_idx)
            :instance_index_difference
        else
            :raw_stack_membership_difference
        end
    return (diagnosis=diagnosis, unmatched...)
end

function _raycore_first_raw_stack_mismatch(reference, candidate; distance_atol::Float32=1.0f-4)
    length(reference.counts) == length(candidate.counts) ||
        return (
            kind=:pixel_count,
            diagnosis=:pixel_count_difference,
            reference_pixels=length(reference.counts),
            candidate_pixels=length(candidate.counts),
        )
    reference.max_hits == candidate.max_hits ||
        return (
            kind=:max_hits,
            diagnosis=:max_hits_difference,
            reference_max_hits=reference.max_hits,
            candidate_max_hits=candidate.max_hits,
        )

    @inbounds for pixel_idx in eachindex(reference.counts)
        ref_count = Int(reference.counts[pixel_idx])
        got_count = Int(candidate.counts[pixel_idx])
        if ref_count != got_count
            details = _raycore_raw_mismatch_details(reference, candidate, pixel_idx, distance_atol)
            return (
                kind=:count,
                details...,
                pixel=pixel_idx,
                reference_count=ref_count,
                candidate_count=got_count,
                reference_overflow=reference.overflow[pixel_idx],
                candidate_overflow=candidate.overflow[pixel_idx],
            )
        end
        if reference.overflow[pixel_idx] != candidate.overflow[pixel_idx]
            return (
                kind=:overflow,
                diagnosis=:overflow_difference,
                pixel=pixel_idx,
                count=ref_count,
                reference_overflow=reference.overflow[pixel_idx],
                candidate_overflow=candidate.overflow[pixel_idx],
            )
        end
        for slot in 1:ref_count
            idx = _raycore_stack_slot(reference, pixel_idx, slot)
            ref_metadata = reference.metadata[idx]
            got_metadata = candidate.metadata[idx]
            if ref_metadata != got_metadata
                details = _raycore_raw_mismatch_details(reference, candidate, pixel_idx, distance_atol)
                return (
                    kind=:metadata,
                    details...,
                    pixel=pixel_idx,
                    slot=slot,
                    reference=ref_metadata,
                    candidate=got_metadata,
                    reference_instance=reference.instance_indices[idx],
                    candidate_instance=candidate.instance_indices[idx],
                    reference_distance=reference.distances[idx],
                    candidate_distance=candidate.distances[idx],
                )
            end
            ref_instance = reference.instance_indices[idx]
            got_instance = candidate.instance_indices[idx]
            if ref_instance != got_instance
                details = _raycore_raw_mismatch_details(reference, candidate, pixel_idx, distance_atol)
                return (
                    kind=:instance_index,
                    details...,
                    pixel=pixel_idx,
                    slot=slot,
                    metadata=ref_metadata,
                    reference=ref_instance,
                    candidate=got_instance,
                    reference_distance=reference.distances[idx],
                    candidate_distance=candidate.distances[idx],
                )
            end
            ref_distance = reference.distances[idx]
            got_distance = candidate.distances[idx]
            if !isapprox(ref_distance, got_distance; atol=distance_atol, rtol=0.0f0)
                details = _raycore_raw_mismatch_details(reference, candidate, pixel_idx, distance_atol)
                return (
                    kind=:distance,
                    details...,
                    pixel=pixel_idx,
                    slot=slot,
                    metadata=ref_metadata,
                    instance_index=ref_instance,
                    reference=ref_distance,
                    candidate=got_distance,
                    distance_atol=distance_atol,
                )
            end
        end
    end
    return nothing
end

function _raycore_raw_stack_comparison(
    reference::RaycoreSceneData,
    candidate::RaycoreSceneData,
    direction,
    options::LightOptions;
    distance_atol::Real=1.0f-4,
)
    ref = _raycore_trace_direction_raw_stacks(reference, direction, options)
    got = _raycore_trace_direction_raw_stacks(candidate, direction, options)
    ref_hits = sum(x -> Int(x), ref.counts; init=0)
    got_hits = sum(x -> Int(x), got.counts; init=0)
    ref_occupied = count(!=(Int32(0)), ref.counts)
    got_occupied = count(!=(Int32(0)), got.counts)
    mismatch = _raycore_first_raw_stack_mismatch(ref, got; distance_atol=Float32(distance_atol))
    return (
        ok=mismatch === nothing,
        diagnosis=mismatch === nothing ? :exact : mismatch.diagnosis,
        reference_hits=ref_hits,
        candidate_hits=got_hits,
        hit_ratio=ref_hits == 0 ? (got_hits == 0 ? 1.0 : Inf) : got_hits / ref_hits,
        reference_occupied=ref_occupied,
        candidate_occupied=got_occupied,
        occupied_ratio=ref_occupied == 0 ? (got_occupied == 0 ? 1.0 : Inf) : got_occupied / ref_occupied,
        reference_overflow=count(ref.overflow),
        candidate_overflow=count(got.overflow),
        reference_max_seen=isempty(ref.counts) ? 0 : maximum(Int.(ref.counts)),
        candidate_max_seen=isempty(got.counts) ? 0 : maximum(Int.(got.counts)),
        mismatch=mismatch,
    )
end

function _raycore_copy_direction_stack_nodes_host!(data::RaycoreSceneData, traced)
    n_pixels = length(traced.overflow)
    max_hits = traced.max_hits
    length(data.stack_counts_host) == n_pixels || resize!(data.stack_counts_host, n_pixels)
    length(data.stack_nodes_host) == n_pixels * max_hits || resize!(data.stack_nodes_host, n_pixels * max_hits)
    copyto!(data.stack_counts_host, traced.counts_dev)
    copyto!(data.stack_nodes_host, traced.nodes_dev)
    return (
        counts=data.stack_counts_host,
        nodes=data.stack_nodes_host,
        overflow=traced.overflow,
        max_hits=max_hits,
    )
end

function _raycore_trace_direction_stacks_direct(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    _raycore_require_all_hits!("Raycore full direction stack tracing")
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    backend = data.backend
    max_hits = data.max_hits_per_pixel
    far = _raycore_projection_far(data, direction, options)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))

    static_tlas = _raycore_kernel_tlas(data)

    counts_dev = data.stack_counts_dev
    nodes_dev = data.stack_nodes_dev
    heights_dev = data.stack_heights_dev
    instance_indices_dev = data.stack_instance_indices_dev
    overflow_dev = data.stack_overflow_dev

    kernel = _raycore_trace_direction_stacks_kernel!(backend, data.workgroupsize)
    kernel(
        counts_dev,
        nodes_dev,
        heights_dev,
        instance_indices_dev,
        overflow_dev,
        data.hit_decoder.node_index_by_instance_metadata_dev,
        data.hit_decoder.metadata_stride,
        static_tlas,
        Float32(plotbox.origin_x),
        Float32(plotbox.origin_y),
        Float32(plotbox.pix_x),
        Float32(plotbox.pix_y),
        plotbox.nx,
        Float32(dir[1]),
        Float32(dir[2]),
        Float32(dir[3]),
        far,
        max_hits,
        data.hit_epsilon;
        ndrange=n_pixels,
    )
    KernelAbstractions.synchronize(backend)
    length(data.stack_counts_host) == n_pixels || resize!(data.stack_counts_host, n_pixels)
    length(data.stack_nodes_host) == n_pixels * max_hits || resize!(data.stack_nodes_host, n_pixels * max_hits)
    length(data.stack_instance_indices_host) == n_pixels * max_hits ||
        resize!(data.stack_instance_indices_host, n_pixels * max_hits)
    length(data.stack_heights_host) == n_pixels * max_hits || resize!(data.stack_heights_host, n_pixels * max_hits)
    length(data.stack_overflow_host) == n_pixels || resize!(data.stack_overflow_host, n_pixels)
    copyto!(data.stack_counts_host, counts_dev)
    copyto!(data.stack_nodes_host, nodes_dev)
    copyto!(data.stack_instance_indices_host, instance_indices_dev)
    copyto!(data.stack_heights_host, heights_dev)
    copyto!(data.stack_overflow_host, overflow_dev)

    return (
        counts=data.stack_counts_host,
        nodes=data.stack_nodes_host,
        instance_indices=data.stack_instance_indices_host,
        heights=data.stack_heights_host,
        overflow=data.stack_overflow_host,
        max_hits=max_hits,
    )
end

function _raycore_projection_far(geometry::InterceptionSceneData, direction; toricity::Bool=false)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))
    pb = geometry.plotbox
    corners = (
        StaticArrays.SVector{3,Float64}(pb.origin_x, pb.origin_y, 0.0),
        StaticArrays.SVector{3,Float64}(pb.origin_x + pb.xdim, pb.origin_y, 0.0),
        StaticArrays.SVector{3,Float64}(pb.origin_x, pb.origin_y + pb.ydim, 0.0),
        StaticArrays.SVector{3,Float64}(pb.origin_x + pb.xdim, pb.origin_y + pb.ydim, 0.0),
    )
    max_distance = 0.0
    @inbounds for v in geometry.vertices
        for c in corners
            max_distance = max(max_distance, dot(v - c, dir))
        end
    end
    zs = (Float64(v[3]) for v in geometry.vertices)
    zmin, zmax = extrema(zs)
    toric_margin = toricity ? (abs(dir[1]) * pb.xdim + abs(dir[2]) * pb.ydim) : 0.0
    margin = sqrt(pb.xdim^2 + pb.ydim^2 + (zmax - zmin)^2) + toric_margin + max(pb.pix_x, pb.pix_y, 1e-3)
    return Float32(max_distance + margin)
end

@inline function _raycore_direction_key(direction)
    return (Float64(direction[1]), Float64(direction[2]), Float64(direction[3]))
end

function _raycore_projection_far(data::RaycoreSceneData, direction, options::LightOptions)
    key = (_raycore_direction_key(direction)..., options.toricity)
    return get!(data.projection_far_cache, key) do
        _raycore_projection_far(data.prepared.geometry, direction; toricity=options.toricity)
    end
end

function _raycore_raster_compat_projection(data::RaycoreSceneData, direction, options::LightOptions)
    prepared = data.prepared
    key = (_raycore_direction_key(direction)..., prepared.upper_hit)
    return get!(data.raster_compat_projection_cache, key) do
        _sort_projection_pixel_stacks!(_prepared_direction_projection(prepared, direction, options))
    end
end

function _raycore_projected_mesh_area_by_node(geometry::InterceptionSceneData, direction)
    projected_mesh_area = zeros(Float64, length(geometry.node_ids))
    vertices = geometry.vertices::Vector{StaticArrays.SVector{3,Float64}}
    faces = geometry.faces
    face2node_index = geometry.face2node_index
    if geometry.plotbox.nx * geometry.plotbox.ny < _AUTO_VECTOR_PIXEL_HITS_MIN_CELLS
        dirx = Float32(direction[1])
        diry = Float32(direction[2])
        dirz = Float32(direction[3])
        dirz == 0.0f0 && return projected_mesh_area

        pb = geometry.plotbox
        ox = Float32(pb.origin_x)
        oy = Float32(pb.origin_y)
        pxs = Float32(pb.pix_x)
        pys = Float32(pb.pix_y)
        @inbounds for fi in eachindex(faces)
            f = faces[fi]
            projected1, _ = _project_vertex_to_ground_pixel(vertices[f[1]], dirx, diry, dirz, ox, oy, pxs, pys, 1.0f0)
            projected2, _ = _project_vertex_to_ground_pixel(vertices[f[2]], dirx, diry, dirz, ox, oy, pxs, pys, 1.0f0)
            projected3, _ = _project_vertex_to_ground_pixel(vertices[f[3]], dirx, diry, dirz, ox, oy, pxs, pys, 1.0f0)
            projected_mesh_area[face2node_index[fi]] += _triangle_area_xy(projected1, projected2, projected3)
        end
        return projected_mesh_area
    end

    dirx = Float32(direction[1])
    diry = Float32(direction[2])
    dirz = Float32(direction[3])
    dirz == 0.0f0 && return projected_mesh_area
    inv_dirz = inv(abs(dirz))

    @inbounds for fi in eachindex(faces)
        f = faces[fi]
        p1 = vertices[f[1]]
        p2 = vertices[f[2]]
        p3 = vertices[f[3]]
        ux = Float32(p2[1] - p1[1])
        uy = Float32(p2[2] - p1[2])
        uz = Float32(p2[3] - p1[3])
        vx = Float32(p3[1] - p1[1])
        vy = Float32(p3[2] - p1[2])
        vz = Float32(p3[3] - p1[3])
        cx = uy * vz - uz * vy
        cy = uz * vx - ux * vz
        cz = ux * vy - uy * vx
        projected_mesh_area[face2node_index[fi]] +=
            Float64(0.5f0 * abs(cx * dirx + cy * diry + cz * dirz) * inv_dirz)
    end
    return projected_mesh_area
end

function _raycore_projected_mesh_area_by_node(data::RaycoreSceneData, direction)
    key = _raycore_direction_key(direction)
    return get!(data.projected_mesh_area_cache, key) do
        _raycore_projected_mesh_area_by_node(data.prepared.geometry, direction)
    end
end

function _raycore_area_ratio_by_node(
    data::RaycoreSceneData,
    mode::Symbol,
    direction,
    projected_pixels_area::Vector{Float64},
)
    key = (mode, _raycore_direction_key(direction)...)
    return get!(data.area_ratio_cache, key) do
        projected_mesh_area = _raycore_projected_mesh_area_by_node(data, direction)
        ratio = ones(Float64, length(projected_pixels_area))
        @inbounds for idx in eachindex(ratio)
            pixels_area = projected_pixels_area[idx]
            ratio[idx] = pixels_area > 0.0 ? projected_mesh_area[idx] / pixels_area : 1.0
        end
        ratio
    end
end

@inline function _raycore_ground_pixel_center(plotbox, pixel_idx::Int)
    zero_idx = pixel_idx - 1
    i = zero_idx % plotbox.nx
    j = zero_idx ÷ plotbox.nx
    x = plotbox.origin_x + (Float64(i) + 0.5) * plotbox.pix_x
    y = plotbox.origin_y + (Float64(j) + 0.5) * plotbox.pix_y
    return StaticArrays.SVector{3,Float64}(x, y, 0.0)
end

function _raycore_sky_to_ground_ray(plotbox, pixel_idx::Int, direction, far::Float32)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))
    ground = _raycore_ground_pixel_center(plotbox, pixel_idx)
    origin = ground + dir * Float64(far)
    return Raycore.Ray(;
        o=_raycore_point3f(origin),
        d=GeometryBasics.Vec3f(Float32(-dir[1]), Float32(-dir[2]), Float32(-dir[3])),
        t_min=0.0f0,
        t_max=2.0f0 * far,
    )
end

function _raycore_trace_direction_stacks(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    return _raycore_trace_direction_stacks_direct(data, direction, options)
end

function _raycore_direction_projection_full_stack(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    stack_type = _pixel_hit_stack_type(options, plotbox)
    pixel_hits = DensePixelHits(stack_type, n_pixels)
    node_hits = zeros(Int, length(geometry.node_ids))
    projected_mesh_area = zeros(Float64, length(geometry.node_ids))
    projected_pixels_area = zeros(Float64, length(geometry.node_ids))

    Float32(direction[3]) <= 0.0f0 && return DenseDirectionProjectionResult(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
    )

    traced = _raycore_trace_direction_stacks(data, direction, options)
    overflow_pixel = findfirst(traced.overflow)
    overflow_pixel === nothing || error(
        "Raycore max_hits_per_pixel=$(traced.max_hits) exceeded for pixel $overflow_pixel. " *
        "Increase RaycoreBackendConfig(max_hits_per_pixel=...).",
    )

    max_hits = traced.max_hits
    @inbounds for pixel_idx in 1:n_pixels
        count = Int(traced.counts[pixel_idx])
        count == 0 && continue
        offset = (pixel_idx - 1) * max_hits
        for slot in 1:count
            node_idx = Int(traced.nodes[offset+slot])
            1 <= node_idx <= length(geometry.node_ids) ||
                error("Raycore hit returned invalid Archimed node index metadata: $node_idx")
            stack = get!(pixel_hits, pixel_idx) do
                _new_hit_stack(stack_type)
            end
            push!(stack, (Float64(traced.heights[offset+slot]), node_idx))
            node_hits[node_idx] += 1
            projected_pixels_area[node_idx] += plotbox.pixel_area
            projected_mesh_area[node_idx] += plotbox.pixel_area
        end
    end

    if options.area_ratio
        projected_mesh_area .= _raycore_projected_mesh_area_by_node(data, direction)
    end

    return DenseDirectionProjectionResult(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
    )
end

function _raycore_direction_projection_top_hit(
    data::RaycoreSceneData,
    direction,
    options::LightOptions,
)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_pixels = plotbox.nx * plotbox.ny
    pixel_hits = DenseUpperPixelHits(n_pixels)
    node_hits = zeros(Int, length(geometry.node_ids))
    projected_mesh_area = zeros(Float64, length(geometry.node_ids))
    projected_pixels_area = zeros(Float64, length(geometry.node_ids))

    Float32(direction[3]) <= 0.0f0 && return DenseDirectionProjectionResult(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
    )

    traced = _raycore_trace_direction_top_hits_direct(data, direction, options)
    @inbounds for pixel_idx in 1:n_pixels
        node_idx = Int(traced.nodes[pixel_idx])
        node_idx == 0 && continue
        1 <= node_idx <= length(geometry.node_ids) ||
            error("Raycore hit returned invalid Archimed node index metadata: $node_idx")
        pixel_hits.heights[pixel_idx] = Float64(traced.heights[pixel_idx])
        pixel_hits.nodes[pixel_idx] = node_idx
        push!(pixel_hits.occupied, pixel_idx)
        node_hits[node_idx] += 1
        projected_pixels_area[node_idx] += plotbox.pixel_area
        projected_mesh_area[node_idx] += plotbox.pixel_area
    end

    if options.area_ratio
        projected_mesh_area .= _raycore_projected_mesh_area_by_node(data, direction)
    end

    return DenseDirectionProjectionResult(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
    )
end

function _check_raycore_top_hit_supported(prepared::PreparedInterceptionData, options::LightOptions)
    return nothing
end

function _raycore_use_raster_compat_projection(data::RaycoreSceneData)
    geometry = data.prepared.geometry
    # Archimed's CPU rasterizer records coverage per projected pixel cell.
    # A single Raycore point ray per cell cannot reproduce that coverage when
    # many faces collapse into fewer cells, especially for coplanar pavement
    # ties in coarse fixtures. Preserve CPU semantics in that regime.
    return geometry.plotbox.nx * geometry.plotbox.ny < length(geometry.faces)
end

function _raycore_geometry_vertical_span(geometry::InterceptionSceneData)
    pb = geometry.plotbox
    isempty(geometry.vertices) && return max(pb.pix_x, pb.pix_y, 1e-3)
    zs = (Float64(v[3]) for v in geometry.vertices)
    zmin, zmax = extrema(zs)
    return max(zmax - zmin, max(pb.pix_x, pb.pix_y, 1e-3))
end

function _raycore_toric_copy_radius_needed(
    geometry::InterceptionSceneData,
    direction,
    vertical_span::Float64=_raycore_geometry_vertical_span(geometry),
)
    dir = _normalize3(StaticArrays.SVector{3,Float64}(direction[1], direction[2], direction[3]))
    dir[3] > 0.0 || return 0
    pb = geometry.plotbox
    x_radius = iszero(pb.xdim) ? typemax(Int) : ceil(Int, abs(dir[1] / dir[3]) * vertical_span / pb.xdim)
    y_radius = iszero(pb.ydim) ? typemax(Int) : ceil(Int, abs(dir[2] / dir[3]) * vertical_span / pb.ydim)
    return max(x_radius, y_radius)
end

function _raycore_use_raster_compat_projection(data::RaycoreSceneData, direction, options::LightOptions)
    _raycore_use_raster_compat_projection(data) && return true
    options.toricity || return false
    # Raycore currently represents toricity with a finite 3x3 copy of the scene.
    # Low-elevation rays can cross more than one plot width/height while moving
    # through the canopy, which the CPU rasterizer handles by wrapping projected
    # pixels. Use the raster projection for those directions to preserve periodic
    # full-stack semantics.
    return _raycore_toric_copy_radius_needed(data.prepared.geometry, direction, data.vertical_span) > 1
end

function _raycore_use_raster_compat_projection(data::RaycoreSceneData, turtle::TurtleGrid, options::LightOptions)
    _raycore_use_raster_compat_projection(data) && return true
    for sector in turtle.sectors
        _raycore_use_raster_compat_projection(data, sector.direction, options) && return true
    end
    return false
end

function _raycore_direction_projection(data::RaycoreSceneData, direction, options::LightOptions)
    if _raycore_use_raster_compat_projection(data, direction, options)
        return _raycore_raster_compat_projection(data, direction, options)
    end
    prepared = data.prepared
    if !prepared.upper_hit
        return _raycore_direction_projection_full_stack(data, direction, options)
    end
    return _raycore_direction_projection_top_hit(data, direction, options)
end

function _compute_first_order_raycore_from_projections(
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    prepared = data.prepared
    geometry = prepared.geometry

    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))

    for (k, sector) in enumerate(turtle.sectors)
        projection = _raycore_direction_projection(data, sector.direction, options)
        visible_area =
            _visible_area_from_projection_dense(
                projection,
                options,
                geometry.plotbox,
                prepared.virtual_node_mask,
                prepared.node_transparency_by_index,
                geometry,
                true,
            )
        _accumulate_projection_hits!(hits_per_node, projection, geometry)

        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if par_flux == 0.0 && nir_flux == 0.0
            continue
        end

        @inbounds for idx in eachindex(visible_area)
            pa = visible_area[idx]
            pa <= 0.0 && continue
            projected_area_per_node[idx] += pa
            incident_power_par[idx] += par_flux * pa
            incident_power_nir[idx] += nir_flux * pa
        end
    end

    if !isempty(prepared.emitter_nodes)
        w = _emitter_transfer_weights_from_raycore(data, turtle, options)
        for ((to, src), ww) in w
            idx = geometry.node_index[to]
            src_idx = geometry.node_index[src]
            incident_power_par[idx] += ww * prepared.emitter_par_power_by_index[src_idx]
            incident_power_nir[idx] += ww * prepared.emitter_nir_power_by_index[src_idx]
        end
    end

    return _first_order_result_from_dense(
        geometry.node_ids,
        projected_area_per_node,
        incident_power_par,
        incident_power_nir,
        hits_per_node,
    )
end

function _compute_first_order_raycore_top_hit(
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    prepared = data.prepared
    geometry = prepared.geometry
    _check_raycore_top_hit_supported(prepared, options)

    if !prepared.upper_hit ||
       !isempty(prepared.emitter_nodes) ||
       geometry.plotbox.nx * geometry.plotbox.ny < _AUTO_VECTOR_PIXEL_HITS_MIN_CELLS
        return _compute_first_order_raycore_from_projections(data, turtle, fluxes, options)
    end

    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))
    pixel_area = geometry.plotbox.pixel_area
    sector_hits = options.area_ratio ? zeros(Int, length(geometry.node_ids)) : Int[]
    sector_projected_pixels_area = options.area_ratio ? zeros(Float64, length(geometry.node_ids)) : Float64[]

    for (k, sector) in enumerate(turtle.sectors)
        if _raycore_use_raster_compat_projection(data, sector.direction, options)
            projection = _raycore_raster_compat_projection(data, sector.direction, options)
            visible_area =
                _visible_area_from_projection_dense(
                    projection,
                    options,
                    geometry.plotbox,
                    prepared.virtual_node_mask,
                    prepared.node_transparency_by_index,
                    geometry,
                    true,
                )
            _accumulate_projection_hits!(hits_per_node, projection, geometry)

            par_flux = fluxes.par[k]
            nir_flux = fluxes.nir[k]
            if par_flux != 0.0 || nir_flux != 0.0
                @inbounds for idx in eachindex(visible_area)
                    pa = visible_area[idx]
                    pa <= 0.0 && continue
                    projected_area_per_node[idx] += pa
                    par_flux != 0.0 && (incident_power_par[idx] += par_flux * pa)
                    nir_flux != 0.0 && (incident_power_nir[idx] += nir_flux * pa)
                end
            end
            continue
        end

        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if !options.area_ratio
            if Float32(sector.direction[3]) > 0.0f0
                nodes = _raycore_trace_direction_top_nodes_direct(data, sector.direction, options)
                @inbounds for pixel_idx in eachindex(nodes)
                    node_idx = Int(nodes[pixel_idx])
                    node_idx == 0 && continue
                    1 <= node_idx <= length(geometry.node_ids) ||
                        error("Raycore hit returned invalid Archimed node index metadata: $node_idx")
                    hits_per_node[node_idx] += 1
                    projected_area_per_node[node_idx] += pixel_area
                    par_flux != 0.0 && (incident_power_par[node_idx] += par_flux * pixel_area)
                    nir_flux != 0.0 && (incident_power_nir[node_idx] += nir_flux * pixel_area)
                end
            end
            continue
        end

        fill!(sector_hits, 0)
        fill!(sector_projected_pixels_area, 0.0)
        if Float32(sector.direction[3]) > 0.0f0
            nodes = _raycore_trace_direction_top_nodes_direct(data, sector.direction, options)
            @inbounds for pixel_idx in eachindex(nodes)
                node_idx = Int(nodes[pixel_idx])
                node_idx == 0 && continue
                1 <= node_idx <= length(geometry.node_ids) ||
                    error("Raycore hit returned invalid Archimed node index metadata: $node_idx")
                sector_hits[node_idx] += 1
                sector_projected_pixels_area[node_idx] += pixel_area
            end
        end

        if par_flux == 0.0 && nir_flux == 0.0
            @inbounds for idx in eachindex(sector_hits)
                hits_per_node[idx] += sector_hits[idx]
            end
            continue
        end

        sector_area_ratio = _raycore_area_ratio_by_node(data, :top_hit, sector.direction, sector_projected_pixels_area)

        @inbounds for idx in eachindex(sector_hits)
            count = sector_hits[idx]
            count == 0 && continue
            hits_per_node[idx] += count
            ratio = sector_area_ratio[idx]
            pa = 0.0
            for _ in 1:count
                pa += pixel_area * ratio
            end
            pa <= 0.0 && continue
            projected_area_per_node[idx] += pa
            incident_power_par[idx] += par_flux * pa
            incident_power_nir[idx] += nir_flux * pa
        end
    end

    if !isempty(prepared.emitter_nodes)
        w = _emitter_transfer_weights_from_raycore(data, turtle, options)
        for ((to, src), ww) in w
            idx = geometry.node_index[to]
            src_idx = geometry.node_index[src]
            incident_power_par[idx] += ww * prepared.emitter_par_power_by_index[src_idx]
            incident_power_nir[idx] += ww * prepared.emitter_nir_power_by_index[src_idx]
        end
    end

    return _first_order_result_from_dense(
        geometry.node_ids,
        projected_area_per_node,
        incident_power_par,
        incident_power_nir,
        hits_per_node,
    )
end

function _interception_output_keys(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    keys_by_node = Dict{Int,Tuple{Int,Int}}()

    pavement_ids = sort(Int[nid for nid in geometry.node_ids if get(geometry.node_group, nid, "") == "pavement"])
    for (i, nid) in enumerate(pavement_ids)
        keys_by_node[nid] = (-1, i + 1)
    end

    for nid in geometry.node_ids
        haskey(keys_by_node, nid) && continue
        object_id = _scene_object_id(scene, nid, 1)
        source_topology_id = _scene_source_topology_id(scene, nid, nid + 1)
        keys_by_node[nid] = (object_id, source_topology_id)
    end
    keys_by_node
end

"""
    compute_first_order(scene, models, turtle, fluxes, options; backend=:raster_cpu)::FirstOrderResult

Compute first-order interception by rasterizing each direction, then integrating projected area,
incident power, and hit counts per geometry node.

`backend` accepts either a symbol (`:raster_cpu`, `:raycore_cpu`) or an
`InterceptionBackend` instance.
"""
function compute_first_order(
    data::RaycoreSceneData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    return _compute_first_order_raycore_top_hit(data, turtle, fluxes, options)
end

function compute_first_order(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions;
    backend=:raster_cpu,
)
    return compute_first_order(scene, models, turtle, fluxes, options, _resolve_interception_backend(backend))
end

function _resolve_interception_backend(backend::InterceptionBackend)
    return backend
end

function _resolve_interception_backend(backend::Symbol)
    backend == :raster_cpu && return RasterCPUBackend()
    backend == :raycore_cpu && return RaycoreInterceptionBackend()
    error("Unsupported interception backend symbol: $backend (supported: :raster_cpu, :raycore_cpu)")
end

function _resolve_interception_backend(backend)
    error(
        "Unsupported interception backend selector type: $(typeof(backend)). " *
        "Use :raster_cpu, :raycore_cpu, RasterCPUBackend(), or RaycoreInterceptionBackend().",
    )
end

function compute_first_order(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    ::RasterCPUBackend,
)
    prepared = _prepare_interception_data(scene, models, options; include_raycore_instancing=true)
    return _compute_first_order(prepared, turtle, fluxes, options)
end

function _compute_first_order(
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    geometry = prepared.geometry

    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))

    for (k, sector) in enumerate(turtle.sectors)
        projection = _prepared_direction_projection(prepared, sector.direction, options)
        visible_area =
            _visible_area_from_projection_dense(
                projection,
                options,
                geometry.plotbox,
                prepared.virtual_node_mask,
                prepared.node_transparency_by_index,
                geometry,
            )
        _accumulate_projection_hits!(hits_per_node, projection, geometry)

        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if par_flux == 0.0 && nir_flux == 0.0
            continue
        end

        @inbounds for idx in eachindex(visible_area)
            pa = visible_area[idx]
            pa <= 0.0 && continue
            projected_area_per_node[idx] += pa
            incident_power_par[idx] += par_flux * pa
            incident_power_nir[idx] += nir_flux * pa
        end
    end

    if !isempty(prepared.emitter_nodes)
        w = _emitter_transfer_weights(
            geometry.vertices,
            geometry.faces,
            geometry.face2node,
            turtle,
            options,
            geometry.plotbox,
            prepared.emitter_nodes,
            prepared.cache_ctx,
        )
        for ((to, src), ww) in w
            idx = geometry.node_index[to]
            src_idx = geometry.node_index[src]
            incident_power_par[idx] += ww * prepared.emitter_par_power_by_index[src_idx]
            incident_power_nir[idx] += ww * prepared.emitter_nir_power_by_index[src_idx]
        end
    end

    return _first_order_result_from_dense(
        geometry.node_ids,
        projected_area_per_node,
        incident_power_par,
        incident_power_nir,
        hits_per_node,
    )
end

function compute_first_order(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    backend::RaycoreInterceptionBackend,
)
    prepared = _prepare_interception_data(scene, models, options; include_raycore_instancing=true)
    prechunk_status =
        _raycore_prechunk_instance_limit_status(prepared, backend.config; toricity=options.toricity)
    if prechunk_status.exceeded
        _raycore_throw_validation_error(
            backend.config,
            :raycore_prechunk_instance_cap,
            :first_order,
            prechunk_status,
        )
    end
    data, data_was_chunked = _raycore_initial_scene_data(prepared, backend.config, options; toricity=options.toricity)
    validation = _raycore_vertical_trace_validation(data, options)
    if !validation.ok
        chunked_data, chunked_validation = _raycore_retry_chunked_scene_data(
            data.prepared,
            backend.config,
            options,
            validation;
            toricity=options.toricity,
            skip_face_chunk_limit=data_was_chunked ? _raycore_face_chunk_limit() : nothing,
        )
        if chunked_data !== nothing
            return _compute_first_order_raycore_top_hit(chunked_data, turtle, fluxes, options)
        end
        validation = chunked_validation
        _raycore_throw_validation_error(backend.config, :raycore_trace_validation, :first_order, validation)
    end
    return _compute_first_order_raycore_top_hit(data, turtle, fluxes, options)
end

function compute_first_order(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    backend::InterceptionBackend,
)
    error("Unsupported interception backend type: $(typeof(backend))")
end
