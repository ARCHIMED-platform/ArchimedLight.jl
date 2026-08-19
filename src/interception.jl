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

struct RasterGPUSceneData{P,B,VX,VY,VZ,F,N,T,M,I,C,S,H,O,G,A,D,K,KC,TC,TF,TO,E}
    prepared::P
    backend::B
    vertex_x_dev::VX
    vertex_y_dev::VY
    vertex_z_dev::VZ
    face_i_dev::F
    face_j_dev::F
    face_k_dev::F
    face2node_index_dev::N
    node_transparency_dev::T
    virtual_node_mask_dev::M
    pavement_node_mask_dev::M
    node_ids_dev::I
    counts_dev::C
    nodes_dev::S
    heights_dev::H
    overflow_dev::O
    node_counts_dev::G
    projected_mesh_area_dev::A
    projected_pixels_area_dev::A
    sector_area_dev::A
    dense_edge_counts_dev::D
    edge_keys_dev::K
    edge_key_counts_dev::KC
    tile_counts_dev::TC
    tile_faces_dev::TF
    tile_unwrapped_i_dev::TF
    tile_unwrapped_j_dev::TF
    tile_overflow_dev::TO
    counts_host::Vector{Int32}
    nodes_host::Vector{UInt32}
    heights_host::Vector{Float32}
    overflow_host::Vector{Bool}
    node_counts_host::Vector{Int32}
    projected_mesh_area_host::Vector{Float32}
    projected_pixels_area_host::Vector{Float32}
    sector_area_host::Vector{Float32}
    dense_edge_counts_host::E
    edge_keys_host::Vector{UInt64}
    edge_key_counts_host::Vector{Int32}
    edge_compact_host::Vector{UInt64}
    tile_counts_host::Vector{Int32}
    tile_overflow_host::Vector{Bool}
    workgroupsize::Int
    max_hits_per_pixel::Int
    tile_size::Int
    tile_face_capacity::Int
    edge_accumulation::Symbol
    dense_edge_limit_bytes::Int
    validate::Bool
end

struct RasterGPUFusedProjectionParams
    dirx::Float32
    diry::Float32
    dirz::Float32
    origin_x::Float32
    origin_y::Float32
    pix_x::Float32
    pix_y::Float32
    pixel_area::Float32
    nx::Int
    ny::Int
    n_tiles_x::Int
    tile_size::Int
    tile_face_capacity::Int
    toricity::Bool
    area_ratio::Bool
    max_hits::Int
    n_nodes::Int
    accumulate_dense_edges::Bool
end
const RASTERGPU_FUSED_DENSE_MAX_HITS = 128
_rastergpu_fused_dense_supported(config::RasterGPUBackendConfig) =
    !(config.backend isa KernelAbstractions.CPU) &&
    config.max_hits_per_pixel <= RASTERGPU_FUSED_DENSE_MAX_HITS

function _rastergpu_use_fused_dense_edges(data::RasterGPUSceneData)
    return !data.prepared.upper_hit &&
           data.dense_edge_counts_dev !== nothing &&
           data.tile_size == 1 &&
           data.max_hits_per_pixel <= RASTERGPU_FUSED_DENSE_MAX_HITS &&
           length(data.nodes_dev) == 1 &&
           length(data.heights_dev) == 1
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

function _is_ignored_node_metadata(group::AbstractString, type_name::AbstractString, ignored::Dict{String,Set{String}})
    isempty(ignored) && return false
    g = _normalize_group_name_local(group)
    t = strip(type_name)
    return haskey(ignored, g) && (t in ignored[g])
end

function _ignored_node_ids(scene::PlantGeom.SceneGeometry, node_ids, ignored::Dict{String,Set{String}})
    out = Set{Int}()
    isempty(ignored) && return out
    for nid in node_ids
        _is_ignored_node(nid, scene, ignored) && push!(out, nid)
    end
    return out
end

function _ignored_node_ids(
    node_ids,
    node_group::Dict{Int,String},
    node_type::Dict{Int,String},
    ignored::Dict{String,Set{String}},
)
    out = Set{Int}()
    isempty(ignored) && return out
    for nid in node_ids
        _is_ignored_node_metadata(get(node_group, nid, ""), get(node_type, nid, ""), ignored) && push!(out, nid)
    end
    return out
end

function _scene_type_from_mtg_node(node, default="")
    for key in (:type, :Type, :functional_type, :functionalType, :organ_type, :organType)
        v = _mtg_node_attr(node, key)
        if v !== nothing
            s = strip(string(v))
            isempty(s) || return s
        end
    end
    node === nothing && return default
    s = string(MultiScaleTreeGraph.symbol(node))
    isempty(s) ? default : s
end

function _scene_group_type_maps(scene::PlantGeom.SceneGeometry, node_ids)
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    wanted = Set{Int}(Int(nid) for nid in node_ids)

    if scene.mtg !== nothing
        inherited_group = Dict{Int,String}()
        MultiScaleTreeGraph.traverse!(scene.mtg) do node
            nid = MultiScaleTreeGraph.node_id(node)
            local_group = nothing
            for key in (:group, :functional_group)
                local_group = _mtg_node_attr(node, key)
                local_group === nothing || break
            end
            parent_group = ""
            if !MultiScaleTreeGraph.isroot(node)
                parent_node = MultiScaleTreeGraph.parent(node)
                parent_node === nothing || (parent_group = get(inherited_group, MultiScaleTreeGraph.node_id(parent_node), ""))
            end
            group = local_group === nothing ? parent_group : string(local_group)
            inherited_group[nid] = group
            if nid in wanted
                node_group[nid] = group
                node_type[nid] = _scene_type_from_mtg_node(node, "")
            end
            return nothing
        end
    end

    for nid in wanted
        haskey(node_group, nid) || (node_group[nid] = _scene_group(scene, nid, ""))
        haskey(node_type, nid) || (node_type[nid] = _scene_type(scene, nid, ""))
    end

    return node_group, node_type
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
    node_ids = unique(face2node)
    node_group, node_type = _scene_group_type_maps(scene, node_ids)
    ignored_nodes = _ignored_node_ids(node_ids, node_group, node_type, ignored)
    for nid in node_ids
        nid in ignored_nodes && continue
        group = strip(get(node_group, nid, ""))
        type_name = strip(get(node_type, nid, ""))
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

function _projection_cache_path(cache_ctx::ProjectionCacheContext, direction, upper_hit::Bool=false)
    scene_hex = string(cache_ctx.scene_key, base=16, pad=16)
    dir_hex = string(_projection_dir_key(direction), base=16, pad=16)
    mode = upper_hit ? "u1" : "u0"
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

    normal = _compute_normal((pix1, pix2, pix3))
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

    normal = _compute_normal((pix1, pix2, pix3))
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

@inline function _rastergpu_project_vertex(
    x::Float32,
    y::Float32,
    z::Float32,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    ox::Float32,
    oy::Float32,
    pxs::Float32,
    pys::Float32,
)
    dz = -z / dirz
    xw = x + dirx * dz
    yw = y + diry * dz
    return xw, yw, (xw - ox) / pxs, (yw - oy) / pys, z
end

@inline function _rastergpu_triangle_area_xy(
    x1::Float32,
    y1::Float32,
    x2::Float32,
    y2::Float32,
    x3::Float32,
    y3::Float32,
)
    return abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) * 0.5f0
end

@inline function _rastergpu_normal(
    x1::Float32,
    y1::Float32,
    z1::Float32,
    x2::Float32,
    y2::Float32,
    z2::Float32,
    x3::Float32,
    y3::Float32,
    z3::Float32,
)
    v1x = x2 - x1
    v1y = y2 - y1
    v1z = z2 - z1
    n1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
    n1 <= 0.0f0 && return 0.0f0, 0.0f0, 0.0f0
    v1x /= n1
    v1y /= n1
    v1z /= n1

    v2x = x3 - x2
    v2y = y3 - y2
    v2z = z3 - z2
    n2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
    n2 <= 0.0f0 && return 0.0f0, 0.0f0, 0.0f0
    v2x /= n2
    v2y /= n2
    v2z /= n2

    nx = (v1y * v2z) - (v1z * v2y)
    ny = (v1z * v2x) - (v1x * v2z)
    nz = (v1x * v2y) - (v1y * v2x)
    nnorm = sqrt(nx * nx + ny * ny + nz * nz)
    nnorm <= 0.0f0 && return 0.0f0, 0.0f0, 0.0f0
    return nx / nnorm, ny / nnorm, nz / nnorm
end

@inline function _rastergpu_edge_bounds(
    i::Int,
    ymin::Int,
    ymax::Int,
    x1::Float32,
    y1::Float32,
    x2::Float32,
    y2::Float32,
)
    pminx = x1
    pminy = y1
    pmaxx = x2
    pmaxy = y2
    if pmaxx < pminx
        pminx = x2
        pminy = y2
        pmaxx = x1
        pmaxy = y1
    end
    dx = pmaxx - pminx
    dx < 1.0f-6 && return ymin, ymax
    i >= ceil(Int, pminx) || return ymin, ymax
    i <= floor(Int, pmaxx) || return ymin, ymax
    slope = (pmaxy - pminy) / dx
    yi = slope * (Float32(i) - pminx) + pminy
    j = trunc(Int, floor(yi) + 0.5f0)
    return min(ymin, j), max(ymax, j)
end

@inline function _rastergpu_wrap_index(i::Int, n::Int)
    ii = i
    while ii < 0
        ii += n
    end
    while ii >= n
        ii -= n
    end
    return ii
end

@inline function _rastergpu_tile_coord(i::Int, tile_size::Int)
    return fld(i, tile_size)
end

KernelAbstractions.@kernel function _rastergpu_clear_direction_kernel!(
    counts,
    overflow,
    node_counts,
    projected_mesh_area,
    projected_pixels_area,
    sector_area,
    n_pixels::Int,
    n_nodes::Int,
)
    idx = @index(Global, Linear)
    @inbounds begin
        if idx <= n_pixels
            counts[idx] = Int32(0)
            overflow[idx] = false
        end
        if idx <= n_nodes
            node_counts[idx] = Int32(0)
            projected_mesh_area[idx] = 0.0f0
            projected_pixels_area[idx] = 0.0f0
            sector_area[idx] = 0.0f0
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_clear_tile_bins_kernel!(
    tile_counts,
    tile_overflow,
    n_tiles::Int,
)
    tile_idx = @index(Global, Linear)
    @inbounds begin
        if tile_idx <= n_tiles
            tile_counts[tile_idx] = Int32(0)
            tile_overflow[tile_idx] = false
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_bin_faces_to_tiles_kernel!(
    tile_counts,
    tile_faces,
    tile_unwrapped_i,
    tile_unwrapped_j,
    tile_overflow,
    projected_mesh_area,
    projected_pixels_area,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    pixel_area::Float32,
    nx::Int,
    ny::Int,
    n_tiles_x::Int,
    n_tiles_y::Int,
    tile_size::Int,
    tile_face_capacity::Int,
    toricity::Bool,
)
    fi = @index(Global, Linear)
    @inbounds begin
        if dirz != 0.0f0
            vi = Int(face_i[fi])
            vj = Int(face_j[fi])
            vk = Int(face_k[fi])
            node_idx = Int(face2node_index[fi])

            x1 = vertex_x[vi]
            y1 = vertex_y[vi]
            z1 = vertex_z[vi]
            x2 = vertex_x[vj]
            y2 = vertex_y[vj]
            z2 = vertex_z[vj]
            x3 = vertex_x[vk]
            y3 = vertex_y[vk]
            z3 = vertex_z[vk]

            px1, py1, pix1x, pix1y, _ =
                _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px2, py2, pix2x, pix2y, _ =
                _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px3, py3, pix3x, pix3y, _ =
                _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

            i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
            i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
            j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
            j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))

            if i_max > i_min && j_max > j_min
                tri_proj_area = _rastergpu_triangle_area_xy(px1, py1, px2, py2, px3, py3)
                nb_hits = Int32(0)
                for i in i_min:(i_max - 1)
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)
                    for _ in ymin_i:(ymax_i - 1)
                        nb_hits += Int32(1)
                    end
                end
                @atomic :monotonic projected_mesh_area[node_idx] += tri_proj_area
                @atomic :monotonic projected_pixels_area[node_idx] += Float32(nb_hits) * pixel_area

                tile_i_min = _rastergpu_tile_coord(i_min, tile_size)
                tile_i_max = _rastergpu_tile_coord(i_max - 1, tile_size)
                tile_j_min = _rastergpu_tile_coord(j_min, tile_size)
                tile_j_max = _rastergpu_tile_coord(j_max - 1, tile_size)
                if !toricity
                    tile_i_min = max(tile_i_min, 0)
                    tile_j_min = max(tile_j_min, 0)
                    tile_i_max = min(tile_i_max, n_tiles_x - 1)
                    tile_j_max = min(tile_j_max, n_tiles_y - 1)
                end

                for tile_j in tile_j_min:tile_j_max
                    for tile_i in tile_i_min:tile_i_max
                        wrapped_i = toricity ? _rastergpu_wrap_index(tile_i, n_tiles_x) : tile_i
                        wrapped_j = toricity ? _rastergpu_wrap_index(tile_j, n_tiles_y) : tile_j
                        if (0 <= wrapped_i) && (wrapped_i < n_tiles_x) &&
                           (0 <= wrapped_j) && (wrapped_j < n_tiles_y)
                            tile_idx = wrapped_i + 1 + wrapped_j * n_tiles_x
                            new_count = @atomic :monotonic tile_counts[tile_idx] += Int32(1)
                            slot = Int(new_count)
                            if slot <= tile_face_capacity
                                candidate_idx = (tile_idx - 1) * tile_face_capacity + slot
                                tile_faces[candidate_idx] = Int32(fi)
                                tile_unwrapped_i[candidate_idx] = Int32(tile_i)
                                tile_unwrapped_j[candidate_idx] = Int32(tile_j)
                            else
                                tile_overflow[tile_idx] = true
                            end
                        end
                    end
                end
            end
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_bin_faces_to_covered_pixel_tiles_kernel!(
    tile_counts,
    tile_faces,
    tile_unwrapped_i,
    tile_unwrapped_j,
    tile_overflow,
    projected_mesh_area,
    projected_pixels_area,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    pixel_area::Float32,
    nx::Int,
    ny::Int,
    tile_face_capacity::Int,
    toricity::Bool,
)
    fi = @index(Global, Linear)
    @inbounds begin
        if dirz != 0.0f0
            vi = Int(face_i[fi])
            vj = Int(face_j[fi])
            vk = Int(face_k[fi])
            node_idx = Int(face2node_index[fi])

            x1 = vertex_x[vi]
            y1 = vertex_y[vi]
            z1 = vertex_z[vi]
            x2 = vertex_x[vj]
            y2 = vertex_y[vj]
            z2 = vertex_z[vj]
            x3 = vertex_x[vk]
            y3 = vertex_y[vk]
            z3 = vertex_z[vk]

            px1, py1, pix1x, pix1y, _ =
                _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px2, py2, pix2x, pix2y, _ =
                _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px3, py3, pix3x, pix3y, _ =
                _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

            i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
            i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
            j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
            j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))

            if i_max > i_min && j_max > j_min
                tri_proj_area = _rastergpu_triangle_area_xy(px1, py1, px2, py2, px3, py3)
                nb_hits = Int32(0)
                for i in i_min:(i_max - 1)
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)
                    for j in ymin_i:(ymax_i - 1)
                        nb_hits += Int32(1)
                        ii = i
                        jj = j
                        if toricity
                            ii = _rastergpu_wrap_index(i, nx)
                            jj = _rastergpu_wrap_index(j, ny)
                        end
                        if toricity || ((0 <= ii) && (ii < nx) && (0 <= jj) && (jj < ny))
                            tile_idx = ii + 1 + jj * nx
                            new_count = @atomic :monotonic tile_counts[tile_idx] += Int32(1)
                            slot = Int(new_count)
                            if slot <= tile_face_capacity
                                candidate_idx = (tile_idx - 1) * tile_face_capacity + slot
                                tile_faces[candidate_idx] = Int32(fi)
                                tile_unwrapped_i[candidate_idx] = Int32(i)
                                tile_unwrapped_j[candidate_idx] = Int32(j)
                            else
                                tile_overflow[tile_idx] = true
                            end
                        end
                    end
                end
                @atomic :monotonic projected_mesh_area[node_idx] += tri_proj_area
                @atomic :monotonic projected_pixels_area[node_idx] += Float32(nb_hits) * pixel_area
            end
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_project_tile_bins_kernel!(
    counts,
    nodes,
    heights,
    overflow,
    node_counts,
    tile_counts,
    tile_faces,
    tile_unwrapped_i,
    tile_unwrapped_j,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    nx::Int,
    ny::Int,
    n_tiles_x::Int,
    tile_size::Int,
    tile_face_capacity::Int,
    toricity::Bool,
    max_hits::Int,
)
    work_idx = @index(Global, Linear)
    @inbounds begin
        tile_area = tile_size * tile_size
        tile_idx = fld(work_idx - 1, tile_area) + 1
        local_idx = (work_idx - 1) - (tile_idx - 1) * tile_area
        local_i = local_idx - fld(local_idx, tile_size) * tile_size
        local_j = fld(local_idx, tile_size)

        candidate_count = Int(tile_counts[tile_idx])
        candidate_count > tile_face_capacity && (candidate_count = tile_face_capacity)
        candidate_base = (tile_idx - 1) * tile_face_capacity
        for candidate in 1:candidate_count
            candidate_idx = candidate_base + candidate
            fi = Int(tile_faces[candidate_idx])
            tile_i = Int(tile_unwrapped_i[candidate_idx])
            tile_j = Int(tile_unwrapped_j[candidate_idx])

            vi = Int(face_i[fi])
            vj = Int(face_j[fi])
            vk = Int(face_k[fi])
            node_idx = Int(face2node_index[fi])

            x1 = vertex_x[vi]
            y1 = vertex_y[vi]
            z1 = vertex_z[vi]
            x2 = vertex_x[vj]
            y2 = vertex_y[vj]
            z2 = vertex_z[vj]
            x3 = vertex_x[vk]
            y3 = vertex_y[vk]
            z3 = vertex_z[vk]

            _, _, pix1x, pix1y, pix1z =
                _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            _, _, pix2x, pix2y, pix2z =
                _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            _, _, pix3x, pix3y, pix3z =
                _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

            i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
            i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
            j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
            j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))
            k_min = min(floor(Int, pix1z), floor(Int, pix2z), floor(Int, pix3z))
            k_max = max(ceil(Int, pix1z), ceil(Int, pix2z), ceil(Int, pix3z))

            if i_max > i_min && j_max > j_min
                normal_x, normal_y, normal_z =
                    _rastergpu_normal(pix1x, pix1y, pix1z, pix2x, pix2y, pix2z, pix3x, pix3y, pix3z)
                slope_x, slope_y =
                    if abs(normal_z) > 1.0f-5
                        normal_x / normal_z, normal_y / normal_z
                    else
                        dirz * normal_x, dirz * normal_y
                    end
                z0 = pix1z + slope_x * (pix1x - Float32(i_min)) + slope_y * (pix1y - Float32(j_min))

                i = tile_i * tile_size + local_i
                j = tile_j * tile_size + local_j
                if (i_min <= i) && (i <= i_max - 1) && (j_min <= j) && (j <= j_max - 1)
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)

                    if (ymin_i <= j) && (j <= ymax_i - 1)
                        ni = i - i_min
                        zi = z0 - slope_x * Float32(ni)
                        nj = j - j_min
                        zpix = zi - slope_y * Float32(nj)
                        zpix = clamp(zpix, Float32(k_min), Float32(k_max))

                        ii = i
                        jj = j
                        if toricity
                            ii = _rastergpu_wrap_index(i, nx)
                            jj = _rastergpu_wrap_index(j, ny)
                        end

                        if toricity || ((0 <= ii) && (ii < nx) && (0 <= jj) && (jj < ny))
                            pixel_idx = ii + 1 + jj * nx
                            new_count = @atomic :monotonic counts[pixel_idx] += Int32(1)
                            slot = Int(new_count)
                            if slot <= max_hits
                                stack_idx = (pixel_idx - 1) * max_hits + slot
                                nodes[stack_idx] = UInt32(node_idx)
                                heights[stack_idx] = zpix
                            else
                                overflow[pixel_idx] = true
                            end
                            @atomic :monotonic node_counts[node_idx] += Int32(1)
                        end
                    end
                end
            end
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_project_tile_bins_top_hit_kernel!(
    node_counts,
    sector_area,
    projected_mesh_area,
    projected_pixels_area,
    tile_counts,
    tile_faces,
    tile_unwrapped_i,
    tile_unwrapped_j,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    pixel_area::Float32,
    nx::Int,
    ny::Int,
    n_tiles_x::Int,
    tile_size::Int,
    tile_face_capacity::Int,
    toricity::Bool,
    area_ratio::Bool,
)
    work_idx = @index(Global, Linear)
    @inbounds begin
        tile_area = tile_size * tile_size
        tile_idx = fld(work_idx - 1, tile_area) + 1
        local_idx = (work_idx - 1) - (tile_idx - 1) * tile_area
        local_i = local_idx - fld(local_idx, tile_size) * tile_size
        local_j = fld(local_idx, tile_size)
        found_top = false
        top_height = -typemax(Float32)
        top_node = 0

        candidate_count = Int(tile_counts[tile_idx])
        candidate_count > tile_face_capacity && (candidate_count = tile_face_capacity)
        candidate_base = (tile_idx - 1) * tile_face_capacity
        for candidate in 1:candidate_count
            candidate_idx = candidate_base + candidate
            fi = Int(tile_faces[candidate_idx])
            tile_i = Int(tile_unwrapped_i[candidate_idx])
            tile_j = Int(tile_unwrapped_j[candidate_idx])

            vi = Int(face_i[fi])
            vj = Int(face_j[fi])
            vk = Int(face_k[fi])
            node_idx = Int(face2node_index[fi])

            x1 = vertex_x[vi]
            y1 = vertex_y[vi]
            z1 = vertex_z[vi]
            x2 = vertex_x[vj]
            y2 = vertex_y[vj]
            z2 = vertex_z[vj]
            x3 = vertex_x[vk]
            y3 = vertex_y[vk]
            z3 = vertex_z[vk]

            _, _, pix1x, pix1y, pix1z =
                _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            _, _, pix2x, pix2y, pix2z =
                _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            _, _, pix3x, pix3y, pix3z =
                _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

            i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
            i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
            j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
            j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))
            k_min = min(floor(Int, pix1z), floor(Int, pix2z), floor(Int, pix3z))
            k_max = max(ceil(Int, pix1z), ceil(Int, pix2z), ceil(Int, pix3z))

            if i_max > i_min && j_max > j_min
                normal_x, normal_y, normal_z =
                    _rastergpu_normal(pix1x, pix1y, pix1z, pix2x, pix2y, pix2z, pix3x, pix3y, pix3z)
                slope_x, slope_y =
                    if abs(normal_z) > 1.0f-5
                        normal_x / normal_z, normal_y / normal_z
                    else
                        dirz * normal_x, dirz * normal_y
                    end
                z0 = pix1z + slope_x * (pix1x - Float32(i_min)) + slope_y * (pix1y - Float32(j_min))

                i = tile_i * tile_size + local_i
                j = tile_j * tile_size + local_j
                if (i_min <= i) && (i <= i_max - 1) && (j_min <= j) && (j <= j_max - 1)
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)

                    if (ymin_i <= j) && (j <= ymax_i - 1)
                        ni = i - i_min
                        zi = z0 - slope_x * Float32(ni)
                        nj = j - j_min
                        zpix = zi - slope_y * Float32(nj)
                        zpix = clamp(zpix, Float32(k_min), Float32(k_max))

                        ii = i
                        jj = j
                        if toricity
                            ii = _rastergpu_wrap_index(i, nx)
                            jj = _rastergpu_wrap_index(j, ny)
                        end

                        if toricity || ((0 <= ii) && (ii < nx) && (0 <= jj) && (jj < ny))
                            if !found_top || zpix > top_height || (zpix == top_height && node_idx < top_node)
                                found_top = true
                                top_height = zpix
                                top_node = node_idx
                            end
                            @atomic :monotonic node_counts[node_idx] += Int32(1)
                        end
                    end
                end
            end
        end

        if found_top
            ratio = 1.0f0
            if area_ratio
                pixels_area = projected_pixels_area[top_node]
                ratio = pixels_area > 0.0f0 ? projected_mesh_area[top_node] / pixels_area : 1.0f0
            end
            @atomic :monotonic sector_area[top_node] += pixel_area * ratio
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_project_direction_kernel!(
    counts,
    nodes,
    heights,
    overflow,
    node_counts,
    projected_mesh_area,
    projected_pixels_area,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    dirx::Float32,
    diry::Float32,
    dirz::Float32,
    origin_x::Float32,
    origin_y::Float32,
    pix_x::Float32,
    pix_y::Float32,
    pixel_area::Float32,
    nx::Int,
    ny::Int,
    toricity::Bool,
    max_hits::Int,
)
    fi = @index(Global, Linear)
    @inbounds begin
        if dirz != 0.0f0
            vi = Int(face_i[fi])
            vj = Int(face_j[fi])
            vk = Int(face_k[fi])
            node_idx = Int(face2node_index[fi])

            x1 = vertex_x[vi]
            y1 = vertex_y[vi]
            z1 = vertex_z[vi]
            x2 = vertex_x[vj]
            y2 = vertex_y[vj]
            z2 = vertex_z[vj]
            x3 = vertex_x[vk]
            y3 = vertex_y[vk]
            z3 = vertex_z[vk]

            px1, py1, pix1x, pix1y, pix1z =
                _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px2, py2, pix2x, pix2y, pix2z =
                _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
            px3, py3, pix3x, pix3y, pix3z =
                _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

            i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
            i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
            j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
            j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))
            k_min = min(floor(Int, pix1z), floor(Int, pix2z), floor(Int, pix3z))
            k_max = max(ceil(Int, pix1z), ceil(Int, pix2z), ceil(Int, pix3z))

            if i_max >= i_min
                normal_x, normal_y, normal_z =
                    _rastergpu_normal(pix1x, pix1y, pix1z, pix2x, pix2y, pix2z, pix3x, pix3y, pix3z)
                slope_x, slope_y =
                    if abs(normal_z) > 1.0f-5
                        normal_x / normal_z, normal_y / normal_z
                    else
                        dirz * normal_x, dirz * normal_y
                    end
                z0 = pix1z + slope_x * (pix1x - Float32(i_min)) + slope_y * (pix1y - Float32(j_min))

                tri_proj_area = _rastergpu_triangle_area_xy(px1, py1, px2, py2, px3, py3)
                nb_hits = Int32(0)
                for i in i_min:(i_max - 1)
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)

                    ni = i - i_min
                    zi = z0 - slope_x * Float32(ni)
                    for j in ymin_i:(ymax_i - 1)
                        nj = j - j_min
                        zpix = zi - slope_y * Float32(nj)
                        zpix = clamp(zpix, Float32(k_min), Float32(k_max))
                        nb_hits += Int32(1)

                        ii = i
                        jj = j
                        if toricity
                            ii = _rastergpu_wrap_index(i, nx)
                            jj = _rastergpu_wrap_index(j, ny)
                        end

                        if toricity || ((0 <= ii) && (ii < nx) && (0 <= jj) && (jj < ny))
                            pixel_idx = ii + 1 + jj * nx
                            new_count = @atomic :monotonic counts[pixel_idx] += Int32(1)
                            slot = Int(new_count)
                            if slot <= max_hits
                                stack_idx = (pixel_idx - 1) * max_hits + slot
                                nodes[stack_idx] = UInt32(node_idx)
                                heights[stack_idx] = zpix
                            else
                                overflow[pixel_idx] = true
                            end
                            @atomic :monotonic node_counts[node_idx] += Int32(1)
                        end
                    end
                end
                @atomic :monotonic projected_mesh_area[node_idx] += tri_proj_area
                @atomic :monotonic projected_pixels_area[node_idx] += Float32(nb_hits) * pixel_area
            end
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_sort_and_reduce_kernel!(
    counts,
    nodes,
    heights,
    overflow,
    sector_area,
    projected_mesh_area,
    projected_pixels_area,
    virtual_node_mask,
    node_transparency,
    pixel_area::Float32,
    max_hits::Int,
    n_nodes::Int,
    area_ratio::Bool,
    upper_hit::Bool,
)
    pixel_idx = @index(Global, Linear)
    @inbounds begin
        n_hits = Int(counts[pixel_idx])
        if n_hits > max_hits
            overflow[pixel_idx] = true
            n_hits = max_hits
            counts[pixel_idx] = Int32(max_hits)
        end
        if n_hits > 0
            stack_base = (pixel_idx - 1) * max_hits

            for h in 2:n_hits
                idx = stack_base + h
                xh = heights[idx]
                xn = nodes[idx]
                j = h - 1
                while j >= 1 &&
                          (heights[stack_base+j] < xh ||
                           (heights[stack_base+j] == xh && nodes[stack_base+j] > xn))
                    heights[stack_base+j+1] = heights[stack_base+j]
                    nodes[stack_base+j+1] = nodes[stack_base+j]
                    j -= 1
                end
                heights[stack_base+j+1] = xh
                nodes[stack_base+j+1] = xn
            end

            if upper_hit
                node_idx = Int(nodes[stack_base+1])
                if 1 <= node_idx <= n_nodes
                    ratio = 1.0f0
                    if area_ratio
                        pixels_area = projected_pixels_area[node_idx]
                        ratio = pixels_area > 0.0f0 ? projected_mesh_area[node_idx] / pixels_area : 1.0f0
                    end
                    @atomic :monotonic sector_area[node_idx] += pixel_area * ratio
                end
            else
                area_active = true
                for h in 1:n_hits
                    node_idx = Int(nodes[stack_base+h])
                    1 <= node_idx <= n_nodes || continue
                    area_active || continue

                    ratio = 1.0f0
                    if area_ratio
                        pixels_area = projected_pixels_area[node_idx]
                        ratio = pixels_area > 0.0f0 ? projected_mesh_area[node_idx] / pixels_area : 1.0f0
                    end

                    if virtual_node_mask[node_idx]
                        @atomic :monotonic sector_area[node_idx] += pixel_area * ratio
                        continue
                    end

                    transparency = node_transparency[node_idx]
                    intercepted_fraction = 1.0f0 - transparency
                    if intercepted_fraction > 0.0f0
                        @atomic :monotonic sector_area[node_idx] += pixel_area * intercepted_fraction * ratio
                    end
                    area_active = transparency > 0.0f0
                end
            end
        end
    end
end

KernelAbstractions.@kernel function _rastergpu_clear_dense_counts_kernel!(dense_counts)
    pair_idx = @index(Global, Linear)
    @inbounds dense_counts[pair_idx] = Int32(0)
end

KernelAbstractions.@kernel function _rastergpu_project_tile_bins_fused_dense_kernel!(
    counts,
    overflow,
    node_counts,
    sector_area,
    dense_edge_counts,
    projected_mesh_area,
    projected_pixels_area,
    tile_counts,
    tile_faces,
    tile_unwrapped_i,
    tile_unwrapped_j,
    vertex_x,
    vertex_y,
    vertex_z,
    face_i,
    face_j,
    face_k,
    face2node_index,
    virtual_node_mask,
    pavement_node_mask,
    node_transparency,
    params::RasterGPUFusedProjectionParams,
)
    work_idx = @index(Global, Linear)
    local_nodes = StaticArrays.MVector{RASTERGPU_FUSED_DENSE_MAX_HITS,UInt32}(undef)
    local_heights = StaticArrays.MVector{RASTERGPU_FUSED_DENSE_MAX_HITS,Float32}(undef)
    @inbounds begin
        dirx = params.dirx
        diry = params.diry
        dirz = params.dirz
        origin_x = params.origin_x
        origin_y = params.origin_y
        pix_x = params.pix_x
        pix_y = params.pix_y
        pixel_area = params.pixel_area
        nx = params.nx
        ny = params.ny
        n_tiles_x = params.n_tiles_x
        tile_size = params.tile_size
        tile_face_capacity = params.tile_face_capacity
        toricity = params.toricity
        area_ratio = params.area_ratio
        max_hits = params.max_hits
        n_nodes = params.n_nodes
        accumulate_dense_edges = params.accumulate_dense_edges
        tile_area = tile_size * tile_size
        tile_idx = fld(work_idx - 1, tile_area) + 1
        local_idx = (work_idx - 1) - (tile_idx - 1) * tile_area
        local_i = local_idx - fld(local_idx, tile_size) * tile_size
        local_j = fld(local_idx, tile_size)
        tile_i_base = (tile_idx - 1) - fld(tile_idx - 1, n_tiles_x) * n_tiles_x
        tile_j_base = fld(tile_idx - 1, n_tiles_x)
        output_i = tile_i_base * tile_size + local_i
        output_j = tile_j_base * tile_size + local_j

        output_pixel_in_bounds = (0 <= output_i) && (output_i < nx) && (0 <= output_j) && (output_j < ny)
        pixel_idx = output_pixel_in_bounds ? output_i + 1 + output_j * nx : 0

        n_hits = 0
        observed_hits = 0
        if pixel_idx != 0 && dirz != 0.0f0
            candidate_count = Int(tile_counts[tile_idx])
            candidate_count > tile_face_capacity && (candidate_count = tile_face_capacity)
            candidate_base = (tile_idx - 1) * tile_face_capacity
            for candidate in 1:candidate_count
                candidate_idx = candidate_base + candidate
                fi = Int(tile_faces[candidate_idx])
                tile_i = Int(tile_unwrapped_i[candidate_idx])
                tile_j = Int(tile_unwrapped_j[candidate_idx])
                i = tile_i * tile_size + local_i
                j = tile_j * tile_size + local_j
                ii = i
                jj = j
                if toricity
                    ii = _rastergpu_wrap_index(i, nx)
                    jj = _rastergpu_wrap_index(j, ny)
                end
                ((0 <= ii) && (ii < nx) && (0 <= jj) && (jj < ny)) || continue
                ii + 1 + jj * nx == pixel_idx || continue

                vi = Int(face_i[fi])
                vj = Int(face_j[fi])
                vk = Int(face_k[fi])
                node_idx = Int(face2node_index[fi])

                x1 = vertex_x[vi]
                y1 = vertex_y[vi]
                z1 = vertex_z[vi]
                x2 = vertex_x[vj]
                y2 = vertex_y[vj]
                z2 = vertex_z[vj]
                x3 = vertex_x[vk]
                y3 = vertex_y[vk]
                z3 = vertex_z[vk]

                _, _, pix1x, pix1y, pix1z =
                    _rastergpu_project_vertex(x1, y1, z1, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
                _, _, pix2x, pix2y, pix2z =
                    _rastergpu_project_vertex(x2, y2, z2, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)
                _, _, pix3x, pix3y, pix3z =
                    _rastergpu_project_vertex(x3, y3, z3, dirx, diry, dirz, origin_x, origin_y, pix_x, pix_y)

                i_min = min(floor(Int, pix1x), floor(Int, pix2x), floor(Int, pix3x))
                i_max = max(ceil(Int, pix1x), ceil(Int, pix2x), ceil(Int, pix3x))
                j_min = min(floor(Int, pix1y), floor(Int, pix2y), floor(Int, pix3y))
                j_max = max(ceil(Int, pix1y), ceil(Int, pix2y), ceil(Int, pix3y))
                k_min = min(floor(Int, pix1z), floor(Int, pix2z), floor(Int, pix3z))
                k_max = max(ceil(Int, pix1z), ceil(Int, pix2z), ceil(Int, pix3z))

                if i_max > i_min && j_max > j_min &&
                   (i_min <= i) && (i <= i_max - 1) && (j_min <= j) && (j <= j_max - 1)
                    normal_x, normal_y, normal_z =
                        _rastergpu_normal(pix1x, pix1y, pix1z, pix2x, pix2y, pix2z, pix3x, pix3y, pix3z)
                    slope_x, slope_y =
                        if abs(normal_z) > 1.0f-5
                            normal_x / normal_z, normal_y / normal_z
                        else
                            dirz * normal_x, dirz * normal_y
                        end
                    z0 = pix1z + slope_x * (pix1x - Float32(i_min)) + slope_y * (pix1y - Float32(j_min))
                    ymin_i = j_max
                    ymax_i = j_min
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix1x, pix1y, pix2x, pix2y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix2x, pix2y, pix3x, pix3y)
                    ymin_i, ymax_i = _rastergpu_edge_bounds(i, ymin_i, ymax_i, pix3x, pix3y, pix1x, pix1y)

                    if (ymin_i <= j) && (j <= ymax_i - 1)
                        ni = i - i_min
                        zi = z0 - slope_x * Float32(ni)
                        nj = j - j_min
                        zpix = zi - slope_y * Float32(nj)
                        zpix = clamp(zpix, Float32(k_min), Float32(k_max))
                        @atomic :monotonic node_counts[node_idx] += Int32(1)
                        observed_hits += 1
                        if n_hits < max_hits
                            n_hits += 1
                            insert_at = n_hits
                            while insert_at > 1 &&
                                  (local_heights[insert_at-1] < zpix ||
                                   (local_heights[insert_at-1] == zpix && local_nodes[insert_at-1] > UInt32(node_idx)))
                                local_heights[insert_at] = local_heights[insert_at-1]
                                local_nodes[insert_at] = local_nodes[insert_at-1]
                                insert_at -= 1
                            end
                            local_heights[insert_at] = zpix
                            local_nodes[insert_at] = UInt32(node_idx)
                        else
                            overflow[pixel_idx] = true
                        end
                    end
                end
            end

            counts[pixel_idx] = Int32(observed_hits)
            if n_hits > 0
                area_active = true
                for h in 1:n_hits
                    node_idx = Int(local_nodes[h])
                    1 <= node_idx <= n_nodes || continue
                    area_active || continue

                    ratio = 1.0f0
                    if area_ratio
                        pixels_area = projected_pixels_area[node_idx]
                        ratio = pixels_area > 0.0f0 ? projected_mesh_area[node_idx] / pixels_area : 1.0f0
                    end

                    if virtual_node_mask[node_idx]
                        @atomic :monotonic sector_area[node_idx] += pixel_area * ratio
                        continue
                    end

                    transparency = node_transparency[node_idx]
                    intercepted_fraction = 1.0f0 - transparency
                    if intercepted_fraction > 0.0f0
                        @atomic :monotonic sector_area[node_idx] += pixel_area * intercepted_fraction * ratio
                    end
                    area_active = transparency > 0.0f0
                end
            end

            if accumulate_dense_edges && n_hits > 1
                nearest_above_idx = 0
                for h in 1:(n_hits - 1)
                    current_idx = Int(local_nodes[h])
                    if !virtual_node_mask[current_idx]
                        nearest_above_idx = current_idx
                    end

                    from_below_idx = 0
                    for k in (h + 1):n_hits
                        candidate_idx = Int(local_nodes[k])
                        if !virtual_node_mask[candidate_idx]
                            from_below_idx = candidate_idx
                            break
                        end
                    end

                    if from_below_idx != 0 &&
                       !(pavement_node_mask[current_idx] && pavement_node_mask[from_below_idx])
                        dense_idx = (current_idx - 1) * n_nodes + from_below_idx
                        @atomic :monotonic dense_edge_counts[dense_idx] += 1
                    end

                    next_idx = Int(local_nodes[h+1])
                    if nearest_above_idx != 0 &&
                       !(pavement_node_mask[next_idx] && pavement_node_mask[nearest_above_idx])
                        dense_idx = (next_idx - 1) * n_nodes + nearest_above_idx
                        @atomic :monotonic dense_edge_counts[dense_idx] += 1
                    end
                end
            end
        end
    end
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
    return _direction_projection_cached(geometry, direction, options, prepared.cache_ctx; upper_hit=prepared.upper_hit)
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
    projection =
        _direction_projection_cached(vertices, faces, face2node, direction, options, plotbox, cache_ctx; upper_hit=upper_hit)
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
    unit_scale = 1.0f0
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
        unit_scale,
        stack_type,
    )
    node_hits_map = _dense_int_node_map(node_ids, node_hits)
    projected_mesh_area_map = _dense_float_node_map(node_ids, projected_mesh_area)
    projected_pixels_area_map = _dense_float_node_map(node_ids, projected_pixels_area)
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
    unit_scale = 1.0f0
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
        unit_scale,
        stack_type,
    )
    use_flat_hit_pool && (pixel_hits = _finalize_flat_pixel_hits(pixel_hits))
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

function _direction_projection_cached(vertices, faces, face2node, direction, options::LightOptions, plotbox, cache_ctx; upper_hit::Bool=false)
    if cache_ctx === nothing
        return _direction_projection(vertices, faces, face2node, direction, options, plotbox; upper_hit=upper_hit)
    end

    path = _projection_cache_path(cache_ctx, direction, upper_hit)
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

    path = _projection_cache_path(cache_ctx, direction, upper_hit)
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
    options::LightOptions,
)
    raw_vertices = GeometryBasics.decompose(GeometryBasics.Point3, scene.merged_mesh)
    vertices = [StaticArrays.SVector{3,Float64}(v[1], v[2], v[3]) for v in raw_vertices]
    all_faces = collect(GeometryBasics.decompose(PlantGeom.Face3, scene.merged_mesh))
    all_face2node = collect(scene.face2node)

    ignored = _ignored_group_types(models)
    all_node_ids = unique(all_face2node)
    all_node_group, all_node_type = _scene_group_type_maps(scene, all_node_ids)
    missing = Set{Tuple{String,String}}()
    ignored_nodes = _ignored_node_ids(all_node_ids, all_node_group, all_node_type, ignored)
    for nid in all_node_ids
        nid in ignored_nodes && continue
        group = strip(get(all_node_group, nid, ""))
        type_name = strip(get(all_node_type, nid, ""))
        _type_model(models, group, type_name) === nothing && push!(missing, (group, type_name))
    end
    if !isempty(missing)
        details = join(["($(repr(g)), $(repr(t)))" for (g, t) in sort!(collect(missing))], ", ")
        error("Missing models for simulated geometry nodes: $details")
    end
    faces, face2node =
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
    isempty(face2node) && error("No intercepting geometry left after applying ignore rules.")

    plotbox = _plotbox(scene, vertices, options.pixel_size)

    node_ids = unique(face2node)
    node_index = Dict{Int,Int}(nid => i for (i, nid) in enumerate(node_ids))
    face2node_index = [node_index[nid] for nid in face2node]
    node_group = Dict{Int,String}(nid => get(all_node_group, nid, "") for nid in node_ids)
    node_group_by_index = [get(node_group, nid, "") for nid in node_ids]
    pavement_node_mask = [group == "pavement" for group in node_group_by_index]
    node_type = Dict{Int,String}(nid => get(all_node_type, nid, "") for nid in node_ids)
    node_type_by_index = [get(node_type, nid, "") for nid in node_ids]
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
)
    geometry = _scene_geometry_for_interception(scene, models, options)
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

function _rastergpu_backend_max_buffer_length(backend)
    occursin("Metal", string(typeof(backend))) || return nothing
    for mod in Base.loaded_modules_array()
        if nameof(mod) == :Metal && isdefined(mod, :device)
            try
                device = Base.invokelatest(getproperty(mod, :device))
                hasproperty(device, :maxBufferLength) || return nothing
                limit = getproperty(device, :maxBufferLength)
                return Int128(limit)
            catch
                return nothing
            end
        end
    end
    return nothing
end

function _rastergpu_backend_maybe_collect(backend)
    occursin("Metal", string(typeof(backend))) || return nothing
    for mod in Base.loaded_modules_array()
        if nameof(mod) == :Metal && isdefined(mod, :maybe_collect)
            try
                Base.invokelatest(getproperty(mod, :maybe_collect); will_block=true)
            catch
            end
            return nothing
        end
    end
    return nothing
end

function _rastergpu_backend_memory_info(backend)
    occursin("Metal", string(typeof(backend))) || return nothing
    for mod in Base.loaded_modules_array()
        if nameof(mod) == :Metal && isdefined(mod, :device)
            try
                device = Base.invokelatest(getproperty(mod, :device))
                working_set =
                    hasproperty(device, :recommendedMaxWorkingSetSize) ?
                    Int128(getproperty(device, :recommendedMaxWorkingSetSize)) :
                    Int128(0)
                allocated =
                    hasproperty(device, :currentAllocatedSize) ?
                    Int128(getproperty(device, :currentAllocatedSize)) :
                    Int128(0)
                free = working_set > 0 ? max(working_set - allocated, Int128(0)) : Int128(0)
                return (working_set=working_set, allocated=allocated, free=free)
            catch
                return nothing
            end
        end
    end
    return nothing
end

_rastergpu_limit_label(limit) = limit === nothing ? "unknown" : string(limit)
_rastergpu_memory_label(x) = x === nothing ? "unknown" : string(x)

function _rastergpu_config_context(
    config::RasterGPUBackendConfig;
    n_pixels=nothing,
    n_tiles=nothing,
    tile_size=nothing,
    tile_face_capacity=nothing,
    max_hits=nothing,
    edge_key_capacity=nothing,
    dense_pairs=nothing,
    stackless_top_hit=nothing,
    fused_dense_edges=nothing,
)
    fields = String[
        "max_hits_per_pixel=$(config.max_hits_per_pixel)",
        "tile_size=$(config.tile_size)",
        "tile_face_capacity=$(config.tile_face_capacity)",
        "top_hit_tile_size=$(config.top_hit_tile_size)",
        "top_hit_tile_face_capacity=$(config.top_hit_tile_face_capacity)",
        "edge_accumulation=$(config.edge_accumulation)",
        "dense_edge_limit_bytes=$(config.dense_edge_limit_bytes)",
    ]
    n_pixels === nothing || push!(fields, "n_pixels=$n_pixels")
    n_tiles === nothing || push!(fields, "n_tiles=$n_tiles")
    tile_size === nothing || push!(fields, "effective_tile_size=$tile_size")
    tile_face_capacity === nothing || push!(fields, "effective_tile_face_capacity=$tile_face_capacity")
    max_hits === nothing || push!(fields, "effective_max_hits_per_pixel=$max_hits")
    edge_key_capacity === nothing || push!(fields, "edge_key_capacity=$edge_key_capacity")
    dense_pairs === nothing || push!(fields, "dense_edge_pairs=$dense_pairs")
    stackless_top_hit === nothing || push!(fields, "stackless_top_hit=$stackless_top_hit")
    fused_dense_edges === nothing || push!(fields, "fused_dense_edges=$fused_dense_edges")
    return join(fields, ", ")
end

function _rastergpu_throw_buffer_size_error(
    buffer_name::AbstractString,
    ::Type{T},
    count::Int128,
    bytes::Int128,
    limit,
    config::RasterGPUBackendConfig;
    reason::AbstractString,
    context::AbstractString,
) where {T}
    error(
        "RasterGPU buffer allocation preflight failed for $buffer_name: $reason. " *
        "requested_count=$count, eltype=$T, requested_bytes=$bytes, " *
        "backend_buffer_limit_bytes=$(_rastergpu_limit_label(limit)). " *
        "Config/context: $context. " *
        "Reduce the controlling RasterGPUBackendConfig fields, for example tile_size, " *
        "tile_face_capacity, top_hit_tile_face_capacity, max_hits_per_pixel, " *
        "edge_accumulation, or dense_edge_limit_bytes.",
    )
end

function _rastergpu_checked_buffer_count(
    buffer_name::AbstractString,
    ::Type{T},
    count,
    config::RasterGPUBackendConfig;
    backend_limit=nothing,
    context::AbstractString="",
) where {T}
    count128 = Int128(count)
    bytes = count128 * Int128(sizeof(T))
    if count128 < 0
        _rastergpu_throw_buffer_size_error(
            buffer_name,
            T,
            count128,
            bytes,
            backend_limit,
            config;
            reason="negative element count",
            context=context,
        )
    end
    if count128 > Int128(typemax(Int))
        _rastergpu_throw_buffer_size_error(
            buffer_name,
            T,
            count128,
            bytes,
            backend_limit,
            config;
            reason="element count exceeds typemax(Int)",
            context=context,
        )
    end
    if bytes > Int128(typemax(Int))
        _rastergpu_throw_buffer_size_error(
            buffer_name,
            T,
            count128,
            bytes,
            backend_limit,
            config;
            reason="requested byte size exceeds typemax(Int)",
            context=context,
        )
    end
    if backend_limit !== nothing && bytes > backend_limit
        _rastergpu_throw_buffer_size_error(
            buffer_name,
            T,
            count128,
            bytes,
            backend_limit,
            config;
            reason="requested byte size exceeds backend buffer limit",
            context=context,
        )
    end
    return Int(count128)
end

function _rastergpu_checked_product(
    name::AbstractString,
    ::Type{T},
    factors,
    config::RasterGPUBackendConfig;
    backend_limit=nothing,
    context::AbstractString="",
) where {T}
    count = Int128(1)
    for factor in factors
        count *= Int128(factor)
    end
    return _rastergpu_checked_buffer_count(
        name,
        T,
        count,
        config;
        backend_limit=backend_limit,
        context=context,
    )
end

_rastergpu_device_allocation_len(n::Int) = max(n, 1)
_rastergpu_device_allocation_len(n::Int128) = max(n, Int128(1))

function _rastergpu_device_buffer_bytes_i128(::Type{T}, count::Integer) where {T}
    return _rastergpu_device_allocation_len(Int128(count)) * Int128(sizeof(T))
end

_rastergpu_allocate(backend, ::Type{T}, n::Integer) where {T} =
    KernelAbstractions.allocate(backend, T, _rastergpu_device_allocation_len(Int(n)); unified=false)

function _rastergpu_device_buffer_bytes(::Type{T}, count::Integer) where {T}
    return _rastergpu_device_buffer_bytes_i128(T, Int(count))
end

function _rastergpu_largest_buffers(buffers, n::Int=6)
    pairs = sort!(collect(buffers); by=last, rev=true)
    isempty(pairs) && return ""
    return join(("$(first(pair))=$(last(pair))" for pair in Iterators.take(pairs, n)), ", ")
end

function _rastergpu_check_total_allocation_bytes(
    buffers,
    config::RasterGPUBackendConfig;
    memory_info=nothing,
    context::AbstractString="",
)
    memory_info === nothing && return sum(last, buffers; init=Int128(0))

    requested = sum(last, buffers; init=Int128(0))
    budget = _rastergpu_total_allocation_budget(memory_info)
    budget === nothing && return requested

    if requested > budget
        error(
            "RasterGPU aggregate device allocation preflight failed: requested_device_bytes=$requested, " *
            "available_budget_bytes=$budget, metal_working_set_bytes=$(memory_info.working_set), " *
            "metal_current_allocated_bytes=$(memory_info.allocated), metal_free_bytes=$(memory_info.free). " *
            "Largest buffers: $(_rastergpu_largest_buffers(buffers)). " *
            "Config/context: $context. " *
            "Reduce RasterGPUBackendConfig(tile_face_capacity, tile_size, max_hits_per_pixel), " *
            "use edge_accumulation=:auto so dense topology can replace oversized sparse edge keys when smaller, " *
            "or use a coarser pixel_size.",
        )
    end
    return requested
end

function _rastergpu_total_allocation_budget(memory_info)
    memory_info === nothing && return nothing
    working_set = memory_info.working_set
    working_set <= 0 && return nothing
    allocated = memory_info.allocated
    free = memory_info.free
    working_budget = max(fld(working_set * Int128(3), Int128(4)) - allocated, Int128(0))
    free_budget = fld(free * Int128(17), Int128(20))
    return min(working_budget, free_budget)
end

function _rastergpu_choose_dense_edge_accumulation(
    config::RasterGPUBackendConfig;
    stackless_top_hit::Bool,
    dense_pairs::Int128,
    dense_bytes::Int128,
    sparse_total_bytes::Int128,
    dense_total_bytes::Int128,
    memory_info=nothing,
)
    stackless_top_hit && return false
    dense_pairs <= Int128(typemax(Int)) || return false

    dense_under_config_limit = dense_bytes <= Int128(config.dense_edge_limit_bytes)
    dense_smaller_than_sparse = dense_total_bytes <= sparse_total_bytes

    if config.edge_accumulation == :dense_atomic
        return dense_under_config_limit
    elseif config.edge_accumulation == :auto
        !dense_smaller_than_sparse && return false
        dense_under_config_limit && return true
        budget = _rastergpu_total_allocation_budget(memory_info)
        budget === nothing && return false
        return dense_total_bytes <= budget
    end

    return false
end

function _rastergpu_scene_data(
    prepared::PreparedInterceptionData,
    config::RasterGPUBackendConfig,
)
    backend = config.backend
    KernelAbstractions.supports_atomics(backend) ||
        error("RasterGPUBackend requires a KernelAbstractions backend with atomic support.")

    geometry = prepared.geometry
    n_vertices = length(geometry.vertices)
    n_faces = length(geometry.faces)
    n_nodes = length(geometry.node_ids)
    buffer_limit = _rastergpu_backend_max_buffer_length(backend)
    max_hits = config.max_hits_per_pixel
    tile_size = prepared.upper_hit ? config.top_hit_tile_size : config.tile_size
    tile_face_capacity = prepared.upper_hit ? config.top_hit_tile_face_capacity : config.tile_face_capacity

    n_pixels = _rastergpu_checked_product(
        "n_pixels",
        UInt8,
        (geometry.plotbox.nx, geometry.plotbox.ny),
        config;
        context=_rastergpu_config_context(config; tile_size=tile_size, tile_face_capacity=tile_face_capacity),
    )
    n_tiles_x = cld(geometry.plotbox.nx, tile_size)
    n_tiles_y = cld(geometry.plotbox.ny, tile_size)
    n_tiles = _rastergpu_checked_product(
        "n_tiles",
        UInt8,
        (n_tiles_x, n_tiles_y),
        config;
        context=_rastergpu_config_context(config; n_pixels=n_pixels, tile_size=tile_size, tile_face_capacity=tile_face_capacity),
    )
    stackless_top_hit = prepared.upper_hit && tile_size == 1
    tile_candidate_len = _rastergpu_checked_product(
        "RasterGPU tile_faces_dev/tile_unwrapped_i_dev/tile_unwrapped_j_dev",
        Int32,
        (n_tiles, tile_face_capacity),
        config;
        backend_limit=buffer_limit,
        context=_rastergpu_config_context(
            config;
            n_pixels=n_pixels,
            n_tiles=n_tiles,
            tile_size=tile_size,
            tile_face_capacity=tile_face_capacity,
            max_hits=max_hits,
            stackless_top_hit=stackless_top_hit,
        ),
    )
    base_context = _rastergpu_config_context(
        config;
        n_pixels=n_pixels,
        n_tiles=n_tiles,
        tile_size=tile_size,
        tile_face_capacity=tile_face_capacity,
        max_hits=max_hits,
        stackless_top_hit=stackless_top_hit,
    )
    _rastergpu_checked_buffer_count(
        "RasterGPU counts_dev",
        Int32,
        n_pixels,
        config;
        backend_limit=buffer_limit,
        context=base_context,
    )
    _rastergpu_checked_buffer_count(
        "RasterGPU overflow_dev",
        Bool,
        n_pixels,
        config;
        backend_limit=buffer_limit,
        context=base_context,
    )
    _rastergpu_checked_buffer_count(
        "RasterGPU tile_counts_dev/tile_overflow_dev",
        Int32,
        n_tiles,
        config;
        backend_limit=buffer_limit,
        context=base_context,
    )
    max_edges = max(2 * (max_hits - 1), 0)
    edge_key_capacity_i128 = stackless_top_hit ? Int128(0) : Int128(n_pixels) * Int128(max_edges)
    dense_pairs = Int128(n_nodes) * Int128(n_nodes)
    dense_bytes = dense_pairs * Int128(sizeof(Int32))
    sparse_edge_bytes =
        _rastergpu_device_buffer_bytes_i128(UInt64, edge_key_capacity_i128) +
        _rastergpu_device_buffer_bytes_i128(Int32, n_pixels)
    full_stack_len_i128 = stackless_top_hit ? Int128(0) : Int128(n_pixels) * Int128(max_hits)
    full_stack_bytes =
        _rastergpu_device_buffer_bytes_i128(UInt32, full_stack_len_i128) +
        _rastergpu_device_buffer_bytes_i128(Float32, full_stack_len_i128)
    stackless_stack_bytes =
        _rastergpu_device_buffer_bytes_i128(UInt32, 1) +
        _rastergpu_device_buffer_bytes_i128(Float32, 1)
    fused_dense_candidate =
        !stackless_top_hit &&
        dense_pairs <= Int128(typemax(Int)) &&
        tile_size == 1 &&
        _rastergpu_fused_dense_supported(config)

    common_device_buffers = Pair{String,Int128}[
        "vertex_x_dev" => _rastergpu_device_buffer_bytes(Float32, n_vertices),
        "vertex_y_dev" => _rastergpu_device_buffer_bytes(Float32, n_vertices),
        "vertex_z_dev" => _rastergpu_device_buffer_bytes(Float32, n_vertices),
        "face_i_dev" => _rastergpu_device_buffer_bytes(Int32, n_faces),
        "face_j_dev" => _rastergpu_device_buffer_bytes(Int32, n_faces),
        "face_k_dev" => _rastergpu_device_buffer_bytes(Int32, n_faces),
        "face2node_index_dev" => _rastergpu_device_buffer_bytes(Int32, n_faces),
        "node_transparency_dev" => _rastergpu_device_buffer_bytes(Float32, n_nodes),
        "virtual_node_mask_dev" => _rastergpu_device_buffer_bytes(Bool, n_nodes),
        "pavement_node_mask_dev" => _rastergpu_device_buffer_bytes(Bool, n_nodes),
        "node_ids_dev" => _rastergpu_device_buffer_bytes(Int, n_nodes),
        "counts_dev" => _rastergpu_device_buffer_bytes(Int32, n_pixels),
        "overflow_dev" => _rastergpu_device_buffer_bytes(Bool, n_pixels),
        "node_counts_dev" => _rastergpu_device_buffer_bytes(Int32, n_nodes),
        "projected_mesh_area_dev" => _rastergpu_device_buffer_bytes(Float32, n_nodes),
        "projected_pixels_area_dev" => _rastergpu_device_buffer_bytes(Float32, n_nodes),
        "sector_area_dev" => _rastergpu_device_buffer_bytes(Float32, n_nodes),
        "tile_counts_dev" => _rastergpu_device_buffer_bytes(Int32, n_tiles),
        "tile_faces_dev" => _rastergpu_device_buffer_bytes(Int32, tile_candidate_len),
        "tile_unwrapped_i_dev" => _rastergpu_device_buffer_bytes(Int32, tile_candidate_len),
        "tile_unwrapped_j_dev" => _rastergpu_device_buffer_bytes(Int32, tile_candidate_len),
        "tile_overflow_dev" => _rastergpu_device_buffer_bytes(Bool, n_tiles),
    ]
    common_device_bytes = sum(last, common_device_buffers; init=Int128(0))
    sparse_total_bytes = common_device_bytes + full_stack_bytes + sparse_edge_bytes
    dense_total_bytes =
        common_device_bytes +
        (fused_dense_candidate ? stackless_stack_bytes : full_stack_bytes) +
        _rastergpu_device_buffer_bytes_i128(Int32, dense_pairs)

    _rastergpu_backend_maybe_collect(backend)
    memory_info = _rastergpu_backend_memory_info(backend)
    dense_enabled = _rastergpu_choose_dense_edge_accumulation(
        config;
        stackless_top_hit=stackless_top_hit,
        dense_pairs=dense_pairs,
        dense_bytes=dense_bytes,
        sparse_total_bytes=sparse_total_bytes,
        dense_total_bytes=dense_total_bytes,
        memory_info=memory_info,
    )
    if !stackless_top_hit && config.edge_accumulation == :dense_atomic && !dense_enabled
        error(
            "RasterGPU edge_accumulation=:dense_atomic requires a dense $n_nodes x $n_nodes edge matrix " *
            "within dense_edge_limit_bytes=$(config.dense_edge_limit_bytes). " *
            "Requested dense_edge_pairs=$dense_pairs, dense_edge_bytes=$dense_bytes. " *
            "Increase RasterGPUBackendConfig(dense_edge_limit_bytes=...) or use edge_accumulation=:auto.",
        )
    end
    if dense_enabled
        _rastergpu_checked_buffer_count(
            "RasterGPU dense_edge_counts_dev/dense_edge_counts_host",
            Int32,
            dense_pairs,
            config;
            backend_limit=buffer_limit,
            context=_rastergpu_config_context(
                config;
                n_pixels=n_pixels,
                n_tiles=n_tiles,
                tile_size=tile_size,
                tile_face_capacity=tile_face_capacity,
                max_hits=max_hits,
                edge_key_capacity=edge_key_capacity_i128,
                dense_pairs=dense_pairs,
                stackless_top_hit=stackless_top_hit,
            ),
        )
    end
    sparse_enabled =
        !stackless_top_hit &&
        !dense_enabled &&
        config.edge_accumulation != :dense_atomic &&
        edge_key_capacity_i128 > 0
    edge_key_capacity =
        sparse_enabled ?
        _rastergpu_checked_product(
            "RasterGPU edge_keys_dev",
            UInt64,
            (n_pixels, max_edges),
            config;
            backend_limit=buffer_limit,
            context=base_context,
        ) :
        0
    fused_dense_enabled = dense_enabled && tile_size == 1 && _rastergpu_fused_dense_supported(config)
    stackless_projection = stackless_top_hit || fused_dense_enabled
    stack_len = stackless_projection ? 0 : _rastergpu_checked_product(
        "RasterGPU nodes_dev/heights_dev",
        UInt32,
        (n_pixels, max_hits),
        config;
        backend_limit=buffer_limit,
        context=_rastergpu_config_context(
            config;
            n_pixels=n_pixels,
            n_tiles=n_tiles,
            tile_size=tile_size,
            tile_face_capacity=tile_face_capacity,
            max_hits=max_hits,
            edge_key_capacity=edge_key_capacity_i128,
            dense_pairs=dense_pairs,
            stackless_top_hit=stackless_top_hit,
            fused_dense_edges=fused_dense_enabled,
        ),
    )
    stack_alloc_len = stackless_projection ? 1 : stack_len

    device_buffers = Pair{String,Int128}[
        common_device_buffers...,
        "nodes_dev" => _rastergpu_device_buffer_bytes(UInt32, stack_alloc_len),
        "heights_dev" => _rastergpu_device_buffer_bytes(Float32, stack_alloc_len),
    ]
    dense_enabled && push!(
        device_buffers,
        "dense_edge_counts_dev" => _rastergpu_device_buffer_bytes(Int32, Int(dense_pairs)),
    )
    sparse_enabled && begin
        push!(device_buffers, "edge_keys_dev" => _rastergpu_device_buffer_bytes(UInt64, edge_key_capacity))
        push!(device_buffers, "edge_key_counts_dev" => _rastergpu_device_buffer_bytes(Int32, n_pixels))
    end
    _rastergpu_check_total_allocation_bytes(
        device_buffers,
        config;
        memory_info=memory_info,
        context=_rastergpu_config_context(
            config;
            n_pixels=n_pixels,
            n_tiles=n_tiles,
            tile_size=tile_size,
            tile_face_capacity=tile_face_capacity,
            max_hits=max_hits,
            edge_key_capacity=edge_key_capacity_i128,
            dense_pairs=dense_pairs,
            stackless_top_hit=stackless_top_hit,
            fused_dense_edges=fused_dense_enabled,
        ),
    )

    vertex_x = Float32[v[1] for v in geometry.vertices]
    vertex_y = Float32[v[2] for v in geometry.vertices]
    vertex_z = Float32[v[3] for v in geometry.vertices]
    face_i = Int32[f[1] for f in geometry.faces]
    face_j = Int32[f[2] for f in geometry.faces]
    face_k = Int32[f[3] for f in geometry.faces]
    face2node_index = Int32.(geometry.face2node_index)

    vertex_x_dev = _rastergpu_allocate(backend, Float32, n_vertices)
    vertex_y_dev = _rastergpu_allocate(backend, Float32, n_vertices)
    vertex_z_dev = _rastergpu_allocate(backend, Float32, n_vertices)
    face_i_dev = _rastergpu_allocate(backend, Int32, n_faces)
    face_j_dev = _rastergpu_allocate(backend, Int32, n_faces)
    face_k_dev = _rastergpu_allocate(backend, Int32, n_faces)
    face2node_index_dev = _rastergpu_allocate(backend, Int32, n_faces)
    node_transparency_dev = _rastergpu_allocate(backend, Float32, n_nodes)
    virtual_node_mask_dev = _rastergpu_allocate(backend, Bool, n_nodes)
    pavement_node_mask_dev = _rastergpu_allocate(backend, Bool, n_nodes)
    node_ids_dev = _rastergpu_allocate(backend, Int, n_nodes)

    n_vertices > 0 && begin
        KernelAbstractions.copyto!(backend, vertex_x_dev, vertex_x)
        KernelAbstractions.copyto!(backend, vertex_y_dev, vertex_y)
        KernelAbstractions.copyto!(backend, vertex_z_dev, vertex_z)
    end
    n_faces > 0 && begin
        KernelAbstractions.copyto!(backend, face_i_dev, face_i)
        KernelAbstractions.copyto!(backend, face_j_dev, face_j)
        KernelAbstractions.copyto!(backend, face_k_dev, face_k)
        KernelAbstractions.copyto!(backend, face2node_index_dev, face2node_index)
    end
    n_nodes > 0 && begin
        KernelAbstractions.copyto!(backend, node_transparency_dev, Float32.(prepared.node_transparency_by_index))
        KernelAbstractions.copyto!(backend, virtual_node_mask_dev, prepared.virtual_node_mask)
        KernelAbstractions.copyto!(backend, pavement_node_mask_dev, geometry.pavement_node_mask)
        KernelAbstractions.copyto!(backend, node_ids_dev, geometry.node_ids)
    end

    counts_dev = _rastergpu_allocate(backend, Int32, n_pixels)
    nodes_dev = _rastergpu_allocate(backend, UInt32, stack_alloc_len)
    heights_dev = _rastergpu_allocate(backend, Float32, stack_alloc_len)
    overflow_dev = _rastergpu_allocate(backend, Bool, n_pixels)
    node_counts_dev = _rastergpu_allocate(backend, Int32, n_nodes)
    projected_mesh_area_dev = _rastergpu_allocate(backend, Float32, n_nodes)
    projected_pixels_area_dev = _rastergpu_allocate(backend, Float32, n_nodes)
    sector_area_dev = _rastergpu_allocate(backend, Float32, n_nodes)
    tile_counts_dev = _rastergpu_allocate(backend, Int32, n_tiles)
    tile_faces_dev = _rastergpu_allocate(backend, Int32, tile_candidate_len)
    tile_unwrapped_i_dev = _rastergpu_allocate(backend, Int32, tile_candidate_len)
    tile_unwrapped_j_dev = _rastergpu_allocate(backend, Int32, tile_candidate_len)
    tile_overflow_dev = _rastergpu_allocate(backend, Bool, n_tiles)

    dense_edge_counts_dev =
        dense_enabled ? _rastergpu_allocate(backend, Int32, Int(dense_pairs)) : nothing
    dense_edge_counts_host = dense_enabled ? zeros(Int32, Int(dense_pairs)) : nothing
    edge_keys_dev = sparse_enabled ? _rastergpu_allocate(backend, UInt64, edge_key_capacity) : nothing
    edge_key_counts_dev = sparse_enabled ? _rastergpu_allocate(backend, Int32, n_pixels) : nothing

    return RasterGPUSceneData(
        prepared,
        backend,
        vertex_x_dev,
        vertex_y_dev,
        vertex_z_dev,
        face_i_dev,
        face_j_dev,
        face_k_dev,
        face2node_index_dev,
        node_transparency_dev,
        virtual_node_mask_dev,
        pavement_node_mask_dev,
        node_ids_dev,
        counts_dev,
        nodes_dev,
        heights_dev,
        overflow_dev,
        node_counts_dev,
        projected_mesh_area_dev,
        projected_pixels_area_dev,
        sector_area_dev,
        dense_edge_counts_dev,
        edge_keys_dev,
        edge_key_counts_dev,
        tile_counts_dev,
        tile_faces_dev,
        tile_unwrapped_i_dev,
        tile_unwrapped_j_dev,
        tile_overflow_dev,
        zeros(Int32, n_pixels),
        zeros(UInt32, stack_len),
        zeros(Float32, stack_len),
        fill(false, n_pixels),
        zeros(Int32, n_nodes),
        zeros(Float32, n_nodes),
        zeros(Float32, n_nodes),
        zeros(Float32, n_nodes),
        dense_edge_counts_host,
        UInt64[],
        Int32[],
        UInt64[],
        zeros(Int32, n_tiles),
        fill(false, n_tiles),
        config.workgroupsize,
        max_hits,
        tile_size,
        tile_face_capacity,
        config.edge_accumulation,
        config.dense_edge_limit_bytes,
        config.validate,
    )
end

function _rastergpu_clear_direction!(data::RasterGPUSceneData)
    geometry = data.prepared.geometry
    n_pixels = geometry.plotbox.nx * geometry.plotbox.ny
    n_nodes = length(geometry.node_ids)
    ndrange = max(n_pixels, n_nodes)
    ndrange == 0 && return data
    kernel = _rastergpu_clear_direction_kernel!(data.backend, data.workgroupsize)
    kernel(
        data.counts_dev,
        data.overflow_dev,
        data.node_counts_dev,
        data.projected_mesh_area_dev,
        data.projected_pixels_area_dev,
        data.sector_area_dev,
        n_pixels,
        n_nodes;
        ndrange=ndrange,
    )
    KernelAbstractions.synchronize(data.backend)
    return data
end

function _rastergpu_clear_tile_bins!(data::RasterGPUSceneData)
    n_tiles = length(data.tile_counts_host)
    n_tiles == 0 && return data
    kernel = _rastergpu_clear_tile_bins_kernel!(data.backend, data.workgroupsize)
    kernel(
        data.tile_counts_dev,
        data.tile_overflow_dev,
        n_tiles;
        ndrange=n_tiles,
    )
    KernelAbstractions.synchronize(data.backend)
    return data
end

function _rastergpu_throw_on_overflow!(data::RasterGPUSceneData)
    copyto!(data.overflow_host, data.overflow_dev)
    pixel = findfirst(data.overflow_host)
    pixel === nothing && return nothing
    copyto!(data.counts_host, data.counts_dev)
    observed_count = Int(data.counts_host[pixel])
    max_count, max_pixel = findmax(data.counts_host)
    error(
        "RasterGPU max_hits_per_pixel=$(data.max_hits_per_pixel) exceeded for pixel $pixel. " *
        "Observed $observed_count hits at that pixel; observed_max_hits_per_pixel=$(Int(max_count)) " *
        "at pixel $max_pixel. Increase RasterGPUBackendConfig(max_hits_per_pixel=...).",
    )
end

function _rastergpu_throw_on_tile_overflow!(data::RasterGPUSceneData)
    copyto!(data.tile_overflow_host, data.tile_overflow_dev)
    tile = findfirst(data.tile_overflow_host)
    tile === nothing && return nothing
    copyto!(data.tile_counts_host, data.tile_counts_dev)
    observed_count = Int(data.tile_counts_host[tile])
    max_count, max_tile = findmax(data.tile_counts_host)
    capacity_field = data.prepared.upper_hit ? "top_hit_tile_face_capacity" : "tile_face_capacity"
    tile_size_field = data.prepared.upper_hit ? "top_hit_tile_size" : "tile_size"
    error(
        "RasterGPU $capacity_field=$(data.tile_face_capacity) exceeded for tile $tile. " *
        "Observed $observed_count candidate faces at that tile; observed_max_tile_face_capacity=$(Int(max_count)) " *
        "at tile $max_tile. " *
        "Increase RasterGPUBackendConfig($capacity_field=...) or use a smaller $tile_size_field.",
    )
end

function _rastergpu_project_direction!(
    data::RasterGPUSceneData,
    direction,
    options::LightOptions;
    accumulate_dense_edges::Bool=false,
)
    geometry = data.prepared.geometry
    plotbox = geometry.plotbox
    n_faces = length(geometry.faces)
    n_pixels = plotbox.nx * plotbox.ny
    n_nodes = length(geometry.node_ids)
    _rastergpu_clear_direction!(data)

    if n_faces > 0 && n_pixels > 0 && n_nodes > 0
        _rastergpu_clear_tile_bins!(data)
        n_tiles_x = cld(plotbox.nx, data.tile_size)
        n_tiles_y = cld(plotbox.ny, data.tile_size)
        n_tiles = n_tiles_x * n_tiles_y
        use_top_hit_kernel =
            data.prepared.upper_hit &&
            (!options.toricity ||
             (plotbox.nx % data.tile_size == 0 && plotbox.ny % data.tile_size == 0))
        use_fused_dense_kernel = _rastergpu_use_fused_dense_edges(data)

        if data.tile_size == 1
            bin_kernel = _rastergpu_bin_faces_to_covered_pixel_tiles_kernel!(data.backend, data.workgroupsize)
            bin_kernel(
                data.tile_counts_dev,
                data.tile_faces_dev,
                data.tile_unwrapped_i_dev,
                data.tile_unwrapped_j_dev,
                data.tile_overflow_dev,
                data.projected_mesh_area_dev,
                data.projected_pixels_area_dev,
                data.vertex_x_dev,
                data.vertex_y_dev,
                data.vertex_z_dev,
                data.face_i_dev,
                data.face_j_dev,
                data.face_k_dev,
                data.face2node_index_dev,
                Float32(direction[1]),
                Float32(direction[2]),
                Float32(direction[3]),
                Float32(plotbox.origin_x),
                Float32(plotbox.origin_y),
                Float32(plotbox.pix_x),
                Float32(plotbox.pix_y),
                Float32(plotbox.pixel_area),
                plotbox.nx,
                plotbox.ny,
                data.tile_face_capacity,
                options.toricity;
                ndrange=n_faces,
            )
        else
            bin_kernel = _rastergpu_bin_faces_to_tiles_kernel!(data.backend, data.workgroupsize)
            bin_kernel(
                data.tile_counts_dev,
                data.tile_faces_dev,
                data.tile_unwrapped_i_dev,
                data.tile_unwrapped_j_dev,
                data.tile_overflow_dev,
                data.projected_mesh_area_dev,
                data.projected_pixels_area_dev,
                data.vertex_x_dev,
                data.vertex_y_dev,
                data.vertex_z_dev,
                data.face_i_dev,
                data.face_j_dev,
                data.face_k_dev,
                data.face2node_index_dev,
                Float32(direction[1]),
                Float32(direction[2]),
                Float32(direction[3]),
                Float32(plotbox.origin_x),
                Float32(plotbox.origin_y),
                Float32(plotbox.pix_x),
                Float32(plotbox.pix_y),
                Float32(plotbox.pixel_area),
                plotbox.nx,
                plotbox.ny,
                n_tiles_x,
                n_tiles_y,
                data.tile_size,
                data.tile_face_capacity,
                options.toricity;
                ndrange=n_faces,
            )
        end
        KernelAbstractions.synchronize(data.backend)
        _rastergpu_throw_on_tile_overflow!(data)

        if use_fused_dense_kernel
            if accumulate_dense_edges
                dense_pairs = length(data.dense_edge_counts_dev)
                if dense_pairs > 0
                    clear_dense_kernel = _rastergpu_clear_dense_counts_kernel!(data.backend, data.workgroupsize)
                    clear_dense_kernel(data.dense_edge_counts_dev; ndrange=dense_pairs)
                    KernelAbstractions.synchronize(data.backend)
                end
            end
            project_kernel = _rastergpu_project_tile_bins_fused_dense_kernel!(data.backend, data.workgroupsize)
            params = RasterGPUFusedProjectionParams(
                Float32(direction[1]),
                Float32(direction[2]),
                Float32(direction[3]),
                Float32(plotbox.origin_x),
                Float32(plotbox.origin_y),
                Float32(plotbox.pix_x),
                Float32(plotbox.pix_y),
                Float32(plotbox.pixel_area),
                plotbox.nx,
                plotbox.ny,
                n_tiles_x,
                data.tile_size,
                data.tile_face_capacity,
                options.toricity,
                options.area_ratio,
                data.max_hits_per_pixel,
                n_nodes,
                accumulate_dense_edges,
            )
            project_kernel(
                data.counts_dev,
                data.overflow_dev,
                data.node_counts_dev,
                data.sector_area_dev,
                data.dense_edge_counts_dev,
                data.projected_mesh_area_dev,
                data.projected_pixels_area_dev,
                data.tile_counts_dev,
                data.tile_faces_dev,
                data.tile_unwrapped_i_dev,
                data.tile_unwrapped_j_dev,
                data.vertex_x_dev,
                data.vertex_y_dev,
                data.vertex_z_dev,
                data.face_i_dev,
                data.face_j_dev,
                data.face_k_dev,
                data.face2node_index_dev,
                data.virtual_node_mask_dev,
                data.pavement_node_mask_dev,
                data.node_transparency_dev,
                params;
                ndrange=n_tiles * data.tile_size * data.tile_size,
            )
            KernelAbstractions.synchronize(data.backend)
            _rastergpu_throw_on_overflow!(data)
        elseif use_top_hit_kernel
            project_kernel = _rastergpu_project_tile_bins_top_hit_kernel!(data.backend, data.workgroupsize)
            project_kernel(
                data.node_counts_dev,
                data.sector_area_dev,
                data.projected_mesh_area_dev,
                data.projected_pixels_area_dev,
                data.tile_counts_dev,
                data.tile_faces_dev,
                data.tile_unwrapped_i_dev,
                data.tile_unwrapped_j_dev,
                data.vertex_x_dev,
                data.vertex_y_dev,
                data.vertex_z_dev,
                data.face_i_dev,
                data.face_j_dev,
                data.face_k_dev,
                data.face2node_index_dev,
                Float32(direction[1]),
                Float32(direction[2]),
                Float32(direction[3]),
                Float32(plotbox.origin_x),
                Float32(plotbox.origin_y),
                Float32(plotbox.pix_x),
                Float32(plotbox.pix_y),
                Float32(plotbox.pixel_area),
                plotbox.nx,
                plotbox.ny,
                n_tiles_x,
                data.tile_size,
                data.tile_face_capacity,
                options.toricity,
                options.area_ratio;
                ndrange=n_tiles * data.tile_size * data.tile_size,
            )
            KernelAbstractions.synchronize(data.backend)
        else
            project_kernel = _rastergpu_project_tile_bins_kernel!(data.backend, data.workgroupsize)
            project_kernel(
                data.counts_dev,
                data.nodes_dev,
                data.heights_dev,
                data.overflow_dev,
                data.node_counts_dev,
                data.tile_counts_dev,
                data.tile_faces_dev,
                data.tile_unwrapped_i_dev,
                data.tile_unwrapped_j_dev,
                data.vertex_x_dev,
                data.vertex_y_dev,
                data.vertex_z_dev,
                data.face_i_dev,
                data.face_j_dev,
                data.face_k_dev,
                data.face2node_index_dev,
                Float32(direction[1]),
                Float32(direction[2]),
                Float32(direction[3]),
                Float32(plotbox.origin_x),
                Float32(plotbox.origin_y),
                Float32(plotbox.pix_x),
                Float32(plotbox.pix_y),
                plotbox.nx,
                plotbox.ny,
                n_tiles_x,
                data.tile_size,
                data.tile_face_capacity,
                options.toricity,
                data.max_hits_per_pixel;
                ndrange=n_tiles * data.tile_size * data.tile_size,
            )
            KernelAbstractions.synchronize(data.backend)
            _rastergpu_throw_on_overflow!(data)

            reduce_kernel = _rastergpu_sort_and_reduce_kernel!(data.backend, data.workgroupsize)
            reduce_kernel(
                data.counts_dev,
                data.nodes_dev,
                data.heights_dev,
                data.overflow_dev,
                data.sector_area_dev,
                data.projected_mesh_area_dev,
                data.projected_pixels_area_dev,
                data.virtual_node_mask_dev,
                data.node_transparency_dev,
                Float32(plotbox.pixel_area),
                data.max_hits_per_pixel,
                n_nodes,
                options.area_ratio,
                data.prepared.upper_hit;
                ndrange=n_pixels,
            )
            KernelAbstractions.synchronize(data.backend)
            _rastergpu_throw_on_overflow!(data)
        end
    end

    return data
end

function _rastergpu_copy_direction_response_host!(data::RasterGPUSceneData)
    copyto!(data.node_counts_host, data.node_counts_dev)
    copyto!(data.sector_area_host, data.sector_area_dev)
    copyto!(data.projected_mesh_area_host, data.projected_mesh_area_dev)
    copyto!(data.projected_pixels_area_host, data.projected_pixels_area_dev)
    return (
        node_counts=data.node_counts_host,
        sector_area=data.sector_area_host,
        projected_mesh_area=data.projected_mesh_area_host,
        projected_pixels_area=data.projected_pixels_area_host,
    )
end

function _rastergpu_direction_response!(data::RasterGPUSceneData, direction, options::LightOptions)
    _rastergpu_project_direction!(data, direction, options)
    return _rastergpu_copy_direction_response_host!(data)
end

function _compute_first_order_rastergpu(
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    prepared = data.prepared
    geometry = prepared.geometry
    isempty(prepared.emitter_nodes) || error(
        "RasterGPUBackend does not yet support emitter transfer; refusing CPU fallback.",
    )

    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))

    for (k, sector) in enumerate(turtle.sectors)
        response = _rastergpu_direction_response!(data, sector.direction, options)
        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        @inbounds for idx in eachindex(geometry.node_ids)
            hits_per_node[idx] += Int(response.node_counts[idx])
            pa = Float64(response.sector_area[idx])
            pa <= 0.0 && continue
            projected_area_per_node[idx] += pa
            par_flux != 0.0 && (incident_power_par[idx] += par_flux * pa)
            nir_flux != 0.0 && (incident_power_nir[idx] += nir_flux * pa)
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
    data::RasterGPUSceneData,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    return _compute_first_order_rastergpu(data, turtle, fluxes, options)
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

`backend` accepts either `:raster_cpu`, `:raster_gpu`, or an
`InterceptionBackend` instance.
"""
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
    backend == :raster_gpu && return RasterGPUBackend()
    error("Unsupported interception backend symbol: $backend (supported: :raster_cpu, :raster_gpu)")
end

function _resolve_interception_backend(backend)
    error(
        "Unsupported interception backend selector type: $(typeof(backend)). " *
        "Use :raster_cpu, :raster_gpu, RasterCPUBackend(), or RasterGPUBackend().",
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
    prepared = _prepare_interception_data(scene, models, options)
    return _compute_first_order(prepared, turtle, fluxes, options)
end

function compute_first_order(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    backend::RasterGPUBackend,
)
    prepared = _prepare_interception_data(scene, models, options)
    data = _rastergpu_scene_data(prepared, backend.config)
    return _compute_first_order_rastergpu(data, turtle, fluxes, options)
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
    backend::InterceptionBackend,
)
    error("Unsupported interception backend type: $(typeof(backend))")
end
