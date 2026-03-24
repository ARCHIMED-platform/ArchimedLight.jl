import GeometryBasics
import StaticArrays
import PlantGeom
import MultiScaleTreeGraph
import LinearAlgebra: norm, cross
import Serialization

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
            stack[j + 1] = stack[j]
            j -= 1
        end
        stack[j + 1] = x
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
        while j >= 1 && heights[start + j - 1] < xh
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
    return @inbounds Float64(stack.parent.heights[_flat_stack_start(stack) + i - 1])
end
@inline function _stack_hit_node(stack::FlatPixelHitStack, i::Int)
    @boundscheck checkbounds(stack, i)
    return @inbounds Int(stack.parent.nodes[_flat_stack_start(stack) + i - 1])
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
        for k in (j + 1):length(stack)
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
    @inbounds for j in 0:(n_hits - 1)
        src_idx = Int(nodes[start + j])
        emitter_node_mask[src_idx] || continue

        to_idx = 0
        for k in (j + 1):(n_hits - 1)
            node_idx = Int(nodes[start + k])
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
    stack,
    options::LightOptions,
    pixel_area::Float64,
    virtual_node_mask::Vector{Bool},
)
    non_virtual_seen = false
    first_non_virtual_idx = 0
    @inbounds for hit in stack
        node_idx = _hit_node(hit)
        if virtual_node_mask[node_idx]
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
    return nothing
end

@inline function _accumulate_visible_area_dense!(
    visible_area::Vector{Float64},
    projection::DenseDirectionProjectionResult,
    stack::FlatPixelHitStack,
    options::LightOptions,
    pixel_area::Float64,
    virtual_node_mask::Vector{Bool},
)
    start = _flat_stack_start(stack)
    nodes = stack.parent.nodes
    n_hits = length(stack)
    first_non_virtual_idx = 0
    @inbounds for h in 0:(n_hits - 1)
        node_idx = Int(nodes[start + h])
        if virtual_node_mask[node_idx]
            ratio = _projection_area_ratio(projection, options, node_idx)
            visible_area[node_idx] += pixel_area * ratio
        else
            first_non_virtual_idx = node_idx
            break
        end
    end
    if first_non_virtual_idx != 0
        ratio = _projection_area_ratio(projection, options, first_non_virtual_idx)
        visible_area[first_non_virtual_idx] += pixel_area * ratio
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
    @inbounds for pos in i:(count - 1)
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

_all_dense_float_node_map(node_ids::Vector{Int}, values::Vector{Float64}) =
    Dict{Int,Float64}(node_ids[i] => values[i] for i in eachindex(node_ids))

_all_dense_int_node_map(node_ids::Vector{Int}, values::Vector{Int}) =
    Dict{Int,Int}(node_ids[i] => values[i] for i in eachindex(node_ids))

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

function _is_ignored_node(node_id::Int, scene::SceneGeometry, ignored::Dict{String,Set{String}})
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

function _validate_scene_models(scene::SceneGeometry, face2node::Vector{Int}, models::LightModels, ignored::Dict{String,Set{String}})
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
    # Java defaults to upper-hit pixel tables unless scattering, virtual sensors,
    # or explicit light emitters require complete interception stacks.
    if options.scattering
        return false
    end
    _has_sensor_models(models) && return false
    isempty(_group_light_emitters(models)) || return false
    return true
end

function _emitter_power_per_node(scene::SceneGeometry, models::LightModels)
    by_group_type = _group_light_emitters(models)
    isempty(by_group_type) && return Dict{Int,Float64}(), Dict{Int,Float64}()

    par = Dict{Int,Float64}()
    nir = Dict{Int,Float64}()
    for ((group, type_name), pwr) in by_group_type
        nids = Int[
            nid for (nid, node) in scene.nodes if node.group == group && node.type == type_name
        ]
        if isempty(nids)
            # Fallback for scenes where type labels are unavailable.
            nids = Int[nid for (nid, node) in scene.nodes if node.group == group]
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
            for k in (j + 1):length(stack)
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

function _extract_scene_xy_bounds(scene::SceneGeometry, vertices)
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

function _plotbox(scene::SceneGeometry, vertices, pixel_size::Float64)
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
            node_idx = geometry.node_index[nid]
            if virtual_node_mask[node_idx]
                if !non_virtual_seen
                    ratio = _projection_area_ratio(projection, options, nid)
                    visible_area[node_idx] += pixel_area * ratio
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
    geometry::InterceptionSceneData,
    stacks_sorted::Bool=false,
)
    visible_area = zeros(Float64, length(geometry.node_ids))
    pixel_area = plotbox.pixel_area
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _accumulate_visible_area_dense!(visible_area, projection, stack, options, pixel_area, virtual_node_mask)
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

function _scene_geometry_for_interception(scene::SceneGeometry, models::LightModels, options::LightOptions)
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
    scene::SceneGeometry,
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
    scene::SceneGeometry,
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

function _interception_output_keys(scene::SceneGeometry, models::LightModels, options::LightOptions)
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

`backend` accepts either a symbol (currently `:raster_cpu`) or an
`InterceptionBackend` instance (currently `RasterCPUBackend()`).
"""
function compute_first_order(
    scene::SceneGeometry,
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
    error("Unsupported interception backend symbol: $backend (supported: :raster_cpu)")
end

function _resolve_interception_backend(backend)
    error(
        "Unsupported interception backend selector type: $(typeof(backend)). " *
        "Use :raster_cpu or RasterCPUBackend().",
    )
end

function compute_first_order(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    ::RasterCPUBackend,
)
    prepared = _prepare_interception_data(scene, models, options)
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

    return FirstOrderResult(
        _all_dense_float_node_map(geometry.node_ids, projected_area_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(geometry.node_ids, incident_power_par),
            _all_dense_float_node_map(geometry.node_ids, incident_power_nir),
        ),
        _all_dense_int_node_map(geometry.node_ids, hits_per_node),
    )
end

function compute_first_order(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
    backend::InterceptionBackend,
)
    error("Unsupported interception backend type: $(typeof(backend))")
end
