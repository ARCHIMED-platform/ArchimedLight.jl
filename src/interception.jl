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

struct DirectionProjectionResult
    pixel_hits::Dict{Int,Vector{Tuple{Float64,Int}}}
    node_hits::Dict{Int,Int}
    projected_mesh_area::Dict{Int,Float64}
    projected_pixels_area::Dict{Int,Float64}
end

struct InterceptionSceneData
    vertices
    faces::Vector{PlantGeom.Face3}
    face2node::Vector{Int}
    node_ids::Vector{Int}
    plotbox
    node_group::Dict{Int,String}
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

    sort!(stack, by=_hit_height, rev=true, alg=Base.Sort.MergeSort)
    top = stack[1]
    top_nid = _hit_node(top)
    top_nid == spec.node_id || return

    deleteat!(stack, 1)
    if isempty(stack)
        delete!(pixel_hits, idx)
    end

    node_hits[top_nid] = max(0, get(node_hits, top_nid, 0) - 1)
    projected_pixels_area[top_nid] = max(0.0, get(projected_pixels_area, top_nid, 0.0) - plotbox.pixel_area)
end

function _virtual_sensor_groups(models::LightModels)
    out = Set{String}()
    for group_model in values(models)
        group = _normalize_group_name_local(group_model.group)
        isempty(group) && continue
        for type_model in values(group_model.types)
            interception = type_model.interception
            interception === nothing && continue
            if lowercase(strip(interception.model)) == "virtualsensor"
                push!(out, group)
                break
            end
        end
    end
    out
end

function _virtual_sensor_node_ids(node_group::Dict{Int,String}, models::LightModels)
    groups = _virtual_sensor_groups(models)
    isempty(groups) && return Set{Int}()
    out = Set{Int}()
    for (nid, g) in node_group
        _normalize_group_name_local(g) in groups && push!(out, nid)
    end
    out
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
    isempty(_virtual_sensor_groups(models)) || return false
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

    pair_counts = Dict{Tuple{Int,Int},Int}()
    total_from = Dict{Int,Int}()

    for sector in turtle.sectors
        sector.source == :sun && continue
        projection =
            _direction_projection_cached(vertices, faces, face2node, sector.direction, options, plotbox, cache_ctx, upper_hit=false)

        for stack in values(projection.pixel_hits)
            length(stack) <= 1 && continue
            # Java uses a stable sort for hit heights; preserve insertion order on ties.
            sort!(stack, by=x -> x[1], rev=true, alg=Base.Sort.MergeSort)

            for i in eachindex(stack)
                src = stack[i][2]
                src in emitter_nodes || continue

                to = 0
                for j in (i+1):length(stack)
                    nid = stack[j][2]
                    nid in emitter_nodes && continue
                    to = nid
                    break
                end
                to == 0 && continue

                pair_counts[(to, src)] = get(pair_counts, (to, src), 0) + 1
                total_from[src] = get(total_from, src, 0) + 1
            end
        end
    end

    weights = Dict{Tuple{Int,Int},Float64}()
    for ((to, src), c) in pair_counts
        n = get(total_from, src, 0)
        n > 0 || continue
        weights[(to, src)] = c / n
    end
    weights
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
    pixel_hits::Dict{Int,Vector{Tuple{Float64,Int}}},
    node_hits::Dict{Int,Int},
    projected_mesh_area::Dict{Int,Float64},
    projected_pixels_area::Dict{Int,Float64},
    node_id::Int,
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
)
    _project_triangle!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        node_id,
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
    )
end

function _project_triangle!(
    pixel_hits::Dict{Int,Vector{Tuple{Float64,Int}}},
    node_hits::Dict{Int,Int},
    projected_mesh_area::Dict{Int,Float64},
    projected_pixels_area::Dict{Int,Float64},
    node_id::Int,
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
)
    v = (p1, p2, p3)
    projected = Vector{StaticArrays.SVector{3,Float64}}(undef, 3)
    pix_points = Vector{StaticArrays.SVector{3,Float64}}(undef, 3)
    u = unit_scale > 0.0f0 ? unit_scale : 1.0f0
    ox = Float32(origin_x) * u
    oy = Float32(origin_y) * u
    pxs = Float32(x_pix) * u
    pys = Float32(y_pix) * u
    dirx = Float32(direction[1])
    diry = Float32(direction[2])
    dirz = Float32(direction[3])

    @inbounds for k in 1:3
        dirz == 0.0f0 && return
        pz = Float32(v[k][3]) * u
        dz = -pz / dirz
        xw = Float32(v[k][1]) * u + dirx * dz
        yw = Float32(v[k][2]) * u + diry * dz
        projected[k] = StaticArrays.SVector{3,Float64}(Float64(xw / u), Float64(yw / u), 0.0)
        x = Float32((xw - ox) / pxs)
        y = Float32((yw - oy) / pys)
        z = pz
        pix_points[k] = StaticArrays.SVector{3,Float64}(Float64(x), Float64(y), Float64(z))
    end

    iMin = floor(Int, pix_points[1][1])
    iMax = ceil(Int, pix_points[1][1])
    jMin = floor(Int, pix_points[1][2])
    jMax = ceil(Int, pix_points[1][2])
    kMin = floor(Int, pix_points[1][3])
    kMax = ceil(Int, pix_points[1][3])

    @inbounds for k in 2:3
        iMin = min(iMin, floor(Int, pix_points[k][1]))
        iMax = max(iMax, ceil(Int, pix_points[k][1]))
        jMin = min(jMin, floor(Int, pix_points[k][2]))
        jMax = max(jMax, ceil(Int, pix_points[k][2]))
        kMin = min(kMin, floor(Int, pix_points[k][3]))
        kMax = max(kMax, ceil(Int, pix_points[k][3]))
    end

    iLength = (iMax - iMin) + 1
    iLength <= 0 && return

    minY = fill(jMax, iLength)
    maxY = fill(jMin, iLength)

    @inbounds for k in 1:3
        a = pix_points[k]
        b = pix_points[k == 3 ? 1 : k + 1]
        _get_border_pixels(a, b, iMin, minY, maxY)
    end

    normal = strict_java_float ? _compute_normal_f32(pix_points) : _compute_normal(pix_points)
    slopeX_f32, slopeY_f32 =
        if abs(normal[3]) > 1e-5
            (Float32(normal[1] / normal[3]), Float32(normal[2] / normal[3]))
        else
            # Java fallback uses signed direction.z when normal.z is almost zero.
            dz = Float32(direction[3])
            (dz * Float32(normal[1]), dz * Float32(normal[2]))
        end

    z0 = Float32(pix_points[1][3])
    z0 += slopeX_f32 * (Float32(pix_points[1][1]) - Float32(iMin))
    z0 += slopeY_f32 * (Float32(pix_points[1][2]) - Float32(jMin))

    tri_proj_area = _triangle_area_xy(projected[1], projected[2], projected[3])
    buffered_hits = Tuple{Int,Float64}[]
    nb_hits = 0
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
                push!(buffered_hits, (idx, zpix_f32))
            end
        end
    end

    for (idx, zpix_f32) in buffered_hits
        h = get!(pixel_hits, idx) do
            Vector{Tuple{Float64,Int}}()
        end
        if upper_hit
            if isempty(h)
                push!(h, (zpix_f32, node_id))
            elseif zpix_f32 > h[1][1]
                h[1] = (zpix_f32, node_id)
            end
        else
            push!(h, (zpix_f32, node_id))
        end
        node_hits[node_id] = get(node_hits, node_id, 0) + 1
    end

    projected_mesh_area[node_id] = get(projected_mesh_area, node_id, 0.0) + tri_proj_area
    projected_pixels_area[node_id] = get(projected_pixels_area, node_id, 0.0) + nb_hits * pixel_area
end

@inline _hit_height(hit) = hit[1]
@inline _hit_node(hit) = hit[2]

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

    ratios = Dict{Int,Float64}()
    for nid in union(keys(projection.projected_mesh_area), keys(projection.projected_pixels_area))
        if !options.area_ratio
            ratios[nid] = 1.0
        else
            ppa = get(projection.projected_pixels_area, nid, 0.0)
            ratios[nid] = ppa > 0.0 ? get(projection.projected_mesh_area, nid, 0.0) / ppa : 1.0
        end
    end

    visible_area = Dict{Int,Float64}()
    for stack in values(projection.pixel_hits)
        isempty(stack) && continue
        # Java uses a stable sort for hit heights; preserve insertion order on ties.
        sort!(stack, by=_hit_height, rev=true, alg=Base.Sort.MergeSort)

        # VirtualSensor nodes are transparent: they receive without occluding other nodes.
        non_virtual_seen = false
        first_non_virtual = 0
        for hit in stack
            nid = _hit_node(hit)
            if nid in virtual_nodes
                if !non_virtual_seen
                    visible_area[nid] = get(visible_area, nid, 0.0) + plotbox.pixel_area * get(ratios, nid, 1.0)
                end
            else
                first_non_virtual = nid
                non_virtual_seen = true
                break
            end
        end
        if first_non_virtual != 0
            visible_area[first_non_virtual] = get(visible_area, first_non_virtual, 0.0) + plotbox.pixel_area * get(ratios, first_non_virtual, 1.0)
        end
    end

    return visible_area, projection.node_hits
end

function _direction_projection(vertices, faces, face2node, direction, options::LightOptions, plotbox; upper_hit::Union{Nothing,Bool}=nothing)
    toricity = _cfg_toricity(options)
    use_upper_hit = upper_hit === nothing ? false : Bool(upper_hit)
    strict_java_float = _strict_java_float(options)
    unit_scale = Float32(_projection_unit_scale(options))

    pixel_hits = Dict{Int,Vector{Tuple{Float64,Int}}}()
    node_hits = Dict{Int,Int}()
    projected_mesh_area = Dict{Int,Float64}()
    projected_pixels_area = Dict{Int,Float64}()

    @inbounds for fi in eachindex(faces)
        f = faces[fi]
        nid = face2node[fi]
        _project_triangle!(
            pixel_hits,
            node_hits,
            projected_mesh_area,
            projected_pixels_area,
            nid,
            vertices[f[1]],
            vertices[f[2]],
            vertices[f[3]],
            direction,
            plotbox.origin_x,
            plotbox.origin_y,
            plotbox.pix_x,
            plotbox.pix_y,
            plotbox.pixel_area,
            plotbox.nx,
            plotbox.ny,
            toricity,
            use_upper_hit,
            strict_java_float,
            unit_scale,
        )
    end

    _apply_debug_drop_leading_hit!(pixel_hits, node_hits, projected_pixels_area, plotbox, options)
    return DirectionProjectionResult(pixel_hits, node_hits, projected_mesh_area, projected_pixels_area)
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
    node_group = Dict{Int,String}(nid => _scene_group(scene, nid, "") for nid in node_ids)

    return InterceptionSceneData(vertices, faces, face2node, unique(node_ids), plotbox, node_group)
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

    geometry = _scene_geometry_for_interception(scene, models, options)
    virtual_nodes = _virtual_sensor_node_ids(geometry.node_group, models)
    upper_hit = _use_upper_hit_pixel_table(models, options)
    cache_ctx = _projection_cache_context(geometry.vertices, geometry.faces, geometry.face2node, geometry.plotbox, options)

    projected_area_per_node = Dict(id => 0.0 for id in geometry.node_ids)
    incident_power = SpectralNodeValues(
        Dict(id => 0.0 for id in geometry.node_ids),
        Dict(id => 0.0 for id in geometry.node_ids),
    )
    hits_per_node = Dict(id => 0 for id in geometry.node_ids)

    for (k, sector) in enumerate(turtle.sectors)
        visible_area, node_hits =
            _rasterize_direction_java(
                geometry.vertices,
                geometry.faces,
                geometry.face2node,
                sector.direction,
                options,
                geometry.plotbox;
                cache_ctx=cache_ctx,
                virtual_nodes=virtual_nodes,
                upper_hit=upper_hit,
            )

        for (nid, h) in node_hits
            hits_per_node[nid] = get(hits_per_node, nid, 0) + h
        end

        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if par_flux == 0.0 && nir_flux == 0.0
            continue
        end

        for (nid, pa) in visible_area
            pa <= 0.0 && continue
            projected_area_per_node[nid] = get(projected_area_per_node, nid, 0.0) + pa
            incident_power.par[nid] = get(incident_power.par, nid, 0.0) + par_flux * pa
            incident_power.nir[nid] = get(incident_power.nir, nid, 0.0) + nir_flux * pa
        end
    end

    emit_par, emit_nir = _emitter_power_per_node(scene, models)
    emitter_nodes = Set(union(keys(emit_par), keys(emit_nir)))
    if !isempty(emitter_nodes)
        w = _emitter_transfer_weights(
            geometry.vertices,
            geometry.faces,
            geometry.face2node,
            turtle,
            options,
            geometry.plotbox,
            emitter_nodes,
            cache_ctx,
        )
        for ((to, src), ww) in w
            incident_power.par[to] = get(incident_power.par, to, 0.0) + ww * get(emit_par, src, 0.0)
            incident_power.nir[to] = get(incident_power.nir, to, 0.0) + ww * get(emit_nir, src, 0.0)
        end
    end

    FirstOrderResult(
        projected_area_per_node,
        incident_power,
        hits_per_node,
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
