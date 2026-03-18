import GeometryBasics
import StaticArrays
import PlantGeom
import MultiScaleTreeGraph
import LinearAlgebra: norm, cross
import Serialization

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

function _cfg_toricity(cfg::LightConfig)
    _as_bool_local(config_value(cfg, "toricity", true), true)
end

function _cfg_debug_drop_leading_hit(cfg::LightConfig)
    spec = config_value(cfg, "debug_drop_leading_hit", nothing)
    spec isa AbstractDict || return nothing
    node_id = try
        Int(get(spec, "node_id", 0))
    catch
        0
    end
    x = try
        Int(get(spec, "x", -1))
    catch
        -1
    end
    y = try
        Int(get(spec, "y", -1))
    catch
        -1
    end
    node_id > 0 || return nothing
    x >= 0 || return nothing
    y >= 0 || return nothing
    return (node_id=node_id, x=x, y=y)
end

function _apply_debug_drop_leading_hit!(pixel_hits, node_hits, projected_pixels_area, plotbox, cfg::LightConfig)
    spec = _cfg_debug_drop_leading_hit(cfg)
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

function _cfg_plot_paving(cfg::LightConfig)
    best = 0
    for spec in model_type_configs(cfg)
        v = try
            Int(get(spec.params, "plot_paving", 0))
        catch
            0
        end
        best = max(best, v)
    end
    return best
end

function _virtual_sensor_groups(cfg::LightConfig)
    out = Set{String}()
    for spec in model_type_configs(cfg)
        group = _normalize_group_name_local(spec.group)
        isempty(group) && continue
        inter = get(spec.params, "Interception", nothing)
        inter isa AbstractDict || continue
        iuse = get(inter, "use", nothing)
        iconf = iuse !== nothing && haskey(inter, string(iuse)) ? inter[string(iuse)] : inter
        iconf isa AbstractDict || continue
        model = lowercase(strip(string(get(iconf, "model", ""))))
        if model == "virtualsensor"
            push!(out, group)
        end
    end
    out
end

function _virtual_sensor_node_ids(node_group::Dict{Int,String}, cfg::LightConfig)
    groups = _virtual_sensor_groups(cfg)
    isempty(groups) && return Set{Int}()
    out = Set{Int}()
    for (nid, g) in node_group
        _normalize_group_name_local(g) in groups && push!(out, nid)
    end
    out
end

function _ignored_group_types(cfg::LightConfig)
    out = Dict{String,Set{String}}()
    for spec in model_type_configs(cfg)
        group = _normalize_group_name_local(spec.group)
        isempty(group) && continue
        inter = get(spec.params, "Interception", nothing)
        inter isa AbstractDict || continue
        iuse = get(inter, "use", nothing)
        iconf = iuse !== nothing && haskey(inter, string(iuse)) ? inter[string(iuse)] : inter
        iconf isa AbstractDict || continue
        model = lowercase(strip(string(get(iconf, "model", ""))))
        model == "ignore" || continue
        tname = strip(string(spec.type))
        isempty(tname) && continue
        s = get!(out, group) do
            Set{String}()
        end
        push!(s, tname)
    end
    out
end

function _is_ignored_node(
    node_id::Int,
    scene::SceneGeometry,
    ignored::Dict{String,Set{String}},
)
    isempty(ignored) && return false
    g = _normalize_group_name_local(get(scene.node_group, node_id, ""))
    t = strip(get(scene.node_type, node_id, ""))
    return haskey(ignored, g) && (t in ignored[g])
end

function _group_light_emitters(cfg::LightConfig)
    out = Dict{Tuple{String,String},NamedTuple{(:par, :nir),Tuple{Float64,Float64}}}()
    for spec in model_type_configs(cfg)
        group = strip(string(spec.group))
        isempty(group) && continue
        type_name = strip(string(spec.type))
        em0 = get(spec.params, "LightEmitter", nothing)
        em0 isa AbstractDict || continue
        em = em0
        model = lowercase(strip(string(get(em, "model", ""))))
        model == "lambertianemitter" || continue

        radiance = try
            Float64(em["radiance"])
        catch
            try
                parse(Float64, string(get(em, "radiance", "0")))
            catch
                0.0
            end
        end
        radiance > 0.0 || continue

        gpar = 0.48
        gnir = 0.52
        g0 = get(em, "gamma", nothing)
        if g0 isa AbstractDict
            gpar = try
                Float64(get(g0, "PAR", gpar))
            catch
                gpar
            end
            gnir = try
                Float64(get(g0, "NIR", gnir))
            catch
                gnir
            end
        end
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
    out
end

function _use_upper_hit_pixel_table(cfg::LightConfig)
    # Default to upper-hit pixel tables unless scattering, virtual sensors,
    # or explicit light emitters require complete interception stacks.
    if get(cfg.general, "scattering", true)
        return false
    end
    isempty(_virtual_sensor_groups(cfg)) || return false
    isempty(_group_light_emitters(cfg)) || return false
    return true
end

function _emitter_power_per_node(scene::SceneGeometry, cfg::LightConfig)
    by_group_type = _group_light_emitters(cfg)
    isempty(by_group_type) && return Dict{Int,Float64}(), Dict{Int,Float64}()

    par = Dict{Int,Float64}()
    nir = Dict{Int,Float64}()
    for ((group, type_name), pwr) in by_group_type
        nids = Int[
            nid for (nid, g) in scene.node_group if g == group && get(scene.node_type, nid, "") == type_name
        ]
        if isempty(nids)
            # Fallback for scenes where type labels are unavailable.
            nids = Int[nid for (nid, g) in scene.node_group if g == group]
        end
        isempty(nids) && continue

        atot = sum(get(scene.total_area_per_node, nid, 0.0) for nid in nids)
        if atot > 0.0
            for nid in nids
                w = get(scene.total_area_per_node, nid, 0.0) / atot
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
    cfg::LightConfig,
    plotbox,
    emitter_nodes::Set{Int},
    cache_ctx,
)
    isempty(emitter_nodes) && return Dict{Tuple{Int,Int},Float64}()

    pair_counts = Dict{Tuple{Int,Int},Int}()
    total_from = Dict{Int,Int}()

    for sector in turtle.sectors
        sector.source == :sun && continue
        pixel_hits, _, _, _ =
            _direction_projection_cached(vertices, faces, face2node, sector.direction, cfg, plotbox, cache_ctx, upper_hit=false)

        for stack in values(pixel_hits)
            length(stack) <= 1 && continue
            # Preserve insertion order on height ties.
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

function _cfg_cache_pixel_table(cfg::LightConfig)
    if haskey(cfg.general, "cache_pixel_table")
        return _as_bool_local(cfg.general["cache_pixel_table"], false)
    end
    false
end

function _projection_cache_dir(cfg::LightConfig)
    base = _config_base_dir(cfg)
    out_rel = cfg.outputs.directory
    isempty(strip(out_rel)) && (out_rel = ".archimedlight_cache")
    out_dir = isabspath(out_rel) ? out_rel : normpath(joinpath(base, out_rel))
    joinpath(out_dir, "pixel_tables_cache")
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

function _projection_scene_key(vertices, faces, face2node, plotbox, cfg::LightConfig)
    h = _FNV64_OFFSET_BASIS
    h = _stable_mix_i64(h, length(vertices))
    h = _stable_mix_i64(h, length(faces))
    h = _stable_mix_i64(h, length(face2node))
    h = _stable_mix_u64(h, _cfg_toricity(cfg) ? UInt64(1) : UInt64(0))
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

function _projection_cache_context(vertices, faces, face2node, plotbox, cfg::LightConfig)
    _cfg_cache_pixel_table(cfg) || return nothing
    return (
        cache_dir=_projection_cache_dir(cfg),
        scene_key=_projection_scene_key(vertices, faces, face2node, plotbox, cfg),
    )
end

function _projection_cache_path(cache_ctx, direction, upper_hit::Bool=false)
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
            return (
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

function _write_projection_cache(path::AbstractString, pixel_hits, node_hits, projected_mesh_area, projected_pixels_area)
    mkpath(dirname(path))
    tmp = path * ".tmp-" * string(getpid()) * "-" * string(time_ns())
    payload = (
        version=1,
        pixel_hits=pixel_hits,
        node_hits=node_hits,
        projected_mesh_area=projected_mesh_area,
        projected_pixels_area=projected_pixels_area,
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
        x0 = p0[1]
        y0 = p0[2]
        x1 = p1[1]
        y1 = p1[2]
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
    end

    xs = [p[1] for p in vertices]
    ys = [p[2] for p in vertices]
    (minimum(xs), minimum(ys), maximum(xs), maximum(ys))
end

function _plotbox(scene::SceneGeometry, vertices, pixel_size::Real)
    pixel_size > zero(pixel_size) || error("pixel_size must be > 0 m")
    pixel_size <= 0.5 || error("pixel_size must be <= 0.5 m (50 cm)")

    x0_raw, y0_raw, x1_raw, y1_raw = _extract_scene_xy_bounds(scene, vertices)
    x0 = min(x0_raw, x1_raw)
    y0 = min(y0_raw, y1_raw)
    x1 = max(x0_raw, x1_raw)
    y1 = max(y0_raw, y1_raw)

    xdim = max(x1 - x0, pixel_size)
    ydim = max(y1 - y0, pixel_size)

    nx = max(1, floor(Int, xdim / pixel_size))
    ny = max(1, floor(Int, ydim / pixel_size))

    pix_x = xdim / nx
    pix_y = ydim / ny
    plot_area = xdim * ydim
    pixel_area = plot_area / (nx * ny)

    (
        origin_x=x0,
        origin_y=y0,
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
    dzden = direction[3]
    iszero(dzden) && return nothing
    pz = point[3]
    dz = -pz / dzden
    x = point[1] + direction[1] * dz
    y = point[2] + direction[2] * dz
    StaticArrays.SVector(x, y, pz)
end

function _triangle_area_xy(p1, p2, p3)
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2])) * 0.5
end

function _compute_normal(points)
    n = length(points)
    if n < 3
        if n == 0
            return StaticArrays.SVector(0.0, 0.0, 0.0)
        end
        z = zero(typeof(points[1][1]))
        return StaticArrays.SVector(z, z, z)
    end
    for k in 1:(n-2)
        p0 = points[k]
        p1 = points[k+1]
        p2 = points[k+2]

        v1x = p1[1] - p0[1]
        v1y = p1[2] - p0[2]
        v1z = p1[3] - p0[3]
        n1 = sqrt((v1x * v1x) + (v1y * v1y) + (v1z * v1z))
        n1 <= 0.0 && continue
        v1x /= n1
        v1y /= n1
        v1z /= n1

        v2x = p2[1] - p1[1]
        v2y = p2[2] - p1[2]
        v2z = p2[3] - p1[3]
        n2 = sqrt((v2x * v2x) + (v2y * v2y) + (v2z * v2z))
        n2 <= 0.0 && continue
        v2x /= n2
        v2y /= n2
        v2z /= n2

        nx = (v1y * v2z) - (v1z * v2y)
        ny = (v1z * v2x) - (v1x * v2z)
        nz = (v1x * v2y) - (v1y * v2x)
        nnorm = sqrt((nx * nx) + (ny * ny) + (nz * nz))
        nnorm <= 0.0 && continue
        return StaticArrays.SVector(nx / nnorm, ny / nnorm, nz / nnorm)
    end
    z = zero(eltype(first(points)))
    StaticArrays.SVector(z, z, z)
end

function _get_border_pixels(p1, p2, i_origin::Int, minY::Vector{Int}, maxY::Vector{Int})
    p_min, p_max = p1[1] < p2[1] ? (p1, p2) : (p2, p1)
    dx = p_max[1] - p_min[1]
    dx < 1e-6 && return

    slope = (p_max[2] - p_min[2]) / dx
    i_min = ceil(Int, p_min[1])
    i_max = floor(Int, p_max[1])

    @inbounds for i in i_min:i_max
        i0 = i - i_origin
        yi = slope * (i - p_min[1]) + p_min[2]
        j = round(Int, yi, RoundNearestTiesUp)
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
    node_hits::Dict{Int,Int},
    projected_mesh_area,
    projected_pixels_area,
    node_id::Int,
    p1,
    p2,
    p3,
    direction,
    origin_x::T,
    origin_y::T,
    x_pix,
    y_pix,
    pixel_area,
    nx::Int,
    ny::Int,
    toricity::Bool,
    upper_hit::Bool,
) where {T}
    v = (p1, p2, p3)
    projected = Vector{typeof(StaticArrays.SVector(origin_x, origin_y, zero(T)))}(undef, 3)
    pix_points = similar(projected)
    ox = origin_x
    oy = origin_y
    pxs = x_pix
    pys = y_pix
    dirx = direction[1]
    diry = direction[2]
    dirz = direction[3]

    @inbounds for k in 1:3
        iszero(dirz) && return
        pz = v[k][3]
        dz = -pz / dirz
        xw = v[k][1] + dirx * dz
        yw = v[k][2] + diry * dz
        projected[k] = StaticArrays.SVector(xw, yw, zero(T))
        x = (xw - ox) / pxs
        y = (yw - oy) / pys
        pix_points[k] = StaticArrays.SVector(x, y, pz)
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

    normal = _compute_normal(pix_points)
    slope_x, slope_y =
        if abs(normal[3]) > 1e-5
            (normal[1] / normal[3], normal[2] / normal[3])
        else
            (direction[3] * normal[1], direction[3] * normal[2])
        end

    z0 = pix_points[1][3]
    z0 += slope_x * (pix_points[1][1] - iMin)
    z0 += slope_y * (pix_points[1][2] - jMin)

    tri_proj_area = _triangle_area_xy(projected[1], projected[2], projected[3])
    buffered_hits = Tuple{Int,T}[]
    nb_hits = 0
    @inbounds for i in iMin:(iMax-1)
        ni = i - iMin
        zi = z0 - slope_x * ni
        ymin_i = minY[ni+1]
        ymax_i = maxY[ni+1]

        for j in ymin_i:(ymax_i-1)
            nj = j - jMin
            zpix = zi - slope_y * nj
            zpix = clamp(zpix, kMin, kMax)
            nb_hits += 1

            ii, jj =
                if toricity
                    (_wrap_index(i, nx), _wrap_index(j, ny))
                else
                    (i, j)
                end

            if toricity || ((0 <= ii < nx) && (0 <= jj < ny))
                idx = ii + 1 + jj * nx
                push!(buffered_hits, (idx, zpix))
            end
        end
    end

    for (idx, zpix) in buffered_hits
        h = get!(pixel_hits, idx) do
            Tuple{T,Int}[]
        end
        if upper_hit
            if isempty(h)
                push!(h, (zpix, node_id))
            elseif zpix > h[1][1]
                h[1] = (zpix, node_id)
            end
        else
            push!(h, (zpix, node_id))
        end
        node_hits[node_id] = get(node_hits, node_id, 0) + 1
    end

    projected_mesh_area[node_id] = get(projected_mesh_area, node_id, zero(T)) + tri_proj_area
    projected_pixels_area[node_id] = get(projected_pixels_area, node_id, zero(T)) + nb_hits * pixel_area
end

@inline _hit_height(hit) = hit[1]
@inline _hit_node(hit) = hit[2]

function _rasterize_direction_java(
    vertices,
    faces,
    face2node,
    direction,
    cfg::LightConfig,
    plotbox;
    cache_ctx=nothing,
    virtual_nodes=Set{Int}(),
    upper_hit::Bool=false,
)
    pixel_hits, node_hits, projected_mesh_area, projected_pixels_area =
        _direction_projection_cached(vertices, faces, face2node, direction, cfg, plotbox, cache_ctx; upper_hit=upper_hit)

    T = isempty(projected_mesh_area) ? typeof(plotbox.pixel_area) : valtype(projected_mesh_area)
    ratios = Dict{Int,T}()
    for nid in union(keys(projected_mesh_area), keys(projected_pixels_area))
        if !get(cfg.general, "area_ratio", true)
            ratios[nid] = one(T)
        else
            ppa = get(projected_pixels_area, nid, zero(T))
            ratios[nid] = ppa > zero(T) ? get(projected_mesh_area, nid, zero(T)) / ppa : one(T)
        end
    end

    visible_area = Dict{Int,T}()
    for stack in values(pixel_hits)
        isempty(stack) && continue
        # Preserve insertion order on height ties.
        sort!(stack, by=_hit_height, rev=true, alg=Base.Sort.MergeSort)

        # VirtualSensor nodes are transparent: they receive without occluding other nodes.
        non_virtual_seen = false
        first_non_virtual = 0
        for hit in stack
            nid = _hit_node(hit)
            if nid in virtual_nodes
                if !non_virtual_seen
                    visible_area[nid] = get(visible_area, nid, zero(T)) + plotbox.pixel_area * get(ratios, nid, one(T))
                end
            else
                first_non_virtual = nid
                non_virtual_seen = true
                break
            end
        end
        if first_non_virtual != 0
            visible_area[first_non_virtual] = get(visible_area, first_non_virtual, zero(T)) + plotbox.pixel_area * get(ratios, first_non_virtual, one(T))
        end
    end

    return visible_area, node_hits
end

function _direction_projection(vertices, faces, face2node, direction, cfg::LightConfig, plotbox; upper_hit::Union{Nothing,Bool}=nothing)
    toricity = _cfg_toricity(cfg)
    use_upper_hit = upper_hit === nothing ? _use_upper_hit_pixel_table(cfg) : Bool(upper_hit)
    T = isempty(vertices) ? typeof(plotbox.pixel_area) : typeof(first(vertices)[1])
    pixel_hits = Dict{Int,Vector{Tuple{T,Int}}}()
    node_hits = Dict{Int,Int}()
    projected_mesh_area = Dict{Int,T}()
    projected_pixels_area = Dict{Int,T}()

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
        )
    end

    _apply_debug_drop_leading_hit!(pixel_hits, node_hits, projected_pixels_area, plotbox, cfg)

    return pixel_hits, node_hits, projected_mesh_area, projected_pixels_area
end

function _direction_projection_cached(vertices, faces, face2node, direction, cfg::LightConfig, plotbox, cache_ctx; upper_hit::Bool=false)
    if cache_ctx === nothing
        return _direction_projection(vertices, faces, face2node, direction, cfg, plotbox; upper_hit=upper_hit)
    end

    path = _projection_cache_path(cache_ctx, direction, upper_hit)
    if isfile(path)
        cached = _read_projection_cache(path)
        cached !== nothing && return cached
        rm(path; force=true)
    end

    pixel_hits, node_hits, projected_mesh_area, projected_pixels_area =
        _direction_projection(vertices, faces, face2node, direction, cfg, plotbox; upper_hit=upper_hit)
    _write_projection_cache(path, pixel_hits, node_hits, projected_mesh_area, projected_pixels_area)
    return pixel_hits, node_hits, projected_mesh_area, projected_pixels_area
end

function _paving_mesh(plotbox, cobble_count::Int, first_node_id::Int)
    T = typeof(plotbox.origin_x)
    x_min = plotbox.origin_x
    y_min = plotbox.origin_y
    x_max = plotbox.origin_x + plotbox.xdim
    y_max = plotbox.origin_y + plotbox.ydim

    plot_x = plotbox.xdim
    plot_y = plotbox.ydim
    plot_area_m2 = plot_x * plot_y
    cobble_area_m2 = plot_area_m2 / max(cobble_count, 1)
    cobble_edge_m = sqrt(cobble_area_m2)

    nx = max(1, floor(Int, plot_x / cobble_edge_m))
    ny = max(1, floor(Int, plot_y / cobble_edge_m))
    cobble_x = plot_x / nx
    cobble_y = plot_y / ny

    min_size = oftype(x_min, 1e-6)
    z = oftype(x_min, 0.005)

    vertices = Vector{typeof(StaticArrays.SVector(x_min, y_min, z))}()
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    node_area = Dict{Int,T}()

    node_id = first_node_id
    x = x_min
    while x < x_max
        x_size = min(cobble_x, x_max - x)
        if x_size > min_size
            y = y_min
            while y < y_max
                y_size = min(cobble_y, y_max - y)
                if y_size > min_size
                    p1 = StaticArrays.SVector(x, y, z)
                    p2 = StaticArrays.SVector(x + x_size, y, z)
                    p3 = StaticArrays.SVector(x + x_size, y + y_size, z)
                    p4 = StaticArrays.SVector(x, y + y_size, z)
                    base = length(vertices)
                    push!(vertices, p1, p2, p3, p4)
                    push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 2, base + 3))
                    push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 3, base + 4))
                    push!(face2node, node_id, node_id)
                    node_area[node_id] = x_size * y_size
                    node_id += 1
                end
                y += cobble_y
            end
        end
        x += cobble_x
    end

    return vertices, faces, face2node, node_area
end

function _scene_geometry_for_interception(scene::SceneGeometry, cfg::LightConfig)
    raw_vertices = GeometryBasics.decompose(GeometryBasics.Point{3,Float64}, scene.merged_mesh)
    vertices = [StaticArrays.SVector(v[1], v[2], v[3]) for v in raw_vertices]
    all_faces = collect(GeometryBasics.decompose(GeometryBasics.TriangleFace{Int}, scene.merged_mesh))
    all_face2node = collect(scene.face2node)

    ignored = _ignored_group_types(cfg)
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    for i in eachindex(all_faces)
        node_id = all_face2node[i]
        _is_ignored_node(node_id, scene, ignored) && continue
        push!(faces, all_faces[i])
        push!(face2node, node_id)
    end
    isempty(face2node) && error("No intercepting geometry left after applying ignore rules.")

    plotbox = _plotbox(scene, vertices, get(cfg.general, "pixel_size", 0.0025))

    node_ids = unique(face2node)
    node_group = Dict{Int,String}(nid => get(scene.node_group, nid, "") for nid in node_ids)
    plot_paving = _cfg_plot_paving(cfg)
    if plot_paving > 0
        first_node = isempty(scene.total_area_per_node) ? 1 : (maximum(keys(scene.total_area_per_node)) + 1)
        paving_vertices, paving_faces, paving_face2node, _ = _paving_mesh(plotbox, plot_paving, first_node)
        v_offset = length(vertices)
        append!(vertices, paving_vertices)
        for f in paving_faces
            push!(faces, GeometryBasics.TriangleFace{Int}(f[1] + v_offset, f[2] + v_offset, f[3] + v_offset))
        end
        append!(face2node, paving_face2node)
        append!(node_ids, unique(paving_face2node))
        for nid in unique(paving_face2node)
            node_group[nid] = "pavement"
        end
    end

    return vertices, faces, face2node, unique(node_ids), plotbox, node_group
end

"""
    compute_first_order(scene, turtle, fluxes, cfg; backend=:raster_cpu)::FirstOrderResult

Compute first-order interception by rasterizing each direction, then integrating projected area,
incident power, and hit counts per geometry node.

`backend` accepts either a symbol (currently `:raster_cpu`) or an
`InterceptionBackend` instance (currently `RasterCPUBackend()`).
"""
function compute_first_order(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    cfg::LightConfig;
    backend=:raster_cpu,
)
    return compute_first_order(scene, turtle, fluxes, cfg, _resolve_interception_backend(backend))
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
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    cfg::LightConfig,
    ::RasterCPUBackend,
)

    vertices, faces, face2node, node_ids, plotbox, node_group = _scene_geometry_for_interception(scene, cfg)
    virtual_nodes = _virtual_sensor_node_ids(node_group, cfg)
    upper_hit = _use_upper_hit_pixel_table(cfg)
    cache_ctx = _projection_cache_context(vertices, faces, face2node, plotbox, cfg)

    projected_area_per_node = Dict(id => 0.0 for id in node_ids)
    incident_par_power_per_node = Dict(id => 0.0 for id in node_ids)
    incident_nir_power_per_node = Dict(id => 0.0 for id in node_ids)
    hits_per_node = Dict(id => 0 for id in node_ids)

    for (k, sector) in enumerate(turtle.sectors)
        visible_area, node_hits =
            _rasterize_direction_java(
                vertices,
                faces,
                face2node,
                sector.direction,
                cfg,
                plotbox;
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
            incident_par_power_per_node[nid] = get(incident_par_power_per_node, nid, 0.0) + par_flux * pa
            incident_nir_power_per_node[nid] = get(incident_nir_power_per_node, nid, 0.0) + nir_flux * pa
        end
    end

    emit_par, emit_nir = _emitter_power_per_node(scene, cfg)
    emitter_nodes = Set(union(keys(emit_par), keys(emit_nir)))
    if !isempty(emitter_nodes)
        w = _emitter_transfer_weights(vertices, faces, face2node, turtle, cfg, plotbox, emitter_nodes, cache_ctx)
        for ((to, src), ww) in w
            incident_par_power_per_node[to] = get(incident_par_power_per_node, to, 0.0) + ww * get(emit_par, src, 0.0)
            incident_nir_power_per_node[to] = get(incident_nir_power_per_node, to, 0.0) + ww * get(emit_nir, src, 0.0)
        end
    end

    FirstOrderResult(
        projected_area_per_node,
        incident_par_power_per_node,
        incident_nir_power_per_node,
        hits_per_node,
    )
end

function compute_first_order(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    cfg::LightConfig,
    backend::InterceptionBackend,
)
    error("Unsupported interception backend type: $(typeof(backend))")
end
