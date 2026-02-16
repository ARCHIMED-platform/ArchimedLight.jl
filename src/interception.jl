import GeometryBasics
import StaticArrays
import PlantGeom
import LinearAlgebra: norm, cross
import Serialization

function _as_bool_local(x, default::Bool)
    x === nothing && return default
    x isa Bool && return x
    x isa Number && return x != 0
    x isa String && return lowercase(strip(x)) in ("1", "true", "yes", "y", "on")
    return default
end

function _cfg_toricity(cfg::LightConfig)
    raw = cfg.raw
    if haskey(raw, "toricity")
        return _as_bool_local(raw["toricity"], true)
    end
    props = get(raw, "prop", nothing)
    if props isa AbstractDict && haskey(props, "toricity")
        return _as_bool_local(props["toricity"], true)
    end
    return true
end

function _cfg_plot_paving(cfg::LightConfig)
    models = get(cfg.raw, "models", nothing)
    models isa AbstractVector || return 0
    base = get(cfg.raw, "__base_dir", dirname(cfg.scene))
    best = 0
    for m in models
        mp = String(m)
        path = isabspath(mp) ? mp : normpath(joinpath(base, mp))
        isfile(path) || continue
        txt = read(path, String)
        for mm in eachmatch(r"plot_paving\s*:\s*([0-9]+)", txt)
            v = try
                parse(Int, mm.captures[1])
            catch
                0
            end
            best = max(best, v)
        end
    end
    best
end

function _cfg_cache_pixel_table(cfg::LightConfig)
    raw = cfg.raw
    if haskey(raw, "cache_pixel_table")
        return _as_bool_local(raw["cache_pixel_table"], false)
    end
    if haskey(raw, "save_on_disk")
        return _as_bool_local(raw["save_on_disk"], false)
    end
    props = get(raw, "prop", nothing)
    if props isa AbstractDict
        haskey(props, "cache_pixel_table") && return _as_bool_local(props["cache_pixel_table"], false)
        haskey(props, "save_on_disk") && return _as_bool_local(props["save_on_disk"], false)
    end
    return false
end

function _projection_cache_dir(cfg::LightConfig)
    base = get(cfg.raw, "__base_dir", dirname(cfg.scene))
    out_rel = haskey(cfg.raw, "output_directory") ? string(cfg.raw["output_directory"]) : ".archimedlight_cache"
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

function _projection_cache_path(cache_ctx, direction)
    scene_hex = string(cache_ctx.scene_key, base=16, pad=16)
    dir_hex = string(_projection_dir_key(direction), base=16, pad=16)
    joinpath(cache_ctx.cache_dir, "proj_" * scene_hex * "_" * dir_hex * ".jls")
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

    attrs =
        if hasproperty(scene.mtg, :attributes)
            getfield(scene.mtg, :attributes)
        else
            nothing
        end

    if attrs isa AbstractDict
        dims =
            if haskey(attrs, :scene_dimensions)
                attrs[:scene_dimensions]
            elseif haskey(attrs, "scene_dimensions")
                attrs["scene_dimensions"]
            else
                nothing
            end

        if dims !== nothing
            try
                p0 = dims[1]
                p1 = dims[2]
                x0 = Float64(p0[1])
                y0 = Float64(p0[2])
                x1 = Float64(p1[1])
                y1 = Float64(p1[2])
                return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
            catch
            end
        end
    end

    xs = Float64[p[1] for p in vertices]
    ys = Float64[p[2] for p in vertices]
    (minimum(xs), minimum(ys), maximum(xs), maximum(ys))
end

function _plotbox(scene::SceneGeometry, vertices, pixel_size::Float64)
    pixel_size > 0.0 || error("pixel_size must be > 0 m")
    # Java fixtures enforce a strict upper bound at 50 cm.
    pixel_size <= 0.5 || error("pixel_size must be <= 0.5 m (50 cm)")

    x0, y0, x1, y1 = _extract_scene_xy_bounds(scene, vertices)
    xdim = max(x1 - x0, pixel_size)
    ydim = max(y1 - y0, pixel_size)

    nx = max(1, floor(Int, xdim / pixel_size))
    ny = max(1, floor(Int, ydim / pixel_size))

    pix_x = xdim / nx
    pix_y = ydim / ny
    pixel_area = (xdim * ydim) / (nx * ny)

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
    abs(direction[3]) <= eps(Float64) && return nothing
    dz = -point[3] / direction[3]
    x = point[1] + direction[1] * dz
    y = point[2] + direction[2] * dz
    # Keep source altitude for stack sorting as in Java.
    StaticArrays.SVector{3,Float64}(x, y, point[3])
end

function _triangle_area_xy(p1, p2, p3)
    abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2])) * 0.5
end

function _compute_normal(points)
    n = length(points)
    n < 3 && return StaticArrays.SVector{3,Float64}(0.0, 0.0, 0.0)
    for k in 1:(n - 2)
        v1 = points[k + 1] - points[k]
        n1 = norm(v1)
        n1 <= 0.0 && continue
        v1n = v1 / n1

        v2 = points[k + 2] - points[k + 1]
        n2 = norm(v2)
        n2 <= 0.0 && continue
        v2n = v2 / n2

        nn = cross(v1n, v2n)
        nnorm = norm(nn)
        nnorm <= 0.0 && continue
        return nn / nnorm
    end
    StaticArrays.SVector{3,Float64}(0.0, 0.0, 0.0)
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
        j = trunc(Int, floor(yi) + 0.5)
        if 0 <= i0 < length(minY)
            idx = i0 + 1
            minY[idx] = min(minY[idx], j)
            maxY[idx] = max(maxY[idx], j)
        end
    end
end

function _wrap_index(i::Int, n::Int)
    mod(i, n)
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
)
    v = (p1, p2, p3)
    projected = Vector{StaticArrays.SVector{3,Float64}}(undef, 3)
    pix_points = Vector{StaticArrays.SVector{3,Float64}}(undef, 3)

    @inbounds for k in 1:3
        pp = _project_point_ground(v[k], direction)
        pp === nothing && return
        projected[k] = StaticArrays.SVector{3,Float64}(pp[1], pp[2], 0.0)
        x = (pp[1] - origin_x) / x_pix
        y = (pp[2] - origin_y) / y_pix
        pix_points[k] = StaticArrays.SVector{3,Float64}(x, y, pp[3])
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
    slopeX, slopeY =
        if abs(normal[3]) > 1e-5
            (normal[1] / normal[3], normal[2] / normal[3])
        else
            (direction[3] * normal[1], direction[3] * normal[2])
        end

    z0 = pix_points[1][3] + slopeX * (pix_points[1][1] - iMin) + slopeY * (pix_points[1][2] - jMin)

    nb_hits = 0
    @inbounds for i in iMin:(iMax - 1)
        ni = i - iMin
        zi = z0 - slopeX * ni
        ymin_i = minY[ni + 1]
        ymax_i = maxY[ni + 1]

        for j in ymin_i:(ymax_i - 1)
            nj = j - jMin
            zpix = zi - slopeY * nj
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
                h = get!(pixel_hits, idx) do
                    Vector{Tuple{Float64,Int}}()
                end
                push!(h, (zpix, node_id))
                node_hits[node_id] = get(node_hits, node_id, 0) + 1
            end
        end
    end

    projected_mesh_area[node_id] = get(projected_mesh_area, node_id, 0.0) + _triangle_area_xy(projected[1], projected[2], projected[3])
    projected_pixels_area[node_id] = get(projected_pixels_area, node_id, 0.0) + nb_hits * pixel_area
end

function _rasterize_direction_java(vertices, faces, face2node, direction, cfg::LightConfig, plotbox; cache_ctx=nothing)
    pixel_hits, node_hits, projected_mesh_area, projected_pixels_area =
        _direction_projection_cached(vertices, faces, face2node, direction, cfg, plotbox, cache_ctx)

    ratios = Dict{Int,Float64}()
    for nid in union(keys(projected_mesh_area), keys(projected_pixels_area))
        if !cfg.area_ratio
            ratios[nid] = 1.0
        else
            ppa = get(projected_pixels_area, nid, 0.0)
            ratios[nid] = ppa > 0.0 ? get(projected_mesh_area, nid, 0.0) / ppa : 1.0
        end
    end

    visible_area = Dict{Int,Float64}()
    for stack in values(pixel_hits)
        isempty(stack) && continue
        sort!(stack, by=x -> x[1], rev=true)

        # Current Julia light-only scope keeps nodes opaque.
        nid = stack[1][2]
        visible_area[nid] = get(visible_area, nid, 0.0) + plotbox.pixel_area * get(ratios, nid, 1.0)
    end

    return visible_area, node_hits
end

function _direction_projection(vertices, faces, face2node, direction, cfg::LightConfig, plotbox)
    toricity = _cfg_toricity(cfg)

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
        )
    end

    return pixel_hits, node_hits, projected_mesh_area, projected_pixels_area
end

function _direction_projection_cached(vertices, faces, face2node, direction, cfg::LightConfig, plotbox, cache_ctx)
    if cache_ctx === nothing
        return _direction_projection(vertices, faces, face2node, direction, cfg, plotbox)
    end

    path = _projection_cache_path(cache_ctx, direction)
    if isfile(path)
        cached = _read_projection_cache(path)
        cached !== nothing && return cached
        rm(path; force=true)
    end

    pixel_hits, node_hits, projected_mesh_area, projected_pixels_area =
        _direction_projection(vertices, faces, face2node, direction, cfg, plotbox)
    _write_projection_cache(path, pixel_hits, node_hits, projected_mesh_area, projected_pixels_area)
    return pixel_hits, node_hits, projected_mesh_area, projected_pixels_area
end

function _paving_mesh(plotbox, cobble_count::Int, first_node_id::Int)
    x_min = plotbox.origin_x
    y_min = plotbox.origin_y
    x_max = plotbox.origin_x + plotbox.xdim
    y_max = plotbox.origin_y + plotbox.ydim

    plot_area = plotbox.xdim * plotbox.ydim
    cobble_area = plot_area / max(cobble_count, 1)
    cobble_edge = sqrt(cobble_area)

    nx = max(1, floor(Int, plotbox.xdim / cobble_edge))
    ny = max(1, floor(Int, plotbox.ydim / cobble_edge))
    cobble_x = plotbox.xdim / nx
    cobble_y = plotbox.ydim / ny

    min_size = 1e-4
    z = 0.005

    vertices = StaticArrays.SVector{3,Float64}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    node_area = Dict{Int,Float64}()

    node_id = first_node_id
    x = x_min
    while x < x_max
        x_size = (x + cobble_x > x_max) ? (x_max - x) : cobble_x
        if x_size > min_size
            y = y_min
            while y < y_max
                y_size = (y + cobble_y > y_max) ? (y_max - y) : cobble_y
                if y_size > min_size
                    p1 = StaticArrays.SVector{3,Float64}(x, y, z)
                    p2 = StaticArrays.SVector{3,Float64}(x + x_size, y, z)
                    p3 = StaticArrays.SVector{3,Float64}(x + x_size, y + y_size, z)
                    p4 = StaticArrays.SVector{3,Float64}(x, y + y_size, z)

                    base = length(vertices)
                    push!(vertices, p1, p2, p3, p4)
                    push!(faces, PlantGeom.Face3(base + 1, base + 2, base + 3))
                    push!(faces, PlantGeom.Face3(base + 1, base + 3, base + 4))
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
    raw_vertices = GeometryBasics.decompose(PlantGeom.Point3, scene.merged_mesh)
    vertices = [StaticArrays.SVector{3,Float64}(v[1], v[2], v[3]) for v in raw_vertices]
    faces = collect(GeometryBasics.decompose(PlantGeom.Face3, scene.merged_mesh))
    face2node = collect(scene.face2node)
    plotbox = _plotbox(scene, vertices, cfg.pixel_size)

    node_ids = collect(keys(scene.total_area_per_node))
    node_group = Dict{Int,String}(k => v for (k, v) in scene.node_group)
    plot_paving = _cfg_plot_paving(cfg)
    if plot_paving > 0
        first_node = isempty(node_ids) ? 1 : (maximum(node_ids) + 1)
        paving_vertices, paving_faces, paving_face2node, _ = _paving_mesh(plotbox, plot_paving, first_node)
        v_offset = length(vertices)
        append!(vertices, paving_vertices)
        for f in paving_faces
            push!(faces, PlantGeom.Face3(f[1] + v_offset, f[2] + v_offset, f[3] + v_offset))
        end
        append!(face2node, paving_face2node)
        append!(node_ids, unique(paving_face2node))
        for nid in unique(paving_face2node)
            node_group[nid] = "pavement"
        end
    end

    return vertices, faces, face2node, unique(node_ids), plotbox, node_group
end

function compute_first_order(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    cfg::LightConfig;
    backend=:raster_cpu,
)
    backend == :raster_cpu || error("Unsupported backend: $backend")

    vertices, faces, face2node, node_ids, plotbox, _ = _scene_geometry_for_interception(scene, cfg)
    cache_ctx = _projection_cache_context(vertices, faces, face2node, plotbox, cfg)

    projected_area_per_node = Dict(id => 0.0 for id in node_ids)
    incident_par_power_per_node = Dict(id => 0.0 for id in node_ids)
    incident_nir_power_per_node = Dict(id => 0.0 for id in node_ids)
    hits_per_node = Dict(id => 0 for id in node_ids)

    for (k, sector) in enumerate(turtle.sectors)
        visible_area, node_hits =
            _rasterize_direction_java(vertices, faces, face2node, sector.direction, cfg, plotbox; cache_ctx=cache_ctx)

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

    FirstOrderResult(
        projected_area_per_node,
        incident_par_power_per_node,
        incident_nir_power_per_node,
        hits_per_node,
    )
end
