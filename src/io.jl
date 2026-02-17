import PlantGeom
import PlantMeteo
import GeometryBasics
import StaticArrays
import YAML
import Tables
import MultiScaleTreeGraph
import LinearAlgebra: norm, cross
import Dates

function _to_string_dict(x)
    out = Dict{String,Any}()
    for (k, v) in x
        ks = string(k)
        out[ks] = v isa AbstractDict ? _to_string_dict(v) : v
    end
    out
end

function _cfg_get(d::Dict{String,Any}, keys::Vector{String}, default)
    for k in keys
        haskey(d, k) && return d[k]
    end
    return default
end

function _as_bool(x, default::Bool)
    x === nothing && return default
    x isa Bool && return x
    x isa Number && return x != 0
    x isa String && return lowercase(strip(x)) in ("1", "true", "yes", "y", "on")
    return default
end

function _as_int(x, default::Int)
    x === nothing && return default
    x isa Integer && return Int(x)
    x isa Number && return round(Int, x)
    x isa String && return parse(Int, x)
    return default
end

function _as_float(x, default::Float64)
    x === nothing && return default
    x isa Number && return Float64(x)
    x isa String && return parse(Float64, x)
    return default
end

function _join_if_relative(base::AbstractString, p::AbstractString)
    isabspath(p) ? p : normpath(joinpath(base, p))
end

"""
    read_light_config(path)::LightConfig

Read a YAML configuration file and normalize ARCHIMED light options into a `LightConfig`.

`pixel_size` is interpreted like Java input files (centimeters) and converted to meters.
"""
function read_light_config(path::AbstractString)
    raw = YAML.load_file(path)
    d = raw isa AbstractDict ? _to_string_dict(raw) : Dict{String,Any}()
    base = dirname(path)
    d["__base_dir"] = base

    scene_rel = _cfg_get(d, ["scene", "scene_file"], "")
    meteo_rel = _cfg_get(d, ["meteo", "meteo_file"], "")
    scene = _join_if_relative(base, string(scene_rel))
    meteo = _join_if_relative(base, string(meteo_rel))

    props = get(d, "prop", Dict{String,Any}())
    propsd = props isa AbstractDict ? _to_string_dict(props) : Dict{String,Any}()

    all_in_turtle = _as_bool(_cfg_get(propsd, ["all_in_turtle"], get(d, "all_in_turtle", nothing)), false)
    turtle_sectors = _as_int(
        _cfg_get(propsd, ["turtle_nb_sectors", "turtle_sectors", "sky_sectors"], _cfg_get(d, ["turtle_sectors", "sky_sectors"], nothing)),
        46
    )
    # Java config uses centimeters; internal computations use meters.
    pixel_size_cm = _as_float(_cfg_get(propsd, ["pixel_size"], get(d, "pixel_size", nothing)), 0.25)
    pixel_size = pixel_size_cm / 100.0
    area_ratio = _as_bool(_cfg_get(propsd, ["area_ratio"], get(d, "area_ratio", nothing)), true)
    scattering = _as_bool(_cfg_get(propsd, ["scattering"], get(d, "scattering", nothing)), true)
    scattering_max_iter = _as_int(
        _cfg_get(propsd, ["scattering_max_iter", "scat_max_iter"], get(d, "scattering_max_iter", nothing)),
        20
    )
    scattering_stop_ratio = _as_float(
        _cfg_get(propsd, ["scattering_stop_ratio", "scat_stop_ratio"], get(d, "scattering_stop_ratio", nothing)),
        0.01
    )
    scattering_coeff_par = _as_float(_cfg_get(propsd, ["scattering_coeff_par"], get(d, "scattering_coeff_par", nothing)), 0.15)
    scattering_coeff_nir = _as_float(_cfg_get(propsd, ["scattering_coeff_nir"], get(d, "scattering_coeff_nir", nothing)), 0.30)
    cache_radiation = _as_bool(_cfg_get(propsd, ["cache_radiation"], get(d, "cache_radiation", nothing)), false)

    return LightConfig(
        scene,
        meteo,
        all_in_turtle,
        turtle_sectors,
        pixel_size,
        area_ratio,
        scattering,
        scattering_max_iter,
        scattering_stop_ratio,
        scattering_coeff_par,
        scattering_coeff_nir,
        cache_radiation,
        d
    )
end

function _triangle_area3d(p1, p2, p3)
    v1 = StaticArrays.SVector{3,Float64}(p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = StaticArrays.SVector{3,Float64}(p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
    0.5 * norm(cross(v1, v2))
end

function _node_attrs(node)
    try
        getfield(node, :attributes)
    catch
        nothing
    end
end

function _node_children(node)
    try
        getfield(node, :children)
    catch
        nothing
    end
end

function _has_attr(attrs, key::Symbol)
    attrs isa AbstractDict || return false
    haskey(attrs, key) || haskey(attrs, String(key))
end

function _dict_attr(attrs, key::Symbol)
    attrs isa AbstractDict || return nothing
    if haskey(attrs, key)
        return attrs[key]
    elseif haskey(attrs, String(key))
        return attrs[String(key)]
    end
    return nothing
end

function _as_int_or(x, default::Int)
    x === nothing && return default
    x isa Integer && return Int(x)
    x isa Number && return round(Int, x)
    if x isa String
        try
            return parse(Int, strip(x))
        catch
            return default
        end
    end
    return default
end

function _node_item_id(attrs, default::Int)
    for key in (:plantID, :plant_id, :item_id, :itemID)
        v = _dict_attr(attrs, key)
        v === nothing && continue
        return _as_int_or(v, default)
    end
    return default
end

function _node_type_name(node, attrs, default::String="")
    if attrs isa AbstractDict
        for key in (:type, :Type, :functional_type, :functionalType, :organ_type, :organType)
            v = _dict_attr(attrs, key)
            v === nothing && continue
            s = strip(string(v))
            isempty(s) || return s
        end
    end
    s = try
        string(MultiScaleTreeGraph.symbol(node))
    catch
        ""
    end
    isempty(s) ? default : s
end

function _node_component_id(node, default::Int)
    idv =
        try
            getfield(node, :id)
        catch
            nothing
        end
    idv === nothing && return default
    # Java component ids match MTG node ids with +1 offset.
    _as_int_or(idv, default) + 1
end

function _opf_geometry_component_ids(path::AbstractString)
    lines = readlines(path)
    stack = Any[]
    out = Int[]
    seen = Set{Int}()

    for ln in lines
        mo = match(r"<([A-Za-z]+)\b[^>]*\bclass=\"[^\"]+\"[^>]*\bid=\"([0-9]+)\"", ln)
        if mo !== nothing
            push!(stack, (tag=mo.captures[1], id=parse(Int, mo.captures[2]), has_shape=false))
        end

        ms = match(r"<shapeIndex>\s*([0-9]+)\s*</shapeIndex>", ln)
        if ms !== nothing && !isempty(stack)
            ctx = stack[end]
            if !ctx.has_shape
                stack[end] = (tag=ctx.tag, id=ctx.id, has_shape=true)
                if !(ctx.id in seen)
                    push!(out, ctx.id)
                    push!(seen, ctx.id)
                end
            end
        end

        mc = match(r"</([A-Za-z]+)>", ln)
        if mc !== nothing && !isempty(stack)
            if mc.captures[1] == stack[end].tag
                pop!(stack)
            end
        end
    end

    out
end

function _ops_component_id_hints(path::AbstractString)
    hints = Dict{Int,Vector{Int}}()
    base = dirname(path)
    lines = readlines(path)
    for line in lines
        s = strip(replace(line, '\r' => ""))
        isempty(s) && continue
        startswith(s, "#") && continue
        startswith(s, "T") && continue
        toks = split(s)
        length(toks) < 3 && continue

        plant_id = try
            parse(Int, toks[2])
        catch
            continue
        end
        f = toks[3]
        full = isabspath(f) ? f : normpath(joinpath(base, f))
        ext = lowercase(splitext(full)[2])
        if ext == ".opf" && isfile(full)
            hints[plant_id] = _opf_geometry_component_ids(full)
        end
    end
    hints
end

function _collect_mesh_nodes!(
    node,
    meshes,
    node_ids,
    node_group,
    node_type,
    node_item_id,
    node_component_id,
    component_id_hints,
    component_id_hint_cursor,
    next_id::Base.RefValue{Int},
    current_group::String="",
    current_item_id::Int=1,
)
    attrs = _node_attrs(node)
    group = current_group
    item_id = _node_item_id(attrs, current_item_id)
    type_name = _node_type_name(node, attrs, "")
    if attrs isa AbstractDict
        if haskey(attrs, :functional_group)
            group = string(attrs[:functional_group])
        elseif haskey(attrs, "functional_group")
            group = string(attrs["functional_group"])
        end
    end

    if _has_attr(attrs, :geometry) && !_has_attr(attrs, :scene_dimensions)
        push!(meshes, PlantGeom.refmesh_to_mesh(node))
        next_id[] += 1
        push!(node_ids, next_id[])
        node_group[next_id[]] = group
        node_type[next_id[]] = type_name
        node_item_id[next_id[]] = item_id

        hinted = nothing
        if haskey(component_id_hints, item_id)
            cursor = get(component_id_hint_cursor, item_id, 1)
            ids = component_id_hints[item_id]
            if 1 <= cursor <= length(ids)
                hinted = ids[cursor]
                component_id_hint_cursor[item_id] = cursor + 1
            end
        end
        node_component_id[next_id[]] = hinted === nothing ? _node_component_id(node, next_id[] + 1) : Int(hinted)
    end

    children = _node_children(node)
    if children !== nothing
        for ch in children
            _collect_mesh_nodes!(
                ch,
                meshes,
                node_ids,
                node_group,
                node_type,
                node_item_id,
                node_component_id,
                component_id_hints,
                component_id_hint_cursor,
                next_id,
                group,
                item_id,
            )
        end
    end
end

function _build_merged_mesh_with_map_local(mtg; component_id_hints::Dict{Int,Vector{Int}}=Dict{Int,Vector{Int}}())
    meshes = Any[]
    mesh_node_ids = Int[]
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    node_item_id = Dict{Int,Int}()
    node_component_id = Dict{Int,Int}()
    component_id_hint_cursor = Dict{Int,Int}(k => 1 for k in keys(component_id_hints))
    next_id = Ref(0)
    _collect_mesh_nodes!(
        mtg,
        meshes,
        mesh_node_ids,
        node_group,
        node_type,
        node_item_id,
        node_component_id,
        component_id_hints,
        component_id_hint_cursor,
        next_id,
    )

    isempty(meshes) && error("No geometry meshes found in scene.")

    merged_mesh = length(meshes) == 1 ? first(meshes) : PlantGeom.merge_simple_meshes(meshes)

    face2node = Int[]
    for (mi, mesh) in enumerate(meshes)
        nfaces = length(GeometryBasics.decompose(PlantGeom.Face3, mesh))
        append!(face2node, fill(mesh_node_ids[mi], nfaces))
    end
    merged_mesh, face2node, node_group, node_type, node_item_id, node_component_id
end

function _build_scene_geometry(mtg, source_path::AbstractString)
    _build_scene_geometry(mtg, source_path, nothing, Dict{Int,Vector{Int}}())
end

function _build_scene_geometry(mtg, source_path::AbstractString, scene_xy_bounds::Union{Nothing,NTuple{4,Float64}})
    _build_scene_geometry(mtg, source_path, scene_xy_bounds, Dict{Int,Vector{Int}}())
end

function _build_scene_geometry(
    mtg,
    source_path::AbstractString,
    scene_xy_bounds::Union{Nothing,NTuple{4,Float64}},
    component_id_hints::Dict{Int,Vector{Int}},
)
    merged_mesh, face2node, node_group, node_type, node_item_id, node_component_id =
        _build_merged_mesh_with_map_local(mtg; component_id_hints=component_id_hints)

    verts = GeometryBasics.decompose(PlantGeom.Point3, merged_mesh)
    faces = GeometryBasics.decompose(PlantGeom.Face3, merged_mesh)
    node_area = Dict{Int,Float64}()
    bary_acc = Dict{Int,NTuple{3,Float64}}()
    for (i, f) in enumerate(faces)
        n = face2node[i]
        p1 = verts[f[1]]
        p2 = verts[f[2]]
        p3 = verts[f[3]]
        area = _triangle_area3d(p1, p2, p3)
        node_area[n] = get(node_area, n, 0.0) + area
        cx = (p1[1] + p2[1] + p3[1]) / 3.0
        cy = (p1[2] + p2[2] + p3[2]) / 3.0
        cz = (p1[3] + p2[3] + p3[3]) / 3.0
        sx, sy, sz = get(bary_acc, n, (0.0, 0.0, 0.0))
        bary_acc[n] = (sx + area * cx, sy + area * cy, sz + area * cz)
    end
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    for (nid, area) in node_area
        if area > 0
            sx, sy, sz = get(bary_acc, nid, (0.0, 0.0, 0.0))
            barycenter_per_node[nid] = (sx / area, sy / area, sz / area)
        else
            barycenter_per_node[nid] = (NaN, NaN, NaN)
        end
    end
    SceneGeometry(
        mtg,
        merged_mesh,
        face2node,
        node_area,
        barycenter_per_node,
        node_group,
        node_type,
        node_item_id,
        node_component_id,
        String(source_path),
        scene_xy_bounds,
    )
end

function _normalize_ops_lines(lines::Vector{String})
    out = String[]
    for line in lines
        s = replace(line, '\r' => "")
        st = strip(s)
        if isempty(st) || startswith(st, "#")
            push!(out, s)
            continue
        end

        if startswith(st, "T ")
            push!(out, join(split(st), " "))
            continue
        end

        toks = split(st)
        if length(toks) >= 10 && occursin(r"\.(opf|gwa)$"i, toks[3])
            # Java reader interprets OPS plant rows as:
            # sceneId plantId file x y z azimuth inclination stemTwist <ignored>,
            # while PlantGeom expects:
            # sceneId plantId file x y z scale inclinationAzimut inclinationAngle rotation
            # Java ignores OPS scale. For parity with current PlantGeom readers:
            # keep OPF at scale=1.0, and keep historical GWA cm->m conversion (0.01).
            scale = endswith(lowercase(toks[3]), ".gwa") ? "0.01" : "1.0"
            row = [toks[1], toks[2], toks[3], toks[4], toks[5], toks[6], scale, toks[7], toks[8], toks[9]]
            push!(out, join(row, '\t'))
        else
            push!(out, s)
        end
    end
    out
end

function _read_ops_relaxed(path::AbstractString)
    lines = readlines(path)
    normalized = _normalize_ops_lines(lines)

    tmp_path = ""
    mtg = mktemp(dirname(path)) do p, io
        tmp_path = p
        write(io, join(normalized, "\n"))
        write(io, "\n")
        flush(io)
        PlantGeom.read_ops(p)
    end

    isfile(tmp_path) && rm(tmp_path; force=true)
    mtg
end

function _ops_scene_xy_bounds(path::AbstractString)
    lines = readlines(path)
    for line in lines
        s = strip(replace(line, '\r' => ""))
        if isempty(s) || startswith(s, "#")
            continue
        end
        if startswith(s, "T")
            toks = split(s)
            length(toks) < 6 && return nothing
            try
                x0 = parse(Float64, toks[2])
                y0 = parse(Float64, toks[3])
                xsize = parse(Float64, toks[5])
                ysize = parse(Float64, toks[6])
                x1 = x0 + xsize
                y1 = y0 + ysize
                return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
            catch
                return nothing
            end
        end
    end
    nothing
end

"""
    read_scene(path; plantgeom_backend=:auto)::SceneGeometry

Read an ARCHIMED scene file (`.ops`, `.opf`, or `.gwa`) and return merged triangle geometry
with per-node metadata used by light interception and scattering.
"""
function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    scene_xy_bounds = nothing
    component_id_hints = Dict{Int,Vector{Int}}()
    mtg =
        if ext == ".ops"
            scene_xy_bounds = _ops_scene_xy_bounds(path)
            component_id_hints = _ops_component_id_hints(path)
            _read_ops_relaxed(path)
        elseif ext == ".opf"
            component_id_hints[1] = _opf_geometry_component_ids(path)
            PlantGeom.read_opf(path)
        elseif ext == ".gwa"
            PlantGeom.read_gwa(path)
        else
            error("Unsupported scene extension: $ext")
        end
    _build_scene_geometry(mtg, path, scene_xy_bounds, component_id_hints)
end

function _rows_to_namedtuples(table)
    Tables.rowtable(table) |> collect
end

function _meta_to_namedtuple(meta)
    if meta isa NamedTuple
        return meta
    elseif meta isa AbstractDict
        pairs = Pair{Symbol,Any}[Symbol(string(k)) => v for (k, v) in meta]
        return (; pairs...)
    else
        return (;)
    end
end

function _namedtuple_with_meta(row::NamedTuple, meta::NamedTuple)
    pairs = Pair{Symbol,Any}[]
    for k in (:latitude, :longitude, :altitude, :use)
        haskey(meta, k) && push!(pairs, k => getfield(meta, k))
    end
    isempty(pairs) ? row : merge(row, (; pairs...))
end

function _forward_fill_row_date(rows::Vector{<:NamedTuple})
    out = Vector{NamedTuple}(undef, length(rows))
    last_date = nothing
    for i in eachindex(rows)
        r = rows[i]
        if :date in propertynames(r)
            d = getproperty(r, :date)
            if d === missing
                out[i] = last_date === nothing ? r : merge(r, (; date=last_date))
            else
                s = lowercase(strip(string(d)))
                if isempty(s) || s == "missing"
                    out[i] = last_date === nothing ? r : merge(r, (; date=last_date))
                else
                    out[i] = r
                    last_date = d
                end
            end
        else
            out[i] = r
        end
    end
    out
end

"""
    read_meteo(path)::MeteoTable

Read a weather/meteo file with `PlantMeteo.read_weather`, preserving parsed metadata
and returning rows as named tuples.
"""
function read_meteo(path::AbstractString)
    weather, meta =
        try
            w = PlantMeteo.read_weather(path)
            m = try
                PlantMeteo.metadata(w)
            catch
                (; file=path)
            end
            (w, m)
        catch
            try
                w = PlantMeteo.read_weather(path; date_format=Dates.DateFormat("yyyy/mm/dd"))
                m = try
                    PlantMeteo.metadata(w)
                catch
                    (; file=path)
                end
                (w, m)
            catch
                data, metadata_ = PlantMeteo.read_weather_(path)
                (data, metadata_)
            end
        end
    meta_nt = merge((; file=path), _meta_to_namedtuple(meta))
    raw_rows = _forward_fill_row_date(_rows_to_namedtuples(weather))
    rows = [_namedtuple_with_meta(r, meta_nt) for r in raw_rows]
    MeteoTable(rows, meta_nt)
end
