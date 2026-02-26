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

function _dict_attr(attrs, key::Symbol)
    hasproperty(attrs, key) ? getproperty(attrs, key) : nothing
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
    for key in (:type, :Type, :functional_type, :functionalType, :organ_type, :organType)
        v = _dict_attr(attrs, key)
        v === nothing && continue
        s = strip(string(v))
        isempty(s) || return s
    end
    s = string(MultiScaleTreeGraph.symbol(node))
    isempty(s) ? default : s
end

function _node_component_id(node, default::Int)
    idv = MultiScaleTreeGraph.index(node)
    # Java component ids match MTG node ids with +1 offset.
    _as_int_or(idv, default) + 1
end

function _opf_geometry_component_ids(path::AbstractString)
    lines = readlines(path)
    stack = Any[]
    out = Int[]
    seen = Set{Int}()

    for ln in lines
        mo = match(r"<([A-Za-z]+)\b([^>]*)>", ln)
        if mo !== nothing
            attrs = mo.captures[2]
            has_class = match(r"(?i)\bclass\s*=\s*\"[^\"]+\"", attrs) !== nothing
            mid = match(r"(?i)\bid\s*=\s*\"([0-9]+)\"", attrs)
            if has_class && mid !== nothing
                push!(stack, (tag=mo.captures[1], id=parse(Int, mid.captures[1]), has_shape=false))
            end
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

function _normalized_abs_path(base::AbstractString, p::AbstractString)
    base_abs = abspath(base)
    isabspath(p) ? normpath(p) : normpath(abspath(joinpath(base_abs, p)))
end

function _ops_component_id_hints_by_item(path::AbstractString)
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
        full = _normalized_abs_path(base, f)
        ext = lowercase(splitext(full)[2])
        if ext == ".opf" && isfile(full)
            ids = _opf_geometry_component_ids(full)
            if haskey(hints, plant_id)
                append!(hints[plant_id], ids)
            else
                hints[plant_id] = copy(ids)
            end
        end
    end
    hints
end

function _component_hint_path(attrs, scene_base_dir::AbstractString)
    fp = _dict_attr(attrs, :filePath)
    fp === nothing && return nothing
    s = strip(string(fp))
    isempty(s) && return nothing
    _normalized_abs_path(scene_base_dir, s)
end

function _collect_mesh_nodes!(
    node,
    meshes,
    node_ids,
    node_group,
    node_type,
    node_item_id,
    node_component_id,
    component_id_hints_by_item,
    component_id_hints_by_object,
    component_id_hint_cursor_by_item,
    per_item_component_counter,
    next_id::Base.RefValue{Int},
    scene_base_dir::AbstractString,
    current_group::String="",
    current_item_id::Int=1,
    current_object_ids::Union{Nothing,Vector{Int}}=nothing,
    current_object_cursor::Union{Nothing,Base.RefValue{Int}}=nothing,
)
    attrs = MultiScaleTreeGraph.attributes(node)
    group = current_group
    item_id = _node_item_id(attrs, current_item_id)
    type_name = _node_type_name(node, attrs, "")
    object_ids = current_object_ids
    object_cursor = current_object_cursor
    fg = _dict_attr(attrs, :functional_group)
    fg !== nothing && (group = string(fg))

    object_path = _component_hint_path(attrs, scene_base_dir)
    if object_path !== nothing && haskey(component_id_hints_by_object, object_path)
        object_ids = component_id_hints_by_object[object_path]
        object_cursor = Ref(1)
    end

    if hasproperty(attrs, :geometry) && !hasproperty(attrs, :scene_dimensions)
        push!(meshes, PlantGeom.refmesh_to_mesh(node))
        next_id[] += 1
        push!(node_ids, next_id[])
        node_group[next_id[]] = group
        node_type[next_id[]] = type_name
        node_item_id[next_id[]] = item_id

        hinted = nothing
        if object_ids !== nothing && object_cursor !== nothing
            cursor = object_cursor[]
            ids = object_ids
            if 1 <= cursor <= length(ids)
                hinted = ids[cursor]
                object_cursor[] = cursor + 1
            end
        elseif haskey(component_id_hints_by_item, item_id)
            cursor = get(component_id_hint_cursor_by_item, item_id, 1)
            ids = component_id_hints_by_item[item_id]
            if 1 <= cursor <= length(ids)
                hinted = ids[cursor]
                component_id_hint_cursor_by_item[item_id] = cursor + 1
            end
        end
        if hinted !== nothing
            node_component_id[next_id[]] = Int(hinted)
        elseif haskey(component_id_hints_by_item, item_id)
            node_component_id[next_id[]] = _node_component_id(node, next_id[] + 1)
        else
            # Java .gwa-like fallback: component ids are local to each item and start at 2.
            k = get(per_item_component_counter, item_id, 0) + 1
            per_item_component_counter[item_id] = k
            node_component_id[next_id[]] = k + 1
        end
    end

    children = MultiScaleTreeGraph.children(node)
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
                component_id_hints_by_item,
                component_id_hints_by_object,
                component_id_hint_cursor_by_item,
                per_item_component_counter,
                next_id,
                scene_base_dir,
                group,
                item_id,
                object_ids,
                object_cursor,
            )
        end
    end
end

function _build_merged_mesh_with_map_local(
    mtg;
    scene_base_dir::AbstractString="",
    component_id_hints_by_item::Dict{Int,Vector{Int}}=Dict{Int,Vector{Int}}(),
    component_id_hints_by_object::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}(),
)
    meshes = Any[]
    mesh_node_ids = Int[]
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    node_item_id = Dict{Int,Int}()
    node_component_id = Dict{Int,Int}()
    component_id_hint_cursor_by_item = Dict{Int,Int}(k => 1 for k in keys(component_id_hints_by_item))
    per_item_component_counter = Dict{Int,Int}()
    next_id = Ref(0)
    _collect_mesh_nodes!(
        mtg,
        meshes,
        mesh_node_ids,
        node_group,
        node_type,
        node_item_id,
        node_component_id,
        component_id_hints_by_item,
        component_id_hints_by_object,
        component_id_hint_cursor_by_item,
        per_item_component_counter,
        next_id,
        scene_base_dir,
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
    component_id_hints_by_item::Dict{Int,Vector{Int}},
    component_id_hints_by_object::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}(),
)
    merged_mesh, face2node, node_group, node_type, node_item_id, node_component_id =
        _build_merged_mesh_with_map_local(
            mtg;
            scene_base_dir=dirname(source_path),
            component_id_hints_by_item=component_id_hints_by_item,
            component_id_hints_by_object=component_id_hints_by_object,
        )

    verts = GeometryBasics.decompose(GeometryBasics.Point3, merged_mesh)
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

function _triangulate_face_indices(ids::Vector{Int})
    out = NTuple{3,Int}[]
    length(ids) < 3 && return out
    for i in 2:(length(ids)-1)
        push!(out, (ids[1], ids[i], ids[i+1]))
    end
    out
end

function _opf_triangulated_text_if_needed(path::AbstractString)
    txt = read(path, String)
    changed = Ref(false)

    new_txt = replace(txt, r"(?s)<faces>.*?</faces>" => (m -> begin
        block = String(m)
        faces = collect(eachmatch(r"(?s)<face\b[^>]*>(.*?)</face>", block))
        isempty(faces) && return String(m)

        tris = NTuple{3,Int}[]
        for f in faces
            ids = Int[parse(Int, mm.match) for mm in eachmatch(r"-?\d+", f.captures[1])]
            length(ids) < 3 && continue
            length(ids) > 3 && (changed[] = true)
            append!(tris, _triangulate_face_indices(ids))
        end
        isempty(tris) && return String(m)

        io = IOBuffer()
        println(io, "<faces>")
        for (i, t) in enumerate(tris)
            println(io, "\t<face Id=\"$(i - 1)\">")
            println(io, "\t\t$(t[1])\t$(t[2])\t$(t[3])")
            println(io, "\t</face>")
            println(io)
        end
        print(io, "</faces>")
        String(take!(io))
    end))

    if occursin(r"(?s)<materialBDD>\s*</materialBDD>", new_txt)
        changed[] = true
        new_txt = replace(
            new_txt,
            r"(?s)<materialBDD>\s*</materialBDD>" => "<materialBDD>\n\t\t<material Id=\"0\">\n\t\t\t<emission>0 0 0 1</emission>\n\t\t\t<ambient>0.2 0.2 0.2 1</ambient>\n\t\t\t<diffuse>0.8 0.8 0.8 1</diffuse>\n\t\t\t<specular>0 0 0 1</specular>\n\t\t\t<shininess>0</shininess>\n\t\t</material>\n\t</materialBDD>",
        )
    end

    return changed[] ? new_txt : nothing
end

function _rewrite_ops_object_paths(lines::Vector{String}, source_ops::AbstractString, tmp_dir::AbstractString)
    base = abspath(dirname(source_ops))
    rewritten = String[]
    opf_cache = Dict{String,String}()
    hints_by_object = Dict{String,Vector{Int}}()

    for line in lines
        s = strip(replace(line, '\r' => ""))
        if isempty(s) || startswith(s, "#") || startswith(s, "T")
            push!(rewritten, line)
            continue
        end

        toks = split(s)
        if length(toks) < 3
            push!(rewritten, line)
            continue
        end

        rel = toks[3]
        ext = lowercase(splitext(rel)[2])
        if ext != ".opf" && ext != ".gwa"
            push!(rewritten, line)
            continue
        end

        src = _normalized_abs_path(base, rel)
        if ext == ".opf" && isfile(src)
            hints = _opf_geometry_component_ids(src)
            dst = get(opf_cache, src, "")
            if isempty(dst)
                tri = _opf_triangulated_text_if_needed(src)
                if tri === nothing
                    dst = src
                else
                    dst = abspath(joinpath(tmp_dir, "$(splitext(basename(src))[1])-triangulated-$(length(opf_cache) + 1).opf"))
                    write(dst, tri)
                end
                opf_cache[src] = dst
            end
            hints_by_object[normpath(dst)] = hints
            toks[3] = dst
        else
            toks[3] = abspath(src)
        end
        push!(rewritten, join(toks, '\t'))
    end

    return rewritten, hints_by_object
end

function _read_ops_relaxed(path::AbstractString)
    lines = readlines(path)
    normalized = _normalize_ops_lines(lines)
    hints_by_object_ref = Ref(Dict{String,Vector{Int}}())
    mtg = mktempdir(dirname(path)) do tmp_dir
        normalized2, hints_by_object = _rewrite_ops_object_paths(normalized, path, tmp_dir)
        hints_by_object_ref[] = hints_by_object
        tmp_ops = joinpath(tmp_dir, "scene-normalized.ops")
        open(tmp_ops, "w") do io
            write(io, join(normalized2, "\n"))
            write(io, "\n")
        end
        PlantGeom.read_ops(tmp_ops)
    end
    mtg, hints_by_object_ref[]
end

function _read_opf_relaxed(path::AbstractString)
    tri = _opf_triangulated_text_if_needed(path)
    if tri === nothing
        return PlantGeom.read_opf(path)
    end
    mktemp(dirname(path)) do p, io
        write(io, tri)
        flush(io)
        PlantGeom.read_opf(p)
    end
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
    component_id_hints_by_item = Dict{Int,Vector{Int}}()
    component_id_hints_by_object = Dict{String,Vector{Int}}()
    mtg =
        if ext == ".ops"
            scene_xy_bounds = _ops_scene_xy_bounds(path)
            component_id_hints_by_item = _ops_component_id_hints_by_item(path)
            mtg0, hints0 = _read_ops_relaxed(path)
            component_id_hints_by_object = hints0
            mtg0
        elseif ext == ".opf"
            component_id_hints_by_item[1] = _opf_geometry_component_ids(path)
            _read_opf_relaxed(path)
        elseif ext == ".gwa"
            PlantGeom.read_gwa(path)
        else
            error("Unsupported scene extension: $ext")
        end
    _build_scene_geometry(mtg, path, scene_xy_bounds, component_id_hints_by_item, component_id_hints_by_object)
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

function _positive_duration_seconds(v; field_name::AbstractString="step_duration")
    PlantMeteo.positive_duration_seconds(v; field_name=field_name)
end

function _normalize_raw_meteo_dates(data)
    hasproperty(data, :date) || return data
    hour_fmt = Dates.DateFormat("HH:MM:SS")
    date_only = (; date=getproperty(data, :date))
    for date_fmt in (
        Dates.DateFormat("yyyy-mm-ddTHH:MM:SS.s"),
        Dates.DateFormat("yyyy/mm/dd"),
        Dates.DateFormat("yyyy-mm-dd"),
    )
        try
            date = PlantMeteo.compute_date(date_only, date_fmt, hour_fmt; forward_fill_date=true)
            return PlantMeteo.set_column(data, :date, date)
        catch
        end
    end
    return data
end

"""
    read_meteo(path)::MeteoTable

Read a weather/meteo file with `PlantMeteo.read_weather`, preserving parsed metadata
and returning rows as named tuples.
"""
function read_meteo(path::AbstractString)
    weather, meta =
        try
            w = PlantMeteo.read_weather(
                path;
                date_formats=(Dates.DateFormat("yyyy/mm/dd"), Dates.DateFormat("yyyy-mm-dd")),
                forward_fill_date=true,
            )
            m = try
                PlantMeteo.metadata(w)
            catch
                (; file=path)
            end
            (w, m)
        catch
            data, metadata_ = PlantMeteo.read_weather_(path)
            data = _normalize_raw_meteo_dates(data)
            (data, metadata_)
        end
    meta_nt = merge((; file=path), _meta_to_namedtuple(meta))
    raw_rows = _rows_to_namedtuples(weather)
    rows = [_namedtuple_with_meta(r, meta_nt) for r in raw_rows]
    MeteoTable(rows, meta_nt)
end
