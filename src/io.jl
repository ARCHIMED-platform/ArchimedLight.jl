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

function _collect_mesh_nodes!(node, meshes, node_ids, node_group, next_id::Base.RefValue{Int}, current_group::String="")
    attrs = _node_attrs(node)
    group = current_group
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
    end

    children = _node_children(node)
    if children !== nothing
        for ch in children
            _collect_mesh_nodes!(ch, meshes, node_ids, node_group, next_id, group)
        end
    end
end

function _build_merged_mesh_with_map_local(mtg)
    meshes = Any[]
    mesh_node_ids = Int[]
    node_group = Dict{Int,String}()
    next_id = Ref(0)
    _collect_mesh_nodes!(mtg, meshes, mesh_node_ids, node_group, next_id)

    isempty(meshes) && error("No geometry meshes found in scene.")

    merged_mesh = length(meshes) == 1 ? first(meshes) : PlantGeom.merge_simple_meshes(meshes)

    face2node = Int[]
    for (mi, mesh) in enumerate(meshes)
        nfaces = length(GeometryBasics.decompose(PlantGeom.Face3, mesh))
        append!(face2node, fill(mesh_node_ids[mi], nfaces))
    end
    merged_mesh, face2node, node_group
end

function _build_scene_geometry(mtg, source_path::AbstractString)
    _build_scene_geometry(mtg, source_path, nothing)
end

function _build_scene_geometry(mtg, source_path::AbstractString, scene_xy_bounds::Union{Nothing,NTuple{4,Float64}})
    merged_mesh, face2node, node_group = _build_merged_mesh_with_map_local(mtg)

    verts = GeometryBasics.decompose(PlantGeom.Point3, merged_mesh)
    faces = GeometryBasics.decompose(PlantGeom.Face3, merged_mesh)
    node_area = Dict{Int,Float64}()
    for (i, f) in enumerate(faces)
        n = face2node[i]
        area = _triangle_area3d(verts[f[1]], verts[f[2]], verts[f[3]])
        node_area[n] = get(node_area, n, 0.0) + area
    end
    SceneGeometry(mtg, merged_mesh, face2node, node_area, node_group, String(source_path), scene_xy_bounds)
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

function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    scene_xy_bounds = nothing
    mtg =
        if ext == ".ops"
            scene_xy_bounds = _ops_scene_xy_bounds(path)
            _read_ops_relaxed(path)
        elseif ext == ".opf"
            PlantGeom.read_opf(path)
        elseif ext == ".gwa"
            PlantGeom.read_gwa(path)
        else
            error("Unsupported scene extension: $ext")
        end
    _build_scene_geometry(mtg, path, scene_xy_bounds)
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
    for k in (:latitude, :longitude, :altitude)
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
