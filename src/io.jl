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

@inline _node_attr(node, key::Symbol) = haskey(node, key) ? node[key] : nothing

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

function _node_item_id(node, default::Int)
    for key in (:plantID, :plant_id, :item_id, :itemID)
        v = _node_attr(node, key)
        v === nothing && continue
        return _as_int_or(v, default)
    end
    return default
end

function _node_type_name(node, default::String="")
    for key in (:type, :Type, :functional_type, :functionalType, :organ_type, :organType)
        v = _node_attr(node, key)
        v === nothing && continue
        s = strip(string(v))
        isempty(s) || return s
    end
    s = string(MultiScaleTreeGraph.symbol(node))
    isempty(s) ? default : s
end

const _OPS_RELAXED_KWARGS = (
    relaxed=true,
    assume_scale_column=false,
    opf_scale=1.0,
    gwa_scale=0.01,
)

@inline function _is_scene_geometry_node(node)
    PlantGeom.has_geometry(node) || return false
    !haskey(node, :scene_dimensions)
end

function _source_topology_id(node)
    sid = _node_attr(node, :source_topology_id)
    sid === nothing && return nothing
    cid = _as_int_or(sid, 0)
    cid > 0 ? cid : nothing
end

function _collect_scene_node_metadata!(
    node,
    node_group,
    node_type,
    node_item_id,
    node_component_id,
    per_item_component_counter,
    current_group::String="",
    current_item_id::Int=1,
)
    group = current_group
    item_id = _node_item_id(node, current_item_id)
    type_name = _node_type_name(node, "")
    fg = _node_attr(node, :functional_group)
    fg !== nothing && (group = string(fg))

    if _is_scene_geometry_node(node)
        nid = Int(MultiScaleTreeGraph.node_id(node))
        node_group[nid] = group
        node_type[nid] = type_name
        node_item_id[nid] = item_id

        sid = _source_topology_id(node)
        if sid !== nothing
            node_component_id[nid] = Int(sid)
        else
            # Java .gwa-like fallback: component ids are local to each item and start at 2.
            k = get(per_item_component_counter, item_id, 0) + 1
            per_item_component_counter[item_id] = k
            node_component_id[nid] = k + 1
        end
    end

    children = MultiScaleTreeGraph.children(node)
    if children !== nothing
        for ch in children
            _collect_scene_node_metadata!(
                ch,
                node_group,
                node_type,
                node_item_id,
                node_component_id,
                per_item_component_counter,
                group,
                item_id,
            )
        end
    end
end

function _build_scene_geometry(mtg, source_path::AbstractString)
    _build_scene_geometry(mtg, source_path, nothing)
end

function _build_scene_geometry(mtg, source_path::AbstractString, scene_xy_bounds::Union{Nothing,NTuple{4,Float64}})
    merged_mesh, face2node = PlantGeom.build_merged_mesh_with_map(mtg; filter_fun=_is_scene_geometry_node)

    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    node_item_id = Dict{Int,Int}()
    node_component_id = Dict{Int,Int}()
    per_item_component_counter = Dict{Int,Int}()
    _collect_scene_node_metadata!(
        mtg,
        node_group,
        node_type,
        node_item_id,
        node_component_id,
        per_item_component_counter,
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

function _read_ops_relaxed(path::AbstractString)
    PlantGeom.read_ops(path; _OPS_RELAXED_KWARGS...)
end

function _scene_xy_bounds_from_mtg(mtg)
    haskey(mtg, :scene_dimensions) || return nothing
    dims = mtg[:scene_dimensions]
    dims === nothing && return nothing
    p0 = dims[1]
    p1 = dims[2]
    return (
        min(Float64(p0[1]), Float64(p1[1])),
        min(Float64(p0[2]), Float64(p1[2])),
        max(Float64(p0[1]), Float64(p1[1])),
        max(Float64(p0[2]), Float64(p1[2])),
    )
end

"""
    read_scene(path; plantgeom_backend=:auto)::SceneGeometry

Read an ARCHIMED scene file (`.ops`, `.opf`, or `.gwa`) and return merged triangle geometry
with per-node metadata used by light interception and scattering.
"""
function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    mtg = if ext == ".ops"
        _read_ops_relaxed(path)
    elseif ext == ".opf"
        PlantGeom.read_opf(path)
    elseif ext == ".gwa"
        PlantGeom.read_gwa(path)
    else
        error("Unsupported scene extension: $ext")
    end
    _build_scene_geometry(mtg, path, _scene_xy_bounds_from_mtg(mtg))
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
