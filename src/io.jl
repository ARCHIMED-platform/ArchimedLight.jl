import PlantGeom
import PlantMeteo
import GeometryBasics
import StaticArrays
import YAML
import Tables
import MultiScaleTreeGraph
import LinearAlgebra: norm, cross
import Dates
import OrderedCollections: OrderedDict

const _NEXT_SCENE_NODE_ID = Ref(1)

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

function _load_yaml_ordered(path::AbstractString)
    YAML.load_file(path; dicttype=OrderedDict{String,Any})
end

function _model_paths_from_raw(d::AbstractDict{String,Any}, scene::AbstractString)
    models = get(d, "models", nothing)
    models isa AbstractVector || return String[]
    base = get(d, "__base_dir", dirname(scene))
    [_join_if_relative(base, string(m)) for m in models]
end

function refresh_light_config!(cfg::LightConfig; reload_models::Bool=false)
    d = cfg.raw
    base = get(d, "__base_dir", dirname(cfg.source_path))
    d["__base_dir"] = base

    haskey(d, "scene") || error("Missing config key `scene`. Expected a top-level key in the YAML config.")
    haskey(d, "meteo") || error("Missing config key `meteo`. Expected a top-level key in the YAML config.")

    cfg.scene = _join_if_relative(base, string(d["scene"]))
    cfg.meteo = _join_if_relative(base, string(d["meteo"]))
    cfg.all_in_turtle = _as_bool(get(d, "all_in_turtle", false), false)
    cfg.turtle_sectors = _as_int(get(d, "sky_sectors", 46), 46)
    cfg.pixel_size = _as_float(get(d, "pixel_size", 0.25), 0.25) / 100.0
    cfg.area_ratio = _as_bool(get(d, "area_ratio", true), true)
    cfg.scattering = _as_bool(get(d, "scattering", true), true)
    cfg.scattering_max_iter = _as_int(get(d, "scattering_max_iter", 20), 20)
    cfg.scattering_stop_ratio = _as_float(get(d, "scattering_stop_ratio", 0.01), 0.01)
    cfg.scattering_coeff_par = _as_float(get(d, "scattering_coeff_par", 0.15), 0.15)
    cfg.scattering_coeff_nir = _as_float(get(d, "scattering_coeff_nir", 0.30), 0.30)
    cfg.cache_radiation = _as_bool(get(d, "cache_radiation", false), false)

    model_paths = _model_paths_from_raw(d, cfg.scene)
    if reload_models || model_paths != cfg.model_paths
        cfg.model_paths = model_paths
        cfg.model_raw = OrderedDict{String,Any}[
            isfile(path) ? _load_yaml_ordered(path) : OrderedDict{String,Any}() for path in model_paths
        ]
    end
    return cfg
end

function _safe_output_relpath(path::AbstractString)
    raw = replace(path, '\\' => '/')
    tokens = filter(!isempty, split(raw, '/'))
    isempty(tokens) && return basename(path)
    rel = joinpath(tokens...)
    if isabspath(rel) || startswith(normpath(rel), "..")
        return basename(rel)
    end
    return rel
end

"""
    read_light_config(path)::LightConfig

Read a YAML configuration file and normalize ARCHIMED light options into a `LightConfig`.

`pixel_size` is interpreted like Java input files (centimeters) and converted to meters.
"""
function read_light_config(path::AbstractString)
    d = _load_yaml_ordered(path)
    d["__base_dir"] = dirname(path)
    cfg = LightConfig(
        "",
        "",
        false,
        46,
        0.0025,
        true,
        true,
        20,
        0.01,
        0.15,
        0.30,
        false,
        d,
        String(path),
        String[],
        OrderedDict{String,Any}[],
    )
    return refresh_light_config!(cfg; reload_models=true)
end

"""
    write_light_inputs(outdir, cfg; scene_rel, config_name="config.yml")::String

Write the ordered YAML inputs back to disk, preserving parameter order from import.
`scene_rel` should point to the scene file to reference from the written config.
"""
function write_light_inputs(
    outdir::AbstractString,
    cfg::LightConfig;
    scene_rel::AbstractString,
    config_name::AbstractString="config.yml",
)
    outroot = String(outdir)
    mkpath(outroot)

    raw = copy(cfg.raw)
    delete!(raw, "__base_dir")

    model_rels = String[]
    model_specs = get(raw, "models", nothing)
    if model_specs isa AbstractVector
        for i in eachindex(cfg.model_raw)
            rel = i <= length(model_specs) ? _safe_output_relpath(string(model_specs[i])) : joinpath("models", "model$(i).yml")
            path = joinpath(outroot, rel)
            mkpath(dirname(path))
            YAML.write_file(path, cfg.model_raw[i])
            push!(model_rels, rel)
        end
        raw["models"] = model_rels
    end

    meteo_rel = basename(cfg.meteo)
    if isfile(cfg.meteo)
        cp(cfg.meteo, joinpath(outroot, meteo_rel); force=true)
    end
    raw["scene"] = _safe_output_relpath(scene_rel)
    raw["meteo"] = meteo_rel

    cfg_path = joinpath(outroot, config_name)
    YAML.write_file(cfg_path, raw)
    return cfg_path
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

function _relabel_scene_node_ids!(root)
    MultiScaleTreeGraph.traverse!(root) do node
        setfield!(node, :id, _NEXT_SCENE_NODE_ID[])
        _NEXT_SCENE_NODE_ID[] += 1
    end
    return root
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
    _relabel_scene_node_ids!(mtg)
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
