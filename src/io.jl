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

const _GENERAL_BOOL_KEYS = Set([
    "all_in_turtle",
    "photosynthesis",
    "scattering",
    "cache_pixel_table",
    "toricity",
    "cache_radiation",
    "area_ratio",
    "nir_interception",
    "nir_scattering",
    "log_debug",
    "debug",
])

const _GENERAL_INT_KEYS = Set([
    "sky_sectors",
    "radiation_timestep",
    "scattering_max_iter",
])

const _GENERAL_FLOAT_KEYS = Set([
    "scene_rotation",
    "scattering_stop_ratio",
    "scattering_coeff_par",
    "scattering_coeff_nir",
])

const _OUTPUT_CONFIG_ORDER = [
    "output_directory",
    "simulation_directory",
    "write_summary",
    "export_ops",
    "component_variables",
    "scene_variables",
    "opf_variables",
    "opf_overwrite_variables",
]

function _normalize_group_name(x)
    s = strip(string(x))
    isempty(s) && return s
    replace(s, r"^#\[[^\]]+\]\s*" => "")
end

function _normalize_ordered_dict(x)
    x isa AbstractDict || return OrderedDict{String,Any}()
    out = OrderedDict{String,Any}()
    for (k, v) in x
        key = string(k)
        out[key] =
            if v isa AbstractDict
                _normalize_ordered_dict(v)
            elseif v isa AbstractVector
                Any[
                    item isa AbstractDict ? _normalize_ordered_dict(item) : item for item in v
                ]
            else
                v
            end
    end
    return out
end

function _normalize_general_value(key::String, value; from_yaml::Bool=false)
    if key == "pixel_size"
        v = _as_float(value, 0.25)
        return from_yaml ? (v / 100.0) : v
    elseif key in _GENERAL_BOOL_KEYS
        return _as_bool(value, false)
    elseif key in _GENERAL_INT_KEYS
        return _as_int(value, 0)
    elseif key in _GENERAL_FLOAT_KEYS
        return _as_float(value, 0.0)
    elseif key == "meteo_range"
        value === nothing && return nothing
        s = strip(string(value))
        return isempty(s) ? nothing : s
    end
    return value isa AbstractDict ? _normalize_ordered_dict(value) : value
end

function _parse_general_config(raw::AbstractDict{String,Any}; from_yaml::Bool=false)
    out = OrderedDict{String,Any}()
    for (k, v) in raw
        key = string(k)
        key in ("scene", "models", "meteo") && continue
        key in _OUTPUT_CONFIG_ORDER && continue
        out[key] = _normalize_general_value(key, v; from_yaml=from_yaml)
    end
    return out
end

function _bool_flags_dict(x)
    x isa AbstractDict || return OrderedDict{String,Bool}()
    out = OrderedDict{String,Bool}()
    for (k, v) in x
        out[string(k)] = _as_bool(v, false)
    end
    return out
end

function _parse_outputs(raw::AbstractDict{String,Any})
    vars = LightOutputVariables(
        _bool_flags_dict(get(raw, "component_variables", OrderedDict{String,Any}())),
        _bool_flags_dict(get(raw, "scene_variables", OrderedDict{String,Any}())),
        _bool_flags_dict(get(raw, "opf_variables", OrderedDict{String,Any}())),
    )
    return LightOutputsConfig(
        string(get(raw, "output_directory", "output")),
        strip(string(get(raw, "simulation_directory", ""))),
        _as_bool(get(raw, "write_summary", false), false),
        get(raw, "export_ops", nothing),
        _as_bool(get(raw, "opf_overwrite_variables", true), true),
        vars,
    )
end

function _normalize_outputs!(cfg::LightConfig)
    out = cfg.outputs
    out.directory = string(out.directory)
    out.simulation_directory = strip(string(out.simulation_directory))
    out.write_summary = _as_bool(out.write_summary, false)
    out.opf_overwrite_variables = _as_bool(out.opf_overwrite_variables, true)
    out.variables.component = _bool_flags_dict(out.variables.component)
    out.variables.scene = _bool_flags_dict(out.variables.scene)
    out.variables.opf = _bool_flags_dict(out.variables.opf)
    return cfg
end

function _model_group_from_file(path::AbstractString)
    raw = _load_yaml_ordered(path)
    group = haskey(raw, "Group") ? _normalize_group_name(raw["Group"]) : ""
    isempty(group) && error("Model file $(path) is missing a non-empty `Group` entry.")
    return group
end

function _source_files_from_yaml(path::AbstractString, raw::AbstractDict{String,Any})
    base = dirname(path)
    haskey(raw, "scene") || error("Missing config key `scene`. Expected a top-level key in the YAML config.")
    haskey(raw, "meteo") || error("Missing config key `meteo`. Expected a top-level key in the YAML config.")
    scene = _join_if_relative(base, string(raw["scene"]))
    meteo = _join_if_relative(base, string(raw["meteo"]))

    model_paths = OrderedDict{String,String}()
    models = get(raw, "models", nothing)
    if models isa AbstractVector
        for entry in models
            abs_path = _join_if_relative(base, string(entry))
            group = _model_group_from_file(abs_path)
            haskey(model_paths, group) && error("Duplicate model group $(repr(group)) in config models list.")
            model_paths[group] = abs_path
        end
    end

    const_path = joinpath(base, "const.yml")
    return LightConfigSourceFiles(
        String(path),
        scene,
        meteo,
        model_paths,
        isfile(const_path) ? const_path : nothing,
        base,
    )
end

function _load_models_from_paths(model_paths::OrderedDict{String,String})
    out = OrderedDict{String,OrderedDict{String,Any}}()
    for (group, path) in model_paths
        raw = _load_yaml_ordered(path)
        types = get(raw, "Type", nothing)
        types isa AbstractDict || (types = OrderedDict{String,Any}())
        out[group] = _normalize_ordered_dict(types)
    end
    return out
end

function _load_constants(const_path::Union{Nothing,String})
    const_path === nothing && return OrderedDict{String,Any}()
    isfile(const_path) || return OrderedDict{String,Any}()
    return _normalize_ordered_dict(_load_yaml_ordered(const_path))
end

function _config_base_dir(cfg::LightConfig)
    return cfg.source_files.base_dir
end

function config_value(cfg::LightConfig, key::String, default=nothing)
    return get(cfg.general, key, default)
end

function output_variable_flags(cfg::LightConfig, kind::Symbol)
    vars = cfg.outputs.variables
    if kind === :component
        return vars.component
    elseif kind === :scene
        return vars.scene
    elseif kind === :opf
        return vars.opf
    end
    error("Unsupported output variable kind: $(kind)")
end

function model_type_configs(cfg::LightConfig)
    out = NamedTuple[]
    for (group, types) in cfg.models
        for (type_name, params) in types
            push!(out, (group=group, type=type_name, params=params))
        end
    end
    return out
end

function refresh_light_config!(cfg::LightConfig; reload_models::Bool=false)
    src = cfg.source_files
    src.base_dir = isempty(src.base_dir) ? dirname(src.config) : src.base_dir
    src.scene = _join_if_relative(src.base_dir, string(src.scene))
    src.meteo = _join_if_relative(src.base_dir, string(src.meteo))

    norm_models = OrderedDict{String,String}()
    for (group, path) in src.models
        norm_models[string(group)] = _join_if_relative(src.base_dir, string(path))
    end
    src.models = norm_models

    normalized_general = OrderedDict{String,Any}()
    for (k, v) in cfg.general
        normalized_general[string(k)] = _normalize_general_value(string(k), v; from_yaml=false)
    end
    cfg.general = normalized_general
    _normalize_outputs!(cfg)

    if reload_models
        cfg.models = _load_models_from_paths(src.models)
    else
        normalized_models = OrderedDict{String,OrderedDict{String,Any}}()
        for (group, types) in cfg.models
            normalized_models[string(group)] = _normalize_ordered_dict(types)
        end
        cfg.models = normalized_models
    end

    cfg.constants = _normalize_ordered_dict(cfg.constants)
    if isempty(cfg.constants)
        cfg.constants = _load_constants(src.constants)
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
    cfg = LightConfig(
        String(path),
        _source_files_from_yaml(path, d),
        _parse_general_config(d; from_yaml=true),
        OrderedDict{String,OrderedDict{String,Any}}(),
        _parse_outputs(d),
        OrderedDict{String,Any}(),
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

    raw = OrderedDict{String,Any}()
    raw["scene"] = _safe_output_relpath(scene_rel)

    model_rels = String[]
    for (group, types) in cfg.models
        src_path = get(cfg.source_files.models, group, joinpath("models", "$(group).yml"))
        rel0 = if isabspath(src_path)
            try
                relpath(src_path, _config_base_dir(cfg))
            catch
                basename(src_path)
            end
        else
            src_path
        end
        rel = _safe_output_relpath(rel0)
        path = joinpath(outroot, rel)
        mkpath(dirname(path))
        YAML.write_file(path, OrderedDict{String,Any}("Group" => group, "Type" => _normalize_ordered_dict(types)))
        push!(model_rels, rel)
    end
    !isempty(model_rels) && (raw["models"] = model_rels)

    meteo_rel = basename(cfg.source_files.meteo)
    if isfile(cfg.source_files.meteo)
        cp(cfg.source_files.meteo, joinpath(outroot, meteo_rel); force=true)
    end
    raw["meteo"] = meteo_rel

    for (k, v) in cfg.general
        raw[k] = k == "pixel_size" ? (Float64(v) * 100.0) : v
    end

    raw["output_directory"] = cfg.outputs.directory
    !isempty(cfg.outputs.simulation_directory) && (raw["simulation_directory"] = cfg.outputs.simulation_directory)
    raw["write_summary"] = cfg.outputs.write_summary
    cfg.outputs.export_ops !== nothing && (raw["export_ops"] = cfg.outputs.export_ops)
    !isempty(cfg.outputs.variables.component) && (raw["component_variables"] = cfg.outputs.variables.component)
    !isempty(cfg.outputs.variables.scene) && (raw["scene_variables"] = cfg.outputs.variables.scene)
    !isempty(cfg.outputs.variables.opf) && (raw["opf_variables"] = cfg.outputs.variables.opf)
    raw["opf_overwrite_variables"] = cfg.outputs.opf_overwrite_variables

    cfg_path = joinpath(outroot, config_name)
    YAML.write_file(cfg_path, raw)

    if !isempty(cfg.constants)
        const_name = cfg.source_files.constants === nothing ? "const.yml" : basename(cfg.source_files.constants)
        YAML.write_file(joinpath(outroot, const_name), cfg.constants)
    end
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
