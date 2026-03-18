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

_ordered_dict_copy(v::AbstractDict) = OrderedDict{String,Any}(string(k) => val for (k, val) in v)
_ordered_dict_copy(::Any) = OrderedDict{String,Any}()

function _ordered_string_dict(v)
    out = OrderedDict{String,Any}()
    v isa AbstractDict || return out
    for (k, val) in v
        out[string(k)] = val
    end
    out
end

function _parse_optical_properties(raw)
    raw isa AbstractDict || return nothing
    values = _ordered_string_dict(raw)
    par = _as_float(get(values, "PAR", 0.0), 0.0)
    nir = _as_float(get(values, "NIR", 0.0), 0.0)
    delete!(values, "PAR")
    delete!(values, "NIR")
    OpticalProperties(par, nir, values)
end

function _serialize_optical_properties(props::OpticalProperties)
    out = copy(props.extras)
    out["PAR"] = props.par
    out["NIR"] = props.nir
    out
end

function _selected_interception_block(raw::OrderedDict{String,Any})
    use_name = get(raw, "use", nothing)
    if use_name !== nothing
        key = string(use_name)
        if haskey(raw, key) && raw[key] isa AbstractDict
            return key, _ordered_string_dict(raw[key])
        end
    end
    return nothing, raw
end

function _parse_interception_model(raw)
    raw isa AbstractDict || return nothing
    raw_dict = _ordered_string_dict(raw)
    use_name, selected = _selected_interception_block(raw_dict)

    variants = OrderedDict{String,OrderedDict{String,Any}}()
    for (k, v) in raw_dict
        v isa AbstractDict || continue
        variants[k] = _ordered_string_dict(v)
    end

    extras = copy(raw_dict)
    delete!(extras, "use")
    for key in keys(variants)
        delete!(extras, key)
    end
    delete!(extras, "model")
    delete!(extras, "transparency")
    delete!(extras, "optical_properties")

    InterceptionModelConfig(
        use_name,
        string(get(selected, "model", "")),
        _as_float(get(selected, "transparency", 0.0), 0.0),
        _parse_optical_properties(get(selected, "optical_properties", nothing)),
        variants,
        extras,
    )
end

function _serialize_interception_model(model::InterceptionModelConfig)
    out = copy(model.extras)
    if model.use !== nothing
        out["use"] = model.use
    end
    for (name, variant) in model.variants
        out[name] = copy(variant)
    end
    target =
        if model.use !== nothing
            get!(out, model.use) do
                OrderedDict{String,Any}()
            end
        else
            out
        end
    target["model"] = model.model
    target["transparency"] = model.transparency
    model.optical_properties === nothing || (target["optical_properties"] = _serialize_optical_properties(model.optical_properties))
    out
end

function _parse_light_emitter(raw)
    raw isa AbstractDict || return nothing
    values = _ordered_string_dict(raw)
    extras = copy(values)
    delete!(extras, "model")
    delete!(extras, "radiance")
    delete!(extras, "gamma")
    gamma = _parse_optical_properties(get(values, "gamma", nothing))
    gamma === nothing && (gamma = OpticalProperties(0.48, 0.52, OrderedDict{String,Any}()))
    LightEmitterConfig(
        string(get(values, "model", "")),
        _as_float(get(values, "radiance", 0.0), 0.0),
        gamma,
        extras,
    )
end

function _serialize_light_emitter(model::LightEmitterConfig)
    out = copy(model.extras)
    out["model"] = model.model
    out["radiance"] = model.radiance
    out["gamma"] = _serialize_optical_properties(model.gamma)
    out
end

function _parse_type_model(raw)
    values = _ordered_string_dict(raw)
    extras = copy(values)
    interception = _parse_interception_model(get(values, "Interception", nothing))
    light_emitter = _parse_light_emitter(get(values, "LightEmitter", nothing))
    plot_paving = _as_int(get(values, "plot_paving", 0), 0)
    delete!(extras, "Interception")
    delete!(extras, "LightEmitter")
    delete!(extras, "plot_paving")
    TypeModelConfig(interception, light_emitter, plot_paving, extras)
end

function _serialize_type_model(model::TypeModelConfig)
    out = copy(model.extras)
    model.interception === nothing || (out["Interception"] = _serialize_interception_model(model.interception))
    model.light_emitter === nothing || (out["LightEmitter"] = _serialize_light_emitter(model.light_emitter))
    model.plot_paving > 0 && (out["plot_paving"] = model.plot_paving)
    out
end

function _parse_group_model(path::AbstractString)
    raw = _load_yaml_ordered(path)
    group = haskey(raw, "Group") ? strip(string(raw["Group"])) : ""
    isempty(group) && error("Model file $(path) is missing `Group`.")
    types_raw = get(raw, "Type", nothing)
    types_raw isa AbstractDict || error("Model file $(path) is missing `Type` definitions.")
    types = OrderedDict{String,TypeModelConfig}()
    for (type_name, spec) in types_raw
        types[string(type_name)] = _parse_type_model(spec)
    end
    extras = copy(raw)
    delete!(extras, "Group")
    delete!(extras, "Type")
    return GroupModelConfig(group, types, extras)
end

function _serialize_group_model(model::GroupModelConfig)
    out = copy(model.extras)
    out["Group"] = model.group
    types = OrderedDict{String,Any}()
    for (type_name, spec) in model.types
        types[type_name] = _serialize_type_model(spec)
    end
    out["Type"] = types
    out
end

function _model_configs_from_paths(paths::Vector{String})
    model_paths = OrderedDict{String,String}()
    models = OrderedDict{String,GroupModelConfig}()
    for path in paths
        model = _parse_group_model(path)
        haskey(models, model.group) && error("Duplicate model group `$(model.group)` in config.")
        model_paths[model.group] = path
        models[model.group] = model
    end
    return model_paths, models
end

function _ordered_bool_dict(v)
    out = OrderedDict{String,Bool}()
    v isa AbstractDict || return out
    for (k, enabled) in v
        out[string(k)] = _as_bool(enabled, false)
    end
    out
end

function _parse_debug_drop_leading_hit(v)
    v isa AbstractDict || return nothing
    node_id = _as_int(get(v, "node_id", nothing), 0)
    x = _as_int(get(v, "x", nothing), -1)
    y = _as_int(get(v, "y", nothing), -1)
    node_id > 0 || return nothing
    x >= 0 || return nothing
    y >= 0 || return nothing
    return (node_id=node_id, x=x, y=y)
end

function _default_general_config()
    LightGeneralConfig(
        false,
        46,
        0.25 / 100.0,
        true,
        true,
        20,
        0.01,
        0.15,
        0.30,
        false,
        false,
        true,
        15.0,
        true,
        true,
        false,
        nothing,
        false,
        false,
        nothing,
    )
end

function _default_outputs_config()
    LightOutputsConfig(
        "output",
        "",
        false,
        nothing,
        OrderedDict{String,Bool}(),
        OrderedDict{String,Bool}(),
        OrderedDict{String,Bool}(),
        true,
    )
end

const _CONFIG_KNOWN_TOP_LEVEL_KEYS = Set([
    "scene",
    "models",
    "meteo",
    "all_in_turtle",
    "sky_sectors",
    "pixel_size",
    "area_ratio",
    "scattering",
    "scattering_max_iter",
    "scattering_stop_ratio",
    "scattering_coeff_par",
    "scattering_coeff_nir",
    "cache_radiation",
    "cache_pixel_table",
    "toricity",
    "radiation_timestep",
    "nir_interception",
    "nir_scattering",
    "java_logged_turtle_dirs",
    "meteo_range",
    "debug",
    "log_debug",
    "debug_drop_leading_hit",
    "output_directory",
    "simulation_directory",
    "write_summary",
    "export_ops",
    "component_variables",
    "scene_variables",
    "opf_variables",
    "opf_overwrite_variables",
])

function _config_extras(raw::OrderedDict{String,Any})
    extras = OrderedDict{String,Any}()
    for (k, v) in raw
        k == "__base_dir" && continue
        k in _CONFIG_KNOWN_TOP_LEVEL_KEYS && continue
        extras[k] = v
    end
    extras
end

function refresh_light_config!(cfg::LightConfig; reload_models::Bool=false)
    base = cfg.paths.base_dir
    cfg.paths.scene = _join_if_relative(base, cfg.paths.scene)
    cfg.paths.meteo = _join_if_relative(base, cfg.paths.meteo)
    for (group, path) in cfg.paths.models
        cfg.paths.models[group] = _join_if_relative(base, path)
    end
    if reload_models
        paths = [path for path in values(cfg.paths.models)]
        cfg.paths.models, cfg.models = _model_configs_from_paths(paths)
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
    raw = _load_yaml_ordered(path)
    raw["__base_dir"] = dirname(path)

    haskey(raw, "scene") || error("Missing config key `scene`. Expected a top-level key in the YAML config.")
    haskey(raw, "meteo") || error("Missing config key `meteo`. Expected a top-level key in the YAML config.")

    base = String(raw["__base_dir"])
    model_paths, models = _model_configs_from_paths(_model_paths_from_raw(raw, path))
    paths = LightConfigPaths(String(path), string(raw["scene"]), string(raw["meteo"]), model_paths, base)

    general = _default_general_config()
    general.all_in_turtle = _as_bool(get(raw, "all_in_turtle", general.all_in_turtle), general.all_in_turtle)
    general.turtle_sectors = _as_int(get(raw, "sky_sectors", general.turtle_sectors), general.turtle_sectors)
    general.pixel_size = _as_float(get(raw, "pixel_size", general.pixel_size * 100.0), general.pixel_size * 100.0) / 100.0
    general.area_ratio = _as_bool(get(raw, "area_ratio", general.area_ratio), general.area_ratio)
    general.scattering = _as_bool(get(raw, "scattering", general.scattering), general.scattering)
    general.scattering_max_iter = _as_int(get(raw, "scattering_max_iter", general.scattering_max_iter), general.scattering_max_iter)
    general.scattering_stop_ratio = _as_float(get(raw, "scattering_stop_ratio", general.scattering_stop_ratio), general.scattering_stop_ratio)
    general.scattering_coeff_par = _as_float(get(raw, "scattering_coeff_par", general.scattering_coeff_par), general.scattering_coeff_par)
    general.scattering_coeff_nir = _as_float(get(raw, "scattering_coeff_nir", general.scattering_coeff_nir), general.scattering_coeff_nir)
    general.cache_radiation = _as_bool(get(raw, "cache_radiation", general.cache_radiation), general.cache_radiation)
    general.cache_pixel_table = _as_bool(get(raw, "cache_pixel_table", general.cache_pixel_table), general.cache_pixel_table)
    general.toricity = _as_bool(get(raw, "toricity", general.toricity), general.toricity)
    general.radiation_timestep_minutes = _as_float(get(raw, "radiation_timestep", general.radiation_timestep_minutes), general.radiation_timestep_minutes)
    general.nir_interception = _as_bool(get(raw, "nir_interception", general.nir_interception), general.nir_interception)
    general.nir_scattering = _as_bool(get(raw, "nir_scattering", general.nir_scattering), general.nir_scattering)
    general.java_logged_turtle_dirs = _as_bool(get(raw, "java_logged_turtle_dirs", general.java_logged_turtle_dirs), general.java_logged_turtle_dirs)
    meteo_range = get(raw, "meteo_range", nothing)
    general.meteo_range = meteo_range === nothing ? nothing : strip(string(meteo_range))
    general.debug = _as_bool(get(raw, "debug", general.debug), general.debug)
    general.log_debug = _as_bool(get(raw, "log_debug", general.log_debug), general.log_debug)
    general.debug_drop_leading_hit = _parse_debug_drop_leading_hit(get(raw, "debug_drop_leading_hit", nothing))

    outputs = _default_outputs_config()
    outputs.output_directory = string(get(raw, "output_directory", outputs.output_directory))
    simdir = get(raw, "simulation_directory", outputs.simulation_directory)
    outputs.simulation_directory = simdir === nothing ? "" : strip(string(simdir))
    outputs.write_summary = _as_bool(get(raw, "write_summary", outputs.write_summary), outputs.write_summary)
    outputs.export_ops = get(raw, "export_ops", outputs.export_ops)
    outputs.component_variables = _ordered_bool_dict(get(raw, "component_variables", nothing))
    outputs.scene_variables = _ordered_bool_dict(get(raw, "scene_variables", nothing))
    outputs.opf_variables = _ordered_bool_dict(get(raw, "opf_variables", nothing))
    outputs.opf_overwrite_variables = _as_bool(get(raw, "opf_overwrite_variables", outputs.opf_overwrite_variables), outputs.opf_overwrite_variables)

    cfg = LightConfig(
        String(path),
        paths,
        general,
        outputs,
        models,
        _config_extras(raw),
    )
    refresh_light_config!(cfg)
end

"""
    write_light_inputs(outdir, cfg; scene_rel, config_name="config.yml")::String

Write the parsed light inputs back to disk.
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

    if !isempty(cfg.paths.models)
        model_rels = String[]
        for (group, model) in cfg.models
            original = get(cfg.paths.models, group, "$(group).yml")
            rel = joinpath("models", basename(original))
            path = joinpath(outroot, rel)
            mkpath(dirname(path))
            YAML.write_file(path, _serialize_group_model(model))
            push!(model_rels, rel)
        end
        raw["models"] = model_rels
    end

    meteo_rel = basename(cfg.paths.meteo)
    if isfile(cfg.paths.meteo)
        cp(cfg.paths.meteo, joinpath(outroot, meteo_rel); force=true)
    end
    raw["meteo"] = meteo_rel

    for (k, v) in cfg.extras
        raw[k] = v
    end

    raw["sky_sectors"] = cfg.general.turtle_sectors
    raw["all_in_turtle"] = cfg.general.all_in_turtle
    raw["radiation_timestep"] = cfg.general.radiation_timestep_minutes
    raw["scattering"] = cfg.general.scattering
    raw["pixel_size"] = cfg.general.pixel_size * 100.0
    raw["cache_pixel_table"] = cfg.general.cache_pixel_table
    raw["toricity"] = cfg.general.toricity
    raw["cache_radiation"] = cfg.general.cache_radiation
    raw["area_ratio"] = cfg.general.area_ratio
    raw["scattering_max_iter"] = cfg.general.scattering_max_iter
    raw["scattering_stop_ratio"] = cfg.general.scattering_stop_ratio
    raw["scattering_coeff_par"] = cfg.general.scattering_coeff_par
    raw["scattering_coeff_nir"] = cfg.general.scattering_coeff_nir
    cfg.general.nir_interception && !haskey(cfg.extras, "nir_interception") || (raw["nir_interception"] = cfg.general.nir_interception)
    cfg.general.nir_scattering && !haskey(cfg.extras, "nir_scattering") || (raw["nir_scattering"] = cfg.general.nir_scattering)
    cfg.general.java_logged_turtle_dirs && (raw["java_logged_turtle_dirs"] = cfg.general.java_logged_turtle_dirs)
    cfg.general.meteo_range === nothing || (raw["meteo_range"] = cfg.general.meteo_range)
    cfg.general.debug && (raw["debug"] = cfg.general.debug)
    cfg.general.log_debug && (raw["log_debug"] = cfg.general.log_debug)
    cfg.general.debug_drop_leading_hit === nothing || (raw["debug_drop_leading_hit"] = OrderedDict(
        "node_id" => cfg.general.debug_drop_leading_hit.node_id,
        "x" => cfg.general.debug_drop_leading_hit.x,
        "y" => cfg.general.debug_drop_leading_hit.y,
    ))

    raw["output_directory"] = cfg.outputs.output_directory
    isempty(cfg.outputs.simulation_directory) || (raw["simulation_directory"] = cfg.outputs.simulation_directory)
    raw["write_summary"] = cfg.outputs.write_summary
    cfg.outputs.export_ops === nothing || (raw["export_ops"] = cfg.outputs.export_ops)
    isempty(cfg.outputs.component_variables) || (raw["component_variables"] = cfg.outputs.component_variables)
    isempty(cfg.outputs.scene_variables) || (raw["scene_variables"] = cfg.outputs.scene_variables)
    isempty(cfg.outputs.opf_variables) || (raw["opf_variables"] = cfg.outputs.opf_variables)
    raw["opf_overwrite_variables"] = cfg.outputs.opf_overwrite_variables

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

function _node_object_id(node, default::Int)
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
    nodes,
    per_object_component_counter,
    current_group::String="",
    current_object_id::Int=1,
)
    group = current_group
    object_id = _node_object_id(node, current_object_id)
    type_name = _node_type_name(node, "")
    fg = _node_attr(node, :functional_group)
    fg !== nothing && (group = string(fg))

    if _is_scene_geometry_node(node)
        nid = Int(MultiScaleTreeGraph.node_id(node))

        sid = _source_topology_id(node)
        if sid !== nothing
            source_topology_id = Int(sid)
        else
            k = get(per_object_component_counter, object_id, 0) + 1
            per_object_component_counter[object_id] = k
            source_topology_id = k + 1
        end
        nodes[nid] = (
            group=group,
            type=type_name,
            source_topology_id=source_topology_id,
            object_id=object_id,
        )
    end

    children = MultiScaleTreeGraph.children(node)
    if children !== nothing
        for ch in children
            _collect_scene_node_metadata!(
                ch,
                nodes,
                per_object_component_counter,
                group,
                object_id,
            )
        end
    end
end

function _build_scene_geometry(mtg, source_path::AbstractString)
    _build_scene_geometry(mtg, source_path, nothing)
end

function _build_scene_geometry(mtg, source_path::AbstractString, scene_xy_bounds::Union{Nothing,NTuple{4,Float64}})
    merged_mesh, face2node = PlantGeom.build_merged_mesh_with_map(mtg; filter_fun=_is_scene_geometry_node)

    node_meta = Dict{Int,NamedTuple{(:group, :type, :source_topology_id, :object_id),Tuple{String,String,Int,Int}}}()
    per_object_component_counter = Dict{Int,Int}()
    _collect_scene_node_metadata!(
        mtg,
        node_meta,
        per_object_component_counter,
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
    nodes = Dict{Int,SceneNodeData{Float64}}()
    for (nid, area) in node_area
        barycenter =
            if area > 0
                sx, sy, sz = get(bary_acc, nid, (0.0, 0.0, 0.0))
                (sx / area, sy / area, sz / area)
            else
                (NaN, NaN, NaN)
            end
        meta = get(node_meta, nid, (group="", type="", source_topology_id=nid, object_id=-1))
        nodes[nid] = SceneNodeData(area, barycenter, meta.group, meta.type, meta.source_topology_id, meta.object_id)
    end
    SceneGeometry(
        mtg,
        merged_mesh,
        face2node,
        nodes,
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
