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

function _as_string(x, default::String)
    x === nothing && return default
    return String(x)
end

function _config_get(raw::AbstractDict, keys::Tuple{Vararg{String}})
    for key in keys
        haskey(raw, key) && return raw[key]
    end
    lower_keys = map(lowercase, keys)
    for (k, v) in raw
        kk = lowercase(string(k))
        kk in lower_keys && return v
    end
    return nothing
end

function _config_output_variable_enabled(raw::AbstractDict, variable::AbstractString)
    wanted = lowercase(variable)
    for block_key in ("component_variables", "opf_variables")
        block = _config_get(raw, (block_key,))
        block isa AbstractDict || continue
        for (k, v) in block
            lowercase(string(k)) == wanted && _as_bool(v, false) && return true
        end
    end
    return false
end

function _join_if_relative(base::AbstractString, p::AbstractString)
    isabspath(p) ? p : normpath(joinpath(base, p))
end

function _load_yaml_ordered(path::AbstractString)
    YAML.load_file(path; dicttype=OrderedDict{String,Any})
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
    has_par = haskey(values, "PAR")
    has_nir = haskey(values, "NIR")
    par = _as_float(get(values, "PAR", 0.0), 0.0)
    nir = _as_float(get(values, "NIR", 0.0), 0.0)
    delete!(values, "PAR")
    delete!(values, "NIR")
    values["__has_par"] = has_par
    values["__has_nir"] = has_nir
    OpticalProperties(par, nir, values)
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
    delete!(extras, "sensor")
    delete!(extras, "transparency")
    delete!(extras, "optical_properties")

    model_name = string(get(selected, "model", ""))
    sensor_flag = _as_bool(get(selected, "sensor", nothing), false) || lowercase(strip(model_name)) == "virtualsensor"

    InterceptionModel(
        use=use_name,
        model=model_name,
        sensor=sensor_flag,
        transparency=_as_float(get(selected, "transparency", 0.0), 0.0),
        optical_properties=_parse_optical_properties(get(selected, "optical_properties", nothing)),
        variants=variants,
        extras=extras,
    )
end

function _parse_emitter_model(raw)
    raw isa AbstractDict || return nothing
    values = _ordered_string_dict(raw)
    extras = copy(values)
    delete!(extras, "model")
    delete!(extras, "radiance")
    delete!(extras, "gamma")
    gamma = _parse_optical_properties(get(values, "gamma", nothing))
    gamma === nothing && (gamma = OpticalProperties(0.48, 0.52))
    EmitterModel(
        model=string(get(values, "model", "")),
        radiance=_as_float(get(values, "radiance", 0.0), 0.0),
        gamma=gamma,
        extras=extras,
    )
end

function _parse_type_model(raw)
    values = _ordered_string_dict(raw)
    extras = copy(values)
    interception = _parse_interception_model(get(values, "Interception", nothing))
    light_emitter = _parse_emitter_model(get(values, "LightEmitter", nothing))
    delete!(extras, "Interception")
    delete!(extras, "LightEmitter")
    TypeModel(; interception=interception, light_emitter=light_emitter, extras=extras)
end

function _parse_group_model(path::AbstractString)
    raw = _load_yaml_ordered(path)
    group = haskey(raw, "Group") ? strip(string(raw["Group"])) : ""
    isempty(group) && error("Model file $(path) is missing `Group`.")
    types_raw = get(raw, "Type", nothing)
    types_raw isa AbstractDict || error("Model file $(path) is missing `Type` definitions.")
    types = OrderedDict{String,TypeModel}()
    for (type_name, spec) in types_raw
        types[string(type_name)] = _parse_type_model(spec)
    end
    extras = copy(raw)
    delete!(extras, "Group")
    delete!(extras, "Type")
    GroupModel(group; types=types, extras=extras)
end

function _model_paths_from_config(path::AbstractString)
    raw = _load_yaml_ordered(path)
    models = get(raw, "models", nothing)
    models isa AbstractVector || return nothing
    base = dirname(path)
    [_join_if_relative(base, string(m)) for m in models]
end

"""
    read_models(path_or_paths)::LightModels

Read one or more ARCHIMED model YAML files and return them as a
[`LightModels`](@ref) collection.

Arguments:

- `path_or_paths`: a model source to read.

`path_or_paths` can be:
- the path to a single model YAML file containing a `Group`
- the path to a config YAML file containing a `models:` list
- a vector of explicit model-file paths
"""
function read_models(path_or_paths)
    model_paths =
        if path_or_paths isa AbstractString
            if endswith(lowercase(path_or_paths), ".yml") || endswith(lowercase(path_or_paths), ".yaml")
                raw = _load_yaml_ordered(path_or_paths)
                if haskey(raw, "Group")
                    [String(path_or_paths)]
                else
                    paths = _model_paths_from_config(path_or_paths)
                    paths === nothing && error("Unsupported models YAML source: $(path_or_paths)")
                    paths
                end
            else
                error("Unsupported models source: $(path_or_paths)")
            end
        elseif path_or_paths isa AbstractVector
            [String(p) for p in path_or_paths]
        else
            error("Unsupported models source type: $(typeof(path_or_paths))")
        end

    groups = OrderedDict{String,GroupModel}()
    for path in model_paths
        group = _parse_group_model(path)
        haskey(groups, group.group) && error("Duplicate model group `$(group.group)`.")
        groups[group.group] = group
    end
    LightModels(groups)
end

"""
    prepare_models(models)

Normalize in-memory model definitions into a [`LightModels`](@ref) object.

Arguments:

- `models`: an in-memory model specification.

Accepted inputs are an existing `LightModels`, a single [`GroupModel`](@ref), a
vector of groups, or an `OrderedDict{String,GroupModel}`.
"""
function prepare_models(models::LightModels)
    models
end

function prepare_models(groups::AbstractVector{<:GroupModel})
    LightModels(groups)
end

function prepare_models(group::GroupModel)
    LightModels([group])
end

function prepare_models(groups::OrderedDict{String,GroupModel})
    LightModels(groups)
end

"""
    read_options(path)::LightOptions

Read runtime light options from a config YAML file.

This parses the ARCHIMED configuration keys such as `sky_sectors`,
`pixel_size`, `toricity`, scattering controls, and meteo-range options into a
[`LightOptions`](@ref) instance.

Arguments:

- `path`: path to an ARCHIMED-style configuration YAML file.
"""
function read_options(path::AbstractString)
    raw = _load_yaml_ordered(path)
    return LightOptions(
        all_in_turtle=_as_bool(get(raw, "all_in_turtle", false), false),
        turtle_sectors=_as_int(get(raw, "sky_sectors", 46), 46),
        pixel_size=_as_float(get(raw, "pixel_size", 0.25), 0.25) / 100.0,
        area_ratio=_as_bool(get(raw, "area_ratio", true), true),
        scattering=_as_bool(get(raw, "scattering", true), true),
        scattering_max_iter=_as_int(get(raw, "scattering_max_iter", 20), 20),
        scattering_stop_ratio=_as_float(get(raw, "scattering_stop_ratio", 0.01), 0.01),
        scattering_coeff_par=_as_float(get(raw, "scattering_coeff_par", 0.15), 0.15),
        scattering_coeff_nir=_as_float(get(raw, "scattering_coeff_nir", 0.30), 0.30),
        cache_radiation=_as_bool(get(raw, "cache_radiation", false), false),
        include_sky_fraction=_config_output_variable_enabled(raw, "sky_fraction"),
        cache_pixel_table=_as_bool(get(raw, "cache_pixel_table", false), false),
        pixel_hit_stack_mode=_as_string(get(raw, "pixel_hit_stack_mode", "auto"), "auto"),
        toricity=_as_bool(get(raw, "toricity", true), true),
        radiation_timestep_minutes=_as_float(get(raw, "radiation_timestep", 15.0), 15.0),
        allow_overlapping_meteo_steps=_as_bool(
            _config_get(raw, ("allow_overlapping_meteo_steps", "allowOverlappingMeteoSteps", "allowOverlappingMeteo")),
            false,
        ),
        nir_interception=_as_bool(get(raw, "nir_interception", true), true),
        nir_scattering=_as_bool(get(raw, "nir_scattering", true), true),
        java_logged_turtle_dirs=_as_bool(get(raw, "java_logged_turtle_dirs", false), false),
        meteo_range=begin
            v = get(raw, "meteo_range", nothing)
            v === nothing ? nothing : strip(string(v))
        end,
        debug=_as_bool(get(raw, "debug", false), false),
        log_debug=_as_bool(get(raw, "log_debug", false), false),
        debug_drop_leading_hit=begin
            v = get(raw, "debug_drop_leading_hit", nothing)
            if v isa AbstractDict
                node_id = _as_int(get(v, "node_id", nothing), 0)
                x = _as_int(get(v, "x", nothing), -1)
                y = _as_int(get(v, "y", nothing), -1)
                node_id > 0 && x >= 0 && y >= 0 ? (node_id=node_id, x=x, y=y) : nothing
            else
                nothing
            end
        end,
    )
end

function _config_ground_spec(models::LightModels)
    count = 0
    group_name = "pavement"
    type_name = "Cobblestone"
    for (group, group_model) in models.groups
        for (type_key, type_model) in group_model.types
            if haskey(type_model.extras, "plot_paving")
                paving_count = _as_int(type_model.extras["plot_paving"], 0)
                if paving_count > count
                    count = paving_count
                    group_name = group
                    type_name = type_key
                end
            end
        end
    end
    return (count=count, group=group_name, type=type_name)
end

function _scene_has_group_type(scene::PlantGeom.SceneGeometry, group::AbstractString, type::AbstractString)
    any(nid -> _scene_group(scene, nid, "") == group && _scene_type(scene, nid, "") == type, keys(scene.nodes))
end

function _materialize_paving!(
    scene::PlantGeom.SceneGeometry,
    count::Int;
    xy_bounds=nothing,
    group::AbstractString="pavement",
    type::AbstractString="Cobblestone",
)
    count > 0 || return scene
    scene.mtg === nothing && error("_materialize_paving! requires an MTG-backed scene.")
    bounds = xy_bounds === nothing ? scene.scene_xy_bounds : xy_bounds
    bounds === nothing && error("Ground bounds are undefined. Pass `xy_bounds=` or use a scene with known bounds.")
    xmin, ymin, xmax, ymax = bounds
    plot_x = Float64(Float32(Float32(xmax) - Float32(xmin)))
    plot_y = Float64(Float32(Float32(ymax) - Float32(ymin)))
    cobble_edge = sqrt((plot_x * plot_y) / max(count, 1))
    nx = max(1, floor(Int, plot_x / cobble_edge))
    ny = max(1, floor(Int, plot_y / cobble_edge))
    PlantGeom.add_ground!(
        scene;
        z=0.005,
        nx=nx,
        ny=ny,
        xy_bounds=(Float64(xmin), Float64(ymin), Float64(xmax), Float64(ymax)),
        group=group,
        type=type,
    )
end

"""
    read_config(path; plot_paving_override=nothing)

Read a complete simulation configuration from `path` and return
`(options, scene, meteo, models)`.

When `plot_paving_override` is provided, it overrides the paving density
declared in the model extras before the ground is materialized into the scene.
"""
function read_config(path::AbstractString; plot_paving_override=nothing)
    raw = _load_yaml_ordered(path)
    base = dirname(path)
    scene_rel = haskey(raw, "scene") ? string(raw["scene"]) : error("Config file $(path) is missing `scene`.")
    meteo_rel = haskey(raw, "meteo") ? string(raw["meteo"]) : error("Config file $(path) is missing `meteo`.")

    options = read_options(path)
    scene = read_scene(_join_if_relative(base, scene_rel))
    meteo = read_meteo(_join_if_relative(base, meteo_rel))
    models = read_models(path)

    ground = _config_ground_spec(models)
    ground_count = plot_paving_override === nothing ? ground.count : _as_int(plot_paving_override, ground.count)
    if ground_count > 0 && !_scene_has_group_type(scene, ground.group, ground.type)
        _materialize_paving!(scene, ground_count; group=ground.group, type=ground.type)
    end

    return options, scene, meteo, models
end

"""
    read_simulation(path; plot_paving_override=nothing, kwargs...)

Read a complete file-based light simulation and return `(sim, meteo)`, where
`sim` is a [`LightSimulation`](@ref).

Arguments:

- `path`: path to the ARCHIMED-style configuration YAML file.

Keywords:

- `plot_paving_override`: optional replacement paving count used when
  materializing model-declared ground geometry.
- `kwargs...`: keyword arguments forwarded to [`LightSimulation`](@ref), such as
  `interception_backend`, `scattering_mode`, `scattering_backend`, and
  `memory_limit_bytes`.
"""
function read_simulation(path::AbstractString; plot_paving_override=nothing, kwargs...)
    options, scene, meteo, models = read_config(path; plot_paving_override=plot_paving_override)
    meteo_eff = prepare_meteo(meteo, options)
    sim_options = options.meteo_range === nothing ? options : LightOptions(options; meteo_range=nothing)
    return LightSimulation(scene, models; options=sim_options, kwargs...), meteo_eff
end

function _triangle_area3d(p1, p2, p3)
    v1 = StaticArrays.SVector{3,Float64}(p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = StaticArrays.SVector{3,Float64}(p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
    0.5 * norm(cross(v1, v2))
end

"""
    read_scene(path)::PlantGeom.SceneGeometry

Read a scene file (`.ops`, `.opf`, or `.gwa`) and return a prepared
`PlantGeom.SceneGeometry`.

The scene is relabelled into a dense node-id space and immediately converted to
the merged-mesh representation expected by the interception pipeline.

Arguments:

- `path`: path to an `.ops`, `.opf`, or `.gwa` scene file.

Keywords:

- `plantgeom_backend`: reserved PlantGeom backend selector. The current default
  is `:auto`.
"""
function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    mtg = if ext == ".ops"
        PlantGeom.read_ops(path; relaxed=true, assume_scale_column=false, opf_scale=1.0, gwa_scale=0.01,)
    elseif ext == ".opf"
        PlantGeom.read_opf(path, attr_type=Dict, attribute_types=Dict("pos" => Float64))
    elseif ext == ".gwa"
        PlantGeom.read_gwa(path)
    else
        error("Unsupported scene extension: $ext")
    end
    PlantGeom.prepare_scene(mtg; source_path=String(path), scene_xy_bounds=PlantGeom._scene_xy_bounds_from_mtg(mtg), relabel_ids=true)
end

function _refresh_ref_mesh_registry!(mtg)
    ref_meshes = PlantGeom.get_ref_meshes(mtg)
    mtg[:ref_meshes] = Dict(i - 1 => mesh_ for (i, mesh_) in enumerate(ref_meshes))
    return mtg[:ref_meshes]
end

function _normalize_source_topology_ids!(mtg)
    MultiScaleTreeGraph.traverse!(mtg) do node
        if !haskey(node, :source_topology_id) || node[:source_topology_id] === nothing
            node[:source_topology_id] = MultiScaleTreeGraph.node_id(node)
        end
        return true
    end
    return mtg
end

"""
    write_scene(path, scene)

Write an MTG-backed `PlantGeom.SceneGeometry` to `path`.

Supported output formats are `.ops`, `.opf`, and `.gwa`. The function refreshes
the reference-mesh registry and normalizes topology ids before export.

Arguments:

- `path`: output path ending in `.ops`, `.opf`, or `.gwa`.
- `scene`: MTG-backed `PlantGeom.SceneGeometry` to write.
"""
function write_scene(path::AbstractString, scene::PlantGeom.SceneGeometry)
    scene.mtg === nothing && error("write_scene requires an MTG-backed scene.")
    ext = lowercase(splitext(path)[2])
    mkpath(dirname(path))
    _normalize_source_topology_ids!(scene.mtg)
    _refresh_ref_mesh_registry!(scene.mtg)
    if ext == ".ops"
        PlantGeom.write_ops(path, scene.mtg; write_objects=true, preserve_file_paths=true)
    elseif ext == ".opf"
        PlantGeom.write_opf(path, scene.mtg)
    elseif ext == ".gwa"
        PlantGeom.write_gwa(path, scene.mtg)
    else
        error("Unsupported scene export extension: $ext")
    end
    return path
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

function _positive_duration_seconds(v; field_name::AbstractString="step_duration")
    PlantMeteo.positive_duration_seconds(v; field_name=field_name)
end

function _table_column(data, name::Symbol)
    Tables.getcolumn(Tables.columns(data), name)
end

function _normalize_raw_meteo_dates(data)
    :date in Tables.columnnames(data) || return data
    hour_fmt = Dates.DateFormat("HH:MM:SS")
    date_only = (; date=data.date)
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

function _table_metadata_namedtuple(data)
    meta =
        try
            PlantMeteo.table_metadata(data)
        catch
            NamedTuple()
        end
    _meta_to_namedtuple(meta)
end

function _with_meteo_metadata(table::PlantMeteo.TimeStepTable, metadata)
    meta = _meta_to_namedtuple(metadata)
    PlantMeteo.TimeStepTable(getfield(table, :names), meta, getfield(table, :ts))
end

function _as_plantmeteo_table(data; metadata=_table_metadata_namedtuple(data))
    data isa PlantMeteo.TimeStepTable && return data
    transformed = Tables.columntable(data)
    column_names = Tables.columnnames(transformed)
    if !(:date in column_names) && :DateTime in column_names
        transformed = PlantMeteo.set_column(transformed, :date, transformed.DateTime)
        column_names = Tables.columnnames(transformed)
    end
    transformed = _normalize_raw_meteo_dates(transformed)
    column_names = Tables.columnnames(transformed)
    if !(:duration in column_names)
        duration = PlantMeteo.compute_duration(transformed, Dates.DateFormat("HH:MM:SS"), nothing)
        transformed = PlantMeteo.set_column(transformed, :duration, duration)
    end
    PlantMeteo.TimeStepTable(transformed, _meta_to_namedtuple(metadata))
end

"""
    read_meteo(path)::PlantMeteo.TimeStepTable
    read_meteo(data)::PlantMeteo.TimeStepTable

Read a meteorological forcing table from `path` and return it as a
`PlantMeteo.TimeStepTable`.

The resulting table keeps available metadata such as latitude, longitude,
altitude, and source file path.

Arguments:

- `path`: path to a meteorological forcing file readable by PlantMeteo.
- `data`: alternatively, a Tables.jl-compatible table or an existing
  `PlantMeteo.TimeStepTable`.
"""
function read_meteo(path::AbstractString)
    try
        weather = PlantMeteo.read_weather(
            path;
            date_formats=(Dates.DateFormat("yyyy/mm/dd"), Dates.DateFormat("yyyy-mm-dd")),
            forward_fill_date=true,
        )
        meta = merge((; file=path), _table_metadata_namedtuple(weather))
        return _with_meteo_metadata(weather, meta)
    catch
        data, metadata_ = PlantMeteo.read_weather_(path)
        data = _normalize_raw_meteo_dates(data)
        meta = merge((; file=path), _meta_to_namedtuple(metadata_))
        return _as_plantmeteo_table(data; metadata=meta)
    end
end

function read_meteo(data)
    Tables.istable(typeof(data)) || error("Unsupported meteo input: expected a path or a Tables.jl-compatible table.")
    data isa PlantMeteo.TimeStepTable && return data
    _as_plantmeteo_table(data)
end
