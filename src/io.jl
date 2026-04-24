import Dates
import GeometryBasics
import LinearAlgebra: cross, norm
import MultiScaleTreeGraph
import OrderedCollections: OrderedDict
import PlantGeom
import PlantMeteo
import StaticArrays
import Tables
import YAML

const _NEXT_SCENE_NODE_ID = Ref(1)

"""
    LightSceneBuilder

Mutable helper used by `light_scene` to assemble an MTG-backed scene before it
is converted to [`SceneGeometry`](@ref).
"""
mutable struct LightSceneBuilder
    mtg
    domain::Union{Nothing,NTuple{4,Float64}}
    source_path::String
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

function _as_string(x, default::String)
    x === nothing && return default
    return String(x)
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
        cache_pixel_table=_as_bool(get(raw, "cache_pixel_table", false), false),
        pixel_hit_stack_mode=_as_string(get(raw, "pixel_hit_stack_mode", "auto"), "auto"),
        toricity=_as_bool(get(raw, "toricity", true), true),
        radiation_timestep_minutes=_as_float(get(raw, "radiation_timestep", 15.0), 15.0),
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

function _scene_has_group_type(scene::SceneGeometry, group::AbstractString, type::AbstractString)
    any(node.group == group && node.type == type for node in values(scene.nodes))
end

function _paving_tile_mesh(vertices, faces_for_tile)
    used = sort!(unique(vcat([[Int(f[1]), Int(f[2]), Int(f[3])] for f in faces_for_tile]...)))
    remap = Dict{Int,Int}(old => new for (new, old) in enumerate(used))
    tile_points = GeometryBasics.Point3f[
        GeometryBasics.Point3f(Float32(vertices[idx][1]), Float32(vertices[idx][2]), Float32(vertices[idx][3])) for idx in used
    ]
    tile_faces = PlantGeom.Face3[
        PlantGeom.Face3(remap[Int(f[1])], remap[Int(f[2])], remap[Int(f[3])]) for f in faces_for_tile
    ]
    return GeometryBasics.Mesh(tile_points, tile_faces)
end

function _materialize_paving!(
    scene::SceneGeometry,
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
    plotbox = (
        origin_x=Float64(xmin),
        origin_y=Float64(ymin),
        xdim=Float64(xmax - xmin),
        ydim=Float64(ymax - ymin),
    )
    first_node = isempty(scene.nodes) ? 1 : (maximum(keys(scene.nodes)) + 1)
    vertices, faces, face2node, _ = _paving_mesh(plotbox, count, first_node)
    faces_by_node = Dict{Int,Vector{PlantGeom.Face3}}()
    for (face, nid) in zip(faces, face2node)
        push!(get!(faces_by_node, nid, PlantGeom.Face3[]), face)
    end

    root_scale = MultiScaleTreeGraph.scale(scene.mtg)
    for nid in sort!(collect(keys(faces_by_node)))
        mesh = _paving_tile_mesh(vertices, faces_by_node[nid])
        ref_mesh = PlantGeom.RefMesh("$(group)_$(nid)", mesh)
        MultiScaleTreeGraph.Node(
            nid,
            scene.mtg,
            MultiScaleTreeGraph.MutableNodeMTG(:+, Symbol(type), nid, root_scale + 1),
            Dict{Symbol,Any}(
                :geometry => PlantGeom.Geometry(ref_mesh=ref_mesh),
                :functional_group => String(group),
                :type => String(type),
                :object_id => -1,
                :source_topology_id => nid,
            ),
        )
    end

    scene.scene_xy_bounds = (Float64(xmin), Float64(ymin), Float64(xmax), Float64(ymax))
    _drop_scene_surface_geometry!(scene.mtg)
    _refresh_ref_mesh_registry!(scene.mtg)
    _refresh_scene!(scene)
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
            _collect_scene_node_metadata!(ch, nodes, per_object_component_counter, group, object_id)
        end
    end
end

function _build_scene_geometry(mtg, source_path::AbstractString, scene_xy_bounds::Union{Nothing,NTuple{4,Float64}})
    merged_mesh, face2node = PlantGeom.build_merged_mesh_with_map(mtg; filter_fun=_is_scene_geometry_node)

    node_meta = Dict{Int,NamedTuple{(:group, :type, :source_topology_id, :object_id),Tuple{String,String,Int,Int}}}()
    per_object_component_counter = Dict{Int,Int}()
    _collect_scene_node_metadata!(mtg, node_meta, per_object_component_counter)

    verts = GeometryBasics.decompose(GeometryBasics.Point3, merged_mesh)
    faces = GeometryBasics.decompose(PlantGeom.Face3, merged_mesh)
    node_area = Dict{Int,Float64}()
    bary_acc = Dict{Int,NTuple{3,Float64}}()
    for (i, f) in enumerate(faces)
        nid = face2node[i]
        p1 = verts[f[1]]
        p2 = verts[f[2]]
        p3 = verts[f[3]]
        area = _triangle_area3d(p1, p2, p3)
        node_area[nid] = get(node_area, nid, 0.0) + area
        cx = (p1[1] + p2[1] + p3[1]) / 3.0
        cy = (p1[2] + p2[2] + p3[2]) / 3.0
        cz = (p1[3] + p2[3] + p3[3]) / 3.0
        sx, sy, sz = get(bary_acc, nid, (0.0, 0.0, 0.0))
        bary_acc[nid] = (sx + area * cx, sy + area * cy, sz + area * cz)
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

    SceneGeometry(mtg, merged_mesh, face2node, nodes, String(source_path), scene_xy_bounds)
end

function _read_ops_relaxed(path::AbstractString)
    PlantGeom.read_ops(path; _OPS_RELAXED_KWARGS...)
end

function _read_opf_relaxed(path::AbstractString)
    PlantGeom.read_opf(path, attr_type=Dict, attribute_types=Dict("pos" => Float64))
end

function _relabel_scene_node_ids!(root)
    MultiScaleTreeGraph.traverse!(root) do node
        setfield!(node, :id, _NEXT_SCENE_NODE_ID[])
        _NEXT_SCENE_NODE_ID[] += 1
    end
    return root
end

function _coerce_domain(domain)
    domain === nothing && return nothing
    length(domain) == 4 || error("Scene domain must be a 4-tuple `(xmin, ymin, xmax, ymax)`.")
    vals = ntuple(i -> Float64(domain[i]), 4)
    vals[3] > vals[1] || error("Scene domain must satisfy xmax > xmin.")
    vals[4] > vals[2] || error("Scene domain must satisfy ymax > ymin.")
    return vals
end

function _set_scene_dimensions!(mtg, bounds::NTuple{4,Float64})
    xmin, ymin, xmax, ymax = bounds
    mtg[:scene_dimensions] = (
        GeometryBasics.Point3f(Float32(xmin), Float32(ymin), 0.0f0),
        GeometryBasics.Point3f(Float32(xmax), Float32(ymax), 0.0f0),
    )
    return mtg
end

function _scene_builder_root(domain::Union{Nothing,NTuple{4,Float64}})
    root = MultiScaleTreeGraph.Node(
        MultiScaleTreeGraph.MutableNodeMTG(:/, :Scene, 1, 0),
        Dict{Symbol,Any}(),
    )
    domain === nothing || _set_scene_dimensions!(root, domain)
    return root
end

function _builder_refresh!(builder::LightSceneBuilder)
    _refresh_ref_mesh_registry!(builder.mtg)
    builder.domain === nothing || _set_scene_dimensions!(builder.mtg, builder.domain)
    return prepare_scene(
        builder.mtg;
        source_path=builder.source_path,
        scene_xy_bounds=builder.domain,
        relabel_ids=true,
    )
end

function _as_tuple3(x, name::AbstractString)
    length(x) == 3 || error("`$name` must be a 3-tuple.")
    ntuple(i -> Float64(x[i]), 3)
end

function _read_scene_object(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext == ".opf" && return _read_opf_relaxed(path)
    ext == ".gwa" && return PlantGeom.read_gwa(path)
    error("add_plant! accepts `.opf` and `.gwa` object paths; use read_scene for `.ops` scenes.")
end

function _annotate_scene_object!(mtg; group::AbstractString, id::Integer, file_path::AbstractString="")
    mtg[:functional_group] = String(group)
    mtg[:object_id] = Int(id)
    mtg[:plantID] = Int(id)
    isempty(file_path) || (mtg[:filePath] = String(file_path))
    MultiScaleTreeGraph.traverse!(mtg) do node
        haskey(node, :functional_group) || (node[:functional_group] = String(group))
        haskey(node, :object_id) || (node[:object_id] = Int(id))
        return true
    end
    return mtg
end

function _translate_scene_object!(mtg, at)
    x, y, z = _as_tuple3(at, "at")
    (x == 0.0 && y == 0.0 && z == 0.0) && return mtg
    translation = PlantGeom.Translation(x, y, z)
    MultiScaleTreeGraph.traverse!(mtg, filter_fun=PlantGeom.has_geometry) do node
        PlantGeom.transform_mesh!(node, translation)
        return true
    end
    return mtg
end

function _add_scene_object!(builder::LightSceneBuilder, obj; group::AbstractString, id::Integer, at=(0.0, 0.0, 0.0), file_path::AbstractString="")
    plant = deepcopy(obj)
    _annotate_scene_object!(plant; group=group, id=id, file_path=file_path)
    _translate_scene_object!(plant, at)
    try
        MultiScaleTreeGraph.addchild!(builder.mtg, plant)
    catch err
        msg = "Could not add plant object to the scene builder. Use MTGs read with PlantGeom defaults (`MutableNodeMTG`, Dict attributes) or pass an `.opf`/`.gwa` path. Original error: $(sprint(showerror, err))"
        throw(ArgumentError(msg))
    end
    return builder
end

"""
    add_plant!(builder, plant; group, id, at=(0, 0, 0))
    add_plant!(builder, path; group, id, at=(0, 0, 0))

Add an in-memory MTG or an `.opf`/`.gwa` object file to a scene builder.
"""
function add_plant!(builder::LightSceneBuilder, plant; group::AbstractString, id::Integer, at=(0.0, 0.0, 0.0))
    _add_scene_object!(builder, plant; group=group, id=id, at=at)
end

function add_plant!(builder::LightSceneBuilder, path::AbstractString; group::AbstractString, id::Integer, at=(0.0, 0.0, 0.0))
    _add_scene_object!(builder, _read_scene_object(path); group=group, id=id, at=at, file_path=path)
end

"""
    light_scene(f; domain, source_path="interactive.scene")

Build a scene interactively and return a prepared [`SceneGeometry`](@ref).
"""
function light_scene(f::Function; domain, source_path::AbstractString="interactive.scene")
    bounds = _coerce_domain(domain)
    builder = LightSceneBuilder(_scene_builder_root(bounds), bounds, String(source_path))
    f(builder)
    return _builder_refresh!(builder)
end

function light_scene(; domain, source_path::AbstractString="interactive.scene")
    bounds = _coerce_domain(domain)
    return _builder_refresh!(LightSceneBuilder(_scene_builder_root(bounds), bounds, String(source_path)))
end

function _scene_xy_bounds_from_mtg(mtg)
    haskey(mtg, :scene_dimensions) || return nothing
    dims = mtg[:scene_dimensions]
    dims === nothing && return nothing
    if dims isa AbstractString
        matches = collect(eachmatch(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", dims))
        length(matches) >= 5 || return nothing
        values = parse.(Float64, getfield.(matches, :match))
        return (
            min(values[1], values[4]),
            min(values[2], values[5]),
            max(values[1], values[4]),
            max(values[2], values[5]),
        )
    end
    p0 = dims[1]
    p1 = dims[2]
    (
        min(Float64(p0[1]), Float64(p1[1])),
        min(Float64(p0[2]), Float64(p1[2])),
        max(Float64(p0[1]), Float64(p1[1])),
        max(Float64(p0[2]), Float64(p1[2])),
    )
end

"""
    prepare_scene(mtg; source_path="interactive.opf", scene_xy_bounds=nothing, relabel_ids=false)::SceneGeometry

Convert an MTG with geometry into the dense scene representation used by the
light solver.

The returned [`SceneGeometry`](@ref) stores the merged mesh, face-to-node map,
per-node metadata, and optional xy bounds used for paving and rasterization.
"""
function prepare_scene(mtg; source_path::AbstractString="interactive.opf", scene_xy_bounds=nothing, relabel_ids::Bool=false)
    relabel_ids && _relabel_scene_node_ids!(mtg)
    bounds = scene_xy_bounds === nothing ? _scene_xy_bounds_from_mtg(mtg) : scene_xy_bounds
    _build_scene_geometry(mtg, source_path, bounds)
end

function _refresh_scene!(scene::SceneGeometry)
    scene.mtg === nothing && error("Scene refresh requires an MTG-backed scene.")
    PlantGeom.bump_scene_version!(scene.mtg)
    refreshed = prepare_scene(scene.mtg; source_path=scene.source_path, scene_xy_bounds=scene.scene_xy_bounds, relabel_ids=false)
    scene.merged_mesh = refreshed.merged_mesh
    scene.face2node = refreshed.face2node
    scene.nodes = refreshed.nodes
    scene.scene_xy_bounds = refreshed.scene_xy_bounds
    return scene
end

"""
    read_scene(path)::SceneGeometry

Read a scene file (`.ops`, `.opf`, or `.gwa`) and return a prepared
[`SceneGeometry`](@ref).

The scene is relabelled into a dense node-id space and immediately converted to
the merged-mesh representation expected by the interception pipeline.
"""
function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    mtg = if ext == ".ops"
        _read_ops_relaxed(path)
    elseif ext == ".opf"
        _read_opf_relaxed(path)
    elseif ext == ".gwa"
        PlantGeom.read_gwa(path)
    else
        error("Unsupported scene extension: $ext")
    end
    prepare_scene(mtg; source_path=String(path), scene_xy_bounds=_scene_xy_bounds_from_mtg(mtg), relabel_ids=true)
end

function _register_ref_mesh!(mtg, ref_mesh)
    if !haskey(mtg, :ref_meshes) || !(mtg[:ref_meshes] isa AbstractDict)
        mtg[:ref_meshes] = Dict{Int,Any}()
    end
    ref_meshes = mtg[:ref_meshes]
    next_id = isempty(ref_meshes) ? 0 : (maximum(Int.(keys(ref_meshes))) + 1)
    ref_meshes[next_id] = ref_mesh
    return next_id
end

function _refresh_ref_mesh_registry!(mtg)
    ref_meshes = PlantGeom.get_ref_meshes(mtg)
    mtg[:ref_meshes] = Dict(i - 1 => mesh_ for (i, mesh_) in enumerate(ref_meshes))
    return mtg[:ref_meshes]
end

function _drop_scene_surface_geometry!(mtg)
    haskey(mtg, :geometry) || return mtg
    mtg[:geometry] = nothing
    return mtg
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
    add_ground!(scene; z=0.0, nx=9, ny=9, xy_bounds=nothing, group="pavement", type="Cobblestone")

Add an explicit rectangular ground paving to an MTG-backed [`SceneGeometry`](@ref).

The paving is discretized into `nx * ny` tiles covering `xy_bounds` (or the
scene xy bounds when omitted), inserted into the scene MTG, and the prepared
scene caches are refreshed.
"""
function add_ground!(
    scene::SceneGeometry{T};
    z::Real=0.0,
    nx::Int=9,
    ny::Int=9,
    xy_bounds=nothing,
    group::AbstractString="pavement",
    type::AbstractString="Cobblestone",
) where {T<:MultiScaleTreeGraph.Node{N}} where N
    scene.mtg === nothing && error("add_ground! requires an MTG-backed scene.")
    bounds = xy_bounds === nothing ? scene.scene_xy_bounds : xy_bounds
    bounds === nothing && error("Ground bounds are undefined. Pass `xy_bounds=` or use a scene with known bounds.")
    xmin, ymin, xmax, ymax = bounds
    x_edges = collect(range(Float64(xmin), Float64(xmax), length=nx + 1))
    y_edges = collect(range(Float64(ymin), Float64(ymax), length=ny + 1))
    root_scale = MultiScaleTreeGraph.scale(scene.mtg)
    next_id = isempty(scene.nodes) ? 1 : (maximum(keys(scene.nodes)) + 1)

    for ix in 1:nx, iy in 1:ny
        x0 = x_edges[ix]
        x1 = x_edges[ix+1]
        y0 = y_edges[iy]
        y1 = y_edges[iy+1]
        points = GeometryBasics.Point3f[
            GeometryBasics.Point3f(Float32(x0), Float32(y0), Float32(z)),
            GeometryBasics.Point3f(Float32(x1), Float32(y0), Float32(z)),
            GeometryBasics.Point3f(Float32(x1), Float32(y1), Float32(z)),
            GeometryBasics.Point3f(Float32(x0), Float32(y1), Float32(z)),
        ]
        faces = PlantGeom.Face3[PlantGeom.Face3(1, 2, 3), PlantGeom.Face3(1, 3, 4)]
        mesh = GeometryBasics.Mesh(points, faces)
        nid = next_id
        next_id += 1
        ref_mesh = PlantGeom.RefMesh("$(group)_$(nid)", mesh)
        MultiScaleTreeGraph.Node(
            nid,
            scene.mtg,
            N(:+, Symbol(type), nid, root_scale + 1),
            Dict{Symbol,Any}(
                :geometry => PlantGeom.Geometry(ref_mesh=ref_mesh),
                :functional_group => String(group),
                :type => String(type),
                :object_id => -1,
                :source_topology_id => nid,
            ),
        )
    end

    scene.scene_xy_bounds = (Float64(xmin), Float64(ymin), Float64(xmax), Float64(ymax))
    _drop_scene_surface_geometry!(scene.mtg)
    _refresh_ref_mesh_registry!(scene.mtg)
    _refresh_scene!(scene)
end

"""
    add_ground!(builder; z=0.0, nx=9, ny=9, group="pavement", type="Cobblestone")

Add explicit ground geometry to a [`LightSceneBuilder`](@ref).
"""
function add_ground!(
    builder::LightSceneBuilder;
    z::Real=0.0,
    nx::Int=9,
    ny::Int=9,
    xy_bounds=nothing,
    group::AbstractString="pavement",
    type::AbstractString="Cobblestone",
)
    bounds = xy_bounds === nothing ? builder.domain : _coerce_domain(xy_bounds)
    bounds === nothing && error("Ground bounds are undefined. Pass `domain=` to light_scene or `xy_bounds=` to add_ground!.")
    scene = _builder_refresh!(builder)
    add_ground!(scene; z=z, nx=nx, ny=ny, xy_bounds=bounds, group=group, type=type)
    builder.mtg = scene.mtg
    builder.domain = scene.scene_xy_bounds
    return builder
end

"""
    write_scene(path, scene)

Write an MTG-backed [`SceneGeometry`](@ref) to `path`.

Supported output formats are `.ops`, `.opf`, and `.gwa`. The function refreshes
the reference-mesh registry and normalizes topology ids before export.
"""
function write_scene(path::AbstractString, scene::SceneGeometry)
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

function _as_plantmeteo_table(data)
    transformed = Tables.columntable(data)
    column_names = Tables.columnnames(data)
    if !(:date in column_names) && :DateTime in column_names
        transformed = PlantMeteo.set_column(transformed, :date, transformed.DateTime)
    end
    transformed = _normalize_raw_meteo_dates(transformed)
    if !(:duration in column_names)
        duration = PlantMeteo.compute_duration(transformed, Dates.DateFormat("HH:MM:SS"), nothing)
        transformed = PlantMeteo.set_column(transformed, :duration, duration)
    end
    PlantMeteo.TimeStepTable(transformed, _table_metadata_namedtuple(data))
end

"""
    read_meteo(path)::MeteoTable

Read a meteorological forcing table from `path` and return it as a
[`MeteoTable`](@ref).

The resulting table stores rows as named tuples and keeps available metadata
such as latitude, longitude, altitude, and source file path.
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

function read_meteo(data)
    Tables.istable(typeof(data)) || error("Unsupported meteo input: expected a path or a Tables.jl-compatible table.")
    data isa MeteoTable && return data
    data isa PlantMeteo.TimeStepTable && return MeteoTable(_meteo_rows(data), _meteo_metadata(data))

    meteo = _as_plantmeteo_table(data)
    MeteoTable(_meteo_rows(meteo), _meteo_metadata(meteo))
end
