import CSV

const _COMPONENT_VARIABLE_ORDER = [
    "step_number",
    "step_duration",
    "item_id",
    "component_id",
    "group",
    "type",
    "area",
    "surface_hits",
    "barycentre_x",
    "barycentre_y",
    "barycentre_z",
    "sky_fraction",
    "Ri_PAR_0_f",
    "Ri_PAR_0_q",
    "Ri_NIR_0_f",
    "Ri_NIR_0_q",
    "Ra_PAR_0_f",
    "Ra_PAR_0_q",
    "Ra_NIR_0_f",
    "Ra_NIR_0_q",
    "Ri_PAR_f",
    "Ri_PAR_q",
    "Ri_NIR_f",
    "Ri_NIR_q",
    "Ra_PAR_f",
    "Ra_PAR_q",
    "Ra_NIR_f",
    "Ra_NIR_q",
    "Ri_TIR_f",
    "Ri_TIR_q",
    "Ra_TIR_f",
    "Ra_TIR_q",
]

const _SCENE_VARIABLE_ORDER = [
    "step_number",
    "step_duration",
    "hour_start",
    "hour_end",
    "RI_SW_f",
    "plot_area",
]

const _SUN_POSITION_LOG_ORDER = [
    "stepNumber",
    "stepStart",
    "stepEnd",
    "azimuthHalf",
    "elevationHalf",
    "azimuthWeighted",
    "elevationWeighted",
]

const _SCATTERING_ITERATION_LOG_ORDER = [
    "step",
    "plantid",
    "nodeid",
    "iter",
    "scat",
]

const _SUMMARY_VARIABLE_ORDER = [
    "step_number",
    "step_duration",
    "group",
    "type",
    "item_id",
    "area",
    "Ri_q",
]

const _NODE_LINKS_STATS_ALLDIRS_ORDER = [
    "plantid",
    "nodeid",
    "links",
    "hits",
]

const _NODE_LINKS_DIR_ORDER = [
    "dir",
    "plantid1",
    "id1",
    "plantid2",
    "id2",
    "n",
]

const _DEFAULT_SIM_COUNTER_DIGITS = 6

const _SCATTERING_REQUIRED_COMPONENT_VARIABLES = Set([
    "Ri_PAR_f",
    "Ri_PAR_q",
    "Ri_NIR_f",
    "Ri_NIR_q",
    "Ra_PAR_f",
    "Ra_PAR_q",
    "Ra_NIR_f",
    "Ra_NIR_q",
])

const _PHOTOSYNTHESIS_REQUIRED_COMPONENT_VARIABLES = Set([
    "An_f",
    "An_q",
    "Gs",
])

const _ENERGY_BALANCE_REQUIRED_COMPONENT_VARIABLES = Set([
    "Ri_TIR_f",
    "Ri_TIR_q",
    "Ra_TIR_f",
    "Ra_TIR_q",
    "LE_f",
    "LE_q",
    "H_f",
    "H_q",
    "Rn_f",
    "Rn_q",
    "Tr_f",
    "Tr_q",
    "T",
])

const _OUT_OF_SCOPE_COMPONENT_VARIABLES = union(
    _PHOTOSYNTHESIS_REQUIRED_COMPONENT_VARIABLES,
    _ENERGY_BALANCE_REQUIRED_COMPONENT_VARIABLES,
)

function _canonical_component_variable_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_COMPONENT_VARIABLE_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_scene_variable_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_SCENE_VARIABLE_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_sun_position_log_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_SUN_POSITION_LOG_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_scattering_iteration_log_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_SCATTERING_ITERATION_LOG_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_summary_variable_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_SUMMARY_VARIABLE_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_node_links_stats_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_NODE_LINKS_STATS_ALLDIRS_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _canonical_node_links_dir_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_NODE_LINKS_DIR_ORDER))
    sort(
        unique(names);
        by=n -> (get(idx, n, typemax(Int)), n),
    )
end

function _default_light_component_variables(cfg::LightConfig)
    cols = [c for c in _COMPONENT_VARIABLE_ORDER if !(c in _OUT_OF_SCOPE_COMPONENT_VARIABLES)]
    if !cfg.general.scattering
        cols = [c for c in cols if !(c in _SCATTERING_REQUIRED_COMPONENT_VARIABLES)]
    end
    return cols
end

function _require_supported_output_configuration(cfg::LightConfig, cols::Vector{String})
    requested = Set(cols)
    bad_scattering = !cfg.general.scattering ? sort!(collect(intersect(requested, _SCATTERING_REQUIRED_COMPONENT_VARIABLES))) : String[]
    bad_scope = sort!(collect(intersect(requested, _OUT_OF_SCOPE_COMPONENT_VARIABLES)))

    isempty(bad_scattering) && isempty(bad_scope) && return nothing

    reasons = String[]
    !isempty(bad_scattering) && push!(reasons, "require scattering=true: $(join(bad_scattering, ", "))")
    if !isempty(bad_scope)
        push!(
            reasons,
            "out of scope in light-only mode (compute with PlantBiophysics): $(join(bad_scope, ", "))",
        )
    end
    error("Unsupported component output variables: $(join(reasons, " | "))")
end

"""
    component_variable_names(cfg)::Vector{String}

Resolve component output variable names from `cfg.outputs.component_variables` when available.
Only variables with truthy values are kept.
"""
function component_variable_names(cfg::LightConfig)
    d = cfg.outputs.component_variables
    isempty(d) && return _default_light_component_variables(cfg)

    vars = String[]
    for (k, v) in d
        _as_bool(v, false) || continue
        push!(vars, string(k))
    end
    cols = isempty(vars) ? _default_light_component_variables(cfg) : _canonical_component_variable_order(vars)
    _require_supported_output_configuration(cfg, cols)
    cols
end

"""
    scene_variable_names(cfg)::Vector{String}

Resolve scene output variable names from `cfg.outputs.scene_variables` when available.
Only variables with truthy values are kept. Defaults match Java `scene_values.csv`.
"""
function scene_variable_names(cfg::LightConfig)
    d = cfg.outputs.scene_variables
    isempty(d) && return copy(_SCENE_VARIABLE_ORDER)

    vars = String[]
    for (k, v) in d
        _as_bool(v, false) || continue
        push!(vars, string(k))
    end
    isempty(vars) ? copy(_SCENE_VARIABLE_ORDER) : _canonical_scene_variable_order(vars)
end

function _config_wants_component_outputs(cfg::LightConfig)
    d = cfg.outputs.component_variables
    for (_, v) in d
        _as_bool(v, false) && return true
    end
    return false
end

function _config_debug_enabled(cfg::LightConfig)
    cfg.general.log_debug || cfg.general.debug
end

function _export_ops_raw_value(cfg::LightConfig)
    cfg.outputs.export_ops
end

function _available_export_step_numbers(start_step_number::Int, nsteps::Int)
    nsteps <= 0 && return Int[]
    collect((start_step_number + 1):(start_step_number + nsteps))
end

function _ops_export_step_numbers(cfg::LightConfig, start_step_number::Int, nsteps::Int)
    raw_value = _export_ops_raw_value(cfg)
    available = _available_export_step_numbers(start_step_number, nsteps)
    isempty(available) && return Int[]

    raw_value === nothing && return Int[]
    raw_value isa Bool && return (raw_value ? [last(available)] : Int[])
    lower = lowercase(strip(string(raw_value)))
    if lower == "false"
        return Int[]
    elseif lower == "true"
        return [last(available)]
    elseif lower == "all"
        return available
    end

    n = try
        parse(Int, lower)
    catch
        error("invalid export_ops parameter: $(raw_value)")
    end
    n in available || error("invalid export_ops parameter: $(raw_value)")
    return [n]
end

function _ops_export_variables(cfg::LightConfig)
    d = cfg.outputs.opf_variables
    isempty(d) && return ["Ri_PAR_0_f", "Ri_NIR_0_f"]
    vars = String[]
    valid = Set(_COMPONENT_VARIABLE_ORDER)
    for (k, v) in d
        _as_bool(v, false) || continue
        ks = string(k)
        ks in valid || continue
        push!(vars, ks)
    end
    isempty(vars) ? ["Ri_PAR_0_f", "Ri_NIR_0_f"] : _canonical_component_variable_order(vars)
end

function _opf_overwrite_variables(cfg::LightConfig)
    cfg.outputs.opf_overwrite_variables
end

function _write_export_ops_outputs(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig,
    meteo_rows,
    outroot::AbstractString,
    start_step_number::Int,
)
    selected_step_numbers = _ops_export_step_numbers(cfg, start_step_number, length(steps))
    isempty(selected_step_numbers) && return Dict{String,String}()

    source_scene = scene.source_path
    ext = lowercase(splitext(source_scene)[2])
    ext == ".ops" || error("export_ops currently requires an OPS scene source, got: $(source_scene)")

    vars = _ops_export_variables(cfg)
    meta = _output_node_metadata(scene, cfg)
    overwrite = _opf_overwrite_variables(cfg)
    step_by_number = Dict((start_step_number + i) => i for i in eachindex(steps))
    out = Dict{String,String}()
    for step_number in selected_step_numbers
        idx = get(step_by_number, step_number, 0)
        idx > 0 || continue
        step = steps[idx]
        row = meteo_rows[idx]
        table = component_values_table(
            scene,
            step,
            cfg;
            meteo_row=row,
            step_number=step_number - 1,
            columns=vars,
            unavailable="NA",
        )
        value_maps = Dict{Symbol,Dict{Int,Any}}(
            Symbol(var) => Dict{Int,Any}(nid => table.rows[i][var] for (i, nid) in enumerate(meta.node_ids)) for var in vars
        )

        scene_copy = deepcopy(scene.mtg)
        scene_node_ids = Set(meta.node_ids)
        MultiScaleTreeGraph.traverse!(scene_copy) do node
            nid = Int(MultiScaleTreeGraph.node_id(node))
            if nid in scene_node_ids
                for (attr, values) in value_maps
                    if overwrite || !haskey(node, attr)
                        node[attr] = get(values, nid, "NA")
                    end
                end
            end
            return true
        end

        scene_stem = splitext(basename(source_scene))[1]
        step_tag = lpad(string(step_number), 5, '0')
        dst_ops = joinpath(outroot, "$(scene_stem)-step$(step_tag).ops")
        PlantGeom.write_ops(dst_ops, scene_copy; write_objects=true, preserve_file_paths=true)
        out["ops_step_$(step_tag)"] = dst_ops
    end
    return out
end

"""
    output_directory(cfg)::String

Resolve absolute output directory from config (`output_directory`, default `"output"`).
"""
function output_directory(cfg::LightConfig)
    out_rel = cfg.outputs.output_directory
    base = cfg.paths.base_dir
    isabspath(out_rel) ? normpath(out_rel) : normpath(joinpath(base, out_rel))
end

function _next_simulation_counter_dirname(base_dir::AbstractString; counter_digits::Int=_DEFAULT_SIM_COUNTER_DIGITS)
    max_counter = 0
    re = Regex("^\\d{$counter_digits}\$")
    if isdir(base_dir)
        for entry in readdir(base_dir)
            occursin(re, entry) || continue
            p = joinpath(base_dir, entry)
            isdir(p) || continue
            v = try
                parse(Int, entry)
            catch
                0
            end
            max_counter = max(max_counter, v)
        end
    end
    return lpad(string(max_counter + 1), counter_digits, '0')
end

"""
    simulation_output_directory(cfg; base_output=nothing, counter_digits=6, create=true, clean_existing=true)::String

Resolve Java-style simulation output directory:
- `<output_directory>/<simulation_directory>` when `simulation_directory` is set in config;
- otherwise `<output_directory>/<counter>` with zero-padded numeric counter (default 6 digits).

When `create=true`, directories are created. If `simulation_directory` is set and the target exists,
`clean_existing=true` removes it first (Java behavior).
"""
function simulation_output_directory(
    cfg::LightConfig;
    base_output::Union{Nothing,AbstractString}=nothing,
    counter_digits::Int=_DEFAULT_SIM_COUNTER_DIGITS,
    create::Bool=true,
    clean_existing::Bool=true,
)
    out_base = base_output === nothing ? output_directory(cfg) : normpath(String(base_output))
    create && mkpath(out_base)

    sim_name = cfg.outputs.simulation_directory
    if isempty(sim_name)
        sim_name = _next_simulation_counter_dirname(out_base; counter_digits=counter_digits)
    end

    out = normpath(joinpath(out_base, sim_name))
    if create && clean_existing && !isempty(cfg.outputs.simulation_directory) && isdir(out)
        rm(out; recursive=true, force=true)
    end
    create && mkpath(out)
    return out
end

function _step_duration_output_local(meteo_row, step_duration_seconds)
    if step_duration_seconds !== nothing
        return Float64(step_duration_seconds)
    end
    meteo_row === nothing && return 1.0
    return _step_duration_seconds_local(meteo_row)
end

function _node_ids_for_output(scene::SceneGeometry)
    ids = collect(keys(scene.total_area_per_node))
    sort!(
        ids;
        by=nid -> (
            get(scene.java_item_id_per_node, nid, 0),
            get(scene.java_component_id_per_node, nid, nid),
            nid,
        ),
    )
    ids
end

function _output_node_metadata(scene::SceneGeometry, cfg::LightConfig)
    vertices, faces, face2node, _, _, node_group_geom = _scene_geometry_for_interception(scene, cfg)
    area_per_node = Dict{Int,Float64}()
    bary_sum_per_node = Dict{Int,NTuple{3,Float64}}()
    for (i, f) in enumerate(faces)
        nid = face2node[i]
        p1 = vertices[f[1]]
        p2 = vertices[f[2]]
        p3 = vertices[f[3]]
        a = _triangle_area3d(p1, p2, p3)
        area_per_node[nid] = get(area_per_node, nid, 0.0) + a
        cx = (p1[1] + p2[1] + p3[1]) / 3.0
        cy = (p1[2] + p2[2] + p3[2]) / 3.0
        cz = (p1[3] + p2[3] + p3[3]) / 3.0
        sx, sy, sz = get(bary_sum_per_node, nid, (0.0, 0.0, 0.0))
        bary_sum_per_node[nid] = (sx + a * cx, sy + a * cy, sz + a * cz)
    end
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    for (nid, a) in area_per_node
        if a > 0.0
            sx, sy, sz = get(bary_sum_per_node, nid, (0.0, 0.0, 0.0))
            barycenter_per_node[nid] = (sx / a, sy / a, sz / a)
        else
            barycenter_per_node[nid] = (NaN, NaN, NaN)
        end
    end

    key_by_node = _interception_java_keys(scene, cfg)
    node_ids = collect(keys(key_by_node))
    sort!(
        node_ids;
        by=nid -> begin
            item_id, component_id = key_by_node[nid]
            (item_id, component_id, nid)
        end,
    )

    item_per_node = Dict{Int,Int}()
    component_per_node = Dict{Int,Int}()
    group_per_node = Dict{Int,String}()
    type_per_node = Dict{Int,String}()
    group_type_hints = _group_type_hints(cfg)
    for nid in node_ids
        item_id, component_id = key_by_node[nid]
        item_per_node[nid] = item_id
        component_per_node[nid] = component_id
        g = haskey(node_group_geom, nid) ? node_group_geom[nid] : get(scene.node_group, nid, "")
        group_per_node[nid] = g
        t = strip(get(scene.node_type, nid, ""))
        if isempty(t)
            if g == "pavement"
                t = "Cobblestone"
            elseif g == "sensors"
                t = "Sensor"
            end
        end
        if (isempty(t) || t == "Mesh") && !isempty(g)
            hints = get(group_type_hints, g, String[])
            length(hints) == 1 && (t = hints[1])
        end
        type_per_node[nid] = t
    end

    return (
        node_ids=node_ids,
        area_per_node=area_per_node,
        barycenter_per_node=barycenter_per_node,
        item_per_node=item_per_node,
        component_per_node=component_per_node,
        group_per_node=group_per_node,
        type_per_node=type_per_node,
    )
end

function _java_id_maps(scene::SceneGeometry, cfg::LightConfig)
    keys_by_node = _interception_java_keys(scene, cfg)
    item_per_node = Dict{Int,Int}()
    comp_per_node = Dict{Int,Int}()
    for (nid, key) in keys_by_node
        item_per_node[nid] = key[1]
        comp_per_node[nid] = key[2]
    end
    return item_per_node, comp_per_node
end

function _group_type_hints(cfg::LightConfig)
    out = Dict{String,Vector{String}}()
    for model in cfg.models
        group = haskey(model, "Group") ? strip(string(model["Group"])) : ""
        isempty(group) && continue
        types = get(model, "Type", nothing)
        types isa AbstractDict || continue
        names = String[]
        for k in keys(types)
            s = strip(string(k))
            isempty(s) || push!(names, s)
        end
        isempty(names) || (out[group] = sort(unique(names)))
    end
    out
end

function _sky_fraction_per_node(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    cfg::LightConfig,
    area_per_node::Dict{Int,Float64},
    node_ids::Vector{Int},
)
    pa_by_sector, _, _, _, _, _ = _build_sector_responses(scene, turtle, cfg)
    sky_count = 0
    visible_sum = Dict{Int,Float64}()
    for (i, sector) in enumerate(turtle.sectors)
        sector.source == :sun && continue
        sky_count += 1
        for (nid, a) in pa_by_sector[i]
            visible_sum[nid] = get(visible_sum, nid, 0.0) + a
        end
    end

    out = Dict{Int,Float64}()
    for nid in node_ids
        area = get(area_per_node, nid, 0.0)
        if area <= 0.0 || sky_count == 0
            out[nid] = 0.0
        else
            out[nid] = get(visible_sum, nid, 0.0) / sky_count / area
        end
    end
    out
end

function _ri_value(
    step::LightStepResult,
    nid::Int,
    band::String,
    interception_only::Bool,
    quantity::Bool,
    area::Float64,
    step_duration::Float64,
)
    b = uppercase(band)
    if b == "PAR"
        if interception_only
            return quantity ? get(step.budget.ri_par_0_q_per_node, nid, 0.0) : get(step.budget.ri_par_0_f_per_node, nid, 0.0)
        end
        return quantity ? get(step.budget.ri_par_q_per_node, nid, 0.0) : get(step.budget.ri_par_f_per_node, nid, 0.0)
    elseif b == "NIR"
        if interception_only
            return quantity ? get(step.budget.ri_nir_0_q_per_node, nid, 0.0) : get(step.budget.ri_nir_0_f_per_node, nid, 0.0)
        end
        return quantity ? get(step.budget.ri_nir_q_per_node, nid, 0.0) : get(step.budget.ri_nir_f_per_node, nid, 0.0)
    end

    qdict = interception_only ? get(step.budget.extra_0_q_per_band, b, Dict{Int,Float64}()) : get(step.budget.extra_q_per_band, b, Dict{Int,Float64}())
    q = get(qdict, nid, 0.0)
    if quantity
        return q
    end
    return q / max(area * step_duration, eps(Float64))
end

function _ra_value(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig,
    nid::Int,
    band::String,
    interception_only::Bool,
    quantity::Bool,
    area::Float64,
    step_duration::Float64,
    absorptance_cache::Dict{String,Dict{Int,Float64}},
)
    b = uppercase(band)
    if b == "PAR"
        if interception_only
            return quantity ? get(step.budget.ra_par_0_q_per_node, nid, 0.0) : get(step.budget.ra_par_0_f_per_node, nid, 0.0)
        end
        return quantity ? get(step.budget.ra_par_q_per_node, nid, 0.0) : get(step.budget.ra_par_f_per_node, nid, 0.0)
    elseif b == "NIR"
        if interception_only
            return quantity ? get(step.budget.ra_nir_0_q_per_node, nid, 0.0) : get(step.budget.ra_nir_0_f_per_node, nid, 0.0)
        end
        return quantity ? get(step.budget.ra_nir_q_per_node, nid, 0.0) : get(step.budget.ra_nir_f_per_node, nid, 0.0)
    end

    abs_band = get!(absorptance_cache, b) do
        _node_absorptance_per_band(scene, cfg, b)
    end
    ri = _ri_value(step, nid, b, interception_only, quantity, area, step_duration)
    return ri * get(abs_band, nid, clamp(1.0 - _default_scattering_factor_local(cfg, b), 0.0, 1.0))
end

function _component_variable_value(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig,
    nid::Int,
    variable::String,
    step_number::Int,
    step_duration::Float64,
    area_per_node::Dict{Int,Float64},
    barycenter_per_node::Dict{Int,NTuple{3,Float64}},
    item_per_node::Dict{Int,Int},
    component_per_node::Dict{Int,Int},
    group_per_node::Dict{Int,String},
    type_per_node::Dict{Int,String},
    absorptance_cache::Dict{String,Dict{Int,Float64}},
    sky_fraction_per_node::Union{Nothing,Dict{Int,Float64}},
    group_type_hints::Dict{String,Vector{String}};
    unavailable::String="NA",
    strict::Bool=false,
)
    area = get(area_per_node, nid, 0.0)
    hits = get(step.first_order.hits_per_node, nid, 0)

    if variable == "step_number"
        return step_number
    elseif variable == "step_duration"
        return step_duration
    elseif variable == "item_id"
        return get(item_per_node, nid, -1)
    elseif variable == "component_id"
        return get(component_per_node, nid, nid)
    elseif variable == "group"
        return get(group_per_node, nid, "")
    elseif variable == "type"
        t = strip(get(type_per_node, nid, ""))
        return if isempty(t)
            g = get(group_per_node, nid, "")
            ts = get(group_type_hints, g, String[])
            length(ts) == 1 ? ts[1] : unavailable
        else
            t
        end
    elseif variable == "area"
        return area
    elseif variable == "surface_hits"
        return area > 0 ? hits / area : 0.0
    elseif variable == "barycentre_x"
        return get(barycenter_per_node, nid, (NaN, NaN, NaN))[1]
    elseif variable == "barycentre_y"
        return get(barycenter_per_node, nid, (NaN, NaN, NaN))[2]
    elseif variable == "barycentre_z"
        return get(barycenter_per_node, nid, (NaN, NaN, NaN))[3]
    elseif variable == "sky_fraction"
        return sky_fraction_per_node === nothing ? unavailable : get(sky_fraction_per_node, nid, 0.0)
    elseif variable == "Ri_TIR_f" || variable == "Ri_TIR_q" || variable == "Ra_TIR_f" || variable == "Ra_TIR_q"
        return unavailable
    end

    m = match(r"^R([ia])_([^_]+)_(0_)?([qf])$", variable)
    if m !== nothing
        kind = String(m.captures[1])
        band = String(m.captures[2])
        interception_only = m.captures[3] !== nothing
        quantity = String(m.captures[4]) == "q"
        if kind == "i"
            return _ri_value(step, nid, band, interception_only, quantity, area, step_duration)
        elseif kind == "a"
            return _ra_value(
                scene,
                step,
                cfg,
                nid,
                band,
                interception_only,
                quantity,
                area,
                step_duration,
                absorptance_cache,
            )
        end
    end

    ms = match(r"^scat_factor_([^_]+)$", variable)
    if ms !== nothing
        band = uppercase(String(ms.captures[1]))
        abs_band = get!(absorptance_cache, band) do
            _node_absorptance_per_band(scene, cfg, band)
        end
        return 1.0 - get(abs_band, nid, clamp(1.0 - _default_scattering_factor_local(cfg, band), 0.0, 1.0))
    end

    if strict
        error("Unsupported component variable: $variable")
    end
    return unavailable
end

"""
    component_values_table(scene, step, cfg; meteo_row=nothing, step_number=0, step_duration_seconds=nothing, columns=nothing, unavailable="NA", strict=false)

Build a Java-like component values table represented as `(columns, rows)`.
`rows` is a vector of dictionaries keyed by column name.
"""
function component_values_table(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=0,
    step_duration_seconds=nothing,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
    strict::Bool=false,
)
    cols = columns === nothing ? component_variable_names(cfg) : String.(columns)
    _require_supported_output_configuration(cfg, cols)
    step_duration = _step_duration_output_local(meteo_row, step_duration_seconds)
    meta = _output_node_metadata(scene, cfg)
    node_ids = meta.node_ids
    area_per_node = meta.area_per_node
    barycenter_per_node = meta.barycenter_per_node
    item_per_node = meta.item_per_node
    component_per_node = meta.component_per_node
    group_per_node = meta.group_per_node
    type_per_node = meta.type_per_node
    rows = Vector{Dict{String,Any}}(undef, length(node_ids))
    absorptance_cache = Dict{String,Dict{Int,Float64}}()
    group_type_hints = _group_type_hints(cfg)
    sky_fraction_per_node =
        ("sky_fraction" in cols) ? _sky_fraction_per_node(scene, step.turtle, cfg, area_per_node, node_ids) : nothing

    for (i, nid) in enumerate(node_ids)
        row = Dict{String,Any}()
        for var in cols
            row[var] = _component_variable_value(
                scene,
                step,
                cfg,
                nid,
                var,
                step_number,
                step_duration,
                area_per_node,
                barycenter_per_node,
                item_per_node,
                component_per_node,
                group_per_node,
                type_per_node,
                absorptance_cache,
                sky_fraction_per_node,
                group_type_hints;
                unavailable=unavailable,
                strict=strict,
            )
        end
        rows[i] = row
    end

    return (columns=cols, rows=rows)
end

function _csv_cell_local(v)
    if v === missing || v === nothing
        return "NA"
    elseif v isa AbstractString
        return replace(v, ';' => ",")
    end
    return string(v)
end

function _write_csv_rows(path::AbstractString, columns::Vector{String}, rows::Vector{<:AbstractDict}; append::Bool=false)
    mkpath(dirname(path))
    mode = append && isfile(path) && filesize(path) > 0 ? "a" : "w"
    open(path, mode) do io
        if mode == "w"
            println(io, join(columns, ';'))
        end
        for row in rows
            vals = [_csv_cell_local(get(row, c, "NA")) for c in columns]
            println(io, join(vals, ';'))
        end
    end
    return path
end

function _parse_sort_columns(sort_columns)
    if sort_columns isa AbstractString
        parts = split(String(sort_columns), r"[;,]")
        cols = String[strip(p) for p in parts if !isempty(strip(p))]
        isempty(cols) && error("sort_csv_file: no sort columns provided.")
        return cols
    elseif sort_columns isa AbstractVector
        cols = String[strip(string(c)) for c in sort_columns if !isempty(strip(string(c)))]
        isempty(cols) && error("sort_csv_file: no sort columns provided.")
        return cols
    else
        error("sort_csv_file: `sort_columns` must be a string or vector of column names.")
    end
end

function _csv_sort_key_value(v)
    if v === missing || v === nothing
        return (2, 0.0, "")
    elseif v isa Number
        return (0, Float64(v), "")
    end
    s = strip(string(v))
    n = tryparse(Float64, s)
    if n === nothing
        return (1, 0.0, lowercase(s))
    end
    return (0, n, "")
end

"""
    sort_csv_file(input_path, output_path; sort_columns)::String

Sort a `;`-separated CSV file by one or more columns and write the sorted file to `output_path`.
`sort_columns` accepts a vector of column names, or a single string with names separated by `;` or `,`.
"""
function sort_csv_file(
    input_path::AbstractString,
    output_path::AbstractString;
    sort_columns,
)
    f = CSV.File(input_path; delim=';', comment="#", normalizenames=false)
    cols = String.(propertynames(f))
    rows = collect(Tables.rowtable(f))

    sort_cols = _parse_sort_columns(sort_columns)
    colset = Set(cols)
    missing_cols = String[c for c in sort_cols if !(c in colset)]
    isempty(missing_cols) || error("sort_csv_file: unknown sort column(s): $(join(missing_cols, ", "))")

    sort_syms = Symbol.(sort_cols)
    sort!(
        rows;
        by=row -> Tuple(_csv_sort_key_value(getproperty(row, s)) for s in sort_syms),
    )

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, join(cols, ';'))
        for row in rows
            vals = [_csv_cell_local(getproperty(row, Symbol(c))) for c in cols]
            println(io, join(vals, ';'))
        end
    end
    return output_path
end

function _row_field_string_local(row, sym::Symbol, default::String="NA")
    row === nothing && return default
    sym in propertynames(row) || return default
    v = getproperty(row, sym)
    v === missing && return default
    s = strip(string(v))
    isempty(s) ? default : s
end

function _scene_plot_area(scene::SceneGeometry, cfg::LightConfig)
    _, _, _, _, plotbox, _ = _scene_geometry_for_interception(scene, cfg)
    return plotbox.xdim * plotbox.ydim
end

function _scene_variable_value(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig,
    variable::String,
    step_number::Int,
    meteo_row,
    plot_area::Float64,
    step_duration_seconds::Float64;
    unavailable::String="NA",
    strict::Bool=false,
)
    if variable == "step_number"
        return step_number
    elseif variable == "step_duration"
        return step_duration_seconds
    elseif variable == "hour_start"
        return _row_field_string_local(meteo_row, :hour_start, unavailable)
    elseif variable == "hour_end"
        return _row_field_string_local(meteo_row, :hour_end, unavailable)
    elseif variable == "date"
        return _row_field_string_local(meteo_row, :date, unavailable)
    elseif variable == "RI_SW_f"
        return step.sky.ri_sw_f
    elseif variable == "plot_area"
        return plot_area
    end

    if strict
        error("Unsupported scene variable: $variable")
    end
    return unavailable
end

function _decimal_hour_to_hm_local(x::Float64)
    h = floor(Int, x)
    m = round(Int, (x - h) * 60)
    if m >= 60
        h += 1
        m -= 60
    end
    return h, m
end

function _java_like_step_datetime_strings(row)
    date = _row_date(row)
    start_h, end_h = _row_step_hours(row)
    start_day = floor(Int, start_h / 24.0)
    end_day = floor(Int, end_h / 24.0)
    sh, sm = _decimal_hour_to_hm_local(start_h - 24.0 * start_day)
    eh, em = _decimal_hour_to_hm_local(end_h - 24.0 * end_day)
    d0 = date + Dates.Day(start_day)
    d1 = date + Dates.Day(end_day)
    s0 = "$(Dates.year(d0))/$(Dates.month(d0))/$(Dates.day(d0))  $(sh):$(lpad(string(sm), 2, '0'))"
    s1 = "$(Dates.year(d1))/$(Dates.month(d1))/$(Dates.day(d1))  $(eh):$(lpad(string(em), 2, '0'))"
    return s0, s1
end

function _sun_position_log_values(row, sky::SkyState)
    date = _row_date(row)
    start_h, end_h = _row_step_hours(row)
    lat_deg = _row_value(row, [:latitude, :lat], 0.0)
    az_half_deg, el_half_deg = _midpoint_sun_position_deg(date, start_h, end_h, deg2rad(lat_deg))
    step_start, step_end = _java_like_step_datetime_strings(row)
    return (
        stepStart=step_start,
        stepEnd=step_end,
        azimuthHalf=deg2rad(az_half_deg),
        elevationHalf=deg2rad(el_half_deg),
        azimuthWeighted=deg2rad(sky.sun_azimuth_deg),
        elevationWeighted=deg2rad(sky.sun_elevation_deg),
    )
end

"""
    scene_values_table(scene, steps, cfg; meteo_rows=nothing, start_step_number=0, columns=nothing, unavailable="NA", strict=false)

Build a Java-like scene values table represented as `(columns, rows)`.
`rows` is a vector of dictionaries keyed by column name.
"""
function scene_values_table(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
    strict::Bool=false,
)
    n = length(steps)
    rows_in = meteo_rows === nothing ? fill(nothing, n) : collect(meteo_rows)
    length(rows_in) == n || error("scene_values_table: `meteo_rows` length must match `steps` length.")
    cols = columns === nothing ? scene_variable_names(cfg) : String.(columns)
    rows = Vector{Dict{String,Any}}(undef, n)
    plot_area = _scene_plot_area(scene, cfg)

    for i in 1:n
        step = steps[i]
        row = rows_in[i]
        step_no = start_step_number + i - 1
        dt = _step_duration_output_local(row, nothing)
        out = Dict{String,Any}()
        for var in cols
            out[var] = _scene_variable_value(
                scene,
                step,
                cfg,
                var,
                step_no,
                row,
                plot_area,
                dt;
                unavailable=unavailable,
                strict=strict,
            )
        end
        rows[i] = out
    end
    return (columns=cols, rows=rows)
end

"""
    write_scene_values_csv(path, scene, steps, cfg; meteo_rows=nothing, start_step_number=0, columns=nothing, unavailable="NA", strict=false)

Write Java-like `scene_values.csv` output for one or many simulation steps.
"""
function write_scene_values_csv(
    path::AbstractString,
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
    strict::Bool=false,
)
    table = scene_values_table(
        scene,
        steps,
        cfg;
        meteo_rows=meteo_rows,
        start_step_number=start_step_number,
        columns=columns,
        unavailable=unavailable,
        strict=strict,
    )
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(table.columns, ';'))
        for row in table.rows
            vals = [_csv_cell_local(get(row, c, unavailable)) for c in table.columns]
            println(io, join(vals, ';'))
        end
    end
    return path
end

"""
    sun_position_log_table(steps, meteo_rows; start_step_number=0, columns=nothing)

Build a Java-like `log-sun-position.csv` table represented as `(columns, rows)`.
Angles are in radians, matching Java logs.
"""
function sun_position_log_table(
    steps::AbstractVector{<:LightStepResult},
    meteo_rows;
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    rows_in = collect(meteo_rows)
    length(rows_in) == length(steps) || error("sun_position_log_table: `meteo_rows` length must match `steps` length.")
    cols = columns === nothing ? copy(_SUN_POSITION_LOG_ORDER) : _canonical_sun_position_log_order(String.(columns))
    rows = Vector{Dict{String,Any}}(undef, length(steps))
    for i in eachindex(steps)
        vals = _sun_position_log_values(rows_in[i], steps[i].sky)
        row = Dict{String,Any}()
        for c in cols
            if c == "stepNumber"
                row[c] = start_step_number + i - 1
            elseif c == "stepStart"
                row[c] = vals.stepStart
            elseif c == "stepEnd"
                row[c] = vals.stepEnd
            elseif c == "azimuthHalf"
                row[c] = vals.azimuthHalf
            elseif c == "elevationHalf"
                row[c] = vals.elevationHalf
            elseif c == "azimuthWeighted"
                row[c] = vals.azimuthWeighted
            elseif c == "elevationWeighted"
                row[c] = vals.elevationWeighted
            else
                row[c] = "NA"
            end
        end
        rows[i] = row
    end
    return (columns=cols, rows=rows)
end

"""
    write_sun_position_log_csv(path, steps, meteo_rows; start_step_number=0, columns=nothing)

Write Java-like `log-sun-position.csv`.
"""
function write_sun_position_log_csv(
    path::AbstractString,
    steps::AbstractVector{<:LightStepResult},
    meteo_rows;
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    table = sun_position_log_table(steps, meteo_rows; start_step_number=start_step_number, columns=columns)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(table.columns, ';'))
        for row in table.rows
            vals = [_csv_cell_local(get(row, c, "NA")) for c in table.columns]
            println(io, join(vals, ';'))
        end
    end
    return path
end

function _scattering_iteration_history_one_band(
    graph::ScatteringTransferGraph,
    initial_power_per_node::Dict{Int,Float64},
    cfg::LightConfig,
    band_key::String,
    default_coeff::Float64,
)
    node_ids = graph.node_ids
    pair_counts = graph.pair_counts
    all_hits = graph.all_hits
    coeff_by_node = _coeff_by_node(graph, band_key, default_coeff)
    current = Dict{Int,Float64}(nid => get(initial_power_per_node, nid, 0.0) for nid in node_ids)
    ref = _sum_dict_values(current)
    thr = cfg.general.scattering_stop_ratio * max(ref, eps(Float64))
    per_iter = Vector{Dict{Int,Float64}}()
    iterations = 0
    converged = false

    for it in 1:cfg.general.scattering_max_iter
        iterations = it
        hit_energy = _dict_zero(node_ids)
        for nid in node_ids
            nh = get(all_hits, nid, 0)
            if nh > 0
                hit_energy[nid] = get(current, nid, 0.0) * get(coeff_by_node, nid, default_coeff) / nh / 2.0
            end
        end

        next = _dict_zero(node_ids)
        for ((to, from), cnt) in pair_counts
            next[to] = get(next, to, 0.0) + cnt * get(hit_energy, from, 0.0)
        end
        push!(per_iter, next)

        total_next = _sum_dict_values(next)
        current = next
        if total_next <= thr
            converged = true
            break
        end
    end
    return per_iter, iterations, converged
end

"""
    scattering_iteration_log_table(scene, step, cfg; meteo_row=nothing, step_number=0, band="PAR", mode=:raycast, backend=nothing, columns=nothing)

Build a Java-like `log-iteration-scat-<band>.csv` table represented as `(columns, rows)`.
"""
function scattering_iteration_log_table(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=0,
    band::AbstractString="PAR",
    mode::Symbol=:raycast,
    backend::Union{Nothing,ScatteringBackend}=nothing,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    cols = columns === nothing ? copy(_SCATTERING_ITERATION_LOG_ORDER) : _canonical_scattering_iteration_log_order(String.(columns))
    b = uppercase(String(band))
    graph = build_scattering_transfer_graph(scene, step.turtle, step.first_order, cfg; mode=mode, backend=backend)
    initial =
        if b == "NIR"
            Dict{Int,Float64}(nid => get(step.first_order.incident_nir_power_per_node, nid, 0.0) for nid in graph.node_ids)
        else
            Dict{Int,Float64}(nid => get(step.first_order.incident_par_power_per_node, nid, 0.0) for nid in graph.node_ids)
        end
    dflt = _default_band_coeff(cfg, b)
    per_iter, _, _ = _scattering_iteration_history_one_band(graph, initial, cfg, b, dflt)
    dt = _step_duration_output_local(meteo_row, nothing)
    w_to_mj = dt / 1e6
    key_by_node = _interception_java_keys(scene, cfg)
    node_ids = collect(keys(key_by_node))
    sort!(
        node_ids;
        by=nid -> begin
            item_id, component_id = key_by_node[nid]
            (item_id, component_id, nid)
        end,
    )
    rows = Dict{String,Any}[]
    for (it, vals) in enumerate(per_iter)
        iter_idx = it - 1
        for nid in node_ids
            row = Dict{String,Any}()
            item_id, component_id = key_by_node[nid]
            for c in cols
                if c == "step"
                    row[c] = step_number
                elseif c == "plantid"
                    row[c] = item_id
                elseif c == "nodeid"
                    row[c] = component_id
                elseif c == "iter"
                    row[c] = iter_idx
                elseif c == "scat"
                    row[c] = get(vals, nid, 0.0) * w_to_mj
                else
                    row[c] = "NA"
                end
            end
            push!(rows, row)
        end
    end
    return (columns=cols, rows=rows)
end

"""
    write_scattering_iteration_log_csv(path, scene, step, cfg; meteo_row=nothing, step_number=0, band="PAR", mode=:raycast, backend=nothing, columns=nothing)

Write Java-like `log-iteration-scat-<band>.csv`.
"""
function write_scattering_iteration_log_csv(
    path::AbstractString,
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=0,
    band::AbstractString="PAR",
    mode::Symbol=:raycast,
    backend::Union{Nothing,ScatteringBackend}=nothing,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    table = scattering_iteration_log_table(
        scene,
        step,
        cfg;
        meteo_row=meteo_row,
        step_number=step_number,
        band=band,
        mode=mode,
        backend=backend,
        columns=columns,
    )
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(table.columns, ';'))
        for row in table.rows
            vals = [_csv_cell_local(get(row, c, "NA")) for c in table.columns]
            println(io, join(vals, ';'))
        end
    end
    return path
end

"""
    node_links_stats_alldirs_table(scene, turtle, cfg; columns=nothing)

Build Java-like `log-nodelinks-stats-alldirs.csv` table represented as `(columns, rows)`.
"""
function node_links_stats_alldirs_table(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    cfg::LightConfig;
    columns::Union{Nothing,AbstractVector}=nothing,
)
    cols = columns === nothing ? copy(_NODE_LINKS_STATS_ALLDIRS_ORDER) : _canonical_node_links_stats_order(String.(columns))
    pair_counts, _, _, _ = _pair_counts_for_scattering(scene, turtle, cfg)
    item_per_node, comp_per_node = _java_id_maps(scene, cfg)

    neighbours_by_node = Dict{Int,Set{Int}}()
    hits_by_node = Dict{Int,Int}()
    for ((to, from), c) in pair_counts
        c > 0 || continue
        if to == from
            s = get!(neighbours_by_node, to) do
                Set{Int}()
            end
            push!(s, from)
            hits_by_node[to] = get(hits_by_node, to, 0) + c
            continue
        end
        s = get!(neighbours_by_node, to) do
            Set{Int}()
        end
        push!(s, from)
        hits_by_node[to] = get(hits_by_node, to, 0) + c
    end

    node_ids = sort(collect(keys(neighbours_by_node)); by=nid -> (get(item_per_node, nid, -1), get(comp_per_node, nid, nid), nid))
    rows = Dict{String,Any}[]
    for nid in node_ids
        row = Dict{String,Any}()
        for c in cols
            if c == "plantid"
                row[c] = get(item_per_node, nid, -1)
            elseif c == "nodeid"
                row[c] = get(comp_per_node, nid, nid)
            elseif c == "links"
                row[c] = length(get(neighbours_by_node, nid, Set{Int}()))
            elseif c == "hits"
                row[c] = get(hits_by_node, nid, 0)
            else
                row[c] = "NA"
            end
        end
        push!(rows, row)
    end
    return (columns=cols, rows=rows)
end

"""
    write_node_links_stats_alldirs_csv(path, scene, turtle, cfg; columns=nothing)

Write Java-like `log-nodelinks-stats-alldirs.csv`.
"""
function write_node_links_stats_alldirs_csv(
    path::AbstractString,
    scene::SceneGeometry,
    turtle::TurtleGrid,
    cfg::LightConfig;
    columns::Union{Nothing,AbstractVector}=nothing,
)
    table = node_links_stats_alldirs_table(scene, turtle, cfg; columns=columns)
    return _write_csv_rows(path, table.columns, table.rows; append=false)
end

"""
    node_links_dir_table(scene, turtle, cfg; direction_index=0, columns=nothing)

Build Java-like `log-nodelinks-dirXX.csv` table represented as `(columns, rows)` for one direction.
"""
function node_links_dir_table(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    cfg::LightConfig;
    direction_index::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    0 <= direction_index < length(turtle.sectors) || error("direction_index out of bounds")
    cols = columns === nothing ? copy(_NODE_LINKS_DIR_ORDER) : _canonical_node_links_dir_order(String.(columns))
    item_per_node, comp_per_node = _java_id_maps(scene, cfg)
    vertices, faces, face2node, _, plotbox, _ = _scene_geometry_for_interception(scene, cfg)
    cache_ctx = _projection_cache_context(vertices, faces, face2node, plotbox, cfg)
    sector = turtle.sectors[direction_index + 1]
    pixel_hits, _, _, _ = _direction_projection_cached(vertices, faces, face2node, sector.direction, cfg, plotbox, cache_ctx)

    counts = Dict{Tuple{Int,Int},Int}()
    for stack in values(pixel_hits)
        length(stack) <= 1 && continue
        sort!(stack, by=x -> x[1], rev=true)
        for h in 1:(length(stack) - 1)
            n1 = stack[h][2]
            n2 = stack[h + 1][2]
            k1 = (get(item_per_node, n1, -1), get(comp_per_node, n1, n1))
            k2 = (get(item_per_node, n2, -1), get(comp_per_node, n2, n2))
            a, b = k1 <= k2 ? (n1, n2) : (n2, n1)
            counts[(a, b)] = get(counts, (a, b), 0) + 1
        end
    end

    keys_sorted = sort(collect(keys(counts)); by=k -> (
        get(item_per_node, k[1], -1),
        get(comp_per_node, k[1], k[1]),
        get(item_per_node, k[2], -1),
        get(comp_per_node, k[2], k[2]),
    ))
    rows = Dict{String,Any}[]
    for (n1, n2) in keys_sorted
        row = Dict{String,Any}()
        for c in cols
            if c == "dir"
                row[c] = direction_index
            elseif c == "plantid1"
                row[c] = get(item_per_node, n1, -1)
            elseif c == "id1"
                row[c] = get(comp_per_node, n1, n1)
            elseif c == "plantid2"
                row[c] = get(item_per_node, n2, -1)
            elseif c == "id2"
                row[c] = get(comp_per_node, n2, n2)
            elseif c == "n"
                row[c] = get(counts, (n1, n2), 0)
            else
                row[c] = "NA"
            end
        end
        push!(rows, row)
    end
    return (columns=cols, rows=rows)
end

"""
    write_node_links_dir_csv(path, scene, turtle, cfg; direction_index=0, columns=nothing)

Write Java-like `log-nodelinks-dirXX.csv` for one direction.
"""
function write_node_links_dir_csv(
    path::AbstractString,
    scene::SceneGeometry,
    turtle::TurtleGrid,
    cfg::LightConfig;
    direction_index::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
)
    table = node_links_dir_table(scene, turtle, cfg; direction_index=direction_index, columns=columns)
    return _write_csv_rows(path, table.columns, table.rows; append=false)
end

"""
    summary_values_table(scene, steps, cfg; meteo_rows=nothing, start_step_number=0, columns=nothing, unavailable="NA")

Build a Java-like `summary.csv` table represented as `(columns, rows)`.
Rows are grouped by `(item_id, group, type)` for each step.
"""
function summary_values_table(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
)
    rows_in = meteo_rows === nothing ? fill(nothing, length(steps)) : collect(meteo_rows)
    length(rows_in) == length(steps) || error("summary_values_table: `meteo_rows` length must match `steps` length.")
    cols = columns === nothing ? copy(_SUMMARY_VARIABLE_ORDER) : _canonical_summary_variable_order(String.(columns))
    rows = Dict{String,Any}[]
    meta = _output_node_metadata(scene, cfg)
    node_ids = meta.node_ids
    area_per_node = meta.area_per_node
    item_per_node = meta.item_per_node
    group_per_node = meta.group_per_node
    type_per_node = meta.type_per_node

    for i in eachindex(steps)
        step = steps[i]
        dt = _step_duration_output_local(rows_in[i], nothing)
        step_no = start_step_number + i - 1
        acc = Dict{Tuple{Int,String,String},Tuple{Float64,Float64}}()
        for nid in node_ids
            item = get(item_per_node, nid, -1)
            group = get(group_per_node, nid, "")
            type = get(type_per_node, nid, unavailable)
            area = get(area_per_node, nid, 0.0)
            ri_q = get(step.budget.ri_par_q_per_node, nid, 0.0) + get(step.budget.ri_nir_q_per_node, nid, 0.0)
            k = (item, group, type)
            a0, r0 = get(acc, k, (0.0, 0.0))
            acc[k] = (a0 + area, r0 + ri_q)
        end

        keys_sorted = sort(collect(keys(acc)); by=k -> (k[1], k[2], k[3]))
        for k in keys_sorted
            area, ri_q = acc[k]
            row = Dict{String,Any}()
            for c in cols
                if c == "step_number"
                    row[c] = step_no
                elseif c == "step_duration"
                    row[c] = dt
                elseif c == "group"
                    row[c] = k[2]
                elseif c == "type"
                    row[c] = k[3]
                elseif c == "item_id"
                    row[c] = k[1]
                elseif c == "area"
                    row[c] = area
                elseif c == "Ri_q"
                    row[c] = ri_q
                else
                    row[c] = unavailable
                end
            end
            push!(rows, row)
        end
    end
    return (columns=cols, rows=rows)
end

"""
    write_summary_csv(path, scene, steps, cfg; meteo_rows=nothing, start_step_number=0, columns=nothing, unavailable="NA")

Write Java-like `summary.csv`.
"""
function write_summary_csv(
    path::AbstractString,
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=0,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
)
    table = summary_values_table(
        scene,
        steps,
        cfg;
        meteo_rows=meteo_rows,
        start_step_number=start_step_number,
        columns=columns,
        unavailable=unavailable,
    )
    return _write_csv_rows(path, table.columns, table.rows; append=false)
end

"""
    write_light_outputs(scene, steps, cfg; meteo_rows=nothing, start_step_number=0, outdir=nothing, write_component=nothing, write_scene=nothing, write_summary=nothing, write_sun_position_log=nothing, write_scattering_log=nothing, scattering_log_bands=["PAR"])

High-level Java-style output writer driven by config defaults and optional debug-log toggles.
Returns a dictionary mapping output kind to written path.
"""
function write_light_outputs(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=0,
    outdir::Union{Nothing,AbstractString}=nothing,
    write_component::Union{Nothing,Bool}=nothing,
    write_scene::Union{Nothing,Bool}=nothing,
    write_summary::Union{Nothing,Bool}=nothing,
    write_sun_position_log::Union{Nothing,Bool}=nothing,
    write_scattering_log::Union{Nothing,Bool}=nothing,
    write_export_ops::Union{Nothing,Bool}=nothing,
    write_inputs::Union{Nothing,Bool}=nothing,
    scattering_log_bands::AbstractVector{<:AbstractString}=["PAR"],
)
    rows_in = meteo_rows === nothing ? fill(nothing, length(steps)) : collect(meteo_rows)
    length(rows_in) == length(steps) || error("write_light_outputs: `meteo_rows` length must match `steps` length.")
    out = Dict{String,String}()
    outroot = outdir === nothing ? simulation_output_directory(cfg) : String(outdir)
    mkpath(outroot)
    out["output_directory"] = outroot

    do_component = isnothing(write_component) ? _config_wants_component_outputs(cfg) : Bool(write_component)
    do_scene = isnothing(write_scene) ? true : Bool(write_scene)
    do_summary = isnothing(write_summary) ? cfg.outputs.write_summary : Bool(write_summary)
    debug_default = _config_debug_enabled(cfg)
    do_sun_log = isnothing(write_sun_position_log) ? debug_default : Bool(write_sun_position_log)
    do_scat_log = isnothing(write_scattering_log) ? (debug_default && cfg.general.scattering) : Bool(write_scattering_log)
    do_export_ops = isnothing(write_export_ops) ? (_export_ops_raw_value(cfg) !== nothing) : Bool(write_export_ops)
    do_inputs = isnothing(write_inputs) ? do_export_ops : Bool(write_inputs)

    if do_component
        path = joinpath(outroot, "component_values.csv")
        cols = component_variable_names(cfg)
        first_write = true
        for i in eachindex(steps)
            t = component_values_table(
                scene,
                steps[i],
                cfg;
                meteo_row=rows_in[i],
                step_number=start_step_number + i - 1,
                columns=cols,
            )
            _write_csv_rows(path, t.columns, t.rows; append=!first_write)
            first_write = false
        end
        out["component_values"] = path
    end

    if do_scene
        path = joinpath(outroot, "scene_values.csv")
        write_scene_values_csv(path, scene, steps, cfg; meteo_rows=rows_in, start_step_number=start_step_number)
        out["scene_values"] = path
    end

    if do_summary
        path = joinpath(outroot, "summary.csv")
        write_summary_csv(path, scene, steps, cfg; meteo_rows=rows_in, start_step_number=start_step_number)
        out["summary"] = path
    end

    if do_sun_log
        path = joinpath(outroot, "log-sun-position.csv")
        write_sun_position_log_csv(path, steps, rows_in; start_step_number=start_step_number)
        out["log_sun_position"] = path
    end

    if do_scat_log && cfg.general.scattering
        for b in scattering_log_bands
            band = uppercase(String(b))
            path = joinpath(outroot, "log-iteration-scat-$(lowercase(band)).csv")
            cols = copy(_SCATTERING_ITERATION_LOG_ORDER)
            first_write = true
            for i in eachindex(steps)
                t = scattering_iteration_log_table(
                    scene,
                    steps[i],
                    cfg;
                    meteo_row=rows_in[i],
                    step_number=start_step_number + i - 1,
                    band=band,
                    columns=cols,
                )
                _write_csv_rows(path, t.columns, t.rows; append=!first_write)
                first_write = false
            end
            out["log_iteration_scat_$(lowercase(band))"] = path
        end
    end

    if do_export_ops
        merge!(out, _write_export_ops_outputs(scene, steps, cfg, rows_in, outroot, start_step_number))
    end

    if do_inputs
        scene_exports = sort(
            [(k, v) for (k, v) in out if startswith(k, "ops_step_")];
            by=first,
        )
        if length(scene_exports) == 1
            _, scene_path = only(scene_exports)
            out["config"] = write_light_inputs(outroot, cfg; scene_rel=relpath(scene_path, outroot))
        elseif length(scene_exports) > 1
            for (key, scene_path) in scene_exports
                tag = replace(key, "ops_step_" => "")
                out["config_step_$(tag)"] = write_light_inputs(
                    outroot,
                    cfg;
                    scene_rel=relpath(scene_path, outroot),
                    config_name="config-step$(tag).yml",
                )
            end
        end
    end

    out
end

function write_light_outputs(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=0,
    kwargs...,
)
    write_light_outputs(scene, [step], cfg; meteo_rows=[meteo_row], start_step_number=step_number, kwargs...)
end

"""
    write_component_values_csv(path, scene, step, cfg; meteo_row=nothing, step_number=0, step_duration_seconds=nothing, columns=nothing, unavailable="NA", strict=false)

Write Java-like `component_values.csv` output for one simulation step.
"""
function write_component_values_csv(
    path::AbstractString,
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=0,
    step_duration_seconds=nothing,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
    strict::Bool=false,
)
    table = component_values_table(
        scene,
        step,
        cfg;
        meteo_row=meteo_row,
        step_number=step_number,
        step_duration_seconds=step_duration_seconds,
        columns=columns,
        unavailable=unavailable,
        strict=strict,
    )
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(table.columns, ';'))
        for row in table.rows
            vals = [_csv_cell_local(get(row, c, unavailable)) for c in table.columns]
            println(io, join(vals, ';'))
        end
    end
    return path
end
