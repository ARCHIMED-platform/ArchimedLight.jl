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

"""
    component_variable_names(cfg)::Vector{String}

Resolve component output variable names from `cfg.raw["component_variables"]` when available.
Only variables with truthy values are kept.
"""
function component_variable_names(cfg::LightConfig)
    d = get(cfg.raw, "component_variables", nothing)
    d isa AbstractDict || return copy(_COMPONENT_VARIABLE_ORDER)

    vars = String[]
    for (k, v) in d
        _as_bool(v, false) || continue
        push!(vars, string(k))
    end
    isempty(vars) ? copy(_COMPONENT_VARIABLE_ORDER) : _canonical_component_variable_order(vars)
end

"""
    scene_variable_names(cfg)::Vector{String}

Resolve scene output variable names from `cfg.raw["scene_variables"]` when available.
Only variables with truthy values are kept. Defaults match Java `scene_values.csv`.
"""
function scene_variable_names(cfg::LightConfig)
    d = get(cfg.raw, "scene_variables", nothing)
    d isa AbstractDict || return copy(_SCENE_VARIABLE_ORDER)

    vars = String[]
    for (k, v) in d
        _as_bool(v, false) || continue
        push!(vars, string(k))
    end
    isempty(vars) ? copy(_SCENE_VARIABLE_ORDER) : _canonical_scene_variable_order(vars)
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

function _group_type_hints(cfg::LightConfig)
    models = get(cfg.raw, "models", nothing)
    models isa AbstractVector || return Dict{String,Vector{String}}()

    base = get(cfg.raw, "__base_dir", dirname(cfg.scene))
    out = Dict{String,Vector{String}}()
    for m in models
        path = isabspath(String(m)) ? String(m) : normpath(joinpath(base, String(m)))
        isfile(path) || continue
        d = try
            YAML.load_file(path)
        catch
            nothing
        end
        d isa AbstractDict || continue
        dd = _to_string_dict(d)
        group = haskey(dd, "Group") ? strip(string(dd["Group"])) : ""
        isempty(group) && continue
        types = get(dd, "Type", nothing)
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
    absorptance_cache::Dict{String,Dict{Int,Float64}},
    group_type_hints::Dict{String,Vector{String}};
    unavailable::String="NA",
    strict::Bool=false,
)
    area = get(scene.total_area_per_node, nid, 0.0)
    hits = get(step.first_order.hits_per_node, nid, 0)

    if variable == "step_number"
        return step_number
    elseif variable == "step_duration"
        return step_duration
    elseif variable == "item_id"
        return get(scene.java_item_id_per_node, nid, -1)
    elseif variable == "component_id"
        return get(scene.java_component_id_per_node, nid, nid)
    elseif variable == "group"
        return get(scene.node_group, nid, "")
    elseif variable == "type"
        t = strip(get(scene.node_type, nid, ""))
        return if isempty(t)
            g = get(scene.node_group, nid, "")
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
        return get(scene.barycenter_per_node, nid, (NaN, NaN, NaN))[1]
    elseif variable == "barycentre_y"
        return get(scene.barycenter_per_node, nid, (NaN, NaN, NaN))[2]
    elseif variable == "barycentre_z"
        return get(scene.barycenter_per_node, nid, (NaN, NaN, NaN))[3]
    elseif variable == "sky_fraction"
        return unavailable
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
    component_values_table(scene, step, cfg; meteo_row=nothing, step_number=1, step_duration_seconds=nothing, columns=nothing, unavailable="NA", strict=false)

Build a Java-like component values table represented as `(columns, rows)`.
`rows` is a vector of dictionaries keyed by column name.
"""
function component_values_table(
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=1,
    step_duration_seconds=nothing,
    columns::Union{Nothing,AbstractVector}=nothing,
    unavailable::String="NA",
    strict::Bool=false,
)
    cols = columns === nothing ? component_variable_names(cfg) : String.(columns)
    step_duration = _step_duration_output_local(meteo_row, step_duration_seconds)
    node_ids = _node_ids_for_output(scene)
    rows = Vector{Dict{String,Any}}(undef, length(node_ids))
    absorptance_cache = Dict{String,Dict{Int,Float64}}()
    group_type_hints = _group_type_hints(cfg)

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
                absorptance_cache,
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

"""
    scene_values_table(scene, steps, cfg; meteo_rows=nothing, start_step_number=1, columns=nothing, unavailable="NA", strict=false)

Build a Java-like scene values table represented as `(columns, rows)`.
`rows` is a vector of dictionaries keyed by column name.
"""
function scene_values_table(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=1,
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
    write_scene_values_csv(path, scene, steps, cfg; meteo_rows=nothing, start_step_number=1, columns=nothing, unavailable="NA", strict=false)

Write Java-like `scene_values.csv` output for one or many simulation steps.
"""
function write_scene_values_csv(
    path::AbstractString,
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult},
    cfg::LightConfig;
    meteo_rows=nothing,
    start_step_number::Int=1,
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
    write_component_values_csv(path, scene, step, cfg; meteo_row=nothing, step_number=1, step_duration_seconds=nothing, columns=nothing, unavailable="NA", strict=false)

Write Java-like `component_values.csv` output for one simulation step.
"""
function write_component_values_csv(
    path::AbstractString,
    scene::SceneGeometry,
    step::LightStepResult,
    cfg::LightConfig;
    meteo_row=nothing,
    step_number::Int=1,
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
