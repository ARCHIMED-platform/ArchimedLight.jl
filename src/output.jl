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

function _canonical_component_variable_order(names::Vector{String})
    idx = Dict{String,Int}(n => i for (i, n) in enumerate(_COMPONENT_VARIABLE_ORDER))
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
    absorptance_cache::Dict{String,Dict{Int,Float64}};
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
        return unavailable
    elseif variable == "area"
        return area
    elseif variable == "surface_hits"
        return area > 0 ? hits / area : 0.0
    elseif variable == "barycentre_x" || variable == "barycentre_y" || variable == "barycentre_z" || variable == "sky_fraction"
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
                absorptance_cache;
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
