const _DEFAULT_BUDGET_ATTRS = Dict{Symbol,Symbol}(
    :incident_par_initial_flux => :Ri_PAR_0_f,
    :incident_nir_initial_flux => :Ri_NIR_0_f,
    :incident_par_flux => :Ri_PAR_f,
    :incident_nir_flux => :Ri_NIR_f,
    :incident_par_initial_energy => :Ri_PAR_0_q,
    :incident_nir_initial_energy => :Ri_NIR_0_q,
    :incident_par_energy => :Ri_PAR_q,
    :incident_nir_energy => :Ri_NIR_q,
    :absorbed_par_initial_flux => :Ra_PAR_0_f,
    :absorbed_nir_initial_flux => :Ra_NIR_0_f,
    :absorbed_par_flux => :Ra_PAR_f,
    :absorbed_nir_flux => :Ra_NIR_f,
    :absorbed_par_initial_energy => :Ra_PAR_0_q,
    :absorbed_nir_initial_energy => :Ra_NIR_0_q,
    :absorbed_par_energy => :Ra_PAR_q,
    :absorbed_nir_energy => :Ra_NIR_q,
)

function _attach_node_values!(
    mtg,
    node_ids,
    attr::Symbol,
    values::AbstractDict{<:Integer};
    fill_value=nothing,
)
    MultiScaleTreeGraph.traverse!(mtg) do node
        nid = Int(MultiScaleTreeGraph.node_id(node))
        if nid in node_ids
            node[attr] = get(values, nid, fill_value)
        end
        return true
    end
    return mtg
end

function _geometry_node_ids(scene::SceneGeometry)
    Set{Int}(keys(scene.nodes))
end

function _budget_node_field(step::LightStepResult, field::Symbol)
    budget = step.budget
    field == :incident_par_initial_flux && return budget.incident_flux.initial.par
    field == :incident_nir_initial_flux && return budget.incident_flux.initial.nir
    field == :incident_par_flux && return budget.incident_flux.total.par
    field == :incident_nir_flux && return budget.incident_flux.total.nir
    field == :incident_par_initial_energy && return budget.incident_energy.initial.par
    field == :incident_nir_initial_energy && return budget.incident_energy.initial.nir
    field == :incident_par_energy && return budget.incident_energy.total.par
    field == :incident_nir_energy && return budget.incident_energy.total.nir
    field == :absorbed_par_initial_flux && return budget.absorbed_flux.initial.par
    field == :absorbed_nir_initial_flux && return budget.absorbed_flux.initial.nir
    field == :absorbed_par_flux && return budget.absorbed_flux.total.par
    field == :absorbed_nir_flux && return budget.absorbed_flux.total.nir
    field == :absorbed_par_initial_energy && return budget.absorbed_energy.initial.par
    field == :absorbed_nir_initial_energy && return budget.absorbed_energy.initial.nir
    field == :absorbed_par_energy && return budget.absorbed_energy.total.par
    field == :absorbed_nir_energy && return budget.absorbed_energy.total.nir
    error("Unknown LightBudget field selector: $field")
end

function _budget_attr_name(field::Symbol, names::AbstractDict{Symbol,Symbol})
    haskey(names, field) && return names[field]
    haskey(_DEFAULT_BUDGET_ATTRS, field) && return _DEFAULT_BUDGET_ATTRS[field]
    error("No default MTG attribute name for `$field`.")
end

function attach_node_values!(
    scene::SceneGeometry,
    attr::Symbol,
    values::AbstractDict{<:Integer};
    fill_value=nothing,
)
    scene.mtg === nothing && error("attach_node_values! requires an MTG-backed scene.")
    _attach_node_values!(scene.mtg, _geometry_node_ids(scene), attr, values; fill_value=fill_value)
end

function attach_light_step!(
    scene::SceneGeometry,
    step::LightStepResult;
    fields::AbstractVector{Symbol}=[:incident_par_flux],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value=nothing,
)
    for field in fields
        attach_node_values!(
            scene,
            _budget_attr_name(field, names),
            _budget_node_field(step, field);
            fill_value=fill_value,
        )
    end
    return scene.mtg
end

function attach_light_series!(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult};
    fields::AbstractVector{Symbol}=[:incident_par_flux],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value::Float64=NaN,
)
    node_ids = scene_node_ids(scene)
    for field in fields
        by_node = Dict{Int,Vector{Float64}}(nid => Float64[] for nid in node_ids)
        for step in steps
            vals = _budget_node_field(step, field)
            for nid in node_ids
                push!(by_node[nid], Float64(get(vals, nid, fill_value)))
            end
        end
        attach_node_values!(scene, _budget_attr_name(field, names), by_node; fill_value=fill(fill_value, length(steps)))
    end
    return scene.mtg
end
