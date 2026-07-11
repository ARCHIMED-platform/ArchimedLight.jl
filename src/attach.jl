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
    :sky_fraction => :sky_fraction,
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

function _geometry_node_ids(scene::PlantGeom.SceneGeometry)
    Set{Int}(keys(scene.nodes))
end

function _node_field(scene::PlantGeom.SceneGeometry, step::LightStepResult, field::Symbol)
    field == :area && return PlantGeom.node_areas(scene)
    return _budget_node_field(step, field)
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
    if field == :sky_fraction
        step.sky_fraction === nothing &&
            error("`sky_fraction` was not stored in this LightStepResult. Re-run with `LightOptions(include_sky_fraction=true)` or request `sky_fraction: true` in `component_variables`/`opf_variables` in the config.")
        return step.sky_fraction
    end
    error("Unknown light field selector: $field")
end

function _budget_attr_name(field::Symbol, names::AbstractDict{Symbol,Symbol})
    haskey(names, field) && return names[field]
    field == :area && return :area
    haskey(_DEFAULT_BUDGET_ATTRS, field) && return _DEFAULT_BUDGET_ATTRS[field]
    error("No default MTG attribute name for `$field`.")
end

"""
    attach_node_values!(scene, attr, values; fill_value=nothing)

Attach a dictionary of per-node values to the MTG stored in `scene`.

Only geometry nodes present in the prepared scene are updated. Missing node ids
receive `fill_value`.

Arguments:

- `scene`: MTG-backed `PlantGeom.SceneGeometry` to mutate.
- `attr`: MTG attribute name to write on each geometry node.
- `values`: dictionary keyed by node id, containing values to attach.

Keywords:

- `fill_value`: value written for geometry nodes absent from `values`.
"""
function attach_node_values!(
    scene::PlantGeom.SceneGeometry,
    attr::Symbol,
    values::AbstractDict{<:Integer};
    fill_value=nothing,
)
    scene.mtg === nothing && error("attach_node_values! requires an MTG-backed scene.")
    _attach_node_values!(scene.mtg, _geometry_node_ids(scene), attr, values; fill_value=fill_value)
end

"""
    attach_light_step!(scene, step; fields=[:incident_par_flux], names=Dict(), fill_value=nothing)

Attach one [`LightStepResult`](@ref) back onto the scene MTG using ARCHIMED-style
attribute names such as `Ri_PAR_f` and `Ra_PAR_q`.

Arguments:

- `scene`: MTG-backed `PlantGeom.SceneGeometry` to mutate.
- `step`: one [`LightStepResult`](@ref) whose node-level values will be
  attached.

Keywords:

- `fields`: budget or metadata selectors to attach.
- `names`: optional mapping from selectors in `fields` to custom MTG attribute
  names.
- `fill_value`: value written for geometry nodes absent from a selected step
  dictionary.

`fields` selects which budget components to export. Each selected field is
attached as one scalar MTG attribute per geometry node. Supported selectors are:

- `:incident_par_initial_flux` => `Ri_PAR_0_f`
- `:incident_nir_initial_flux` => `Ri_NIR_0_f`
- `:incident_par_flux` => `Ri_PAR_f`
- `:incident_nir_flux` => `Ri_NIR_f`
- `:incident_par_initial_energy` => `Ri_PAR_0_q`
- `:incident_nir_initial_energy` => `Ri_NIR_0_q`
- `:incident_par_energy` => `Ri_PAR_q`
- `:incident_nir_energy` => `Ri_NIR_q`
- `:absorbed_par_initial_flux` => `Ra_PAR_0_f`
- `:absorbed_nir_initial_flux` => `Ra_NIR_0_f`
- `:absorbed_par_flux` => `Ra_PAR_f`
- `:absorbed_nir_flux` => `Ra_NIR_f`
- `:absorbed_par_initial_energy` => `Ra_PAR_0_q`
- `:absorbed_nir_initial_energy` => `Ra_NIR_0_q`
- `:absorbed_par_energy` => `Ra_PAR_q`
- `:absorbed_nir_energy` => `Ra_NIR_q`
- `:sky_fraction` => `sky_fraction`
- `:area` => `area`

Selector naming follows the budget hierarchy:

- `incident` means intercepted radiation, corresponding to historical `Ri`
- `absorbed` means absorbed radiation, corresponding to historical `Ra`
- `initial` means first-order only, corresponding to historical `_0_`
- no `initial` means first-order plus scattering
- `flux` means `W m^-2`, corresponding to historical `_f`
- `energy` means `J` per component and per step, corresponding to historical `_q`
- `area` means the prepared object surface area in `m^2`

`names` is an optional dictionary that remaps those exported fields to custom MTG
attribute names. For example:

```julia
attach_light_step!(
    scene,
    step;
    fields=[:incident_par_flux, :absorbed_par_energy],
    names=Dict(
        :incident_par_flux => :my_par_flux,
        :absorbed_par_energy => :my_par_energy,
    ),
)
```

With that override, the values are attached on `:my_par_flux` and
`:my_par_energy` instead of the default ARCHIMED names.

To expose the sky-view fraction or PlantBiophysics-specific names, you can mix
selectors and overrides:

```julia
attach_light_step!(
    scene,
    step;
    fields=[:area, :absorbed_par_flux, :absorbed_nir_flux, :sky_fraction],
    names=Dict(:absorbed_nir_flux => :Ra_SW_f),
)
```

`sky_fraction` is only available when the step was produced with
`LightOptions(include_sky_fraction=true)`, or from a config that requests
`sky_fraction: true` in `component_variables` or `opf_variables`.
"""
function attach_light_step!(
    scene::PlantGeom.SceneGeometry,
    step::LightStepResult;
    fields::AbstractVector{Symbol}=[:incident_par_flux],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value=nothing,
)
    for field in fields
        attach_node_values!(
            scene,
            _budget_attr_name(field, names),
            _node_field(scene, step, field);
            fill_value=fill_value,
        )
    end
    return scene.mtg
end

"""
    attach_light_series!(scene, steps; fields=[:incident_par_flux], names=Dict(), fill_value=NaN)

Attach a time series of [`LightStepResult`](@ref) values to the scene MTG.

Arguments:

- `scene`: MTG-backed `PlantGeom.SceneGeometry` to mutate.
- `steps`: ordered vector of [`LightStepResult`](@ref) values to attach.

Keywords:

- `fields`: budget or metadata selectors to attach.
- `names`: optional mapping from selectors in `fields` to custom MTG attribute
  names.
- `fill_value`: value stored for geometry nodes absent from a selected step
  dictionary.

For each selected field, every geometry node receives a vector ordered like
`steps`, which is convenient for downstream plotting or coupled simulations.
For example, with the default field, each geometry node gets
`node[:Ri_PAR_f]::Vector{Float64}` with one value per light step.

`fields` accepts the same selectors as [`attach_light_step!`](@ref):

- `:area` => `area`
- `:incident_par_initial_flux` => `Ri_PAR_0_f`
- `:incident_nir_initial_flux` => `Ri_NIR_0_f`
- `:incident_par_flux` => `Ri_PAR_f`
- `:incident_nir_flux` => `Ri_NIR_f`
- `:incident_par_initial_energy` => `Ri_PAR_0_q`
- `:incident_nir_initial_energy` => `Ri_NIR_0_q`
- `:incident_par_energy` => `Ri_PAR_q`
- `:incident_nir_energy` => `Ri_NIR_q`
- `:absorbed_par_initial_flux` => `Ra_PAR_0_f`
- `:absorbed_nir_initial_flux` => `Ra_NIR_0_f`
- `:absorbed_par_flux` => `Ra_PAR_f`
- `:absorbed_nir_flux` => `Ra_NIR_f`
- `:absorbed_par_initial_energy` => `Ra_PAR_0_q`
- `:absorbed_nir_initial_energy` => `Ra_NIR_0_q`
- `:absorbed_par_energy` => `Ra_PAR_q`
- `:absorbed_nir_energy` => `Ra_NIR_q`
- `:sky_fraction` => `sky_fraction`

For `:area`, series attachment stores the same area value repeated once per
step so the result has the same vector shape as the light fields.

Use `names` to override attached attribute names. This is useful for downstream
packages that expect different names:

```julia
attach_light_series!(
    scene,
    steps;
    fields=[:area, :absorbed_par_flux, :absorbed_nir_flux, :sky_fraction],
    names=Dict(:absorbed_nir_flux => :Ra_SW_f),
)
```

`fill_value` is used for geometry nodes that are missing from a given step
dictionary. The default is `NaN`, so missing values remain visible in the
attached vectors.

`sky_fraction` is only available when each step was produced with
`LightOptions(include_sky_fraction=true)`, or from a config that requests
`sky_fraction: true` in `component_variables` or `opf_variables`.
"""
function attach_light_series!(
    scene::PlantGeom.SceneGeometry,
    steps::AbstractVector{<:LightStepResult};
    fields::AbstractVector{Symbol}=[:incident_par_flux],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value::Float64=NaN,
)
    node_ids = PlantGeom.scene_node_ids(scene)
    for field in fields
        by_node = Dict{Int,Vector{Float64}}(nid => Float64[] for nid in node_ids)
        for step in steps
            vals = _node_field(scene, step, field)
            for nid in node_ids
                push!(by_node[nid], Float64(get(vals, nid, fill_value)))
            end
        end
        attach_node_values!(scene, _budget_attr_name(field, names), by_node; fill_value=fill(fill_value, length(steps)))
    end
    return scene.mtg
end
