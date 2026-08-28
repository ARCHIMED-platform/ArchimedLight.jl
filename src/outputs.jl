"""
    CompiledComponentAggregation(metadata)

Precompute the stable component-to-owner mapping used to aggregate one
[`LightStepResult`](@ref) from radiative scene components to botanical
`PlantGeom.SourceOwnerKey`s.

Component energy quantities (`*_q`) are summed. Component flux densities
(`*_f`) are averaged with `radiative_area` as weight. The mapping is valid only
for results sharing the exact [`LightComponentMetadata`](@ref) snapshot used at
construction; rebuild it after [`update_scene!`](@ref).
"""
struct CompiledComponentAggregation
    metadata::LightComponentMetadata
    owner_keys::_ReadOnlyVector{PlantGeom.SourceOwnerKey,Vector{PlantGeom.SourceOwnerKey}}
    component_to_owner::_ReadOnlyVector{Int,Vector{Int}}
    owner_radiative_area::_ReadOnlyVector{Float64,Vector{Float64}}

    function CompiledComponentAggregation(
        metadata::LightComponentMetadata,
        owner_keys::Vector{PlantGeom.SourceOwnerKey},
        component_to_owner::Vector{Int},
        owner_radiative_area::Vector{Float64},
    )
        ncomponents = length(metadata.node_id)
        nowners = length(owner_keys)
        length(component_to_owner) == ncomponents || throw(DimensionMismatch(
            "component_to_owner has length $(length(component_to_owner)); " *
            "expected $ncomponents.",
        ))
        length(owner_radiative_area) == nowners || throw(DimensionMismatch(
            "owner_radiative_area has length $(length(owner_radiative_area)); " *
            "expected $nowners.",
        ))
        length(unique(owner_keys)) == nowners || throw(ArgumentError(
            "Compiled component owner keys must be unique.",
        ))

        recomputed_area = zeros(Float64, nowners)
        @inbounds for component in eachindex(component_to_owner)
            destination = component_to_owner[component]
            1 <= destination <= nowners || throw(ArgumentError(
                "Component $component maps to invalid owner position $destination; " *
                "expected a position in 1:$nowners.",
            ))
            owner_keys[destination] == metadata.source_owner[component] || throw(ArgumentError(
                "Component $component does not map to its declared source owner.",
            ))
            recomputed_area[destination] += metadata.radiative_area[component]
        end
        recomputed_area == owner_radiative_area || throw(ArgumentError(
            "owner_radiative_area is inconsistent with component ownership and areas.",
        ))
        new(
            metadata,
            _ReadOnlyVector(copy(owner_keys)),
            _ReadOnlyVector(copy(component_to_owner)),
            _ReadOnlyVector(copy(owner_radiative_area)),
        )
    end
end

function CompiledComponentAggregation(metadata::LightComponentMetadata)
    owner_keys = unique(copy(metadata.source_owner))
    sort!(owner_keys; by=key -> (key.source_instance_id, key.source_node_id))
    owner_position = Dict(key => i for (i, key) in pairs(owner_keys))
    component_to_owner = Vector{Int}(undef, length(metadata.source_owner))
    owner_area = zeros(Float64, length(owner_keys))

    for i in eachindex(metadata.source_owner)
        position = owner_position[metadata.source_owner[i]]
        area = metadata.radiative_area[i]
        @inbounds component_to_owner[i] = position
        @inbounds owner_area[position] += area
    end
    for (i, area) in pairs(owner_area)
        isfinite(area) && area > 0.0 || throw(ArgumentError(
            "Source owner $(owner_keys[i]) has zero or invalid total radiative area $area; " *
            "area-weighted flux aggregation is undefined.",
        ))
    end

    return CompiledComponentAggregation(
        metadata,
        owner_keys,
        component_to_owner,
        owner_area,
    )
end

CompiledComponentAggregation(step::LightStepResult) =
    CompiledComponentAggregation(_component_metadata(step))

@inline function _component_metadata(step::LightStepResult)
    metadata = step.component_metadata
    metadata === nothing && throw(ArgumentError(
        "This LightStepResult has no retained component ownership metadata. " *
        "Use a result returned by `run_light` with a current PlantGeom scene.",
    ))
    return metadata
end

function _new_light_metric_columns(n::Int)
    return (
        Ri_PAR_0_f=zeros(Float64, n),
        Ri_NIR_0_f=zeros(Float64, n),
        Ri_PAR_f=zeros(Float64, n),
        Ri_NIR_f=zeros(Float64, n),
        Ri_PAR_0_q=zeros(Float64, n),
        Ri_NIR_0_q=zeros(Float64, n),
        Ri_PAR_q=zeros(Float64, n),
        Ri_NIR_q=zeros(Float64, n),
        Ra_PAR_0_f=zeros(Float64, n),
        Ra_NIR_0_f=zeros(Float64, n),
        Ra_PAR_f=zeros(Float64, n),
        Ra_NIR_f=zeros(Float64, n),
        Ra_PAR_0_q=zeros(Float64, n),
        Ra_NIR_0_q=zeros(Float64, n),
        Ra_PAR_q=zeros(Float64, n),
        Ra_NIR_q=zeros(Float64, n),
        Ra_SW_f=zeros(Float64, n),
        aPPFD=zeros(Float64, n),
    )
end

@inline function _validated_par_energy_to_photon(value::Real)
    factor = Float64(value)
    isfinite(factor) && factor > 0.0 || throw(ArgumentError(
        "par_energy_to_photon must be a finite positive conversion factor; got $value.",
    ))
    return factor
end

@inline function _require_light_metric_column(table, ::Val{name}, n::Int) where {name}
    hasproperty(table, name) || throw(ArgumentError(
        "The component-values buffer is missing declared column `$name`.",
    ))
    column = getproperty(table, name)
    column isa Vector{Float64} || throw(ArgumentError(
        "The component-values buffer column `$name` must be a Vector{Float64} " *
        "returned by `component_values`; " *
        "got $(typeof(column)).",
    ))
    length(column) == n || throw(DimensionMismatch(
        "The component-values buffer column `$name` has length $(length(column)); expected $n.",
    ))
    return column
end


function _validate_light_metric_buffer(table, n::Int)
    # Keep this explicit so the normal refill path is fully inferred and does
    # not allocate a dynamic collection of column names on every light step.
    _require_light_metric_column(table, Val(:Ri_PAR_0_f), n)
    _require_light_metric_column(table, Val(:Ri_NIR_0_f), n)
    _require_light_metric_column(table, Val(:Ri_PAR_f), n)
    _require_light_metric_column(table, Val(:Ri_NIR_f), n)
    _require_light_metric_column(table, Val(:Ri_PAR_0_q), n)
    _require_light_metric_column(table, Val(:Ri_NIR_0_q), n)
    _require_light_metric_column(table, Val(:Ri_PAR_q), n)
    _require_light_metric_column(table, Val(:Ri_NIR_q), n)
    _require_light_metric_column(table, Val(:Ra_PAR_0_f), n)
    _require_light_metric_column(table, Val(:Ra_NIR_0_f), n)
    _require_light_metric_column(table, Val(:Ra_PAR_f), n)
    _require_light_metric_column(table, Val(:Ra_NIR_f), n)
    _require_light_metric_column(table, Val(:Ra_PAR_0_q), n)
    _require_light_metric_column(table, Val(:Ra_NIR_0_q), n)
    _require_light_metric_column(table, Val(:Ra_PAR_q), n)
    _require_light_metric_column(table, Val(:Ra_NIR_q), n)
    _require_light_metric_column(table, Val(:Ra_SW_f), n)
    _require_light_metric_column(table, Val(:aPPFD), n)
    return table
end

@inline function _light_budget_maps(step::LightStepResult)
    budget = step.budget
    return (
        Ri_PAR_0_f=budget.incident_flux.initial.par,
        Ri_NIR_0_f=budget.incident_flux.initial.nir,
        Ri_PAR_f=budget.incident_flux.total.par,
        Ri_NIR_f=budget.incident_flux.total.nir,
        Ri_PAR_0_q=budget.incident_energy.initial.par,
        Ri_NIR_0_q=budget.incident_energy.initial.nir,
        Ri_PAR_q=budget.incident_energy.total.par,
        Ri_NIR_q=budget.incident_energy.total.nir,
        Ra_PAR_0_f=budget.absorbed_flux.initial.par,
        Ra_NIR_0_f=budget.absorbed_flux.initial.nir,
        Ra_PAR_f=budget.absorbed_flux.total.par,
        Ra_NIR_f=budget.absorbed_flux.total.nir,
        Ra_PAR_0_q=budget.absorbed_energy.initial.par,
        Ra_NIR_0_q=budget.absorbed_energy.initial.nir,
        Ra_PAR_q=budget.absorbed_energy.total.par,
        Ra_NIR_q=budget.absorbed_energy.total.nir,
    )
end

@inline function _copy_component_metric!(destination, source, node_ids)
    @inbounds for i in eachindex(node_ids)
        destination[i] = get(source, node_ids[i], 0.0)
    end
    return destination
end

function _fill_component_metric_columns!(columns, step::LightStepResult, metadata)
    maps = _light_budget_maps(step)
    ids = metadata.node_id
    _copy_component_metric!(columns.Ri_PAR_0_f, maps.Ri_PAR_0_f, ids)
    _copy_component_metric!(columns.Ri_NIR_0_f, maps.Ri_NIR_0_f, ids)
    _copy_component_metric!(columns.Ri_PAR_f, maps.Ri_PAR_f, ids)
    _copy_component_metric!(columns.Ri_NIR_f, maps.Ri_NIR_f, ids)
    _copy_component_metric!(columns.Ri_PAR_0_q, maps.Ri_PAR_0_q, ids)
    _copy_component_metric!(columns.Ri_NIR_0_q, maps.Ri_NIR_0_q, ids)
    _copy_component_metric!(columns.Ri_PAR_q, maps.Ri_PAR_q, ids)
    _copy_component_metric!(columns.Ri_NIR_q, maps.Ri_NIR_q, ids)
    _copy_component_metric!(columns.Ra_PAR_0_f, maps.Ra_PAR_0_f, ids)
    _copy_component_metric!(columns.Ra_NIR_0_f, maps.Ra_NIR_0_f, ids)
    _copy_component_metric!(columns.Ra_PAR_f, maps.Ra_PAR_f, ids)
    _copy_component_metric!(columns.Ra_NIR_f, maps.Ra_NIR_f, ids)
    _copy_component_metric!(columns.Ra_PAR_0_q, maps.Ra_PAR_0_q, ids)
    _copy_component_metric!(columns.Ra_NIR_0_q, maps.Ra_NIR_0_q, ids)
    _copy_component_metric!(columns.Ra_PAR_q, maps.Ra_PAR_q, ids)
    _copy_component_metric!(columns.Ra_NIR_q, maps.Ra_NIR_q, ids)
    return columns
end

@inline function _sum_owner_metric!(destination, source, node_ids, component_to_owner)
    fill!(destination, 0.0)
    @inbounds for i in eachindex(node_ids)
        destination[component_to_owner[i]] += get(source, node_ids[i], 0.0)
    end
    return destination
end

@inline function _mean_owner_metric!(
    destination,
    source,
    node_ids,
    component_to_owner,
    component_area,
    owner_area,
)
    fill!(destination, 0.0)
    @inbounds for i in eachindex(node_ids)
        destination[component_to_owner[i]] += get(source, node_ids[i], 0.0) * component_area[i]
    end
    @inbounds for i in eachindex(destination)
        destination[i] /= owner_area[i]
    end
    return destination
end

function _fill_aggregated_metric_columns!(
    columns,
    step::LightStepResult,
    metadata::LightComponentMetadata,
    component_to_destination,
    destination_area,
)
    maps = _light_budget_maps(step)
    ids = metadata.node_id
    positions = component_to_destination

    # Energy variables are extensive and therefore sum over components.
    _sum_owner_metric!(columns.Ri_PAR_0_q, maps.Ri_PAR_0_q, ids, positions)
    _sum_owner_metric!(columns.Ri_NIR_0_q, maps.Ri_NIR_0_q, ids, positions)
    _sum_owner_metric!(columns.Ri_PAR_q, maps.Ri_PAR_q, ids, positions)
    _sum_owner_metric!(columns.Ri_NIR_q, maps.Ri_NIR_q, ids, positions)
    _sum_owner_metric!(columns.Ra_PAR_0_q, maps.Ra_PAR_0_q, ids, positions)
    _sum_owner_metric!(columns.Ra_NIR_0_q, maps.Ra_NIR_0_q, ids, positions)
    _sum_owner_metric!(columns.Ra_PAR_q, maps.Ra_PAR_q, ids, positions)
    _sum_owner_metric!(columns.Ra_NIR_q, maps.Ra_NIR_q, ids, positions)

    # Flux densities are intensive and use the exact radiative denominator
    # employed by the light solver for each component.
    area = metadata.radiative_area
    _mean_owner_metric!(columns.Ri_PAR_0_f, maps.Ri_PAR_0_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ri_NIR_0_f, maps.Ri_NIR_0_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ri_PAR_f, maps.Ri_PAR_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ri_NIR_f, maps.Ri_NIR_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ra_PAR_0_f, maps.Ra_PAR_0_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ra_NIR_0_f, maps.Ra_NIR_0_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ra_PAR_f, maps.Ra_PAR_f, ids, positions, area, destination_area)
    _mean_owner_metric!(columns.Ra_NIR_f, maps.Ra_NIR_f, ids, positions, area, destination_area)
    return columns
end

function _fill_owner_metric_columns!(columns, step::LightStepResult, aggregation)
    return _fill_aggregated_metric_columns!(
        columns,
        step,
        aggregation.metadata,
        aggregation.component_to_owner,
        aggregation.owner_radiative_area,
    )
end

function _derive_light_coupling_columns!(columns, par_energy_to_photon::Float64)
    @inbounds for i in eachindex(columns.Ra_SW_f)
        columns.Ra_SW_f[i] = columns.Ra_PAR_f[i] + columns.Ra_NIR_f[i]
        columns.aPPFD[i] = columns.Ra_PAR_f[i] * par_energy_to_photon
    end
    return columns
end

"""
    component_values(step; level=:component, par_energy_to_photon=4.57)
    component_values(step, aggregation; par_energy_to_photon=4.57)

Return a typed Tables.jl-compatible column table for one light step.

At `level=:component`, rows are the filtered geometric components and carry
`node_id`, `source_owner`, `radiative_mesh_area`, all initial and total PAR/NIR
incident and absorbed flux/energy variables, `Ra_SW_f`, and `aPPFD`.

At `level=:source_owner`, components sharing a botanical owner are coalesced.
Energy quantities (`*_q`) are summed and flux densities (`*_f`) are weighted by
radiative area. Pass a reusable [`CompiledComponentAggregation`](@ref) as the
second argument to avoid recompiling this identity mapping.

`par_energy_to_photon` converts absorbed PAR energy flux (W m⁻²) to absorbed
PPFD (μmol m⁻² s⁻¹). The default `4.57` is the conventional broadband PAR
factor; supply a spectrum-specific value when available.
"""
function component_values(
    step::LightStepResult;
    level::Symbol=:component,
    par_energy_to_photon::Real=4.57,
)
    metadata = _component_metadata(step)
    conversion = _validated_par_energy_to_photon(par_energy_to_photon)
    if level === :component
        columns = _new_light_metric_columns(length(metadata.node_id))
        _fill_component_metric_columns!(columns, step, metadata)
        _derive_light_coupling_columns!(columns, conversion)
        return merge(
            (
                node_id=metadata.node_id,
                source_owner=metadata.source_owner,
                radiative_mesh_area=metadata.radiative_area,
            ),
            columns,
        )
    elseif level === :source_owner
        return component_values(
            step,
            CompiledComponentAggregation(metadata);
            par_energy_to_photon=par_energy_to_photon,
        )
    end
    throw(ArgumentError("level must be :component or :source_owner; got $(repr(level))."))
end

function component_values(
    step::LightStepResult,
    aggregation::CompiledComponentAggregation;
    par_energy_to_photon::Real=4.57,
)
    metadata = _component_metadata(step)
    metadata === aggregation.metadata || throw(ArgumentError(
        "The compiled component aggregation belongs to another scene generation. " *
        "Rebuild it from this result's component metadata.",
    ))
    table = merge(
        (
            source_owner=aggregation.owner_keys,
            radiative_mesh_area=aggregation.owner_radiative_area,
        ),
        _new_light_metric_columns(length(aggregation.owner_keys)),
    )
    return component_values!(
        table,
        step,
        aggregation;
        par_energy_to_photon=par_energy_to_photon,
    )
end

"""
    component_values!(table, step, aggregation; par_energy_to_photon=4.57)

Refill a preallocated source-owner table returned by
[`component_values`](@ref). The identity and result columns are reused in
place, so scene-level coupling code can compile the owner mapping once and
avoid allocating a result table at every light step.

`table.source_owner` and `table.radiative_mesh_area` must be the exact columns
owned by `aggregation`; this prevents a buffer from another scene generation
from being reused accidentally.
"""
function component_values!(
    table,
    step::LightStepResult,
    aggregation::CompiledComponentAggregation;
    par_energy_to_photon::Real=4.57,
)
    metadata = _component_metadata(step)
    metadata === aggregation.metadata || throw(ArgumentError(
        "The compiled component aggregation belongs to another scene generation. " *
        "Rebuild it from this result's component metadata.",
    ))
    factor = _validated_par_energy_to_photon(par_energy_to_photon)
    table.source_owner === aggregation.owner_keys || throw(ArgumentError(
        "The result buffer belongs to another component aggregation.",
    ))
    table.radiative_mesh_area === aggregation.owner_radiative_area || throw(ArgumentError(
        "The result buffer radiative-area column belongs to another component aggregation.",
    ))
    _validate_light_metric_buffer(table, length(aggregation.owner_keys))
    _fill_owner_metric_columns!(table, step, aggregation)
    _derive_light_coupling_columns!(table, factor)
    return table
end

function _component_output_node_ids(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    keys_by_node = _interception_output_keys(scene, models, options)
    ids = collect(geometry.node_ids)
    sort!(
        ids;
        by=nid -> (
            get(keys_by_node, nid, (_scene_object_id(scene, nid, 1), _scene_source_topology_id(scene, nid, nid)))[1],
            get(keys_by_node, nid, (_scene_object_id(scene, nid, 1), _scene_source_topology_id(scene, nid, nid)))[2],
            nid,
        ),
    )
    return ids
end

function _component_sky_fraction_per_node(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    step::LightStepResult,
    options::LightOptions,
    node_ids,
)
    step.sky_fraction !== nothing && return step.sky_fraction

    responses = _build_sector_responses(scene, models, step.turtle, options)
    sky_count = 0
    visible_sum = Dict{Int,Float64}()
    for (i, sector) in enumerate(step.turtle.sectors)
        sector.source == :sun && continue
        sky_count += 1
        for (nid, area) in responses.projected_area_per_sector[i]
            visible_sum[nid] = get(visible_sum, nid, 0.0) + area
        end
    end

    out = Dict{Int,Float64}()
    for nid in node_ids
        area = _scene_area(scene, nid, 0.0)
        if area <= 0.0 || sky_count == 0
            out[nid] = 0.0
        else
            out[nid] = get(visible_sum, nid, 0.0) / sky_count / area
        end
    end
    return out
end

function _component_display_type(scene::PlantGeom.SceneGeometry, models::LightModels, nid::Int)
    _scene_display_type(scene, models, nid)
end

function _component_rows_for_step(
    sim::LightSimulation,
    step::LightStepResult,
    step_number::Int,
    node_ids,
    keys_by_node;
    include_scattering_columns::Bool,
)
    scene = sim.scene
    models = sim.models
    options = sim.options
    sky_fraction = _component_sky_fraction_per_node(scene, models, step, options, node_ids)
    rows = OrderedDict{String,Any}[]
    for nid in node_ids
        barycenter = _scene_barycenter(scene, nid, (NaN, NaN, NaN))
        object_id = _scene_object_id(scene, nid, -1)
        source_topology_id = _scene_source_topology_id(scene, nid, nid)
        item_id, component_id = get(keys_by_node, nid, (object_id, source_topology_id))
        row = OrderedDict{String,Any}(
            "step_number" => step_number,
            "node_id" => nid,
            "source_topology_id" => source_topology_id,
            "object_id" => object_id,
            "item_id" => item_id,
            "component_id" => component_id,
            "group" => _scene_group(scene, nid, ""),
            "type" => _component_display_type(scene, models, nid),
            "area" => _scene_area(scene, nid, 0.0),
            "barycentre_z" => barycenter[3],
            "sky_fraction" => get(sky_fraction, nid, 0.0),
            "Ri_PAR_0_f" => get(step.budget.incident_flux.initial.par, nid, 0.0),
            "Ri_NIR_0_f" => get(step.budget.incident_flux.initial.nir, nid, 0.0),
            "Ri_PAR_0_q" => get(step.budget.incident_energy.initial.par, nid, 0.0),
            "Ri_NIR_0_q" => get(step.budget.incident_energy.initial.nir, nid, 0.0),
            "Ra_PAR_0_q" => get(step.budget.absorbed_energy.initial.par, nid, 0.0),
            "Ra_NIR_0_q" => get(step.budget.absorbed_energy.initial.nir, nid, 0.0),
        )
        if include_scattering_columns
            row["Ri_PAR_f"] = get(step.budget.incident_flux.total.par, nid, 0.0)
            row["Ri_NIR_f"] = get(step.budget.incident_flux.total.nir, nid, 0.0)
            row["Ri_PAR_q"] = get(step.budget.incident_energy.total.par, nid, 0.0)
            row["Ri_NIR_q"] = get(step.budget.incident_energy.total.nir, nid, 0.0)
            row["Ra_PAR_q"] = get(step.budget.absorbed_energy.total.par, nid, 0.0)
            row["Ra_NIR_q"] = get(step.budget.absorbed_energy.total.nir, nid, 0.0)
        end
        push!(rows, row)
    end
    return rows
end

_component_steps(step::LightStepResult) = LightStepResult[step]
_component_steps(series) = collect(series)

"""
    write_component_values(path::AbstractString, sim::LightSimulation, series; step_index_base::Integer=1)

Write an ARCHIMED-style `component_values.csv` file from already-computed
[`LightStepResult`](@ref) values.

`series` may be a single `LightStepResult` or a collection of them. The function
serializes the supplied results only; it does not run the simulation. Step
numbers are 1-based by default. Use `step_index_base=0` only for historical
harness compatibility.

Arguments:

- `path`: destination CSV file path.
- `sim`: [`LightSimulation`](@ref) that supplies scene, model, and option
  context for component rows.
- `series`: one [`LightStepResult`](@ref) or a collection of results to
  serialize.

Keywords:

- `step_index_base`: first step number written to the CSV. The default is `1`.
"""
function write_component_values(
    path::AbstractString,
    sim::LightSimulation,
    series;
    step_index_base::Integer=1,
)
    steps = _component_steps(series)
    all(step -> step isa LightStepResult, steps) ||
        throw(ArgumentError("series must be a LightStepResult or a collection of LightStepResult values"))

    parent = dirname(String(path))
    if !isempty(parent) && parent != "."
        mkpath(parent)
    end

    node_ids = _component_output_node_ids(sim.scene, sim.models, sim.options)
    keys_by_node = _interception_output_keys(sim.scene, sim.models, sim.options)
    include_scattering_columns = any(step -> step.scattering !== nothing, steps)

    rows = OrderedDict{String,Any}[]
    for (i, step) in enumerate(steps)
        append!(
            rows,
            _component_rows_for_step(
                sim,
                step,
                Int(step_index_base) + i - 1,
                node_ids,
                keys_by_node;
                include_scattering_columns=include_scattering_columns,
            ),
        )
    end
    CSV.write(path, rows; delim=';')
    return path
end
