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
