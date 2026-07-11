const _LIGHT_QUERY_BASE_COLUMNS = (
    :step_number,
    :node_id,
    :source_topology_id,
    :object_id,
    :item_id,
    :component_id,
    :group,
    :type,
    :symbol,
    :scale,
)

_light_query_scene(scene::PlantGeom.SceneGeometry) = scene
_light_query_scene(sim::LightSimulation) = sim.scene

function _light_query_candidate_ids(scene::PlantGeom.SceneGeometry)
    sort!(Int[nid for nid in keys(scene.nodes)])
end

function _light_query_candidate_ids(sim::LightSimulation)
    geometry = _scene_geometry_for_interception(sim.scene, sim.models, sim.options)
    sort!(collect(Int, geometry.node_ids))
end

function _light_query_output_keys(scene::PlantGeom.SceneGeometry, node_ids)
    _interception_output_keys_for_node_ids(scene, node_ids)
end

_light_query_output_keys(sim::LightSimulation, node_ids) =
    _interception_output_keys(sim.scene, sim.models, sim.options)

_light_query_type(scene::PlantGeom.SceneGeometry, nid::Int) = _scene_type(scene, nid, "")
_light_query_type(sim::LightSimulation, nid::Int) =
    _scene_display_type(sim.scene, sim.models, nid)

function _light_query_values(value)
    value isa AbstractString && return (value,)
    value isa Symbol && return (value,)
    value isa Pair && return (value,)
    value isa Tuple && return value
    value isa AbstractArray && return value
    value isa AbstractSet && return value
    return (value,)
end

_light_query_matches(actual, expected) = any(candidate -> isequal(actual, candidate), _light_query_values(expected))

function _light_query_symbol_matches(actual, expected)
    any(candidate -> isequal(actual, Symbol(candidate)), _light_query_values(expected))
end

function _light_query_string_matches(actual, expected)
    any(candidate -> isequal(string(actual), string(candidate)), _light_query_values(expected))
end

function _light_query_node(scene::PlantGeom.SceneGeometry, nid::Int)
    scene.mtg === nothing && throw(
        ArgumentError("Scene metadata queries require an MTG-backed scene."),
    )
    node = MultiScaleTreeGraph.get_node(scene.mtg, nid)
    node === nothing && throw(
        ArgumentError("Geometric scene node id $nid is absent from the scene MTG."),
    )
    node
end

function _light_query_attribute_pairs(attributes)
    attributes === nothing && return Pair{Symbol,Any}[]
    (attributes isa NamedTuple || attributes isa AbstractDict) || throw(
        ArgumentError("`attributes` must be a named tuple or dictionary."),
    )
    Pair{Symbol,Any}[Symbol(key) => value for (key, value) in pairs(attributes)]
end

function _light_query_validate_attributes(scene::PlantGeom.SceneGeometry, attribute_pairs)
    isempty(attribute_pairs) && return nothing
    scene.mtg === nothing && throw(
        ArgumentError("Attribute filters require an MTG-backed scene."),
    )
    available = Set(Symbol.(MultiScaleTreeGraph.get_attributes(scene.mtg)))
    missing_names = sort!(Symbol[key for (key, _) in attribute_pairs if !(key in available)])
    isempty(missing_names) || throw(
        ArgumentError(
            "Unknown scene attribute(s): $(join(string.(missing_names), ", ")). " *
            "Available attributes include: $(join(string.(sort!(collect(available))), ", ")).",
        ),
    )
    return nothing
end

function _light_query_attribute(
    scene::PlantGeom.SceneGeometry,
    node,
    nid::Int,
    key::Symbol,
    inherit_attributes::Bool,
)
    inherit_attributes && return _inherited_attr(scene, nid, (key,), nothing)
    MultiScaleTreeGraph.attribute(node, key; default=nothing)
end

function _light_query_requested_ids(node_ids, candidates::Vector{Int})
    node_ids === nothing && return candidates
    requested = sort!(unique!(Int[Int(value) for value in _light_query_values(node_ids)]))
    candidate_set = Set(candidates)
    unknown = Int[nid for nid in requested if !(nid in candidate_set)]
    isempty(unknown) || throw(
        ArgumentError(
            "Unknown or non-geometric scene node id(s): $(join(unknown, ", ")).",
        ),
    )
    requested
end

"""
    light_node_ids(scene_or_sim; filters...)

Return the sorted runtime node ids of geometric scene components matching all
requested filters. Pass the returned ids back through `node_ids=` to reuse a
selection across light metrics or timesteps.

Standard filters are `node_ids`, `source_topology_id`, `group` (or its alias
`species`), `object_id`, `symbol`, `scale`, and `type`. Scalar filters use
equality; collections use membership. `attributes` accepts a named tuple or
dictionary of exact attribute matches. Arbitrary attributes are node-local
unless `inherit_attributes=true`. The optional `where` predicate receives the
corresponding MTG node.

`group`/`species` and `object_id` use inherited scene metadata, matching the
way OPS object roots identify their descendant geometry.
"""
function light_node_ids(
    scene_or_sim;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes::Bool=false,
    where=nothing,
)
    (scene_or_sim isa PlantGeom.SceneGeometry || scene_or_sim isa LightSimulation) || throw(
        ArgumentError("The first argument must be a PlantGeom.SceneGeometry or LightSimulation."),
    )
    group !== nothing && species !== nothing && throw(
        ArgumentError("`group` and `species` are aliases; pass only one of them."),
    )
    scene = _light_query_scene(scene_or_sim)
    candidates = _light_query_requested_ids(
        node_ids,
        _light_query_candidate_ids(scene_or_sim),
    )
    attribute_pairs = _light_query_attribute_pairs(attributes)
    _light_query_validate_attributes(scene, attribute_pairs)

    needs_mtg = source_topology_id !== nothing || group !== nothing || species !== nothing || object_id !== nothing ||
                symbol !== nothing || scale !== nothing || type !== nothing ||
                !isempty(attribute_pairs) || where !== nothing
    needs_mtg && scene.mtg === nothing && throw(
        ArgumentError("Metadata filters require an MTG-backed scene."),
    )
    where !== nothing && !isempty(candidates) &&
        !applicable(where, _light_query_node(scene, first(candidates))) && throw(
        ArgumentError("`where` must be a callable predicate accepting an MTG node."),
    )

    selected_group = group === nothing ? species : group
    selected = Int[]
    for nid in candidates
        node = needs_mtg ? _light_query_node(scene, nid) : nothing
        selected_group !== nothing &&
            !_light_query_string_matches(_scene_group(scene, nid, ""), selected_group) && continue
        source_topology_id !== nothing &&
            !_light_query_matches(
                _scene_source_topology_id(scene, nid, nid),
                source_topology_id,
            ) && continue
        object_id !== nothing &&
            !_light_query_matches(_scene_object_id(scene, nid, -1), object_id) && continue
        symbol !== nothing &&
            !_light_query_symbol_matches(MultiScaleTreeGraph.symbol(node), symbol) && continue
        scale !== nothing &&
            !_light_query_matches(MultiScaleTreeGraph.scale(node), scale) && continue
        type !== nothing &&
            !_light_query_string_matches(_light_query_type(scene_or_sim, nid), type) && continue

        matches_attributes = all(attribute_pairs) do (key, expected)
            actual = _light_query_attribute(scene, node, nid, key, inherit_attributes)
            _light_query_matches(actual, expected)
        end === true
        matches_attributes || continue
        if where !== nothing
            keep = where(node)
            keep isa Bool || throw(
                ArgumentError("`where` must return Bool, got $(typeof(keep))."),
            )
            keep || continue
        end
        push!(selected, nid)
    end
    selected
end

function _light_metadata_attribute_table(metadata::LightNodeMetadata, inherit_attributes::Bool)
    inherit_attributes ? metadata.inherited_attributes : metadata.attributes
end

function _light_metadata_selected_indices(
    metadata::LightNodeMetadata;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes::Bool=false,
    where=nothing,
)
    group !== nothing && species !== nothing && throw(
        ArgumentError("`group` and `species` are aliases; pass only one of them."),
    )
    where === nothing || throw(ArgumentError(
        "`where` requires a live scene MTG and cannot be applied to a retained metadata snapshot. " *
        "Capture the required scalar attribute with `node_metadata_attributes` instead.",
    ))
    attribute_pairs = _light_query_attribute_pairs(attributes)
    attribute_table = _light_metadata_attribute_table(metadata, inherit_attributes)
    available_attributes = Set(propertynames(attribute_table))
    missing_attributes = Symbol[key for (key, _) in attribute_pairs if !(key in available_attributes)]
    isempty(missing_attributes) || throw(ArgumentError(
        "Attribute(s) were not retained in this light result: " *
        join(string.(missing_attributes), ", ") *
        ". Add them to `LightOptions(node_metadata_attributes=...)` before running the step.",
    ))

    candidate_indices = collect(eachindex(metadata.node_id))
    if node_ids !== nothing
        requested = Set(Int(value) for value in _light_query_values(node_ids))
        available = Set(metadata.node_id)
        unknown = sort!(collect(setdiff(requested, available)))
        isempty(unknown) || throw(ArgumentError(
            "Unknown or non-geometric scene node id(s) in retained metadata: $(join(unknown, ", ")).",
        ))
        filter!(i -> metadata.node_id[i] in requested, candidate_indices)
    end

    selected_group = group === nothing ? species : group
    selected = Int[]
    for i in candidate_indices
        source_topology_id !== nothing &&
            !_light_query_matches(metadata.source_topology_id[i], source_topology_id) && continue
        selected_group !== nothing &&
            !_light_query_string_matches(metadata.group[i], selected_group) && continue
        object_id !== nothing && !_light_query_matches(metadata.object_id[i], object_id) && continue
        symbol !== nothing && !_light_query_symbol_matches(metadata.symbol[i], symbol) && continue
        scale !== nothing && !_light_query_matches(metadata.scale[i], scale) && continue
        type !== nothing && !_light_query_string_matches(metadata.type[i], type) && continue
        matches_attributes = all(attribute_pairs) do (key, expected)
            _light_query_matches(getproperty(attribute_table, key)[i], expected)
        end === true
        matches_attributes || continue
        push!(selected, i)
    end
    selected
end

function light_node_ids(
    step::LightStepResult;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes::Bool=false,
    where=nothing,
)
    metadata = step.node_metadata
    metadata === nothing && throw(ArgumentError(
        "This LightStepResult has no retained node metadata. Re-run with " *
        "LightOptions(store_node_metadata=true), or query it with its original scene.",
    ))
    indices = _light_metadata_selected_indices(
        metadata;
        node_ids=node_ids,
        source_topology_id=source_topology_id,
        group=group,
        species=species,
        object_id=object_id,
        symbol=symbol,
        scale=scale,
        type=type,
        attributes=attributes,
        inherit_attributes=inherit_attributes,
        where=where,
    )
    metadata.node_id[indices]
end

function _light_query_table(names::Vector{Symbol}, columns::Vector{<:AbstractVector})
    NamedTuple{Tuple(names)}(Tuple(columns))
end

function _light_query_narrow(values::Vector{Any})
    isempty(values) && return Any[]
    T = foldl((acc, value) -> typejoin(acc, typeof(value)), values; init=Union{})
    T === Union{} && return Any[]
    T[values...]
end

function _light_query_materialize(table, sink)
    sink === nothing && return table
    applicable(sink, table) && return sink(table)
    materializer = try
        Tables.materializer(sink)
    catch
        try
            Tables.materializer(typeof(sink))
        catch
            nothing
        end
    end
    materializer !== nothing && return materializer(table)
    throw(
        ArgumentError(
            "Unsupported `sink`. Pass a Tables.jl materializer such as DataFrame or a callable accepting a table.",
        ),
    )
end

function _light_query_steps(step::LightStepResult)
    LightStepResult[step], true
end

function _light_query_steps(steps::AbstractVector{<:LightStepResult})
    collect(steps), false
end

function _light_query_steps(data)
    throw(ArgumentError("Light data must be a LightStepResult or a vector of LightStepResult values."))
end

function _light_query_detailed_table(
    scene_or_sim,
    steps,
    metric::Symbol,
    selected_ids::Vector{Int},
    attribute_pairs,
    inherit_attributes::Bool,
)
    scene = _light_query_scene(scene_or_sim)
    scene.mtg === nothing && throw(
        ArgumentError("Scene-aware light value tables require an MTG-backed scene."),
    )
    output_keys = _light_query_output_keys(scene_or_sim, selected_ids)
    attribute_names = Symbol[first(pair) for pair in attribute_pairs]
    reserved = intersect(Set(attribute_names), Set((_LIGHT_QUERY_BASE_COLUMNS..., :value)))
    isempty(reserved) || throw(
        ArgumentError(
            "Attribute output name(s) conflict with standard query columns: " *
            join(string.(sort!(collect(reserved))), ", "),
        ),
    )

    names = Symbol[_LIGHT_QUERY_BASE_COLUMNS...]
    append!(names, attribute_names)
    push!(names, :value)
    columns = AbstractVector[
        Int[], Int[], Int[], Int[], Int[], Int[], String[], String[], Symbol[], Int[],
    ]
    for _ in attribute_names
        push!(columns, Any[])
    end
    push!(columns, Float64[])

    for (step_number, step) in enumerate(steps)
        values = light_metric_values(step, metric)
        absent = Int[nid for nid in selected_ids if !haskey(values, nid)]
        isempty(absent) || throw(
            ArgumentError(
                "Selected scene node id(s) are absent from the light result: $(join(absent, ", ")). " *
                "Ensure the result was computed from this scene and simulation configuration.",
            ),
        )
        for nid in selected_ids
            node = _light_query_node(scene, nid)
            object_id = _scene_object_id(scene, nid, -1)
            source_topology_id = _scene_source_topology_id(scene, nid, nid)
            item_id, component_id = get(output_keys, nid, (object_id, source_topology_id))
            base_values = (
                step_number,
                nid,
                source_topology_id,
                object_id,
                item_id,
                component_id,
                _scene_group(scene, nid, ""),
                _light_query_type(scene_or_sim, nid),
                Symbol(MultiScaleTreeGraph.symbol(node)),
                Int(MultiScaleTreeGraph.scale(node)),
            )
            for (column, value) in zip(columns, base_values)
                push!(column, value)
            end
            offset = length(_LIGHT_QUERY_BASE_COLUMNS)
            for (i, (key, _)) in enumerate(attribute_pairs)
                push!(
                    columns[offset + i],
                    _light_query_attribute(scene, node, nid, key, inherit_attributes),
                )
            end
            push!(last(columns), Float64(values[nid]))
        end
    end
    for i in eachindex(attribute_names)
        column_index = length(_LIGHT_QUERY_BASE_COLUMNS) + i
        columns[column_index] = _light_query_narrow(columns[column_index])
    end
    _light_query_table(names, columns)
end

function _light_query_snapshot_table(
    steps,
    metric::Symbol,
    attribute_pairs;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes::Bool=false,
    where=nothing,
)
    attribute_names = Symbol[first(pair) for pair in attribute_pairs]
    reserved = intersect(Set(attribute_names), Set((_LIGHT_QUERY_BASE_COLUMNS..., :value)))
    isempty(reserved) || throw(ArgumentError(
        "Attribute output name(s) conflict with standard query columns: " *
        join(string.(sort!(collect(reserved))), ", "),
    ))

    names = Symbol[_LIGHT_QUERY_BASE_COLUMNS...]
    append!(names, attribute_names)
    push!(names, :value)
    columns = AbstractVector[
        Int[], Int[], Int[], Int[], Int[], Int[], String[], String[], Symbol[], Int[],
    ]
    for _ in attribute_names
        push!(columns, Any[])
    end
    push!(columns, Float64[])

    for (step_number, step) in enumerate(steps)
        metadata = step.node_metadata
        metadata === nothing && throw(ArgumentError(
            "Step $step_number has no retained node metadata. Dynamic-series queries require " *
            "LightOptions(store_node_metadata=true) for every step.",
        ))
        selected_indices = _light_metadata_selected_indices(
            metadata;
            node_ids=node_ids,
            source_topology_id=source_topology_id,
            group=group,
            species=species,
            object_id=object_id,
            symbol=symbol,
            scale=scale,
            type=type,
            attributes=attributes,
            inherit_attributes=inherit_attributes,
            where=where,
        )
        attribute_table = _light_metadata_attribute_table(metadata, inherit_attributes)
        values = light_metric_values(step, metric)
        for i in selected_indices
            nid = metadata.node_id[i]
            haskey(values, nid) || throw(ArgumentError(
                "Retained scene node id $nid is absent from light result step $step_number.",
            ))
            base_values = (
                step_number,
                nid,
                metadata.source_topology_id[i],
                metadata.object_id[i],
                metadata.item_id[i],
                metadata.component_id[i],
                metadata.group[i],
                metadata.type[i],
                metadata.symbol[i],
                metadata.scale[i],
            )
            for (column, value) in zip(columns, base_values)
                push!(column, value)
            end
            offset = length(_LIGHT_QUERY_BASE_COLUMNS)
            for (j, (key, _)) in enumerate(attribute_pairs)
                push!(columns[offset + j], getproperty(attribute_table, key)[i])
            end
            push!(last(columns), Float64(values[nid]))
        end
    end
    for i in eachindex(attribute_names)
        column_index = length(_LIGHT_QUERY_BASE_COLUMNS) + i
        columns[column_index] = _light_query_narrow(columns[column_index])
    end
    _light_query_table(names, columns)
end

function _light_query_by_columns(by)
    by === nothing && return Symbol[]
    values = _light_query_values(by)
    Symbol[Symbol(value) for value in values]
end

function _light_query_reduce(reducer, values, label)
    try
        reducer(values)
    catch err
        isempty(values) || rethrow()
        throw(
            ArgumentError(
                "Reducer $(repr(reducer)) has no result for the empty selection at $label.",
            ),
        )
    end
end

function _light_query_grouped_table(table, reducer, by_columns, is_single::Bool)
    available = Set(propertynames(table))
    invalid = Symbol[column for column in by_columns if !(column in available) || column == :value]
    isempty(invalid) || throw(
        ArgumentError(
            "Unknown or non-groupable `by` column(s): $(join(string.(invalid), ", ")).",
        ),
    )
    group_columns = is_single ? copy(by_columns) : unique!(Symbol[:step_number; by_columns])
    groups = Dict{Tuple,Vector{Float64}}()
    order = Tuple[]
    for row in Tables.rows(table)
        key = Tuple(getproperty(row, column) for column in group_columns)
        if !haskey(groups, key)
            groups[key] = Float64[]
            push!(order, key)
        end
        push!(groups[key], Float64(row.value))
    end

    columns = AbstractVector[Any[] for _ in group_columns]
    reduced = Any[]
    for key in order
        for (i, value) in enumerate(key)
            push!(columns[i], value)
        end
        push!(reduced, _light_query_reduce(reducer, groups[key], "group $(repr(key))"))
    end
    for i in eachindex(columns)
        columns[i] isa Vector{Any} && (columns[i] = _light_query_narrow(columns[i]))
    end
    push!(columns, _light_query_narrow(reduced))
    _light_query_table(Symbol[group_columns; :value], columns)
end

"""
    light_metric_values(scene_or_sim, step_or_series, metric; filters..., reduce=nothing, by=nothing, sink=nothing)

Return a filtered, Tables.jl-compatible column table for one light metric.
`metric` accepts package selectors such as `:absorbed_par_energy` and historical
ARCHIMED names such as `:Ra_PAR_q`.

The node filters are the same as [`light_node_ids`](@ref). Without `reduce`, the
result contains `step_number`, scene and output node identifiers, group/type and
MTG metadata, requested attribute columns, and `value`. A series is returned in
long form with one-based step numbers.

With `reduce=sum`, a single step returns a scalar and a series returns one value
per timestep. Pass `by` to obtain a grouped table instead; series grouping always
retains `step_number`. `sink` optionally materializes table results with a
Tables.jl-compatible sink such as DataFrames.DataFrame.

Summing component energy metrics is physically meaningful. Component fluxes are
area-normalized, so summing them generally is not a scene-scale energy total.
"""
function light_metric_values(
    scene_or_sim,
    data,
    metric::Symbol;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes::Bool=false,
    where=nothing,
    reduce=nothing,
    by=nothing,
    sink=nothing,
)
    by !== nothing && reduce === nothing && throw(
        ArgumentError("`by` requires a reducer such as `reduce=sum`."),
    )
    steps, is_single = _light_query_steps(data)
    attribute_pairs = _light_query_attribute_pairs(attributes)
    all_snapshots = !isempty(steps) && all(step -> step.node_metadata !== nothing, steps)
    snapshots_have_attributes = all_snapshots && all(steps) do step
        metadata = step.node_metadata::LightNodeMetadata
        table = _light_metadata_attribute_table(metadata, inherit_attributes)
        all(
            pair -> any(name -> isequal(name, first(pair)), propertynames(table)),
            attribute_pairs,
        ) === true
    end === true
    if !is_single && all_snapshots && where !== nothing
        throw(ArgumentError(
            "`where` cannot query a series through retained metadata. Capture the predicate's " *
            "scalar inputs with `node_metadata_attributes` and use `attributes` instead.",
        ))
    end
    if !is_single && all_snapshots && !snapshots_have_attributes
        missing_names = join(string.(first.(attribute_pairs)), ", ")
        throw(ArgumentError(
            "Series query attribute(s) were not retained: $missing_names. Add them to " *
            "LightOptions(node_metadata_attributes=...) before running the series.",
        ))
    end
    use_snapshots = all_snapshots && where === nothing && snapshots_have_attributes
    table = if use_snapshots
        _light_query_snapshot_table(
            steps,
            metric,
            attribute_pairs;
            node_ids=node_ids,
            source_topology_id=source_topology_id,
            group=group,
            species=species,
            object_id=object_id,
            symbol=symbol,
            scale=scale,
            type=type,
            attributes=attributes,
            inherit_attributes=inherit_attributes,
            where=where,
        )
    else
        selected_ids = light_node_ids(
            scene_or_sim;
            node_ids=node_ids,
            source_topology_id=source_topology_id,
            group=group,
            species=species,
            object_id=object_id,
            symbol=symbol,
            scale=scale,
            type=type,
            attributes=attributes,
            inherit_attributes=inherit_attributes,
            where=where,
        )
        _light_query_detailed_table(
            scene_or_sim,
            steps,
            metric,
            selected_ids,
            attribute_pairs,
            inherit_attributes,
        )
    end

    reduce === nothing && return _light_query_materialize(table, sink)
    by_columns = _light_query_by_columns(by)
    if isempty(by_columns)
        sink === nothing || throw(
            ArgumentError("`sink` is only supported when the query returns a table."),
        )
        if is_single
            return _light_query_reduce(reduce, table.value, "step 1")
        end
        return [
            _light_query_reduce(
                reduce,
                Float64[table.value[i] for i in eachindex(table.value) if table.step_number[i] == step_number],
                "step $step_number",
            ) for step_number in eachindex(steps)
        ]
    end
    grouped = _light_query_grouped_table(table, reduce, by_columns, is_single)
    _light_query_materialize(grouped, sink)
end
