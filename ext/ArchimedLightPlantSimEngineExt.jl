module ArchimedLightPlantSimEngineExt

using ArchimedLight
import PlantGeom
import PlantSimEngine

PlantSimEngine.@process "light_interception" verbose = false

const _COUPLING_OUTPUT_NAMES = (
    :Ri_PAR_f,
    :Ri_NIR_f,
    :Ra_PAR_f,
    :Ra_NIR_f,
    :Ra_SW_f,
    :aPPFD,
    :radiative_mesh_area,
)

const _FULL_OUTPUT_NAMES = (
    :Ri_PAR_0_f,
    :Ri_NIR_0_f,
    :Ri_PAR_f,
    :Ri_NIR_f,
    :Ri_PAR_0_q,
    :Ri_NIR_0_q,
    :Ri_PAR_q,
    :Ri_NIR_q,
    :Ra_PAR_0_f,
    :Ra_NIR_0_f,
    :Ra_PAR_f,
    :Ra_NIR_f,
    :Ra_PAR_0_q,
    :Ra_NIR_0_q,
    :Ra_PAR_q,
    :Ra_NIR_q,
    :Ra_SW_f,
    :aPPFD,
    :radiative_mesh_area,
)

@inline _output_names(::Val{:coupling}) = _COUPLING_OUTPUT_NAMES
@inline _output_names(::Val{:full}) = _FULL_OUTPUT_NAMES

Base.@constprop :aggressive function _output_schema(schema::Symbol)
    schema === :coupling && return Val(:coupling)
    schema === :full && return Val(:full)
    throw(ArgumentError(
        "`output_schema` must be `:coupling` or `:full`; got $(repr(schema)).",
    ))
end

const _COUPLING_OUTPUTS = (
    Ri_PAR_f=PlantSimEngine.Default(0.0),
    Ri_NIR_f=PlantSimEngine.Default(0.0),
    Ra_PAR_f=PlantSimEngine.Default(0.0),
    Ra_NIR_f=PlantSimEngine.Default(0.0),
    Ra_SW_f=PlantSimEngine.Default(0.0),
    aPPFD=PlantSimEngine.Default(0.0),
    radiative_mesh_area=PlantSimEngine.Default(0.0),
)

const _FULL_OUTPUTS = (
    Ri_PAR_0_f=PlantSimEngine.Default(0.0),
    Ri_NIR_0_f=PlantSimEngine.Default(0.0),
    Ri_PAR_f=PlantSimEngine.Default(0.0),
    Ri_NIR_f=PlantSimEngine.Default(0.0),
    Ri_PAR_0_q=PlantSimEngine.Default(0.0),
    Ri_NIR_0_q=PlantSimEngine.Default(0.0),
    Ri_PAR_q=PlantSimEngine.Default(0.0),
    Ri_NIR_q=PlantSimEngine.Default(0.0),
    Ra_PAR_0_f=PlantSimEngine.Default(0.0),
    Ra_NIR_0_f=PlantSimEngine.Default(0.0),
    Ra_PAR_f=PlantSimEngine.Default(0.0),
    Ra_NIR_f=PlantSimEngine.Default(0.0),
    Ra_PAR_0_q=PlantSimEngine.Default(0.0),
    Ra_NIR_0_q=PlantSimEngine.Default(0.0),
    Ra_PAR_q=PlantSimEngine.Default(0.0),
    Ra_NIR_q=PlantSimEngine.Default(0.0),
    Ra_SW_f=PlantSimEngine.Default(0.0),
    aPPFD=PlantSimEngine.Default(0.0),
    radiative_mesh_area=PlantSimEngine.Default(0.0),
)

@inline _output_declarations(::Val{:coupling}) = _COUPLING_OUTPUTS
@inline _output_declarations(::Val{:full}) = _FULL_OUTPUTS

Base.@constprop :aggressive function ArchimedLight.archimed_light_outputs(
    schema::Symbol=:coupling,
)
    schema === :coupling && return _COUPLING_OUTPUTS
    schema === :full && return _FULL_OUTPUTS
    throw(ArgumentError(
        "`schema` must be `:coupling` or `:full`; got $(repr(schema)).",
    ))
end

const _LightOutputColumns = NamedTuple{
    _FULL_OUTPUT_NAMES,
    NTuple{length(_FULL_OUTPUT_NAMES),Vector{Float64}},
}

struct _PreparedDistributedLight
    scene::PlantGeom.SceneGeometry
    scene_revision::Int
    result_metadata::ArchimedLight.LightComponentMetadata
    selected_metadata::ArchimedLight.LightComponentMetadata
    owner_map::AbstractDict{
        PlantGeom.SourceOwnerKey,
        PlantSimEngine.ObjectId,
    }
    runtime_model::WeakRef
    execution_id::PlantSimEngine.ObjectId
    model_revision::Int
    environment_revision::Int
    target_ids::Vector{PlantSimEngine.ObjectId}
    component_to_target::Vector{Int}
    target_radiative_area::Vector{Float64}
    columns::_LightOutputColumns
end

mutable struct _ArchimedLightModel{S,SP,SR,OR,Schema,Group} <:
               AbstractLight_InterceptionModel
    simulation::S
    scene_provider::SP
    source_roots::SR
    object_resolver::OR
    par_energy_to_photon::Float64
    runtime_cache::Union{Nothing,_PreparedDistributedLight}
end

function ArchimedLight.ArchimedLightModel(
    simulation::ArchimedLight.LightSimulation;
    output_schema::Symbol=:coupling,
    output_group::Symbol=:organs,
    scene_provider=nothing,
    source_roots=nothing,
    object_resolver=nothing,
    par_energy_to_photon::Real=4.57,
)
    schema = _output_schema(output_schema)
    source_roots !== nothing && object_resolver !== nothing && throw(ArgumentError(
        "`source_roots` and `object_resolver` are mutually exclusive.",
    ))
    factor = ArchimedLight._validated_par_energy_to_photon(par_energy_to_photon)
    return _archimed_light_model(
        simulation,
        scene_provider,
        source_roots,
        object_resolver,
        factor,
        schema,
        Val(output_group),
    )
end

function _archimed_light_model(
    simulation::S,
    scene_provider::SP,
    source_roots::SR,
    object_resolver::OR,
    factor::Float64,
    ::Val{Schema},
    ::Val{Group},
) where {S,SP,SR,OR,Schema,Group}
    return _ArchimedLightModel{
        S,
        SP,
        SR,
        OR,
        Schema,
        Group,
    }(
        simulation,
        scene_provider,
        source_roots,
        object_resolver,
        factor,
        nothing,
    )
end

PlantSimEngine.inputs_(::_ArchimedLightModel) = NamedTuple()
PlantSimEngine.outputs_(::_ArchimedLightModel) = NamedTuple()
PlantSimEngine.environment_inputs_(::_ArchimedLightModel) = (
    sun_azimuth_deg=0.0,
    sun_elevation_deg=0.0,
    Ri_PAR_f=0.0,
    Ri_NIR_f=0.0,
    direct_fraction=0.0,
)

@inline function _require_finite(name::Symbol, value)
    value isa Real || throw(ArgumentError(
        "ArchimedLight environment `$name` must be a real number; got $(typeof(value)).",
    ))
    result = Float64(value)
    isfinite(result) || throw(ArgumentError(
        "ArchimedLight environment `$name` must be finite; got $value.",
    ))
    return result
end

@inline function _require_closed_interval(name::Symbol, value, lower, upper)
    result = _require_finite(name, value)
    lower <= result <= upper || throw(ArgumentError(
        "ArchimedLight environment `$name` must be in [$lower, $upper]; got $value.",
    ))
    return result
end

@inline function _require_nonnegative(name::Symbol, value)
    result = _require_finite(name, value)
    result >= 0.0 || throw(ArgumentError(
        "ArchimedLight environment `$name` must be non-negative; got $value.",
    ))
    return result
end

function _validated_forcing(environment)
    hasproperty(environment, :duration) || throw(ArgumentError(
        "ArchimedLight requires PlantSimEngine's sampled environment `duration`.",
    ))
    azimuth = _require_finite(:sun_azimuth_deg, environment.sun_azimuth_deg)
    0.0 <= azimuth < 360.0 || throw(ArgumentError(
        "ArchimedLight environment `sun_azimuth_deg` must be in [0, 360); " *
        "got $(environment.sun_azimuth_deg).",
    ))
    elevation = _require_closed_interval(
        :sun_elevation_deg,
        environment.sun_elevation_deg,
        -90.0,
        90.0,
    )
    ri_par = _require_nonnegative(:Ri_PAR_f, environment.Ri_PAR_f)
    ri_nir = _require_nonnegative(:Ri_NIR_f, environment.Ri_NIR_f)
    direct = _require_closed_interval(
        :direct_fraction,
        environment.direct_fraction,
        0.0,
        1.0,
    )
    duration = try
        ArchimedLight._duration_seconds_strict(
            environment.duration;
            field_name="duration",
        )
    catch error
        error isa InterruptException && rethrow()
        throw(ArgumentError(
            "Invalid ArchimedLight environment `duration`: " *
            sprint(showerror, error),
        ))
    end
    sky = ArchimedLight.SkyState(
        azimuth,
        elevation,
        ri_par,
        ri_nir,
        direct,
        1.0 - direct,
    )
    return sky, duration
end

@inline _object_value(id::PlantSimEngine.ObjectId) = id.value

function _copy_target_ids(targets)
    ids = PlantSimEngine.object_ids(targets)
    copied = Vector{PlantSimEngine.ObjectId}(undef, length(ids))
    copyto!(copied, ids)
    return copied
end

function _validate_target_schema(targets, ::Val{Schema}) where {Schema}
    expected = _output_names(Val(Schema))
    actual = propertynames(targets.columns)
    actual == expected || throw(ArgumentError(
        "ArchimedLightModel(output_schema=$(repr(Schema))) requires the matching " *
        "`archimed_light_outputs($(repr(Schema)))` declaration. Expected variables " *
        "$expected, got $actual.",
    ))
    return nothing
end

function _source_root_entries(source_roots::Pair)
    return (source_roots,)
end

_source_root_entries(source_roots::Union{AbstractDict,NamedTuple,Base.Pairs}) =
    pairs(source_roots)
_source_root_entries(source_roots) = source_roots

_materialize_source_roots(source_roots) = collect(_source_root_entries(source_roots))

function _source_instance_ids(source_roots)
    ids = Set{Int}()
    for entry in _source_root_entries(source_roots)
        entry isa Pair || throw(ArgumentError(
            "`source_roots` must be a mapping or iterable of " *
            "`instance_id => exact_mtg_root` pairs.",
        ))
        raw_id = first(entry)
        raw_id isa Integer || throw(ArgumentError(
            "Source-root instance ids must be integers; got $(typeof(raw_id)).",
        ))
        id = Int(raw_id)
        id > 0 || throw(ArgumentError(
            "Source-root instance ids must be positive; got $id.",
        ))
        id in ids && throw(ArgumentError(
            "`source_roots` contains source instance $id more than once.",
        ))
        push!(ids, id)
    end
    return ids
end

function _resolve_source_roots(source_roots, context)
    source_roots === nothing && return nothing
    if applicable(source_roots, context)
        roots = source_roots(context)
        roots === nothing && throw(ArgumentError(
            "`source_roots(context)` returned `nothing`; expected a source-root mapping.",
        ))
        return _materialize_source_roots(roots)
    end
    return _materialize_source_roots(source_roots)
end

function _resolve_scene_from_provider(provider, context)
    provider === nothing && throw(ArgumentError(
        "PlantSimEngine topology/environment bindings changed without a coordinated " *
        "ArchimedLight scene refresh. Pass `scene_provider=context -> scene` to " *
        "ArchimedLightModel, or call `update_scene!` with the rebuilt scene before " *
        "continuing the simulation.",
    ))
    applicable(provider, context) || throw(ArgumentError(
        "`scene_provider` must be callable with the PlantSimEngine RunContext; " *
        "got $(typeof(provider)).",
    ))
    scene = provider(context)
    scene isa PlantGeom.SceneGeometry || throw(ArgumentError(
        "`scene_provider(context)` must return PlantGeom.SceneGeometry; " *
        "got $(typeof(scene)).",
    ))
    return scene
end

function _owner_keys_for_mapping(metadata, roots)
    owners = metadata.source_owner
    if roots === nothing
        return collect(owners)
    end
    instance_ids = _source_instance_ids(roots)
    selected = PlantGeom.SourceOwnerKey[]
    sizehint!(selected, length(owners))
    for owner in owners
        owner.source_instance_id in instance_ids && push!(selected, owner)
    end
    isempty(selected) && throw(ArgumentError(
        "No ArchimedLight component belongs to an instance declared by `source_roots`.",
    ))
    return selected
end

function _selected_component_mapping(metadata, owner_map, target_ids)
    target_index = Dict{PlantSimEngine.ObjectId,Int}(
        id => i for (i, id) in pairs(target_ids)
    )
    length(target_index) == length(target_ids) || throw(ArgumentError(
        "The ArchimedLight distributed-output selector resolved duplicate ObjectIds.",
    ))

    node_ids = Int[]
    owner_keys = PlantGeom.SourceOwnerKey[]
    component_areas = Float64[]
    component_to_target = Int[]
    sizehint!(node_ids, length(metadata.node_id))
    sizehint!(owner_keys, length(metadata.node_id))
    sizehint!(component_areas, length(metadata.node_id))
    sizehint!(component_to_target, length(metadata.node_id))
    target_area = zeros(Float64, length(target_ids))

    for i in eachindex(metadata.node_id)
        owner = metadata.source_owner[i]
        haskey(owner_map, owner) || continue
        object_id = owner_map[owner]
        target_position = get(target_index, object_id, 0)
        target_position == 0 && continue
        area = metadata.radiative_area[i]
        push!(node_ids, metadata.node_id[i])
        push!(owner_keys, owner)
        push!(component_areas, area)
        push!(component_to_target, target_position)
        @inbounds target_area[target_position] += area
    end

    missing_positions = findall(iszero, target_area)
    isempty(missing_positions) || throw(ArgumentError(
        "ArchimedLight cannot provide exact distributed-output coverage. " *
        "No selected radiative component maps to target ObjectId(s): " *
        join((_object_value(target_ids[i]) for i in missing_positions), ", ") * ".",
    ))
    selected_metadata = ArchimedLight.LightComponentMetadata(
        node_ids,
        owner_keys,
        component_areas,
    )
    return selected_metadata, component_to_target, target_area
end

function _new_output_columns(target_area)
    return merge(
        ArchimedLight._new_light_metric_columns(length(target_area)),
        (radiative_mesh_area=target_area,),
    )
end

function _prepare_distributed_light(
    model::_ArchimedLightModel,
    context,
    targets,
    scene,
    execution_id,
    model_revision::Int,
    environment_revision::Int,
)
    light_cache = ArchimedLight._ensure_light_cache!(model.simulation)
    metadata = light_cache.component_metadata
    metadata === nothing && throw(ArgumentError(
        "ArchimedLight cannot publish identity-keyed organ outputs because the " *
        "prepared scene has incomplete source-owner metadata.",
    ))

    roots = _resolve_source_roots(model.source_roots, context)
    owner_keys = _owner_keys_for_mapping(metadata, roots)
    owner_map = PlantGeom.compile_source_owner_map(
        scene,
        context;
        owner_keys=owner_keys,
        source_roots=roots,
        object_resolver=model.object_resolver,
    )
    target_ids = _copy_target_ids(targets)
    selected_metadata, component_to_target, target_area =
        _selected_component_mapping(metadata, owner_map, target_ids)
    columns = _new_output_columns(target_area)
    return _PreparedDistributedLight(
        scene,
        PlantGeom.scene_version(scene),
        metadata,
        selected_metadata,
        owner_map,
        WeakRef(PlantSimEngine.runtime_model(context)),
        execution_id,
        model_revision,
        environment_revision,
        target_ids,
        component_to_target,
        target_area,
        columns,
    )
end

function _run_distributed_light!(
    model::_ArchimedLightModel{S,SP,SR,OR,Schema,Group},
    ::Nothing,
    environment,
    context,
) where {S,SP,SR,OR,Schema,Group}
    targets = PlantSimEngine.output_targets(context, Val(Group))
    _validate_target_schema(targets, Val(Schema))
    sky, duration = _validated_forcing(environment)
    runtime = PlantSimEngine.runtime_model(context)
    model_revision = PlantSimEngine.Advanced.model_revision(runtime)
    environment_revision = PlantSimEngine.Advanced.environment_revision(runtime)
    execution_id = context.object_id
    cache = _prepare_distributed_light(
        model,
        context,
        targets,
        model.simulation.scene,
        execution_id,
        model_revision,
        environment_revision,
    )
    model.runtime_cache = cache
    return _run_light_and_publish!(model, cache, targets, sky, duration)
end

function _prepared_distributed_light_cached!(
    model::_ArchimedLightModel,
    context,
    targets,
    cache::_PreparedDistributedLight,
    execution_id,
    model_revision::Int,
    environment_revision::Int,
)
    runtime = PlantSimEngine.runtime_model(context)
    if cache.runtime_model.value !== runtime
        throw(ArgumentError(
            "One ArchimedLightModel kernel instance cannot be shared by distinct " *
            "PlantSimEngine CompositeModels. Construct one kernel per Scene " *
            "application and runtime model.",
        ))
    end
    if cache.execution_id != execution_id
        throw(ArgumentError(
            "One ArchimedLightModel kernel instance cannot execute on multiple Scene " *
            "objects. It first ran on ObjectId $(_object_value(cache.execution_id)) " *
            "and is now being reused on $(_object_value(execution_id)). Construct one " *
            "kernel per Scene application.",
        ))
    end

    scene = model.simulation.scene
    coordinated_scene_change = false
    if scene !== cache.scene
        ArchimedLight.update_scene!(model.simulation, scene)
        coordinated_scene_change = true
    elseif PlantGeom.scene_version(scene) != cache.scene_revision
        # A topology-version bump does not prove that the merged radiative
        # geometry was rebuilt. Require the provider to return (or explicitly
        # affirm) the SceneGeometry prepared for this generation.
        scene = _resolve_scene_from_provider(model.scene_provider, context)
        ArchimedLight.update_scene!(model.simulation, scene)
        coordinated_scene_change = true
    end
    revisions_changed =
        model_revision != cache.model_revision ||
        environment_revision != cache.environment_revision
    if revisions_changed && !coordinated_scene_change
        scene = _resolve_scene_from_provider(model.scene_provider, context)
        # A provider may explicitly affirm that an unrelated runtime revision
        # did not change the scene. Preserve the expensive radiative cache in
        # that case while still rebuilding owner/target bindings below.
        if scene !== cache.scene ||
           PlantGeom.scene_version(scene) != cache.scene_revision
            ArchimedLight.update_scene!(model.simulation, scene)
        end
        coordinated_scene_change = true
    end

    if coordinated_scene_change
        cache = _prepare_distributed_light(
            model,
            context,
            targets,
            scene,
            execution_id,
            model_revision,
            environment_revision,
        )
        model.runtime_cache = cache
        return cache
    end

    # A direct update_models!/update_options! call invalidates ArchimedLight's
    # numerical cache without changing PSE or scene revisions.
    light_cache = ArchimedLight._ensure_light_cache!(model.simulation)
    metadata_changed = light_cache.component_metadata !== cache.result_metadata
    owner_map_current = PlantGeom.source_owner_map_iscurrent(cache.owner_map, scene, context)
    if metadata_changed || !owner_map_current
        if !owner_map_current &&
           model_revision == cache.model_revision &&
           PlantGeom.scene_version(scene) == cache.scene_revision
            throw(ArgumentError(
                "The cached ArchimedLight source-owner mapping became stale without " *
                "a coordinated scene or PlantSimEngine revision change.",
            ))
        end
        cache = _prepare_distributed_light(
            model,
            context,
            targets,
            scene,
            execution_id,
            model_revision,
            environment_revision,
        )
        model.runtime_cache = cache
    end
    return cache
end

@noinline function _refresh_distributed_light_and_run!(
    model::_ArchimedLightModel,
    cache::_PreparedDistributedLight,
    context,
    targets,
    sky,
    duration,
    execution_id,
    model_revision::Int,
    environment_revision::Int,
)
    refreshed = _prepared_distributed_light_cached!(
        model,
        context,
        targets,
        cache,
        execution_id,
        model_revision,
        environment_revision,
    )
    return _run_light_and_publish!(
        model,
        refreshed,
        targets,
        sky,
        duration,
    )
end

function _run_distributed_light!(
    model::_ArchimedLightModel{S,SP,SR,OR,Schema,Group},
    cache::_PreparedDistributedLight,
    environment,
    context,
) where {S,SP,SR,OR,Schema,Group}
    targets = PlantSimEngine.output_targets(context, Val(Group))
    _validate_target_schema(targets, Val(Schema))
    sky, duration = _validated_forcing(environment)
    runtime = PlantSimEngine.runtime_model(context)
    execution_id = context.object_id
    if cache.runtime_model.value !== runtime
        throw(ArgumentError(
            "One ArchimedLightModel kernel instance cannot be shared by distinct " *
            "PlantSimEngine CompositeModels. Construct one kernel per Scene " *
            "application and runtime model.",
        ))
    end
    if cache.execution_id != execution_id
        throw(ArgumentError(
            "One ArchimedLightModel kernel instance cannot execute on multiple Scene " *
            "objects. It first ran on ObjectId $(_object_value(cache.execution_id)) " *
            "and is now being reused on $(_object_value(execution_id)). Construct one " *
            "kernel per Scene application.",
        ))
    end

    model_revision = PlantSimEngine.Advanced.model_revision(runtime)
    environment_revision = PlantSimEngine.Advanced.environment_revision(runtime)
    scene = model.simulation.scene
    revision_changed =
        model_revision != cache.model_revision ||
        environment_revision != cache.environment_revision
    scene_changed =
        scene !== cache.scene ||
        PlantGeom.scene_version(scene) != cache.scene_revision
    if scene_changed || revision_changed
        return _refresh_distributed_light_and_run!(
            model,
            cache,
            context,
            targets,
            sky,
            duration,
            execution_id,
            model_revision,
            environment_revision,
        )
    end

    # Direct ArchimedLight model/option updates invalidate this cache without a
    # PlantSimEngine revision. Inspect the already prepared cache by identity;
    # `run_light` performs the ordinary lazy ensure below.
    light_cache = model.simulation.cache
    metadata_changed =
        light_cache === nothing ||
        light_cache.component_metadata !== cache.result_metadata
    owner_map_current = PlantGeom.source_owner_map_iscurrent(
        cache.owner_map,
        scene,
        context,
    )
    if metadata_changed || !owner_map_current
        if !owner_map_current
            throw(ArgumentError(
                "The cached ArchimedLight source-owner mapping became stale without " *
                "a coordinated scene or PlantSimEngine revision change.",
            ))
        end
        return _refresh_distributed_light_and_run!(
            model,
            cache,
            context,
            targets,
            sky,
            duration,
            execution_id,
            model_revision,
            environment_revision,
        )
    end

    return _run_light_and_publish!(model, cache, targets, sky, duration)
end

function _publish_light_step!(
    model::_ArchimedLightModel,
    cache::_PreparedDistributedLight,
    targets,
    step,
)
    step.component_metadata === cache.result_metadata || throw(ArgumentError(
        "ArchimedLight returned a result from another scene generation; no organ " *
        "status was modified.",
    ))
    ArchimedLight._fill_aggregated_metric_columns!(
        cache.columns,
        step,
        cache.selected_metadata,
        cache.component_to_target,
        cache.target_radiative_area,
    )
    ArchimedLight._derive_light_coupling_columns!(
        cache.columns,
        model.par_energy_to_photon,
    )
    PlantSimEngine.assign_outputs!(targets, cache.target_ids, cache.columns)
    return nothing
end

function _run_light_and_publish!(
    model::_ArchimedLightModel,
    cache::_PreparedDistributedLight,
    targets,
    sky,
    duration,
)
    # All forcing, identity mapping, and exact target coverage checks occur
    # before the numerical solver. Status publication remains the final step.
    step = ArchimedLight.run_light(
        model.simulation,
        sky;
        step_duration_seconds=duration,
    )
    _publish_light_step!(model, cache, targets, step)
    return nothing
end

function PlantSimEngine.run!(
    model::_ArchimedLightModel{S,SP,SR,OR,Schema,Group},
    status,
    environment,
    constants,
    context,
) where {S,SP,SR,OR,Schema,Group}
    return _run_distributed_light!(
        model,
        model.runtime_cache,
        environment,
        context,
    )
end

end
