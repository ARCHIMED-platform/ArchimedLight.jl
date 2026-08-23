const _PSE_BENCHMARK_COMPONENT_COUNTS = (1, 100, 1_000)
const _PSE_BENCHMARK_DURATION_SECONDS = 60.0
const _PSE_BENCHMARK_EXTENSION = Base.get_extension(
    ArchimedLight,
    :ArchimedLightPlantSimEngineExt,
)

_PSE_BENCHMARK_EXTENSION === nothing && error(
    "The ArchimedLight PlantSimEngine extension must be loaded for its benchmarks.",
)

function _pse_benchmark_triangle()
    edge = 0.20f0
    return GeometryBasics.Mesh(
        GeometryBasics.Point3f[
            GeometryBasics.Point3f(0, 0, 0),
            GeometryBasics.Point3f(edge, 0, 0),
            GeometryBasics.Point3f(0, edge, 0),
        ],
        GeometryBasics.TriangleFace{Int}[
            GeometryBasics.TriangleFace{Int}(1, 2, 3),
        ],
    )
end

function _pse_benchmark_scene(ncomponents::Int)
    ncomponents > 0 || throw(ArgumentError("ncomponents must be positive"))
    side = ceil(Int, sqrt(ncomponents))
    spacing = 0.35
    extent = (side - 1) * spacing + 0.20
    mesh = _pse_benchmark_triangle()
    return PlantGeom.make_scene(
        domain=(0.0, 0.0, extent, extent),
        source_path="benchmark_plantsimengine_$ncomponents",
    ) do builder
        for component in 1:ncomponents
            row, column = fldmod(component - 1, side)
            PlantGeom.add_object!(
                builder,
                mesh;
                group="plants",
                type="Leaf",
                id=component,
                at=(column * spacing, row * spacing, 1.0),
            )
        end
    end
end

function _pse_benchmark_models()
    return ArchimedLight.models_for(
        "plants" => (
            "Leaf" => ArchimedLight.translucent(par=0.15, nir=0.30),
        ),
    )
end

function _pse_benchmark_options()
    return ArchimedLight.LightOptions(
        turtle_sectors=1,
        all_in_turtle=false,
        scattering=false,
        toricity=false,
        pixel_size=0.05,
        cache_radiation=true,
    )
end

function _pse_benchmark_environment()
    return (
        duration=Dates.Minute(1),
        sun_azimuth_deg=180.0,
        sun_elevation_deg=90.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        direct_fraction=1.0,
    )
end

function _pse_benchmark_sky(environment)
    return ArchimedLight.SkyState(
        environment.sun_azimuth_deg,
        environment.sun_elevation_deg,
        environment.Ri_PAR_f,
        environment.Ri_NIR_f,
        environment.direct_fraction,
        1.0 - environment.direct_fraction,
    )
end

_pse_benchmark_leaf_id(index::Int) = Symbol(:leaf_, index)

function _pse_benchmark_resolver(destination_ids)
    return owner -> begin
        source = owner.source_instance_id
        1 <= source <= length(destination_ids) || throw(ArgumentError(
            "No benchmark destination for source instance $source.",
        ))
        return @inbounds destination_ids[source]
    end
end

function _pse_benchmark_objects(destination_ids)
    objects = PlantSimEngine.Object[
        PlantSimEngine.Object(:scene; scale=:Scene),
    ]
    sizehint!(objects, length(destination_ids) + 1)
    for id in destination_ids
        push!(
            objects,
            PlantSimEngine.Object(id; scale=:Leaf, parent=:scene),
        )
    end
    return objects
end

function _pse_benchmark_execution_target(simulation, kernel)
    matches = Any[]
    for batch in simulation.execution_plan.batches
        hasproperty(batch, :targets) || continue
        for target in batch.targets
            target.model === kernel && push!(matches, target)
        end
    end
    return only(matches)
end

function _pse_benchmark_model_case(
    scene,
    models,
    options,
    environment,
    destination_ids,
    schema::Symbol;
    resolver_ids=destination_ids,
    scene_provider=nothing,
)
    light = ArchimedLight.LightSimulation(scene, models; options=options)
    kernel = ArchimedLight.ArchimedLightModel(
        light;
        output_schema=schema,
        scene_provider=scene_provider,
        object_resolver=_pse_benchmark_resolver(resolver_ids),
    )
    application = PlantSimEngine.ModelSpec(
        kernel;
        name=Symbol(:archimed_light_, schema),
        on=PlantSimEngine.One(scale=:Scene),
        outputs_to=(
            organs=PlantSimEngine.OutputTo(
                PlantSimEngine.Many(
                    scale=:Leaf,
                    within=PlantSimEngine.SceneScope(),
                );
                vars=ArchimedLight.archimed_light_outputs(schema),
            ),
        ),
    )
    runtime = PlantSimEngine.CompositeModel(
        _pse_benchmark_objects(destination_ids)...;
        applications=(application,),
        environment=environment,
    )
    simulation = PlantSimEngine.run!(runtime; outputs=:none)
    target = _pse_benchmark_execution_target(simulation, kernel)
    return (
        light=light,
        kernel=kernel,
        runtime=runtime,
        simulation=simulation,
        target=target,
        context=target.context,
        environment=environment,
    )
end

function _pse_benchmark_fixture(ncomponents::Int)
    scene = _pse_benchmark_scene(ncomponents)
    models = _pse_benchmark_models()
    options = _pse_benchmark_options()
    environment = _pse_benchmark_environment()
    sky = _pse_benchmark_sky(environment)
    destination_ids = [_pse_benchmark_leaf_id(i) for i in 1:ncomponents]
    direct_light = ArchimedLight.LightSimulation(scene, models; options=options)
    coupling = _pse_benchmark_model_case(
        scene,
        models,
        options,
        environment,
        destination_ids,
        :coupling,
    )
    full = _pse_benchmark_model_case(
        scene,
        models,
        options,
        environment,
        destination_ids,
        :full,
    )
    return (
        ncomponents=ncomponents,
        sky=sky,
        direct_light=direct_light,
        coupling=coupling,
        full=full,
    )
end

@noinline function _pse_benchmark_direct_step!(fixture)
    ArchimedLight.run_light(
        fixture.direct_light,
        fixture.sky;
        step_duration_seconds=_PSE_BENCHMARK_DURATION_SECONDS,
    )
    return nothing
end

@noinline function _pse_benchmark_kernel_step!(case)
    PlantSimEngine.run!(
        case.kernel,
        case.target.status,
        case.environment,
        case.simulation.constants,
        case.context,
    )
    return nothing
end

@noinline function _pse_benchmark_continue!(case)
    PlantSimEngine.continue!(case.simulation; steps=1)
    return nothing
end

@noinline function _pse_benchmark_publish!(case, step)
    cache = case.kernel.runtime_cache
    cache === nothing && error("The benchmark kernel has no prepared cache.")
    targets = PlantSimEngine.output_targets(case.context, Val(:organs))
    _PSE_BENCHMARK_EXTENSION._publish_light_step!(
        case.kernel,
        cache,
        targets,
        step,
    )
    return nothing
end

function _pse_benchmark_adapter_allocation_gate(fixture)
    sky = fixture.sky
    coupling = fixture.coupling
    full = fixture.full

    _pse_benchmark_direct_step!(fixture)
    _pse_benchmark_kernel_step!(coupling)
    _pse_benchmark_kernel_step!(full)
    direct_bytes = @allocated _pse_benchmark_direct_step!(fixture)
    coupling_bytes = @allocated _pse_benchmark_kernel_step!(coupling)
    full_bytes = @allocated _pse_benchmark_kernel_step!(full)

    coupling_step = ArchimedLight.run_light(
        coupling.light,
        sky;
        step_duration_seconds=_PSE_BENCHMARK_DURATION_SECONDS,
    )
    full_step = ArchimedLight.run_light(
        full.light,
        sky;
        step_duration_seconds=_PSE_BENCHMARK_DURATION_SECONDS,
    )
    _pse_benchmark_publish!(coupling, coupling_step)
    _pse_benchmark_publish!(full, full_step)
    coupling_publish_bytes = @allocated _pse_benchmark_publish!(
        coupling,
        coupling_step,
    )
    full_publish_bytes = @allocated _pse_benchmark_publish!(full, full_step)

    coupling_bytes == direct_bytes || error(
        "The coupling adapter allocated $(coupling_bytes - direct_bytes) extra " *
        "bytes for $(fixture.ncomponents) components.",
    )
    full_bytes == direct_bytes || error(
        "The full adapter allocated $(full_bytes - direct_bytes) extra bytes " *
        "for $(fixture.ncomponents) components.",
    )
    coupling_publish_bytes == 0 || error(
        "Coupling publication allocated $coupling_publish_bytes bytes for " *
        "$(fixture.ncomponents) components.",
    )
    full_publish_bytes == 0 || error(
        "Full publication allocated $full_publish_bytes bytes for " *
        "$(fixture.ncomponents) components.",
    )
    return (
        direct_bytes=direct_bytes,
        coupling_adapter_bytes=coupling_bytes - direct_bytes,
        full_adapter_bytes=full_bytes - direct_bytes,
        coupling_publish_bytes=coupling_publish_bytes,
        full_publish_bytes=full_publish_bytes,
    )
end

function _pse_benchmark_lifecycle_fixture(ncomponents::Int=100)
    base_scene = _pse_benchmark_scene(ncomponents)
    grown_scene = _pse_benchmark_scene(ncomponents + 1)
    models = _pse_benchmark_models()
    options = _pse_benchmark_options()
    environment = _pse_benchmark_environment()
    current_scene = Ref(base_scene)
    provider_calls = Ref(0)
    provider = _ -> begin
        provider_calls[] += 1
        return current_scene[]
    end
    all_ids = [_pse_benchmark_leaf_id(i) for i in 1:(ncomponents + 1)]
    case = _pse_benchmark_model_case(
        base_scene,
        models,
        options,
        environment,
        all_ids[1:ncomponents],
        :coupling;
        resolver_ids=all_ids,
        scene_provider=provider,
    )
    return (
        case=case,
        current_scene=current_scene,
        grown_scene=grown_scene,
        provider_calls=provider_calls,
        new_id=last(all_ids),
    )
end

@noinline function _pse_benchmark_lifecycle_step!(fixture)
    calls_before = fixture.provider_calls[]
    fixture.current_scene[] = fixture.grown_scene
    PlantSimEngine.register_object!(
        fixture.case.runtime,
        PlantSimEngine.Object(fixture.new_id; scale=:Leaf);
        parent=:scene,
    )
    PlantSimEngine.continue!(fixture.case.simulation; steps=1)
    fixture.provider_calls[] == calls_before + 1 || error(
        "The lifecycle scene provider did not run exactly once.",
    )
    return nothing
end

const PSE_BENCHMARK_FIXTURES = Dict(
    count => _pse_benchmark_fixture(count)
    for count in _PSE_BENCHMARK_COMPONENT_COUNTS
)

# Executable steady-state gates: the adapter may add work, but no allocation
# beyond the underlying light solve. Publication itself must remain at zero.
const PSE_ADAPTER_ALLOCATION_GATE = Dict(
    count => _pse_benchmark_adapter_allocation_gate(
        PSE_BENCHMARK_FIXTURES[count],
    )
    for count in _PSE_BENCHMARK_COMPONENT_COUNTS
)

SUITE["PlantSimEngine"] = BenchmarkGroup()
SUITE["PlantSimEngine"]["steady"] = BenchmarkGroup()
for count in _PSE_BENCHMARK_COMPONENT_COUNTS
    fixture = PSE_BENCHMARK_FIXTURES[count]
    group = SUITE["PlantSimEngine"]["steady"]["$count components"] =
        BenchmarkGroup()
    group["direct run_light"] =
        @benchmarkable _pse_benchmark_direct_step!($fixture) evals = 1
    group["coupling kernel"] =
        @benchmarkable _pse_benchmark_kernel_step!($(fixture.coupling)) evals = 1
    group["coupling continue"] =
        @benchmarkable _pse_benchmark_continue!($(fixture.coupling)) evals = 1
    group["full kernel"] =
        @benchmarkable _pse_benchmark_kernel_step!($(fixture.full)) evals = 1
    group["full continue"] =
        @benchmarkable _pse_benchmark_continue!($(fixture.full)) evals = 1
end

SUITE["PlantSimEngine"]["lifecycle"] = BenchmarkGroup()
SUITE["PlantSimEngine"]["lifecycle"]["100 to 101 components"] =
    @benchmarkable _pse_benchmark_lifecycle_step!(fixture) setup = (
        fixture = _pse_benchmark_lifecycle_fixture(100)
    ) samples = 10 evals = 1 seconds = 2
