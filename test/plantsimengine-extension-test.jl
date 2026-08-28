@testmodule ArchimedLightPlantSimEngineTestSupport begin
    using ArchimedLight
    using Dates
    using GeometryBasics
    using PlantGeom
    using PlantSimEngine

    const COUPLING_NAMES = (
        :Ri_PAR_f,
        :Ri_NIR_f,
        :Ra_PAR_f,
        :Ra_NIR_f,
        :Ra_SW_f,
        :aPPFD,
        :radiative_mesh_area,
    )

    const FULL_NAMES = (
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

    PlantSimEngine.@process "archimed_extension_test_consumer" verbose = false

    struct ArchimedExtensionTestConsumerModel <:
           AbstractArchimed_Extension_Test_ConsumerModel end

    PlantSimEngine.inputs_(::ArchimedExtensionTestConsumerModel) =
        (Ra_PAR_f=PlantSimEngine.Required(Float64),)
    PlantSimEngine.outputs_(::ArchimedExtensionTestConsumerModel) = (seen=0.0,)

    function PlantSimEngine.run!(
        ::ArchimedExtensionTestConsumerModel,
        status,
        environment,
        constants,
        context,
    )
        status.seen = status.Ra_PAR_f
        return nothing
    end

    function triangle_mesh()
        return GeometryBasics.Mesh(
            GeometryBasics.Point3f[
                GeometryBasics.Point3f(0, 0, 0),
                GeometryBasics.Point3f(1, 0, 0),
                GeometryBasics.Point3f(0, 1, 0),
            ],
            GeometryBasics.TriangleFace{Int}[
                GeometryBasics.TriangleFace{Int}(1, 2, 3),
            ],
        )
    end

    function light_fixture(component_count::Int=1; ground::Bool=false)
        component_count in (1, 2) || error("The test fixture supports one or two leaf components.")
        mesh = triangle_mesh()
        scene = PlantGeom.make_scene(domain=(-1.0, -1.0, 5.0, 3.0)) do builder
            PlantGeom.add_object!(
                builder,
                mesh;
                group="plants",
                type="Leaf",
                id=11,
                at=(0.0, 0.0, 1.0),
            )
            if component_count == 2
                PlantGeom.add_object!(
                    builder,
                    mesh;
                    group="plants",
                    type="Leaf",
                    id=12,
                    at=(2.0, 0.0, 1.0),
                    scale=1.5,
                    rotate=(x=60.0,),
                    deg=true,
                )
            end
            ground && PlantGeom.add_ground!(
                builder;
                nx=1,
                ny=1,
                group="soil",
                type="Ground",
                z=0.0,
            )
        end

        model_specs = Pair{String,Any}[
            "plants" => (
                "Leaf" => ArchimedLight.translucent(
                    par=0.15,
                    nir=0.30,
                ),
            ),
        ]
        ground && push!(
            model_specs,
            "soil" => (
                "Ground" => ArchimedLight.translucent(
                    par=0.20,
                    nir=0.40,
                ),
            ),
        )
        models = ArchimedLight.models_for(model_specs...)
        options = ArchimedLight.LightOptions(
            turtle_sectors=1,
            all_in_turtle=false,
            scattering=false,
            toricity=false,
            pixel_size=0.05,
            cache_radiation=true,
        )
        return scene, models, options
    end

    forcing(; duration=Dates.Minute(1)) = (
        duration=duration,
        sun_azimuth_deg=180.0,
        sun_elevation_deg=90.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        direct_fraction=1.0,
    )

    function sky(row=forcing())
        return ArchimedLight.SkyState(
            Float64(row.sun_azimuth_deg),
            Float64(row.sun_elevation_deg),
            Float64(row.Ri_PAR_f),
            Float64(row.Ri_NIR_f),
            Float64(row.direct_fraction),
            1.0 - Float64(row.direct_fraction),
        )
    end

    function light_application(kernel; schema::Symbol=:coupling)
        return PlantSimEngine.ModelSpec(
            kernel;
            name=:archimed_light,
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
    end

    function object_runtime(
        kernel,
        leaf_ids;
        schema::Symbol=:coupling,
        environment=forcing(),
        consumer::Bool=false,
        extra_objects=PlantSimEngine.Object[],
    )
        objects = PlantSimEngine.Object[
            PlantSimEngine.Object(:scene; scale=:Scene),
        ]
        for leaf_id in leaf_ids
            push!(
                objects,
                PlantSimEngine.Object(
                    leaf_id;
                    scale=:Leaf,
                    parent=:scene,
                ),
            )
        end
        append!(objects, extra_objects)
        writer = light_application(kernel; schema=schema)
        applications = if consumer
            (
                PlantSimEngine.ModelSpec(
                    ArchimedExtensionTestConsumerModel();
                    name=:leaf_consumer,
                    on=PlantSimEngine.Many(scale=:Leaf),
                ),
                writer,
            )
        else
            (writer,)
        end
        return PlantSimEngine.CompositeModel(
            objects...;
            applications=applications,
            environment=environment,
        )
    end

    function captured_error(f)
        try
            f()
            return nothing
        catch error
            error isa InterruptException && rethrow()
            return error
        end
    end

    status_values(status, names) =
        NamedTuple{names}(map(name -> getproperty(status, name), names))

    function sentinel_status(offset::Real)
        values = ntuple(i -> Float64(offset + i), length(COUPLING_NAMES))
        return PlantSimEngine.Status(NamedTuple{COUPLING_NAMES}(values))
    end
end

@testitem "PlantSimEngine extension schemas and canonical forcing" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using Dates
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    trait_kernel = ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        object_resolver=_ -> :leaf,
    )
    contracts = variable_contracts(trait_kernel)
    @test keys(contracts) == (:aPPFD, :Ra_SW_f, :radiative_mesh_area)
    @test contracts.aPPFD == VariableContract(
        unit=:micromol_photon,
        basis=:radiative_mesh_area,
        temporal=:second,
        aggregation=:rate,
        extent=:intensive,
    )
    @test contracts.Ra_SW_f == VariableContract(
        unit=:joule,
        basis=:radiative_mesh_area,
        temporal=:second,
        aggregation=:rate,
        extent=:intensive,
    )
    @test contracts.radiative_mesh_area == VariableContract(
        unit=:square_metre_radiative,
        basis=:organ,
        temporal=nothing,
        aggregation=:total,
        extent=:extensive,
    )
    @test all(
        name ∉ keys(contracts)
        for name in (:Ri_PAR_f, :Ri_NIR_f, :Ra_PAR_f, :Ra_NIR_f)
    )
    @test propertynames(PlantSimEngine.environment_inputs_(trait_kernel)) == (
        :sun_azimuth_deg,
        :sun_elevation_deg,
        :Ri_PAR_f,
        :Ri_NIR_f,
        :direct_fraction,
    )
    @test :duration ∉ propertynames(PlantSimEngine.environment_inputs_(trait_kernel))
    @test propertynames(archimed_light_outputs()) == H.COUPLING_NAMES
    @test propertynames(archimed_light_outputs(:full)) == H.FULL_NAMES
    @test all(value -> value isa Default, values(archimed_light_outputs()))
    @test_throws ArgumentError archimed_light_outputs(:unknown)
    @test_throws ArgumentError ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        output_schema=:unknown,
    )

    function run_full(duration)
        simulation = LightSimulation(scene, models; options=options)
        kernel = ArchimedLightModel(
            simulation;
            output_schema=:full,
            object_resolver=_ -> :leaf,
        )
        runtime = H.object_runtime(
            kernel,
            (:leaf,);
            schema=:full,
            environment=H.forcing(; duration=duration),
        )
        return run!(runtime; outputs=:none)
    end

    one_minute = run_full(Minute(1))
    ninety_seconds = run_full(Minute(1) + Second(30))
    two_minutes = run_full(120.0)
    one_state = final_state(one_minute, :leaf)
    ninety_state = final_state(ninety_seconds, :leaf)
    two_state = final_state(two_minutes, :leaf)
    @test one_state.Ra_PAR_f > 0.0
    @test ninety_state.Ra_PAR_f ≈ one_state.Ra_PAR_f
    @test two_state.Ra_PAR_f ≈ one_state.Ra_PAR_f
    @test ninety_state.Ra_PAR_q ≈ 1.5 * one_state.Ra_PAR_q
    @test two_state.Ra_PAR_q ≈ 2.0 * one_state.Ra_PAR_q
    @test two_state.Ri_PAR_q ≈ 2.0 * one_state.Ri_PAR_q

    for invalid_environment in (
        merge(H.forcing(), (sun_azimuth_deg=360.0,)),
        merge(H.forcing(), (sun_elevation_deg=91.0,)),
        merge(H.forcing(), (Ri_PAR_f=-1.0,)),
        merge(H.forcing(), (Ri_NIR_f=Inf,)),
        merge(H.forcing(), (direct_fraction=1.01,)),
        merge(H.forcing(), (duration=0.0,)),
    )
        invalid_kernel = ArchimedLightModel(
            LightSimulation(scene, models; options=options);
            object_resolver=_ -> :leaf,
        )
        invalid_model = H.object_runtime(
            invalid_kernel,
            (:leaf,);
            environment=invalid_environment,
        )
        invalid_error = H.captured_error() do
            run!(invalid_model; outputs=:none)
        end
        @test invalid_error !== nothing
    end

    missing_nir_environment = (
        duration=Minute(1),
        sun_azimuth_deg=180.0,
        sun_elevation_deg=90.0,
        Ri_PAR_f=100.0,
        direct_fraction=1.0,
    )
    missing_kernel = ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        object_resolver=_ -> :leaf,
    )
    missing_model = H.object_runtime(
        missing_kernel,
        (:leaf,);
        environment=missing_nir_environment,
    )
    missing_error = H.captured_error() do
        run!(missing_model; outputs=:none)
    end
    @test missing_error !== nothing
    @test occursin("ri_nir_f", lowercase(sprint(showerror, missing_error)))

    invalid_duration = ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        object_resolver=_ -> :leaf,
    )
    invalid_runtime = H.object_runtime(
        invalid_duration,
        (:leaf,);
        environment=H.forcing(; duration=Month(1)),
    )
    error = H.captured_error() do
        run!(invalid_runtime; outputs=:none)
    end
    @test error !== nothing
    @test occursin("non-fixed period", lowercase(sprint(showerror, error)))

    mismatch_kernel = ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        output_schema=:full,
        object_resolver=_ -> :leaf,
    )
    mismatch_runtime = H.object_runtime(
        mismatch_kernel,
        (:leaf,);
        schema=:coupling,
    )
    mismatch_error = H.captured_error() do
        run!(mismatch_runtime; outputs=:none)
    end
    @test mismatch_error isa ArgumentError
    @test occursin("matching", lowercase(sprint(showerror, mismatch_error)))
end

@testitem "Exact scene MTG identity resolves PlantGeom source owners" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    light = LightSimulation(scene, models; options=options)
    kernel = ArchimedLightModel(light)
    runtime = CompositeModel(
        scene.mtg;
        applications=(H.light_application(kernel),),
        environment=H.forcing(),
    )
    simulation = run!(runtime; outputs=:none)

    leaf_id = only(object_ids(runtime; scale=:Leaf))
    leaf = final_state(simulation, leaf_id)
    @test leaf.Ri_PAR_f > 0.0
    @test leaf.Ra_PAR_f > 0.0
    @test leaf.radiative_mesh_area > 0.0
end

@testitem "ArchimedLight is the distributed writer and assigns by identity" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture(2)
    light = LightSimulation(scene, models; options=options)
    expected_step = run_light(
        light,
        H.sky();
        step_duration_seconds=60.0,
    )
    expected = component_values(expected_step)
    @test length(expected.source_owner) == 2

    destinations = Dict(
        expected.source_owner[1] => :leaf_b,
        expected.source_owner[2] => :leaf_a,
    )
    resolver = owner -> destinations[owner]
    kernel = ArchimedLightModel(light; object_resolver=resolver)
    runtime = H.object_runtime(
        kernel,
        (:leaf_b, :leaf_a);
        consumer=true,
    )
    compiled = Advanced.refresh_bindings!(runtime)

    schedule = Dict(
        row.application_id => row.execution_index
        for row in Diagnostics.explain_schedule(compiled)
    )
    @test schedule[:archimed_light] < schedule[:leaf_consumer]
    consumer_bindings = [
        row for row in Diagnostics.explain_bindings(compiled)
        if row.application_id == :leaf_consumer && row.input == :Ra_PAR_f
    ]
    @test length(consumer_bindings) == 2
    @test all(
        row.source_application_ids == [:archimed_light]
        for row in consumer_bindings
    )
    writer_rows = [
        row for row in Diagnostics.explain_writers(compiled)
        if row.variable == :Ra_PAR_f && row.object_id in (:leaf_a, :leaf_b)
    ]
    @test length(writer_rows) == 2
    @test all(:archimed_light in row.application_ids for row in writer_rows)

    simulation = run!(runtime; outputs=:none)
    for (row, owner) in pairs(expected.source_owner)
        object_id = destinations[owner]
        status = final_state(simulation, object_id)
        for name in H.COUPLING_NAMES
            @test getproperty(status, name) ≈ getproperty(expected, name)[row]
        end
        @test status.seen ≈ status.Ra_PAR_f
    end
end

@testitem "Ground is simulated but excluded from organ publication" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture(1; ground=true)
    light = LightSimulation(scene, models; options=options)
    expected = component_values(run_light(
        light,
        H.sky();
        step_duration_seconds=60.0,
    ))
    @test length(unique(expected.source_owner)) > 1
    leaf_row = only(
        i for i in eachindex(expected.source_owner)
        if expected.source_owner[i].source_instance_id == 1
    )

    resolver = owner -> owner.source_instance_id == 1 ? :leaf : :ground
    kernel = ArchimedLightModel(light; object_resolver=resolver)
    ground_object = PlantSimEngine.Object(
        :ground;
        scale=:Ground,
        parent=:scene,
    )
    runtime = H.object_runtime(
        kernel,
        (:leaf,);
        extra_objects=PlantSimEngine.Object[ground_object],
    )
    simulation = run!(runtime; outputs=:none)
    leaf = final_state(simulation, :leaf)
    for name in H.COUPLING_NAMES
        @test getproperty(leaf, name) ≈ getproperty(expected, name)[leaf_row]
    end
    ground_status = model_object(runtime, :ground).status
    @test ground_status === nothing || !hasproperty(ground_status, :Ra_PAR_f)
end

@testitem "Missing geometrized target fails before any status write" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    kernel = ArchimedLightModel(
        LightSimulation(scene, models; options=options);
        object_resolver=_ -> :leaf_present,
    )
    present = PlantSimEngine.Object(
        :leaf_present;
        scale=:Leaf,
        parent=:scene,
        status=H.sentinel_status(10.0),
    )
    missing = PlantSimEngine.Object(
        :leaf_missing;
        scale=:Leaf,
        parent=:scene,
        status=H.sentinel_status(20.0),
    )
    runtime = PlantSimEngine.CompositeModel(
        PlantSimEngine.Object(:scene; scale=:Scene),
        present,
        missing;
        applications=(H.light_application(kernel),),
        environment=H.forcing(),
    )
    Advanced.refresh_bindings!(runtime)
    before_present = H.status_values(
        model_object(runtime, :leaf_present).status,
        H.COUPLING_NAMES,
    )
    before_missing = H.status_values(
        model_object(runtime, :leaf_missing).status,
        H.COUPLING_NAMES,
    )

    error = H.captured_error() do
        run!(runtime; outputs=:none)
    end
    @test error isa ArgumentError
    message = lowercase(sprint(showerror, error))
    @test occursin("exact distributed-output coverage", message)
    @test occursin("leaf_missing", message)
    @test H.status_values(
        model_object(runtime, :leaf_present).status,
        H.COUPLING_NAMES,
    ) == before_present
    @test H.status_values(
        model_object(runtime, :leaf_missing).status,
        H.COUPLING_NAMES,
    ) == before_missing
end

@testitem "Lifecycle provider refreshes once and retention follows membership" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    first_scene, first_models, options = H.light_fixture(1)
    grown_scene, _, _ = H.light_fixture(2)
    current_scene = Ref(first_scene)
    provider_calls = Ref(0)
    provider = context -> begin
        provider_calls[] += 1
        current_scene[]
    end
    resolver = owner -> owner.source_instance_id == 1 ? :leaf_a : :leaf_b
    light = LightSimulation(first_scene, first_models; options=options)
    kernel = ArchimedLightModel(
        light;
        scene_provider=provider,
        object_resolver=resolver,
    )
    runtime = H.object_runtime(kernel, (:leaf_a,))

    simulation = run!(runtime; steps=2, outputs=:all)
    @test provider_calls[] == 0
    @test length(outputs(simulation)[
        (:archimed_light, ObjectId(:leaf_a), :Ra_PAR_f)
    ]) == 2

    current_scene[] = grown_scene
    register_object!(
        runtime,
        PlantSimEngine.Object(:leaf_b; scale=:Leaf);
        parent=:scene,
    )
    continue!(simulation)

    @test provider_calls[] == 1
    @test light.scene === grown_scene
    @test final_state(simulation, :leaf_b).Ra_PAR_f >= 0.0
    leaf_a_history = outputs(simulation)[
        (:archimed_light, ObjectId(:leaf_a), :Ra_PAR_f)
    ]
    leaf_b_history = outputs(simulation)[
        (:archimed_light, ObjectId(:leaf_b), :Ra_PAR_f)
    ]
    @test first.(leaf_a_history) == [1.0, 2.0, 3.0]
    @test first.(leaf_b_history) == [3.0]
end

@testitem "A same-object scene version bump requires a provider before solving" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantGeom
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    light = LightSimulation(scene, models; options=options)
    kernel = ArchimedLightModel(light; object_resolver=_ -> :leaf)
    runtime = H.object_runtime(kernel, (:leaf,))
    simulation = run!(runtime; outputs=:none)

    prepared_before = kernel.runtime_cache
    light_cache_before = light.cache
    solver_tick_before = light_cache_before.tick
    status_before = H.status_values(
        model_object(runtime, :leaf).status,
        H.COUPLING_NAMES,
    )
    scene_revision_before = scene_version(scene)
    bump_scene_version!(scene.mtg)
    @test scene_version(scene) == scene_revision_before + 1

    error = H.captured_error() do
        continue!(simulation)
    end
    @test error isa ArgumentError
    message = lowercase(sprint(showerror, error))
    @test occursin("scene_provider", message)
    @test occursin("coordinated", message)
    @test kernel.runtime_cache === prepared_before
    @test light.cache === light_cache_before
    @test light.cache.tick == solver_tick_before
    @test H.status_values(
        model_object(runtime, :leaf).status,
        H.COUPLING_NAMES,
    ) == status_before
end

@testitem "An unrelated runtime revision recompiles bindings without dropping the light cache" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantGeom
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    provider_calls = Ref(0)
    provider = context -> begin
        provider_calls[] += 1
        scene
    end
    light = LightSimulation(scene, models; options=options)
    kernel = ArchimedLightModel(
        light;
        scene_provider=provider,
        object_resolver=_ -> :leaf,
    )
    runtime = H.object_runtime(kernel, (:leaf,))
    simulation = run!(runtime; outputs=:none)

    prepared_before = kernel.runtime_cache
    light_cache_before = light.cache
    solver_tick_before = light_cache_before.tick
    model_revision_before = Advanced.model_revision(runtime)
    scene_revision_before = scene_version(scene)
    @test provider_calls[] == 0

    register_object!(
        runtime,
        PlantSimEngine.Object(:observer; scale=:Observer);
        parent=:scene,
    )
    @test Advanced.model_revision(runtime) > model_revision_before
    @test scene_version(scene) == scene_revision_before
    continue!(simulation)

    prepared_after = kernel.runtime_cache
    @test provider_calls[] == 1
    @test prepared_after !== prepared_before
    @test prepared_after.scene === scene
    @test prepared_after.scene_revision == scene_revision_before
    @test prepared_after.model_revision == Advanced.model_revision(runtime)
    @test light.cache === light_cache_before
    @test light.cache.tick > solver_tick_before
    @test final_state(simulation, :leaf).Ra_PAR_f > 0.0
end

@testitem "One kernel cannot be shared by distinct runtimes with the same Scene ObjectId" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture()
    light = LightSimulation(scene, models; options=options)
    kernel = ArchimedLightModel(light; object_resolver=_ -> :leaf)
    first_runtime = H.object_runtime(kernel, (:leaf,))
    second_runtime = H.object_runtime(kernel, (:leaf,))
    Advanced.refresh_bindings!(first_runtime)
    Advanced.refresh_bindings!(second_runtime)
    @test only(object_ids(first_runtime; scale=:Scene)) ==
          only(object_ids(second_runtime; scale=:Scene))

    run!(first_runtime; outputs=:none)
    prepared_before = kernel.runtime_cache
    solver_tick_before = light.cache.tick
    second_status_before = H.status_values(
        model_object(second_runtime, :leaf).status,
        H.COUPLING_NAMES,
    )

    error = H.captured_error() do
        run!(second_runtime; outputs=:none)
    end
    @test error isa ArgumentError
    message = lowercase(sprint(showerror, error))
    @test occursin("distinct", message)
    @test occursin("compositemodel", message)
    @test kernel.runtime_cache === prepared_before
    @test light.cache.tick == solver_tick_before
    @test H.status_values(
        model_object(second_runtime, :leaf).status,
        H.COUPLING_NAMES,
    ) == second_status_before
end

@testitem "Several scene owners coalesce onto one botanical target" tags = [:plantsimengine, :fast] setup = [ArchimedLightPlantSimEngineTestSupport] begin
    using ArchimedLight
    using PlantSimEngine

    const H = ArchimedLightPlantSimEngineTestSupport

    scene, models, options = H.light_fixture(2)
    light = LightSimulation(scene, models; options=options)
    component_owner_values = component_values(run_light(
        light,
        H.sky();
        step_duration_seconds=60.0,
    ))
    @test length(component_owner_values.source_owner) == 2
    areas = component_owner_values.radiative_mesh_area
    total_area = sum(areas)

    kernel = ArchimedLightModel(light; object_resolver=_ -> :compound_leaf)
    runtime = H.object_runtime(kernel, (:compound_leaf,))
    simulation = run!(runtime; outputs=:none)
    status = final_state(simulation, :compound_leaf)

    @test status.radiative_mesh_area ≈ total_area
    for name in (
        :Ri_PAR_f,
        :Ri_NIR_f,
        :Ra_PAR_f,
        :Ra_NIR_f,
        :Ra_SW_f,
        :aPPFD,
    )
        expected = sum(
            areas[i] * getproperty(component_owner_values, name)[i]
            for i in eachindex(areas)
        ) / total_area
        @test getproperty(status, name) ≈ expected
    end
end
