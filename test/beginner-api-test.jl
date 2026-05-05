@testitem "Beginner API model helpers and validation" tags = [:beginner_api, :fast] begin
    using ArchimedLight

    scene = light_scene(domain=(0.0, 0.0, 1.0, 1.0)) do s
        add_plant!(s, joinpath(@__DIR__, "..", "example_1", "scene", "opf", "simple_OPF_shapes.opf"); group="simple_plant", id=1)
    end

    models = models_for(
        "simple_plant" => (
            "Leaf" => translucent(par=0.15, nir=0.90),
            "Internode" => translucent(par=0.20, nir=0.50),
        ),
    )
    report = check_models(scene, models)
    @test isempty(report.errors)

    wildcard = models_for("simple_plant" => ("*" => translucent(par=0.15, nir=0.90),))
    @test isempty(check_models(scene, wildcard).errors)

    missing = models_for("simple_plant" => ("Leaf" => translucent(par=0.15, nir=0.90),))
    @test !isempty(check_models(scene, missing).errors)

    sensor_models = models_for("sensor" => ("plate" => virtual_sensor(),))
    @test sensor_models["sensor"].types["plate"].interception.sensor

    emitter_models = models_for("lamp" => ("bulb" => emitter(radiance=10.0),))
    @test emitter_models["lamp"].types["bulb"].light_emitter.radiance == 10.0
end

@testitem "Beginner API scene builder object placement" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using GeometryBasics

    mesh = GeometryBasics.Mesh(
        GeometryBasics.Point3f[
            GeometryBasics.Point3f(0, 0, 0),
            GeometryBasics.Point3f(1, 0, 0),
            GeometryBasics.Point3f(0, 1, 0),
        ],
        GeometryBasics.TriangleFace{Int}[GeometryBasics.TriangleFace{Int}(1, 2, 3)],
    )

    scene = light_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=4, at=(1.0, 2.0, 3.0), scale=2.0)
    end
    nid = only(scene_node_ids(scene))
    @test scene.nodes[nid].group == "sensor"
    @test scene.nodes[nid].type == "panel"
    @test scene.nodes[nid].object_id == 4
    @test node_areas(scene)[nid] ≈ 2.0
    @test node_barycenters(scene)[nid] == (5 / 3, 8 / 3, 3.0)

    rotated = light_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=5, rotate=(z=90.0,), deg=true)
    end
    rid = only(scene_node_ids(rotated))
    @test collect(node_barycenters(rotated)[rid]) ≈ [-1 / 3, 1 / 3, 0.0]

    xy_order = light_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=6, rotate=(x=90.0, y=90.0), deg=true)
    end
    yx_order = light_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=7, rotate=(y=90.0, x=90.0), deg=true)
    end
    xy_id = only(scene_node_ids(xy_order))
    yx_id = only(scene_node_ids(yx_order))
    @test !(collect(node_barycenters(xy_order)[xy_id]) ≈ collect(node_barycenters(yx_order)[yx_id]))

    @test_throws ErrorException light_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
        add_object!(s, "panel.obj"; group="sensor", type="panel", id=6)
    end
end

@testitem "Beginner API run_light parity and cache lifecycle" tags = [:beginner_api, :fast] begin
    using ArchimedLight

    config = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config)
    row = first(prepare_meteo(meteo, options).rows)

    old_step = ArchimedLight.run_light_step(scene, models, row, options)
    sim = LightSimulation(scene, models; options=options)
    new_step = run_light(sim, row)
    @test old_step.budget.incident_energy.total.par == new_step.budget.incident_energy.total.par

    cached = LightSimulation(scene, models; options=LightOptions(options; cache_radiation=true))
    series = run_light(cached, meteo)
    @test length(series) == length(prepare_meteo(meteo, options).rows)
    @test cache_summary(cached).cached_turtle_count >= 1

    scene2 = light_scene(domain=(0.0, 0.0, 1.0, 1.0)) do s
        add_plant!(s, joinpath(@__DIR__, "..", "example_1", "scene", "opf", "simple_OPF_shapes.opf"); group="simple_plant", id=1)
    end
    update_scene!(cached, scene2)
    @test cache_summary(cached).mode == :unprepared
    @test cached.cache === nothing
end

@testitem "Beginner API meteo validation" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using Dates

    good = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            latitude=15.0,
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]
    @test isempty(check_meteo(good).errors)

    missing_lat = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]
    @test any(contains("latitude"), check_meteo(missing_lat).errors)

    explicit_sun = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            sun_azimuth=180.0,
            sun_elevation=60.0,
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]
    @test isempty(check_meteo(explicit_sun).errors)
end
