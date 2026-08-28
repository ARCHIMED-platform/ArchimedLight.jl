@testitem "Beginner API model helpers and validation" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using PlantGeom

    scene = make_scene(domain=(0.0, 0.0, 1.0, 1.0)) do s
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
    summary = summarize_scene(scene; models=models)
    @test summary.node_count > 0
    @test summary.face_count > 0
    @test isempty(summary.missing_models)

    wildcard = models_for("simple_plant" => ("*" => translucent(par=0.15, nir=0.90),))
    @test isempty(check_models(scene, wildcard).errors)

    missing = models_for("simple_plant" => ("Leaf" => translucent(par=0.15, nir=0.90),))
    missing_report = check_models(scene, missing)
    @test !isempty(missing_report.errors)
    @test occursin("models_for", missing_report.errors[1])
    @test !isempty(summarize_scene(scene; models=missing).missing_models)

    sensor_models = models_for("sensor" => ("plate" => virtual_sensor(),))
    @test sensor_models["sensor"].types["plate"].interception.sensor

    emitter_models = models_for("lamp" => ("bulb" => emitter(radiance=10.0),))
    @test emitter_models["lamp"].types["bulb"].light_emitter.radiance == 10.0
end

@testitem "Beginner API scene builder object placement" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using GeometryBasics
    using PlantGeom

    mesh = GeometryBasics.Mesh(
        GeometryBasics.Point3f[
            GeometryBasics.Point3f(0, 0, 0),
            GeometryBasics.Point3f(1, 0, 0),
            GeometryBasics.Point3f(0, 1, 0),
        ],
        GeometryBasics.TriangleFace{Int}[GeometryBasics.TriangleFace{Int}(1, 2, 3)],
    )

    scene = make_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=4, at=(1.0, 2.0, 3.0), scale=2.0)
    end
    nid = only(PlantGeom.scene_node_ids(scene))
    sensor_summary = only(summarize_scene(scene).group_types)
    @test sensor_summary.group == "sensor"
    @test sensor_summary.type == "panel"
    @test sensor_summary.object_ids == [4]
    @test PlantGeom.node_areas(scene)[nid] ≈ 2.0
    @test PlantGeom.node_barycenters(scene)[nid] == (5 / 3, 8 / 3, 3.0)

    rotated = make_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=5, rotate=(z=90.0,), deg=true)
    end
    rid = only(PlantGeom.scene_node_ids(rotated))
    @test collect(PlantGeom.node_barycenters(rotated)[rid]) ≈ [-1 / 3, 1 / 3, 0.0]

    xy_order = make_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=6, rotate=(x=90.0, y=90.0), deg=true)
    end
    yx_order = make_scene(domain=(-5.0, -5.0, 5.0, 5.0)) do s
        add_object!(s, mesh; group="sensor", type="panel", id=7, rotate=(y=90.0, x=90.0), deg=true)
    end
    xy_id = only(PlantGeom.scene_node_ids(xy_order))
    yx_id = only(PlantGeom.scene_node_ids(yx_order))
    @test !(collect(PlantGeom.node_barycenters(xy_order)[xy_id]) ≈ collect(PlantGeom.node_barycenters(yx_order)[yx_id]))

    @test_throws ErrorException make_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
        add_object!(s, "panel.obj"; group="sensor", type="panel", id=6)
    end
end

@testitem "Beginner API run_light parity and cache lifecycle" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using PlantGeom

    config = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config)
    row = first(prepare_meteo(meteo, options))

    old_step = ArchimedLight.run_light_step(scene, models, row, options)
    sim = LightSimulation(scene, models; options=options)
    new_step = run_light(sim, row)
    @test old_step.budget.incident_energy.total.par == new_step.budget.incident_energy.total.par

    cached = LightSimulation(scene, models; options=LightOptions(options; cache_radiation=true))
    @test cache_summary(cached).component_metadata_count == 0
    @test cache_summary(cached).component_metadata_bytes == 0
    series = run_light(cached, meteo)
    @test length(series) == length(prepare_meteo(meteo, options))
    @test cache_summary(cached).cached_turtle_count >= 1

    scene2 = make_scene(domain=(0.0, 0.0, 1.0, 1.0)) do s
        add_plant!(s, joinpath(@__DIR__, "..", "example_1", "scene", "opf", "simple_OPF_shapes.opf"); group="simple_plant", id=1)
    end
    update_scene!(cached, scene2)
    @test cache_summary(cached).mode == :unprepared
    @test cached.cache === nothing
end

@testitem "Beginner API write_component_values" tags = [:beginner_api, :fast] begin
    using ArchimedLight
    using CSV
    using Tables

    config = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config)
    sim = LightSimulation(scene, models; options=options)
    series = run_light(sim, meteo)
    step = first(series)
    node_count = length(ArchimedLight._component_output_node_ids(sim.scene, sim.models, sim.options))

    single_path = joinpath(mktempdir(), "nested", "component_values.csv")
    @test write_component_values(single_path, sim, step) == single_path
    @test isfile(single_path)
    @test occursin(';', readline(single_path))
    single_rows = collect(Tables.rowtable(CSV.File(single_path; delim=';', normalizenames=false)))
    @test length(single_rows) == node_count
    @test unique(Int(row.step_number) for row in single_rows) == [1]

    expected_columns = Set([
        "step_number",
        "node_id",
        "source_topology_id",
        "object_id",
        "item_id",
        "component_id",
        "group",
        "type",
        "area",
        "barycentre_z",
        "sky_fraction",
        "Ri_PAR_0_f",
        "Ri_NIR_0_f",
        "Ri_PAR_0_q",
        "Ri_NIR_0_q",
        "Ra_PAR_0_q",
        "Ra_NIR_0_q",
    ])
    @test expected_columns ⊆ Set(string.(propertynames(first(single_rows))))

    series_path = joinpath(mktempdir(), "component_values.csv")
    write_component_values(series_path, sim, series)
    series_rows = collect(Tables.rowtable(CSV.File(series_path; delim=';', normalizenames=false)))
    @test length(series_rows) == node_count * length(series)
    @test sort!(unique(Int(row.step_number) for row in series_rows)) == collect(1:length(series))

    zero_path = joinpath(mktempdir(), "component_values.csv")
    write_component_values(zero_path, sim, series; step_index_base=0)
    zero_rows = collect(Tables.rowtable(CSV.File(zero_path; delim=';', normalizenames=false)))
    @test sort!(unique(Int(row.step_number) for row in zero_rows)) == collect(0:(length(series)-1))

    scattering_options = LightOptions(options; scattering=true)
    scattering_sim = LightSimulation(scene, models; options=scattering_options)
    scattering_step = run_light(scattering_sim, first(prepare_meteo(meteo, scattering_options)))
    @test scattering_step.scattering !== nothing
    scattering_path = joinpath(mktempdir(), "component_values.csv")
    write_component_values(scattering_path, scattering_sim, scattering_step)
    scattering_rows = collect(Tables.rowtable(CSV.File(scattering_path; delim=';', normalizenames=false)))
    scattering_columns = Set(string.(propertynames(first(scattering_rows))))
    @test Set(["Ri_PAR_f", "Ri_NIR_f", "Ri_PAR_q", "Ri_NIR_q", "Ra_PAR_q", "Ra_NIR_q"]) ⊆ scattering_columns
end

@testitem "Beginner API run_light accepts SkyState" tags = [:beginner_api, :fast] begin
    using ArchimedLight

    config = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml")
    options, scene, _, models = ArchimedLight.read_config(config)
    sky = SkyState(135.0, 35.0, 200.0, 180.0, 0.5, 0.5)

    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
    old_budget = ArchimedLight.integrate_light(scene, models, first, scat, options; step_duration_seconds=1800.0)

    sim = LightSimulation(scene, models; options=options)
    step = run_light(sim, sky; step_duration_seconds=1800.0)
    @test step.budget.incident_energy.total.par == old_budget.incident_energy.total.par
    @test_throws ErrorException run_light(sim, sky)

    skies = [
        sky,
        SkyState(180.0, 45.0, 120.0, 130.0, 0.0, 1.0),
    ]
    shared_duration = run_light(sim, skies; step_duration_seconds=1800.0)
    @test shared_duration isa Vector{LightStepResult}
    @test length(shared_duration) == length(skies)
    @test shared_duration[1].budget.incident_energy.total.par ==
          step.budget.incident_energy.total.par

    durations = [900.0, 3600.0]
    variable_duration = run_light(sim, skies; step_duration_seconds=durations)
    for i in eachindex(skies)
        individual = run_light(sim, skies[i]; step_duration_seconds=durations[i])
        @test variable_duration[i].budget.incident_energy.total.par ==
              individual.budget.incident_energy.total.par
    end

    @test isempty(run_light(sim, SkyState[]; step_duration_seconds=1800.0))
    @test_throws ErrorException run_light(sim, skies)
    @test_throws ArgumentError run_light(sim, skies; step_duration_seconds=[1800.0])
    @test_throws ErrorException run_light(sim, skies; step_duration_seconds=[1800.0, 0.0])
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
    meteo_summary = summarize_meteo(good)
    @test meteo_summary.row_count == 1
    @test meteo_summary.duration_seconds == 3600.0
    @test meteo_summary.solar_geometry == "reconstructed from date/time and latitude"
    @test "RI_PAR_f" in meteo_summary.radiation_inputs
    @test "RI_NIR_f" in meteo_summary.radiation_inputs

    missing_lat = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]
    lat_report = check_meteo(missing_lat)
    @test any(contains("latitude"), lat_report.errors)
    @test any(contains("Available columns"), lat_report.errors)
    @test summarize_meteo(missing_lat).solar_geometry == "missing latitude or explicit sun position"

    explicit_sun = [
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            sun_azimuth=180.0,
            sun_elevation=60.0,
            direct_fraction=0.5,
            RI_PAR_f=120.0,
            RI_NIR_f=80.0,
        ),
    ]
    @test isempty(check_meteo(explicit_sun).errors)
end
