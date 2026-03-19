@testset "Core smoke" begin
    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    selected = ArchimedLight.prepare_meteo(fixture.meteo, fixture.options)
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)

    @test !isempty(selected.rows)
    @test !isempty(series)
    @test length(series) == length(selected.rows)
    @test !isempty(series[1].budget.incident_energy.total.par)
    @test length(series[1].turtle.sectors) == 17

    options2, scene2, meteo2, models2 = ArchimedLight.read_config(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml"))
    @test options2.turtle_sectors == fixture.options.turtle_sectors
    @test length(meteo2.rows) == length(fixture.meteo.rows)
    @test collect(keys(models2.groups)) == collect(keys(fixture.models.groups))
    @test !isempty([nid for (nid, node) in scene2.nodes if node.group == "pavement"])
    @test haskey(scene2.mtg, :geometry)
    @test scene2.mtg[:geometry] === nothing

    _, scene3, _, _ = ArchimedLight.read_config(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml");
        plot_paving_override=25,
    )
    @test length([nid for (nid, node) in scene3.nodes if node.group == "pavement"]) == 25
    @test scene3.mtg[:geometry] === nothing

    raw_scene = ArchimedLight.read_scene(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "scene", "simple.ops"),
    )
    @test haskey(raw_scene.mtg, :geometry)
    @test raw_scene.mtg[:geometry] !== nothing
    ArchimedLight.add_ground!(raw_scene; nx=2, ny=2)
    @test raw_scene.mtg[:geometry] === nothing
end
