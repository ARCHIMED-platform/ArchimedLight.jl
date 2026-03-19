@testset "Core smoke" begin
    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    selected = ArchimedLight.prepare_meteo(fixture.meteo, fixture.options)
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)

    @test !isempty(selected.rows)
    @test !isempty(series)
    @test length(series) == length(selected.rows)
    @test !isempty(series[1].budget.incident_energy.total.par)
    @test length(series[1].turtle.sectors) == 17
end
