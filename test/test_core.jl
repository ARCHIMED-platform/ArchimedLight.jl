using Test

@testset "Core smoke" begin
    case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input")
    cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
    scene = ArchimedLight.read_scene(cfg.source_files.scene)
    meteo = ArchimedLight.read_meteo(cfg.source_files.meteo)
    selected = ArchimedLight.prepare_meteo(meteo, cfg)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)

    @test !isempty(selected.rows)
    @test !isempty(series)
    @test length(series) == length(selected.rows)
    @test !isempty(series[1].budget.incident.par.energy_per_node)
    @test length(series[1].turtle.sectors) == 17
end
