using Test

const _TEST_PROFILE = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_PROFILE", "all"))
const _RUN_CORE_TESTS = _TEST_PROFILE in ("all", "core")

if _RUN_CORE_TESTS
    @testset "Core smoke" begin
        case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input")
        cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)

        @test !isempty(selected.rows)
        @test !isempty(series)
        @test length(series) == length(selected.rows)
        @test !isempty(series[1].budget.ri_par_q_per_node)
        @test length(series[1].turtle.sectors) == 17
    end
end
