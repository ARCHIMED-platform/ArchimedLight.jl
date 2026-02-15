@testset "Phase 0 artifacts" begin
    note = joinpath(dirname(@__DIR__), "docs", "phase0_reverse_engineering.md")
    @test isfile(note)
end

@testset "Phase 1 parity harness" begin
    fixtures = light_parity_fixtures()
    @test length(fixtures) >= 10
    for fx in fixtures
        @test isfile(fx.config_path)
        @test isfile(fx.scene_path)
        @test isfile(fx.meteo_path)
    end
end

@testset "Pipeline smoke test" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-customband")
    cfg = ArchimedLight.read_light_config(joinpath(root, "config.yml"))
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    @test !isempty(meteo.rows)

    step = ArchimedLight.run_light_step(scene, first(meteo.rows), cfg)
    @test length(step.turtle.sectors) >= 1
    @test !isempty(step.budget.ri_par_q_per_node)

    cfg_cache = ArchimedLight.LightConfig(
        cfg.scene,
        cfg.meteo,
        cfg.all_in_turtle,
        cfg.turtle_sectors,
        cfg.pixel_size,
        cfg.area_ratio,
        cfg.scattering,
        cfg.scattering_max_iter,
        cfg.scattering_stop_ratio,
        cfg.scattering_coeff_par,
        cfg.scattering_coeff_nir,
        true,
        cfg.raw
    )
    subset = ArchimedLight.MeteoTable(meteo.rows[1:min(2, end)], meteo.metadata)
    series = ArchimedLight.run_light_series(scene, subset, cfg_cache)
    @test length(series) == length(subset.rows)
end

@testset "OPS GWA scene load" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-hitcount")
    cfg = ArchimedLight.read_light_config(joinpath(root, "config.yml"))
    scene = ArchimedLight.read_scene(cfg.scene)
    @test !isempty(scene.face2node)
end

@testset "Parity reports" begin
    relerr(a, b) = abs(a - b) / max(abs(b), eps(Float64))
    function max_abs_float_dict_diff(a::Dict{Int,Float64}, b::Dict{Int,Float64})
        ids = union(keys(a), keys(b))
        m = 0.0
        for id in ids
            m = max(m, abs(get(a, id, 0.0) - get(b, id, 0.0)))
        end
        m
    end
    function max_abs_int_dict_diff(a::Dict{Int,Int}, b::Dict{Int,Int})
        ids = union(keys(a), keys(b))
        m = 0
        for id in ids
            m = max(m, abs(get(a, id, 0) - get(b, id, 0)))
        end
        m
    end

    fixtures = Dict(f.name => f for f in light_parity_fixtures())

    f_hit = fixtures["test-hitcount2"]
    r_false = fixture_parity_report(f_hit; all_in_turtle=false)
    r_true = fixture_parity_report(f_hit; all_in_turtle=true)
    @test r_false.expected_hitcount_total !== nothing
    @test r_true.expected_hitcount_total !== nothing
    @test r_false.snapshot.total_hits > 0
    @test r_true.snapshot.total_hits > 0
    @test relerr(r_false.snapshot.total_hits, r_false.expected_hitcount_total) < 0.01
    @test relerr(r_true.snapshot.total_hits, r_true.expected_hitcount_total) < 0.01

    f_wsun = fixtures["test-weighted-sun"]
    r_wsun = fixture_parity_report(f_wsun)
    @test r_wsun.expected_sun_azimuth !== nothing
    @test r_wsun.expected_sun_elevation !== nothing
    @test isfinite(r_wsun.snapshot.sun_azimuth)
    @test isfinite(r_wsun.snapshot.sun_elevation)
    @test abs(r_wsun.snapshot.sun_azimuth - r_wsun.expected_sun_azimuth) < 0.2
    @test abs(r_wsun.snapshot.sun_elevation - r_wsun.expected_sun_elevation) < 0.2

    f_scat = fixtures["test-scattering-one-plate"]
    r_scat = fixture_parity_report(f_scat)
    @test r_scat.expected_scattering_total !== nothing
    @test relerr(r_scat.julia_scattering_total, r_scat.expected_scattering_total) < 0.05

    f_scat2 = fixtures["test-scattering-two-plates"]
    r_scat2 = fixture_parity_report(f_scat2)
    @test r_scat2.expected_scattering_total !== nothing
    @test relerr(r_scat2.julia_scattering_total, r_scat2.expected_scattering_total) < 0.01

    f_links = fixtures["test-links-pixeltable"]
    cfg_links = ArchimedLight.read_light_config(f_links.config_path)
    scene_links = ArchimedLight.read_scene(cfg_links.scene)
    meteo_links = ArchimedLight.read_meteo(cfg_links.meteo)
    row_links = first(meteo_links.rows)
    sky_links = ArchimedLight.compute_sky(row_links, cfg_links)
    turtle_links = ArchimedLight.build_turtle(cfg_links, sky_links)
    fluxes_links = ArchimedLight.compute_directional_fluxes(sky_links, turtle_links, cfg_links)
    first_links = ArchimedLight.compute_first_order(scene_links, turtle_links, fluxes_links, cfg_links)
    scat_pix = ArchimedLight.compute_scattering(scene_links, turtle_links, first_links, cfg_links; mode=:raycast)
    scat_lnk = ArchimedLight.compute_scattering(scene_links, turtle_links, first_links, cfg_links; mode=:links)
    bud_pix = ArchimedLight.integrate_light(first_links, scat_pix, cfg_links)
    bud_lnk = ArchimedLight.integrate_light(first_links, scat_lnk, cfg_links)
    @test max_abs_float_dict_diff(bud_pix.ri_par_q_per_node, bud_lnk.ri_par_q_per_node) < 1e-12
    @test max_abs_float_dict_diff(bud_pix.ri_nir_q_per_node, bud_lnk.ri_nir_q_per_node) < 1e-12

    function run_fixture_series(cfg_path::String)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        subset = ArchimedLight.MeteoTable(meteo.rows[1:min(2, end)], meteo.metadata)
        ArchimedLight.run_light_series(scene, subset, cfg)
    end

    for cached_name in ("test-cached-radiation2", "test-cached-radiation3")
        fx = fixtures[cached_name]
        cfg1 = joinpath(fixture_path(fx), "config.yml")
        cfg2 = joinpath(fixture_path(fx), "config2.yml")
        @test isfile(cfg1)
        @test isfile(cfg2)

        s1 = run_fixture_series(cfg1)
        s2 = run_fixture_series(cfg2)
        @test length(s1) == length(s2)

        for i in eachindex(s1)
            @test max_abs_float_dict_diff(s1[i].first_order.projected_area_per_node, s2[i].first_order.projected_area_per_node) == 0.0
            @test max_abs_int_dict_diff(s1[i].first_order.hits_per_node, s2[i].first_order.hits_per_node) == 0
            @test max_abs_float_dict_diff(s1[i].budget.ri_par_q_per_node, s2[i].budget.ri_par_q_per_node) == 0.0
            @test max_abs_float_dict_diff(s1[i].budget.ri_nir_q_per_node, s2[i].budget.ri_nir_q_per_node) == 0.0
        end
    end
end
