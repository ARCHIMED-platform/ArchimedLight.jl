using Statistics

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

    row = first(meteo.rows)
    step = ArchimedLight.run_light_step(scene, row, cfg)
    @test length(step.turtle.sectors) >= 1
    @test !isempty(step.budget.ri_par_q_per_node)
    step_bk = ArchimedLight.run_light_step(
        scene,
        row,
        cfg;
        interception_backend=ArchimedLight.RasterCPUBackend(),
        scattering_backend=ArchimedLight.RaycastScatteringBackend(),
    )
    step_lk = ArchimedLight.run_light_step(
        scene,
        row,
        cfg;
        interception_backend=ArchimedLight.RasterCPUBackend(),
        scattering_backend=ArchimedLight.LinksScatteringBackend(),
    )

    # Function-first composability check: manual staged execution must match run_light_step.
    sky = ArchimedLight.compute_sky(row, cfg)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    max_abs_float_dict_diff(a::Dict{Int,Float64}, b::Dict{Int,Float64}) = maximum(abs(get(a, id, 0.0) - get(b, id, 0.0)) for id in union(keys(a), keys(b)); init=0.0)
    max_abs_int_dict_diff(a::Dict{Int,Int}, b::Dict{Int,Int}) = maximum(abs(get(a, id, 0) - get(b, id, 0)) for id in union(keys(a), keys(b)); init=0)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
    first_order_cpu_obj = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg; backend=ArchimedLight.RasterCPUBackend())
    scat =
        if cfg.scattering
            graph = ArchimedLight.build_scattering_transfer_graph(scene, turtle, first_order, cfg)
            graph_cpu = ArchimedLight.build_scattering_transfer_graph(
                scene,
                turtle,
                first_order,
                cfg;
                backend=ArchimedLight.RaycastScatteringBackend(),
            )
            scat_graph = ArchimedLight.compute_scattering(graph, first_order, cfg)
            scat_graph_cpu = ArchimedLight.compute_scattering(
                graph_cpu,
                first_order,
                cfg;
                backend=ArchimedLight.RaycastScatteringBackend(),
            )
            scat_scene = ArchimedLight.compute_scattering(scene, turtle, first_order, cfg)
            scat_scene_links = ArchimedLight.compute_scattering(
                scene,
                turtle,
                first_order,
                cfg;
                backend=ArchimedLight.LinksScatteringBackend(),
            )
            @test max_abs_float_dict_diff(scat_graph.added_par_power_per_node, scat_scene.added_par_power_per_node) == 0.0
            @test max_abs_float_dict_diff(scat_graph.added_nir_power_per_node, scat_scene.added_nir_power_per_node) == 0.0
            @test max_abs_float_dict_diff(scat_graph_cpu.added_par_power_per_node, scat_scene.added_par_power_per_node) == 0.0
            @test max_abs_float_dict_diff(scat_graph_cpu.added_nir_power_per_node, scat_scene.added_nir_power_per_node) == 0.0
            @test max_abs_float_dict_diff(scat_scene_links.added_par_power_per_node, scat_scene.added_par_power_per_node) == 0.0
            @test max_abs_float_dict_diff(scat_scene_links.added_nir_power_per_node, scat_scene.added_nir_power_per_node) == 0.0
            @test_throws ErrorException ArchimedLight.compute_scattering(scene, turtle, first_order, cfg; mode=:gpu)
            scat_scene
        else
            nothing
        end
    budget = ArchimedLight.integrate_light(
        first_order,
        scat,
        cfg;
        step_duration_seconds=_step_duration_seconds(row),
        component_area_per_node=scene.total_area_per_node,
    )
    dt_seconds = _step_duration_seconds(row)

    @test abs(step.sky.ri_sw_f - sky.ri_sw_f) < 1e-12
    @test abs(step.sky.direct_fraction - sky.direct_fraction) < 1e-12
    @test length(step.turtle.sectors) == length(turtle.sectors)
    @test max_abs_float_dict_diff(step.first_order.projected_area_per_node, step_bk.first_order.projected_area_per_node) == 0.0
    @test max_abs_int_dict_diff(step.first_order.hits_per_node, step_bk.first_order.hits_per_node) == 0
    @test max_abs_float_dict_diff(step.budget.ri_par_q_per_node, step_bk.budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step.budget.ri_nir_q_per_node, step_bk.budget.ri_nir_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step.budget.ri_par_q_per_node, step_lk.budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step.budget.ri_nir_q_per_node, step_lk.budget.ri_nir_q_per_node) == 0.0
    @test step.fluxes.sector_ids == fluxes.sector_ids
    @test maximum(abs.(step.fluxes.par .- fluxes.par)) == 0.0
    @test maximum(abs.(step.fluxes.nir .- fluxes.nir)) == 0.0
    @test max_abs_float_dict_diff(step.first_order.projected_area_per_node, first_order.projected_area_per_node) == 0.0
    @test max_abs_int_dict_diff(step.first_order.hits_per_node, first_order.hits_per_node) == 0
    @test max_abs_float_dict_diff(first_order.projected_area_per_node, first_order_cpu_obj.projected_area_per_node) == 0.0
    @test max_abs_int_dict_diff(first_order.hits_per_node, first_order_cpu_obj.hits_per_node) == 0
    @test max_abs_float_dict_diff(step.budget.ri_par_q_per_node, budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step.budget.ri_nir_q_per_node, budget.ri_nir_q_per_node) == 0.0
    for nid in keys(step.budget.ri_par_q_per_node)
        area = get(scene.total_area_per_node, nid, 0.0)
        area > 0 || continue
        @test isapprox(step.budget.ri_par_0_q_per_node[nid], step.budget.ri_par_0_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ri_par_q_per_node[nid], step.budget.ri_par_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ri_nir_0_q_per_node[nid], step.budget.ri_nir_0_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ri_nir_q_per_node[nid], step.budget.ri_nir_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ra_par_0_q_per_node[nid], step.budget.ra_par_0_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ra_par_q_per_node[nid], step.budget.ra_par_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ra_nir_0_q_per_node[nid], step.budget.ra_nir_0_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test isapprox(step.budget.ra_nir_q_per_node[nid], step.budget.ra_nir_f_per_node[nid] * area * dt_seconds; atol=1e-12, rtol=1e-12)
        @test step.budget.ra_par_q_per_node[nid] <= step.budget.ri_par_q_per_node[nid] + 1e-12
        @test step.budget.ra_nir_q_per_node[nid] <= step.budget.ri_nir_q_per_node[nid] + 1e-12
    end

    out_cols = [
        "step_number",
        "step_duration",
        "item_id",
        "component_id",
        "group",
        "type",
        "area",
        "barycentre_x",
        "barycentre_y",
        "barycentre_z",
        "sky_fraction",
        "Ri_PAR_0_f",
        "Ri_PAR_0_q",
        "Ra_PAR_0_f",
        "Ra_PAR_0_q",
        "Ri_NIR_f",
        "Ri_NIR_q",
        "Ra_NIR_f",
        "Ra_NIR_q",
        "Ri_custom_q",
    ]
    tbl = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, step_number=1, columns=out_cols)
    @test tbl.columns == out_cols
    @test length(tbl.rows) == length(scene.total_area_per_node)
    @test all(haskey(r, "Ri_PAR_0_q") for r in tbl.rows)
    @test all(haskey(r, "Ra_PAR_0_q") for r in tbl.rows)
    @test all(haskey(r, "Ri_custom_q") for r in tbl.rows)
    @test all(get(r, "type", "NA") != "NA" for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_x", NaN))) for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_y", NaN))) for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_z", NaN))) for r in tbl.rows)
    @test all(Float64(get(r, "sky_fraction", -1.0)) >= 0.0 for r in tbl.rows)

    mktempdir() do tmp
        out_csv = joinpath(tmp, "component_values.csv")
        ArchimedLight.write_component_values_csv(out_csv, scene, step, cfg; meteo_row=row, step_number=1, columns=out_cols)
        @test isfile(out_csv)
        rows_csv = read_java_csv(out_csv)
        @test length(rows_csv) == length(scene.total_area_per_node)
        @test :Ri_PAR_0_q in propertynames(first(rows_csv))
        @test :Ra_PAR_0_q in propertynames(first(rows_csv))
    end

    @test_throws ErrorException ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg; backend=:gpu)
    if cfg.scattering
        @test step.scattering !== nothing
        @test scat !== nothing
        @test max_abs_float_dict_diff(step.scattering.added_par_power_per_node, scat.added_par_power_per_node) == 0.0
        @test max_abs_float_dict_diff(step.scattering.added_nir_power_per_node, scat.added_nir_power_per_node) == 0.0
    else
        @test step.scattering === nothing
    end

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
    series_bk = ArchimedLight.run_light_series(
        scene,
        subset,
        cfg_cache;
        interception_backend=ArchimedLight.RasterCPUBackend(),
        scattering_backend=ArchimedLight.RaycastScatteringBackend(),
    )
    series_lk = ArchimedLight.run_light_series(
        scene,
        subset,
        cfg_cache;
        interception_backend=ArchimedLight.RasterCPUBackend(),
        scattering_backend=ArchimedLight.LinksScatteringBackend(),
    )
    @test length(series) == length(subset.rows)
    @test length(series_bk) == length(series)
    @test length(series_lk) == length(series)
    for i in eachindex(series)
        @test max_abs_float_dict_diff(series[i].first_order.projected_area_per_node, series_bk[i].first_order.projected_area_per_node) == 0.0
        @test max_abs_int_dict_diff(series[i].first_order.hits_per_node, series_bk[i].first_order.hits_per_node) == 0
        @test max_abs_float_dict_diff(series[i].budget.ri_par_q_per_node, series_bk[i].budget.ri_par_q_per_node) == 0.0
        @test max_abs_float_dict_diff(series[i].budget.ri_nir_q_per_node, series_bk[i].budget.ri_nir_q_per_node) == 0.0
        @test max_abs_float_dict_diff(series[i].budget.ri_par_q_per_node, series_lk[i].budget.ri_par_q_per_node) == 0.0
        @test max_abs_float_dict_diff(series[i].budget.ri_nir_q_per_node, series_lk[i].budget.ri_nir_q_per_node) == 0.0
    end

    scene_cols = ["step_number", "step_duration", "hour_start", "hour_end", "RI_SW_f", "plot_area"]
    scene_tbl = ArchimedLight.scene_values_table(scene, series, cfg; meteo_rows=subset.rows, columns=scene_cols)
    @test scene_tbl.columns == scene_cols
    @test length(scene_tbl.rows) == length(series)
    @test all(haskey(r, "RI_SW_f") for r in scene_tbl.rows)
    @test all(Float64(r["plot_area"]) > 0 for r in scene_tbl.rows)

    mktempdir() do tmp
        out_csv = joinpath(tmp, "scene_values.csv")
        ArchimedLight.write_scene_values_csv(out_csv, scene, series, cfg; meteo_rows=subset.rows, columns=scene_cols)
        @test isfile(out_csv)
        rows_csv = read_java_csv(out_csv)
        @test length(rows_csv) == length(series)
        @test :RI_SW_f in propertynames(first(rows_csv))
        got = Float64[getproperty(r, :RI_SW_f) for r in rows_csv]
        exp = Float64[s.sky.ri_sw_f for s in series]
        @test maximum(abs.(got .- exp)) < 1e-12
    end

    sun_tbl = ArchimedLight.sun_position_log_table(series, subset.rows; start_step_number=0)
    @test length(sun_tbl.rows) == length(series)
    @test "azimuthWeighted" in sun_tbl.columns
    @test "elevationWeighted" in sun_tbl.columns
    for i in eachindex(series)
        r = sun_tbl.rows[i]
        @test isapprox(Float64(r["azimuthWeighted"]), deg2rad(series[i].sky.sun_azimuth_deg); atol=1e-12, rtol=1e-12)
        @test isapprox(Float64(r["elevationWeighted"]), deg2rad(series[i].sky.sun_elevation_deg); atol=1e-12, rtol=1e-12)
    end

    mktempdir() do tmp
        sun_csv = joinpath(tmp, "log-sun-position.csv")
        ArchimedLight.write_sun_position_log_csv(sun_csv, series, subset.rows; start_step_number=0)
        @test isfile(sun_csv)
        sun_rows = read_java_csv(sun_csv)
        @test length(sun_rows) == length(series)
        @test :azimuthWeighted in propertynames(first(sun_rows))
    end

    if cfg.scattering
        scat_tbl = ArchimedLight.scattering_iteration_log_table(scene, step, cfg; meteo_row=row, step_number=0, band="PAR")
        @test !isempty(scat_tbl.rows)
        @test scat_tbl.columns == ["step", "plantid", "nodeid", "iter", "scat"]
        scat_sum = sum(Float64(r["scat"]) for r in scat_tbl.rows)
        expected_sum = sum(get(step.scattering.added_par_power_per_node, nid, 0.0) for nid in keys(scene.total_area_per_node)) * _step_duration_seconds(row) / 1e6
        @test isapprox(scat_sum, expected_sum; atol=1e-12, rtol=1e-12)

        mktempdir() do tmp
            scat_csv = joinpath(tmp, "log-iteration-scat-par.csv")
            ArchimedLight.write_scattering_iteration_log_csv(scat_csv, scene, step, cfg; meteo_row=row, step_number=0, band="PAR")
            @test isfile(scat_csv)
            scat_rows = read_java_csv(scat_csv)
            @test !isempty(scat_rows)
            @test :scat in propertynames(first(scat_rows))
        end
    end

    function with_raw_overrides(cfg::ArchimedLight.LightConfig, raw_override::Dict{String,Any})
        ArchimedLight.LightConfig(
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
            cfg.cache_radiation,
            raw_override,
        )
    end

    mktempdir() do tmp
        out_base = joinpath(tmp, "output")
        raw_auto = copy(cfg.raw)
        raw_auto["output_directory"] = out_base
        delete!(raw_auto, "simulation_directory")
        cfg_auto = with_raw_overrides(cfg, raw_auto)

        sim1 = ArchimedLight.simulation_output_directory(cfg_auto)
        sim2 = ArchimedLight.simulation_output_directory(cfg_auto)
        @test basename(sim1) == "000001"
        @test basename(sim2) == "000002"
        @test dirname(sim1) == out_base
        @test dirname(sim2) == out_base

        raw_named = copy(cfg.raw)
        raw_named["output_directory"] = out_base
        raw_named["simulation_directory"] = "simdir"
        cfg_named = with_raw_overrides(cfg, raw_named)
        sim_named = ArchimedLight.simulation_output_directory(cfg_named)
        dummy = joinpath(sim_named, "dummy.txt")
        write(dummy, "x")
        @test isfile(dummy)
        sim_named2 = ArchimedLight.simulation_output_directory(cfg_named)
        @test sim_named2 == sim_named
        @test !isfile(dummy)
    end

    mktempdir() do tmp
        out_base = joinpath(tmp, "output")
        raw_out = copy(cfg.raw)
        raw_out["output_directory"] = out_base
        delete!(raw_out, "simulation_directory")
        cfg_out = with_raw_overrides(cfg, raw_out)
        out_paths = ArchimedLight.write_light_outputs(
            scene,
            series,
            cfg_out;
            meteo_rows=subset.rows,
            start_step_number=0,
            write_component=true,
            write_scene=true,
            write_summary=true,
            write_sun_position_log=true,
            write_scattering_log=cfg_out.scattering,
            scattering_log_bands=["PAR"],
        )
        @test haskey(out_paths, "output_directory")
        @test basename(out_paths["output_directory"]) == "000001"
        @test haskey(out_paths, "component_values")
        @test haskey(out_paths, "scene_values")
        @test haskey(out_paths, "summary")
        @test haskey(out_paths, "log_sun_position")
        @test isfile(out_paths["component_values"])
        @test isfile(out_paths["scene_values"])
        @test isfile(out_paths["summary"])
        @test isfile(out_paths["log_sun_position"])
        if cfg_out.scattering
            @test haskey(out_paths, "log_iteration_scat_par")
            @test isfile(out_paths["log_iteration_scat_par"])
        end

        comp_rows = read_java_csv(out_paths["component_values"])
        @test !isempty(comp_rows)
        @test maximum(Int[getproperty(r, :step_number) for r in comp_rows]) == length(series) - 1
        scene_rows = read_java_csv(out_paths["scene_values"])
        @test length(scene_rows) == length(series)
        summary_rows = read_java_csv(out_paths["summary"])
        @test !isempty(summary_rows)
    end
end

@testset "OPS GWA scene load" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-hitcount")
    cfg = ArchimedLight.read_light_config(joinpath(root, "config.yml"))
    scene = ArchimedLight.read_scene(cfg.scene)
    @test !isempty(scene.face2node)

    # Legacy OPS rows are interpreted with extension-dependent scale normalization:
    # .gwa keeps historical cm->m conversion, .opf stays at scale=1.0.
    cfg_gwa = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-hitcount2", "config.yml"))
    scene_gwa = ArchimedLight.read_scene(cfg_gwa.scene)
    @test isapprox(sum(values(scene_gwa.total_area_per_node)), 1.0; atol=1e-9, rtol=1e-9)

    cfg_opf = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-customband", "config.yml"))
    scene_opf = ArchimedLight.read_scene(cfg_opf.scene)
    @test isapprox(sum(values(scene_opf.total_area_per_node)), 0.007472579789211572; atol=1e-12, rtol=1e-9)
end

@testset "Turtle geometry parity" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    sky_ref = ArchimedLight.SkyState(0.0, 45.0, 1.0, 1.0, 0.0, 1.0)

    function with_sectors(cfg::ArchimedLight.LightConfig, n::Int)
        ArchimedLight.LightConfig(
            cfg.scene,
            cfg.meteo,
            true,
            n,
            cfg.pixel_size,
            cfg.area_ratio,
            cfg.scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            cfg.cache_radiation,
            cfg.raw,
        )
    end

    for n in (1, 6, 16, 46, 136, 406)
        t = ArchimedLight.build_turtle(with_sectors(cfg_ref, n), sky_ref)
        @test length(t.sectors) == n
    end

    t6 = ArchimedLight.build_turtle(with_sectors(cfg_ref, 6), sky_ref)
    got6 = [s.direction for s in t6.sectors]
    exp6 = [
        (0.0, 0.0, -1.0),
        (1.5637987e-7, -0.89438856, -0.44729084),
        (-0.85061413, -0.27638108, -0.44729084),
        (-0.5257086, 0.7235755, -0.4472909),
        (0.5257085, 0.7235755, -0.44729084),
        (0.85061413, -0.2763814, -0.4472909),
    ]
    @test length(got6) == length(exp6)
    for i in eachindex(exp6)
        @test abs(got6[i][1] - exp6[i][1]) < 1e-6
        @test abs(got6[i][2] - exp6[i][2]) < 1e-6
        @test abs(got6[i][3] - exp6[i][3]) < 1e-6
    end
end

@testset "Pixel size validation parity" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-pixelsize")
    cfg = ArchimedLight.read_light_config(joinpath(root, "config.yml"))
    scene = ArchimedLight.read_scene(cfg.scene)
    row = first(ArchimedLight.read_meteo(cfg.meteo).rows)

    function with_pixel_size(cfg::ArchimedLight.LightConfig, pixel_size_m::Float64)
        ArchimedLight.LightConfig(
            cfg.scene,
            cfg.meteo,
            cfg.all_in_turtle,
            cfg.turtle_sectors,
            pixel_size_m,
            cfg.area_ratio,
            cfg.scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            cfg.cache_radiation,
            cfg.raw,
        )
    end

    cfg_ok = with_pixel_size(cfg, 49.9999 / 100.0)
    @test begin
        step = ArchimedLight.run_light_step(scene, row, cfg_ok)
        !isempty(step.budget.ri_par_q_per_node)
    end

    cfg_bad = with_pixel_size(cfg, 50.0001 / 100.0)
    @test_throws ErrorException ArchimedLight.run_light_step(scene, row, cfg_bad)
end

@testset "Clearness/global parity" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-clearness-global")
    cfg = ArchimedLight.read_light_config(joinpath(root, "clearness", "config.yml"))
    meteo_c = ArchimedLight.read_meteo(joinpath(root, "clearness", "meteo.csv"))
    meteo_g = ArchimedLight.read_meteo(joinpath(root, "global", "meteo.csv"))
    @test length(meteo_c.rows) == length(meteo_g.rows)
    @test !isempty(meteo_c.rows)

    # Java meteo conventions leave date empty after the first line.
    # We forward-fill so day-of-year dependent sky computations remain stable.
    @test getproperty(meteo_c.rows[1], :date) !== missing
    @test getproperty(meteo_c.rows[2], :date) == getproperty(meteo_c.rows[1], :date)

    for i in eachindex(meteo_c.rows)
        sky_c = ArchimedLight.compute_sky(meteo_c.rows[i], cfg)
        sky_g = ArchimedLight.compute_sky(meteo_g.rows[i], cfg)

        row_g = meteo_g.rows[i]
        global_ref =
            if :RI_SW_f in propertynames(row_g)
                Float64(getproperty(row_g, :RI_SW_f))
            else
                Float64(getproperty(row_g, :global))
            end
        @test abs(sky_c.ri_sw_f - global_ref) < 2e-4
        @test abs(sky_g.ri_sw_f - global_ref) < 2e-9
        @test abs(sky_c.direct_fraction - sky_g.direct_fraction) < 0.06
    end
end

@testset "Meteo use parity" begin
    root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-clearness-global2")
    cfg = ArchimedLight.read_light_config(joinpath(root, "config.yml"))

    function sky_from_file(n::Int)
        m = ArchimedLight.read_meteo(joinpath(root, "meteo$(n).csv"))
        ArchimedLight.compute_sky(first(m.rows), cfg)
    end

    for n in (1, 2, 3, 4, 5, 6, 7, 17, 18, 19, 20, 21, 23, 24, 25)
        s = sky_from_file(n)
        @test isfinite(s.ri_par_f + s.ri_nir_f)
        @test isfinite(s.direct_fraction)
    end

    for n in (8, 9, 10, 11, 12, 16, 22)
        @test_throws ErrorException sky_from_file(n)
    end

    # Java behavior: with both clearness and RI_SW_f columns present, use lines validate
    # consistency but the provided RI_SW_f remains the irradiance source.
    s5 = sky_from_file(5)
    s6 = sky_from_file(6)
    s7 = sky_from_file(7)
    @test abs(s5.ri_sw_f - 100.0) < 1e-12
    @test abs(s6.ri_sw_f - 100.0) < 1e-12
    @test abs(s7.ri_sw_f - 100.0) < 1e-12

    # Java ClearnessGlobalRelation precedence checks on valid variants.
    s1 = sky_from_file(1)   # clearness only
    s2 = sky_from_file(2)   # global only
    s17 = sky_from_file(17) # PAR only
    s18 = sky_from_file(18) # NIR only
    s19 = sky_from_file(19) # PAR + NIR
    s20 = sky_from_file(20) # PAR + global
    s21 = sky_from_file(21) # NIR + global
    s23 = sky_from_file(23) # PAR + NIR + global (+use)
    s24 = sky_from_file(24)
    s25 = sky_from_file(25)

    @test abs(s1.ri_sw_f - 528.1693424202258) < 1e-4
    @test abs(s2.ri_sw_f - 0.67) < 1e-12

    @test abs(s17.ri_sw_f - (10 / 0.48)) < 1e-9
    @test abs(s18.ri_sw_f - (10 / 0.52)) < 1e-9
    @test abs(s19.ri_sw_f - 20.0) < 1e-12

    @test abs(s20.ri_sw_f - 10.0) < 1e-12
    @test abs(s21.ri_sw_f - 10.0) < 1e-12
    @test abs(s23.ri_sw_f - 10.0) < 1e-12
    @test abs(s24.ri_sw_f - 10.0) < 1e-12
    @test abs(s25.ri_sw_f - 10.0) < 1e-12

    # When RI_SW_f is defined and only one waveband is provided, missing band comes from RI_SW_f.
    @test abs(s20.ri_nir_f - (10.0 * 0.52)) < 1e-12
    @test abs(s21.ri_par_f - (10.0 * 0.48)) < 1e-12
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
    to_int(v) = v isa Number ? Int(round(v)) : parse(Int, string(v))
    mean_std(vals::Vector{Float64}) = (mean(vals), std(vals))
    function with_cache_pixel_table(cfg::ArchimedLight.LightConfig, enabled::Bool)
        raw = copy(cfg.raw)
        raw["cache_pixel_table"] = enabled
        raw["save_on_disk"] = enabled
        ArchimedLight.LightConfig(
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
            cfg.cache_radiation,
            raw,
        )
    end
    function expected_component_hits(path::String)
        rows = read_java_csv(path)
        vals = Int[]
        for r in rows
            names = propertynames(r)
            (:area in names && :surface_hits in names) || continue
            area = Float64(getproperty(r, :area))
            hits = Float64(getproperty(r, :surface_hits))
            push!(vals, Int(floor(area * hits)))
        end
        sort(vals)
    end
    function expected_component_metrics(path::String)
        rows = read_java_csv(path)
        out = Dict{Tuple{Int,Int},NamedTuple{(:area,:hits,:irr0),Tuple{Float64,Int,Float64}}}()
        for r in rows
            item = to_int(getproperty(r, :item_id))
            comp = to_int(getproperty(r, :component_id))
            area = Float64(getproperty(r, :area))
            hits = Int(floor(area * Float64(getproperty(r, :surface_hits))))
            irr0 = Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR))
            out[(item, comp)] = (area=area, hits=hits, irr0=irr0)
        end
        out
    end
    function observed_component_metrics(scene::ArchimedLight.SceneGeometry, cfg::ArchimedLight.LightConfig, step::ArchimedLight.LightStepResult)
        key_by_node = ArchimedLight._interception_java_keys(scene, cfg)
        out = Dict{Tuple{Int,Int},NamedTuple{(:hits,:power0),Tuple{Int,Float64}}}()
        for (nid, key) in key_by_node
            p0 = get(step.first_order.incident_par_power_per_node, nid, 0.0) + get(step.first_order.incident_nir_power_per_node, nid, 0.0)
            out[key] = (hits=get(step.first_order.hits_per_node, nid, 0), power0=p0)
        end
        out
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

    for all_in_turtle in (false, true)
        cfg = ArchimedLight.read_light_config(f_hit.config_path)
        cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)

        expected_path = _expected_component_values_path(f_hit; all_in_turtle=all_in_turtle)
        @test expected_path !== nothing
        expected = expected_component_metrics(expected_path)
        observed = observed_component_metrics(scene, cfg, step)

        @test Set(keys(observed)) == Set(keys(expected))
        for key in keys(expected)
            exp = expected[key]
            got = observed[key]
            @test abs(got.hits - exp.hits) / max(exp.hits, 1) < 0.003
        end

        # Item component IDs are now mapped explicitly; plant component irradiance is in tight parity.
        plant_key = (1, 2)
        @test haskey(expected, plant_key)
        @test haskey(observed, plant_key)
        expp = expected[plant_key]
        gotp = observed[plant_key]
        got_irr = gotp.power0 / max(expp.area, eps(Float64))
        @test abs(got_irr - expp.irr0) / max(abs(expp.irr0), eps(Float64)) < 1e-6
    end

    f_hit3 = fixtures["test-hitcount3"]
    for all_in_turtle in (false, true)
        cfg = ArchimedLight.read_light_config(f_hit3.config_path)
        cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        expected_path = _expected_component_values_path(f_hit3; all_in_turtle=all_in_turtle)
        @test expected_path !== nothing
        expected = expected_component_metrics(expected_path)
        observed = observed_component_metrics(scene, cfg, step)

        missing = setdiff(Set(keys(expected)), Set(keys(observed)))
        extra = setdiff(Set(keys(observed)), Set(keys(expected)))
        @test isempty(missing)
        @test isempty(extra)

        hit_rel_hi = Float64[]
        irr_rel_hi = Float64[]
        irr_abs_hi = Float64[]
        for key in keys(expected)
            exp = expected[key]
            got = observed[key]
            rel_hit = abs(got.hits - exp.hits) / max(exp.hits, 1)
            if exp.hits >= 1000
                push!(hit_rel_hi, rel_hit)
            end

            got_irr = got.power0 / max(exp.area, eps(Float64))
            abs_err_irr = abs(got_irr - exp.irr0)
            if abs(exp.irr0) >= 1.0
                push!(irr_abs_hi, abs_err_irr)
                push!(irr_rel_hi, abs_err_irr / max(abs(exp.irr0), eps(Float64)))
            end
        end

        @test !isempty(hit_rel_hi)
        @test !isempty(irr_abs_hi)
        if all_in_turtle
            @test maximum(hit_rel_hi) < 0.003
            @test maximum(irr_rel_hi) < 0.4
        else
            @test maximum(hit_rel_hi) < 0.05
            @test maximum(irr_abs_hi) < 650.0
        end
    end

    for (all_in_turtle, max_abs_err, max_rel_err) in ((false, 900, 0.002), (true, 500, 0.0012))
        step = run_fixture_once(f_hit; all_in_turtle=all_in_turtle)
        expected_path = _expected_component_values_path(f_hit; all_in_turtle=all_in_turtle)
        @test expected_path !== nothing
        expected_hits = expected_component_hits(expected_path)
        got_hits = sort(collect(values(step.first_order.hits_per_node)))
        @test length(got_hits) == length(expected_hits)
        abs_err = abs.(got_hits .- expected_hits)
        rel_errs = abs_err ./ max.(abs.(Float64.(expected_hits)), 1.0)
        @test maximum(abs_err) <= max_abs_err
        @test maximum(rel_errs) <= max_rel_err
    end

    f_wsun = fixtures["test-weighted-sun"]
    r_wsun = fixture_parity_report(f_wsun)
    @test r_wsun.expected_sun_azimuth !== nothing
    @test r_wsun.expected_sun_elevation !== nothing
    @test isfinite(r_wsun.snapshot.sun_azimuth)
    @test isfinite(r_wsun.snapshot.sun_elevation)
    @test abs(r_wsun.snapshot.sun_azimuth - r_wsun.expected_sun_azimuth) < 0.05
    @test abs(r_wsun.snapshot.sun_elevation - r_wsun.expected_sun_elevation) < 0.05

    cfg_wsun = ArchimedLight.read_light_config(f_wsun.config_path)
    meteo_wsun = ArchimedLight.read_meteo(cfg_wsun.meteo)
    sun_log_wsun = _expected_sun_log_path(f_wsun)
    @test sun_log_wsun !== nothing
    rows_wsun = read_java_csv(sun_log_wsun)
    @test length(rows_wsun) == length(meteo_wsun.rows)
    max_az_err = 0.0
    max_el_err = 0.0
    for i in eachindex(meteo_wsun.rows)
        sky_i = ArchimedLight.compute_sky(meteo_wsun.rows[i], cfg_wsun)
        az_exp = _to_degrees_if_radians(Float64(getproperty(rows_wsun[i], :azimuthWeighted)))
        el_exp = _to_degrees_if_radians(Float64(getproperty(rows_wsun[i], :elevationWeighted)))
        da = abs(sky_i.sun_azimuth_deg - az_exp)
        da = min(da, 360.0 - da)
        de = abs(sky_i.sun_elevation_deg - el_exp)
        max_az_err = max(max_az_err, da)
        max_el_err = max(max_el_err, de)
    end
    @test max_az_err < 0.02
    @test max_el_err < 0.05

    for (fx_name, max_abs_err, max_mean_rel_err) in (
        ("test-weighted-sun", 1.0, 0.01),
        ("test-save_on_disk1", 1e-3, 1e-3),
        ("test-save_on_disk6", 1e-3, 1e-3),
    )
        fx = fixtures[fx_name]
        expected_path = _expected_scene_values_path(fx)
        @test expected_path !== nothing
        expected_rows = read_java_csv(expected_path)
        expected_sw = Float64[getproperty(r, :RI_SW_f) for r in expected_rows]

        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        got_sw = Float64[s.sky.ri_sw_f for s in series]
        @test length(got_sw) == length(expected_sw)

        abs_err = abs.(got_sw .- expected_sw)
        rel_err = abs_err ./ max.(abs.(expected_sw), eps(Float64))
        @test maximum(abs_err) < max_abs_err
        @test mean(rel_err) < max_mean_rel_err
    end

    f_scat = fixtures["test-scattering-one-plate"]
    r_scat = fixture_parity_report(f_scat)
    @test r_scat.expected_scattering_total !== nothing
    @test relerr(r_scat.julia_scattering_total, r_scat.expected_scattering_total) < 0.05

    f_scat2 = fixtures["test-scattering-two-plates"]
    r_scat2 = fixture_parity_report(f_scat2)
    @test r_scat2.expected_scattering_total !== nothing
    @test relerr(r_scat2.julia_scattering_total, r_scat2.expected_scattering_total) < 0.01

    for links_fx_name in ("test-links-pixeltable", "test-links-pixeltable2")
        f_links = fixtures[links_fx_name]
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
    end

    area_ratio_metrics = Dict{String,Tuple{Float64,Float64,Float64,Float64}}()
    for area_name in ("test-area_ratio", "test-area_ratio2", "test-area_ratio3", "test-area_ratio4")
        fx = fixtures[area_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        dt_seconds = _step_duration_seconds(row)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        @test !isempty(step.budget.ri_par_0_q_per_node)

        base_max = maximum(keys(scene.total_area_per_node))
        got = Float64[]
        for nid in keys(step.budget.ri_par_0_q_per_node)
            nid > base_max || continue
            push!(got, step.budget.ri_par_0_q_per_node[nid] + step.budget.ri_nir_0_q_per_node[nid])
        end
        @test length(got) == 100
        got_mean, got_std = mean_std(got)

        expected_path = _expected_component_values_path(fx)
        @test expected_path !== nothing
        expected_rows = read_java_csv(expected_path)
        expected = Float64[
            Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR)) * Float64(getproperty(r, :area)) * dt_seconds for
            r in expected_rows if to_int(getproperty(r, :item_id)) == -1
        ]
        @test length(expected) == 100
        exp_mean, exp_std = mean_std(expected)

        # Java intent: pavement intercepted energy is effectively uniform with area-ratio correction.
        @test got_std / max(abs(got_mean), eps(Float64)) < 1e-8
        @test exp_std / max(abs(exp_mean), eps(Float64)) < 1e-5
        # Java sky conversion and units are matched closely.
        @test relerr(got_mean, exp_mean) < 1e-5

        area_ratio_metrics[area_name] = (got_mean, got_std, exp_mean, exp_std)
    end

    # Java runs this fixture family under four cache/all_in_turtle combinations.
    # Means should stay stable across those runtime switches.
    got_means = [v[1] for v in values(area_ratio_metrics)]
    @test maximum(got_means) - minimum(got_means) < 1e-9

    f_links_stats = fixtures["test-links"]
    cfg_ls = ArchimedLight.read_light_config(f_links_stats.config_path)
    scene_ls = ArchimedLight.read_scene(cfg_ls.scene)
    meteo_ls = ArchimedLight.read_meteo(cfg_ls.meteo)
    row_ls = first(meteo_ls.rows)
    sky_ls = ArchimedLight.compute_sky(row_ls, cfg_ls)
    turtle_ls = ArchimedLight.build_turtle(cfg_ls, sky_ls)
    flux_ls = ArchimedLight.compute_directional_fluxes(sky_ls, turtle_ls, cfg_ls)
    first_ls = ArchimedLight.compute_first_order(scene_ls, turtle_ls, flux_ls, cfg_ls)
    pair_counts_ls, sun_hits_ls, node_ids_ls, _ = ArchimedLight._pair_counts_for_scattering(scene_ls, turtle_ls, cfg_ls)
    all_hits_ls = ArchimedLight._all_dir_hits_for_scattering(first_ls, sun_hits_ls, cfg_ls, node_ids_ls)

    got_hits = sort([v for v in values(all_hits_ls) if v > 0])
    neighbours_by_node = Dict{Int,Set{Int}}()
    for ((to, from), c) in pair_counts_ls
        c > 0 || continue
        s = get!(neighbours_by_node, to) do
            Set{Int}()
        end
        push!(s, from)
    end
    got_links = sort([length(v) for v in values(neighbours_by_node) if !isempty(v)])

    expected_links_path = joinpath(fixture_path(f_links_stats), "expected", "log-nodelinks-stats-alldirs.csv")
    @test isfile(expected_links_path)
    expected_rows = read_java_csv(expected_links_path)
    expected_hits = sort(Int[to_int(getproperty(r, :hits)) for r in expected_rows])
    expected_links = sort(Int[to_int(getproperty(r, :links)) for r in expected_rows])
    @test got_hits == expected_hits
    @test got_links == expected_links

    function link_metrics_from_fixture_config(cfg_path::String)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        row = first(meteo.rows)
        sky = ArchimedLight.compute_sky(row, cfg)
        turtle = ArchimedLight.build_turtle(cfg, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
        first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
        pair_counts, sun_hits, node_ids, _ = ArchimedLight._pair_counts_for_scattering(scene, turtle, cfg)
        all_hits = ArchimedLight._all_dir_hits_for_scattering(first_order, sun_hits, cfg, node_ids)

        neighbours_by_node = Dict{Int,Set{Int}}()
        hits_by_node = Dict{Int,Int}()
        undirected_edges = Dict{Tuple{Int,Int},Int}()
        for ((to, from), c) in pair_counts
            c > 0 || continue
            s = get!(neighbours_by_node, to) do
                Set{Int}()
            end
            push!(s, from)
            hits_by_node[to] = get(hits_by_node, to, 0) + c
            e = to < from ? (to, from) : (from, to)
            undirected_edges[e] = max(get(undirected_edges, e, 0), c)
        end

        return (
            links=sort([length(v) for v in values(neighbours_by_node) if !isempty(v)]),
            hits=sort([v for v in values(hits_by_node) if v > 0]),
            all_hits=sort([v for v in values(all_hits) if v > 0]),
            edge_counts=sort(collect(values(undirected_edges))),
        )
    end

    for links_name in ("test-links2", "test-links3", "test-links4")
        fx = fixtures[links_name]
        metrics = link_metrics_from_fixture_config(fx.config_path)
        expected_stats = read_java_csv(joinpath(fixture_path(fx), "expected", "log-nodelinks-stats-alldirs.csv"))
        expected_dir = read_java_csv(joinpath(fixture_path(fx), "expected", "log-nodelinks-dir00.csv"))
        expected_links = sort(Int[to_int(getproperty(r, :links)) for r in expected_stats])
        expected_hits = sort(Int[to_int(getproperty(r, :hits)) for r in expected_stats])
        expected_edge_counts = sort(Int[to_int(getproperty(r, :n)) for r in expected_dir])

        @test metrics.links == expected_links
        @test metrics.hits == expected_hits
        @test metrics.edge_counts == expected_edge_counts
    end

    # Java fixture `test-links-stats` compares node-link counts across two pixel sizes.
    # It does not compare our internal all-direction hit count proxy.
    links_stats_cfg = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml")
    cfg_links_stats = ArchimedLight.read_light_config(links_stats_cfg)
    cfg_links_stats_017 = ArchimedLight.LightConfig(
        cfg_links_stats.scene,
        cfg_links_stats.meteo,
        cfg_links_stats.all_in_turtle,
        cfg_links_stats.turtle_sectors,
        0.17 / 100.0,
        cfg_links_stats.area_ratio,
        cfg_links_stats.scattering,
        cfg_links_stats.scattering_max_iter,
        cfg_links_stats.scattering_stop_ratio,
        cfg_links_stats.scattering_coeff_par,
        cfg_links_stats.scattering_coeff_nir,
        cfg_links_stats.cache_radiation,
        cfg_links_stats.raw,
    )
    cfg_links_stats_003 = ArchimedLight.LightConfig(
        cfg_links_stats.scene,
        cfg_links_stats.meteo,
        cfg_links_stats.all_in_turtle,
        cfg_links_stats.turtle_sectors,
        0.03 / 100.0,
        cfg_links_stats.area_ratio,
        cfg_links_stats.scattering,
        cfg_links_stats.scattering_max_iter,
        cfg_links_stats.scattering_stop_ratio,
        cfg_links_stats.scattering_coeff_par,
        cfg_links_stats.scattering_coeff_nir,
        cfg_links_stats.cache_radiation,
        cfg_links_stats.raw,
    )
    function link_metrics_from_cfg(cfg::ArchimedLight.LightConfig)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        row = first(meteo.rows)
        sky = ArchimedLight.compute_sky(row, cfg)
        turtle = ArchimedLight.build_turtle(cfg, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
        first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
        pair_counts, sun_hits, node_ids, _ = ArchimedLight._pair_counts_for_scattering(scene, turtle, cfg)
        _ = ArchimedLight._all_dir_hits_for_scattering(first_order, sun_hits, cfg, node_ids)
        neighbours_by_node = Dict{Int,Set{Int}}()
        for ((to, from), c) in pair_counts
            c > 0 || continue
            s = get!(neighbours_by_node, to) do
                Set{Int}()
            end
            push!(s, from)
        end

        vertices, faces, face2node, _, plotbox, _ = ArchimedLight._scene_geometry_for_interception(scene, cfg)
        perdir_links = Int[]
        for sector in turtle.sectors
            pixel_hits, _, _, _ = ArchimedLight._direction_projection(vertices, faces, face2node, sector.direction, cfg, plotbox)
            dir_neigh = Dict{Int,Set{Int}}()
            for stack in values(pixel_hits)
                length(stack) <= 1 && continue
                sort!(stack, by=x -> x[1], rev=true)
                for h in 1:(length(stack) - 1)
                    above = stack[h][2]
                    below = stack[h + 1][2]
                    s1 = get!(dir_neigh, above) do
                        Set{Int}()
                    end
                    s2 = get!(dir_neigh, below) do
                        Set{Int}()
                    end
                    push!(s1, below)
                    push!(s2, above)
                end
            end
            append!(perdir_links, [length(v) for v in values(dir_neigh) if !isempty(v)])
        end

        (
            links=sort([length(v) for v in values(neighbours_by_node) if !isempty(v)]),
            perdir_links=sort(perdir_links),
        )
    end
    m_links_017 = link_metrics_from_cfg(cfg_links_stats_017)
    m_links_003 = link_metrics_from_cfg(cfg_links_stats_003)
    @test m_links_017.links == m_links_003.links
    @test m_links_017.perdir_links == m_links_003.perdir_links

    f_div = fixtures["test-scattering-divergence"]
    cfg_div = ArchimedLight.read_light_config(f_div.config_path)
    scene_div = ArchimedLight.read_scene(cfg_div.scene)
    meteo_div = ArchimedLight.read_meteo(cfg_div.meteo)
    step_div = ArchimedLight.run_light_step(scene_div, first(meteo_div.rows), cfg_div)
    @test step_div.scattering !== nothing
    @test step_div.scattering.converged
    @test step_div.scattering.iterations <= cfg_div.scattering_max_iter
    @test all(isfinite, values(step_div.scattering.added_par_power_per_node))
    @test all(isfinite, values(step_div.scattering.added_nir_power_per_node))

    f_custom = fixtures["test-customband"]
    cfg_custom = ArchimedLight.read_light_config(f_custom.config_path)
    scene_custom = ArchimedLight.read_scene(cfg_custom.scene)
    meteo_custom = ArchimedLight.read_meteo(cfg_custom.meteo)
    row_custom = first(meteo_custom.rows)
    step_custom = ArchimedLight.run_light_step(scene_custom, row_custom, cfg_custom)
    @test haskey(step_custom.extra_band_irradiance, "CUSTOM")
    @test haskey(step_custom.budget.extra_0_q_per_band, "CUSTOM")
    @test haskey(step_custom.budget.extra_q_per_band, "CUSTOM")

    _, _, _, _, plotbox_custom, _ = ArchimedLight._scene_geometry_for_interception(scene_custom, cfg_custom)
    plot_area = plotbox_custom.xdim * plotbox_custom.ydim
    step_seconds = _step_duration_seconds(row_custom)
    custom_irr = step_custom.extra_band_irradiance["CUSTOM"]

    scene_incid_par = step_custom.sky.ri_par_f * plot_area * step_seconds
    scene_incid_custom = custom_irr * plot_area * step_seconds
    scene_par_0 = sum(values(step_custom.budget.ri_par_0_q_per_node))
    scene_par_n = sum(values(step_custom.budget.ri_par_q_per_node))
    scene_custom_0 = sum(values(step_custom.budget.extra_0_q_per_band["CUSTOM"]))
    scene_custom_n = sum(values(step_custom.budget.extra_q_per_band["CUSTOM"]))

    r_par_0 = scene_par_0 / max(scene_incid_par, eps(Float64))
    r_custom_0 = scene_custom_0 / max(scene_incid_custom, eps(Float64))
    r_par_n = scene_par_n / max(scene_incid_par, eps(Float64))
    r_custom_n = scene_custom_n / max(scene_incid_custom, eps(Float64))
    @test relerr(r_par_0, r_custom_0) < 4e-4
    @test relerr(r_par_n, r_custom_n) < 4e-4

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

    f_disk = fixtures["test-save_on_disk1"]
    cfg_disk_base = ArchimedLight.read_light_config(f_disk.config_path)
    cfg_disk_on = with_cache_pixel_table(cfg_disk_base, true)
    cfg_disk_off = with_cache_pixel_table(cfg_disk_base, false)
    cache_dir = ArchimedLight._projection_cache_dir(cfg_disk_on)
    isdir(cache_dir) && rm(cache_dir; recursive=true, force=true)

    scene_disk = ArchimedLight.read_scene(cfg_disk_on.scene)
    row_disk = first(ArchimedLight.read_meteo(cfg_disk_on.meteo).rows)

    step_disk_1 = ArchimedLight.run_light_step(scene_disk, row_disk, cfg_disk_on)
    @test isdir(cache_dir)
    files_1 = sort(filter(f -> endswith(f, ".jls"), readdir(cache_dir; join=true)))
    @test !isempty(files_1)
    mtimes_1 = Dict(f => stat(f).mtime for f in files_1)

    step_disk_2 = ArchimedLight.run_light_step(scene_disk, row_disk, cfg_disk_on)
    files_2 = sort(filter(f -> endswith(f, ".jls"), readdir(cache_dir; join=true)))
    @test basename.(files_2) == basename.(files_1)
    mtimes_2 = Dict(f => stat(f).mtime for f in files_2)
    @test all(mtimes_2[f] == mtimes_1[f] for f in files_1)

    step_mem = ArchimedLight.run_light_step(scene_disk, row_disk, cfg_disk_off)
    @test max_abs_float_dict_diff(step_disk_2.first_order.projected_area_per_node, step_mem.first_order.projected_area_per_node) == 0.0
    @test max_abs_int_dict_diff(step_disk_2.first_order.hits_per_node, step_mem.first_order.hits_per_node) == 0
    @test max_abs_float_dict_diff(step_disk_2.budget.ri_par_q_per_node, step_mem.budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step_disk_2.budget.ri_nir_q_per_node, step_mem.budget.ri_nir_q_per_node) == 0.0
    rm(cache_dir; recursive=true, force=true)

    f_scat_cache = fixtures["test-scattering-one-plate"]
    cfg_scat_base = ArchimedLight.read_light_config(f_scat_cache.config_path)
    cfg_scat_on = with_cache_pixel_table(cfg_scat_base, true)
    cfg_scat_off = with_cache_pixel_table(cfg_scat_base, false)
    scat_cache_dir = ArchimedLight._projection_cache_dir(cfg_scat_on)
    isdir(scat_cache_dir) && rm(scat_cache_dir; recursive=true, force=true)

    scene_scat = ArchimedLight.read_scene(cfg_scat_on.scene)
    row_scat = first(ArchimedLight.read_meteo(cfg_scat_on.meteo).rows)
    sky_scat = ArchimedLight.compute_sky(row_scat, cfg_scat_on)
    turtle_scat = ArchimedLight.build_turtle(cfg_scat_on, sky_scat)
    flux_scat = ArchimedLight.compute_directional_fluxes(sky_scat, turtle_scat, cfg_scat_on)
    first_scat_on = ArchimedLight.compute_first_order(scene_scat, turtle_scat, flux_scat, cfg_scat_on)
    @test isdir(scat_cache_dir)
    scat_files_1 = sort(filter(f -> endswith(f, ".jls"), readdir(scat_cache_dir; join=true)))
    @test !isempty(scat_files_1)
    scat_mtimes_1 = Dict(f => stat(f).mtime for f in scat_files_1)
    scat_on = ArchimedLight.compute_scattering(scene_scat, turtle_scat, first_scat_on, cfg_scat_on)
    scat_files_2 = sort(filter(f -> endswith(f, ".jls"), readdir(scat_cache_dir; join=true)))
    @test basename.(scat_files_2) == basename.(scat_files_1)
    scat_mtimes_2 = Dict(f => stat(f).mtime for f in scat_files_2)
    @test all(scat_mtimes_2[f] == scat_mtimes_1[f] for f in scat_files_1)

    first_scat_off = ArchimedLight.compute_first_order(scene_scat, turtle_scat, flux_scat, cfg_scat_off)
    scat_off = ArchimedLight.compute_scattering(scene_scat, turtle_scat, first_scat_off, cfg_scat_off)
    @test max_abs_float_dict_diff(first_scat_on.projected_area_per_node, first_scat_off.projected_area_per_node) == 0.0
    @test max_abs_int_dict_diff(first_scat_on.hits_per_node, first_scat_off.hits_per_node) == 0
    @test max_abs_float_dict_diff(scat_on.added_par_power_per_node, scat_off.added_par_power_per_node) == 0.0
    @test max_abs_float_dict_diff(scat_on.added_nir_power_per_node, scat_off.added_nir_power_per_node) == 0.0
    rm(scat_cache_dir; recursive=true, force=true)
end
