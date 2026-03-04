using Statistics
using LinearAlgebra: norm
using Rotations: RotZ, AngleAxis, RotMatrix
using StaticArrays: SVector
using GeometryBasics
import PlantGeom

const _TEST_PROFILE = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_PROFILE", "all"))
const _RUN_CORE_TESTS = _TEST_PROFILE in ("all", "core")
const _RUN_PARITY_TESTS = _TEST_PROFILE in ("all", "parity")
const _RUN_GUARD_TESTS = _TEST_PROFILE == "guard"
const _RUN_SKY_MATRIX_TESTS = _TEST_PROFILE == "sky_matrix"

if _RUN_CORE_TESTS
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
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, cfg)
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
    expected_component_count = length(ArchimedLight._interception_java_keys(scene, cfg))
    @test tbl.columns == out_cols
    @test length(tbl.rows) == expected_component_count
    @test all(haskey(r, "Ri_PAR_0_q") for r in tbl.rows)
    @test all(haskey(r, "Ra_PAR_0_q") for r in tbl.rows)
    @test all(haskey(r, "Ri_custom_q") for r in tbl.rows)
    @test all(get(r, "type", "NA") != "NA" for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_x", NaN))) for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_y", NaN))) for r in tbl.rows)
    @test all(isfinite(Float64(get(r, "barycentre_z", NaN))) for r in tbl.rows)
    @test all(Float64(get(r, "sky_fraction", -1.0)) >= 0.0 for r in tbl.rows)

    # Java parity: default step numbering is 0-based.
    tbl_default_step = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, columns=["step_number"])
    @test all(Int(r["step_number"]) == 0 for r in tbl_default_step.rows)

    mktempdir() do tmp
        out_csv = joinpath(tmp, "component_values.csv")
        ArchimedLight.write_component_values_csv(out_csv, scene, step, cfg; meteo_row=row, step_number=1, columns=out_cols)
        @test isfile(out_csv)
        rows_csv = read_java_csv(out_csv)
        @test length(rows_csv) == expected_component_count
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
    scene_tbl_default_step = ArchimedLight.scene_values_table(scene, series, cfg; meteo_rows=subset.rows, columns=["step_number"])
    @test [Int(r["step_number"]) for r in scene_tbl_default_step.rows] == collect(0:(length(series) - 1))

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

    raw_scat = copy(cfg.raw)
    raw_scat["component_variables"] = Dict("Ri_PAR_f" => true)
    cfg_scat_off = ArchimedLight.LightConfig(
        cfg.scene,
        cfg.meteo,
        cfg.all_in_turtle,
        cfg.turtle_sectors,
        cfg.pixel_size,
        cfg.area_ratio,
        false,
        cfg.scattering_max_iter,
        cfg.scattering_stop_ratio,
        cfg.scattering_coeff_par,
        cfg.scattering_coeff_nir,
        cfg.cache_radiation,
        raw_scat,
    )
    @test_throws ErrorException ArchimedLight.component_variable_names(cfg_scat_off)

    raw_photo = copy(cfg.raw)
    raw_photo["photosynthesis"] = true
    raw_photo["component_variables"] = Dict("An_f" => true)
    cfg_photo = with_raw_overrides(cfg, raw_photo)
    @test_throws ErrorException ArchimedLight.component_variable_names(cfg_photo)

    raw_enb = copy(cfg.raw)
    raw_enb["photosynthesis"] = true
    raw_enb["energy_balance"] = true
    raw_enb["component_variables"] = Dict("Ri_TIR_f" => true)
    cfg_enb = with_raw_overrides(cfg, raw_enb)
    @test_throws ErrorException ArchimedLight.component_variable_names(cfg_enb)

    raw_ok = copy(cfg.raw)
    raw_ok["component_variables"] = Dict("Ri_PAR_0_f" => true, "Ra_PAR_0_f" => true)
    cfg_ok = with_raw_overrides(cfg, raw_ok)
    @test ArchimedLight.component_variable_names(cfg_ok) == ["Ri_PAR_0_f", "Ra_PAR_0_f"]

    raw_default = copy(cfg.raw)
    delete!(raw_default, "component_variables")
    cfg_default = with_raw_overrides(cfg, raw_default)
    cols_default = ArchimedLight.component_variable_names(cfg_default)
    @test !("Ri_TIR_f" in cols_default)
    @test !("Ri_TIR_q" in cols_default)
    @test !("Ra_TIR_f" in cols_default)
    @test !("Ra_TIR_q" in cols_default)

    @test_throws ErrorException ArchimedLight.component_values_table(scene, step, cfg_enb; meteo_row=row, columns=["Ri_TIR_f"])
    @test_throws ErrorException ArchimedLight.component_values_table(scene, step, cfg_photo; meteo_row=row, columns=["An_f"])

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

@testset "Projection tie ownership at toric borders" begin
    p1 = SVector(-0.2, -0.2, 1.0)
    p2 = SVector(1.8, -0.2, 1.0)
    p3 = SVector(0.8, 1.8, 1.0)
    direction = SVector(0.0, 0.0, 1.0)

    pixel_hits = Dict{Int,Vector{Tuple{Float64,Int}}}()
    node_hits = Dict{Int,Int}()
    projected_mesh_area = Dict{Int,Float64}()
    projected_pixels_area = Dict{Int,Float64}()

    ArchimedLight._project_triangle!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        101,
        p1,
        p2,
        p3,
        direction,
        1.5,
        0.0,
        1.0,
        1.0,
        1.0,
        4,
        4,
        false,
        true,
    )
    @test isempty(pixel_hits)

    empty!(pixel_hits)
    empty!(node_hits)
    empty!(projected_mesh_area)
    empty!(projected_pixels_area)

    ArchimedLight._project_triangle!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        101,
        p1,
        p2,
        p3,
        direction,
        1.5,
        0.0,
        1.0,
        1.0,
        1.0,
        4,
        4,
        true,
        true,
    )

    ArchimedLight._project_triangle!(
        pixel_hits,
        node_hits,
        projected_mesh_area,
        projected_pixels_area,
        202,
        p1,
        p2,
        p3,
        direction,
        1.5,
        0.0,
        1.0,
        1.0,
        1.0,
        4,
        4,
        true,
        true,
    )

    @test haskey(pixel_hits, 4)
    @test pixel_hits[4][1][2] == 101
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

    @testset "OPS inclination transform parity" begin
        opf = joinpath(dirname(@__DIR__), "example", "scene", "opf", "simple_OPF_shapes.opf")
        @test isfile(opf)

        mktempdir() do tmp
            function _write_ops(path::AbstractString, inclination_angle_rad::Float64)
                open(path, "w") do io
                    println(io, "T 0 0 0 1 1 flat")
                    println(io, "1 1 $(opf) 1.25 -2.0 0.5 0.37 $(inclination_angle_rad) 0 0")
                    println(io, "-1 1 1")
                end
            end

            ops_base = joinpath(tmp, "base.ops")
            ops_inclined = joinpath(tmp, "inclined.ops")
            _write_ops(ops_base, 0.0)
            _write_ops(ops_inclined, 0.24)

            scene_base = ArchimedLight.read_scene(ops_base)
            scene_inclined = ArchimedLight.read_scene(ops_inclined)

            @test Set(keys(scene_base.total_area_per_node)) == Set(keys(scene_inclined.total_area_per_node))
            @test scene_base.java_item_id_per_node == scene_inclined.java_item_id_per_node
            @test scene_base.java_component_id_per_node == scene_inclined.java_component_id_per_node

            pivot = SVector{3,Float64}(1.25, -2.0, 0.5)
            axis = RotZ(0.37) * SVector{3,Float64}(0.0, 1.0, 0.0)
            axis_u = axis / norm(axis)
            rot = RotMatrix(AngleAxis(0.24, axis_u[1], axis_u[2], axis_u[3]))

            for nid in keys(scene_base.total_area_per_node)
                @test isapprox(
                    scene_inclined.total_area_per_node[nid],
                    scene_base.total_area_per_node[nid];
                    atol=1e-12,
                    rtol=1e-12,
                )

                b0 = scene_base.barycenter_per_node[nid]
                b1 = scene_inclined.barycenter_per_node[nid]
                expected = rot * (SVector{3,Float64}(b0[1], b0[2], b0[3]) - pivot) + pivot
                @test isapprox(b1[1], expected[1]; atol=1e-9, rtol=1e-9)
                @test isapprox(b1[2], expected[2]; atol=1e-9, rtol=1e-9)
                @test isapprox(b1[3], expected[3]; atol=1e-9, rtol=1e-9)
            end
        end
    end
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
        (0.0, 0.0, 1.0),
        (1.5637987e-7, 0.89438856, 0.44729084),
        (0.85061413, 0.27638108, 0.44729084),
        (0.5257086, -0.7235755, 0.4472909),
        (-0.5257085, -0.7235755, 0.44729084),
        (-0.85061413, 0.2763814, 0.4472909),
    ]
    @test length(got6) == length(exp6)
    for i in eachindex(exp6)
        @test abs(got6[i][1] - exp6[i][1]) < 1e-6
        @test abs(got6[i][2] - exp6[i][2]) < 1e-6
        @test abs(got6[i][3] - exp6[i][3]) < 1e-6
    end

    for n in (1, 6, 16, 46)
        t = ArchimedLight.build_turtle(with_sectors(cfg_ref, n), sky_ref)
        for sector in t.sectors
            @test abs(norm(sector.direction) - 1.0) < 1e-6
            @test sector.direction[3] >= -1e-8
        end
    end
end

@testset "Sky sector energy unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))

    function with_turtle(cfg::ArchimedLight.LightConfig; sectors::Int=cfg.turtle_sectors, all_in_turtle::Bool=cfg.all_in_turtle, scattering::Bool=cfg.scattering)
        raw = copy(cfg.raw)
        raw["all_in_turtle"] = all_in_turtle
        raw["sky_sectors"] = sectors
        raw["scattering"] = scattering
        ArchimedLight.LightConfig(
            cfg.scene,
            cfg.meteo,
            all_in_turtle,
            sectors,
            cfg.pixel_size,
            cfg.area_ratio,
            scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            cfg.cache_radiation,
            raw,
        )
    end

    sun_direction(azimuth_deg::Float64, elevation_deg::Float64) = begin
        az = deg2rad(azimuth_deg)
        el = deg2rad(elevation_deg)
        SVector{3,Float64}(cos(el) * sin(az), cos(el) * cos(az), sin(el))
    end

    sky_split = ArchimedLight.SkyState(215.0, 52.0, 240.0, 260.0, 0.4, 0.6)
    cfg_split = with_turtle(cfg_ref; sectors=6, all_in_turtle=false, scattering=false)
    turtle_split = ArchimedLight.build_turtle(cfg_split, sky_split)
    flux_split = ArchimedLight.compute_directional_fluxes(sky_split, turtle_split, cfg_split)
    sky_ids = findall(s -> s.source == :sky, turtle_split.sectors)
    sun_ids = findall(s -> s.source == :sun, turtle_split.sectors)
    @test length(turtle_split.sectors) == 7
    @test length(sun_ids) == 1
    @test isapprox(sum(flux_split.par), sky_split.ri_par_f; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(flux_split.nir), sky_split.ri_nir_f; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(flux_split.par[sky_ids]), sky_split.ri_par_f * sky_split.diffuse_fraction; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(flux_split.nir[sky_ids]), sky_split.ri_nir_f * sky_split.diffuse_fraction; atol=1e-12, rtol=1e-12)
    @test isapprox(flux_split.par[first(sun_ids)], sky_split.ri_par_f * sky_split.direct_fraction; atol=1e-12, rtol=1e-12)
    @test isapprox(flux_split.nir[first(sun_ids)], sky_split.ri_nir_f * sky_split.direct_fraction; atol=1e-12, rtol=1e-12)
    @test maximum(abs.(turtle_split.sectors[first(sun_ids)].direction .- sun_direction(sky_split.sun_azimuth_deg, sky_split.sun_elevation_deg))) < 1e-12

    sky_diffuse = ArchimedLight.SkyState(120.0, 35.0, 180.0, 120.0, 0.0, 1.0)
    cfg_diffuse = with_turtle(cfg_ref; sectors=16, all_in_turtle=true, scattering=false)
    turtle_diffuse = ArchimedLight.build_turtle(cfg_diffuse, sky_diffuse)
    flux_diffuse = ArchimedLight.compute_directional_fluxes(sky_diffuse, turtle_diffuse, cfg_diffuse)
    @test length(turtle_diffuse.sectors) == 16
    @test all(s -> s.source == :sky, turtle_diffuse.sectors)
    @test isapprox(sum(flux_diffuse.par), sky_diffuse.ri_par_f; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(flux_diffuse.nir), sky_diffuse.ri_nir_f; atol=1e-12, rtol=1e-12)
    @test all(>=(0.0), flux_diffuse.par)
    @test all(>=(0.0), flux_diffuse.nir)
    par_weights = flux_diffuse.par ./ sum(flux_diffuse.par)
    nir_weights = flux_diffuse.nir ./ sum(flux_diffuse.nir)
    @test maximum(abs.(par_weights .- nir_weights)) < 1e-12
end

function _synthetic_plate_cfg(
    cfg::ArchimedLight.LightConfig;
    sectors::Int=1,
    all_in_turtle::Bool=false,
    scattering::Bool=false,
    pixel_size::Float64=cfg.pixel_size,
    models::Vector{String}=String[],
)
    raw = copy(cfg.raw)
    raw["all_in_turtle"] = all_in_turtle
    raw["sky_sectors"] = sectors
    raw["scattering"] = scattering
    raw["pixel_size"] = pixel_size * 100.0
    raw["models"] = models
    ArchimedLight.LightConfig(
        cfg.scene,
        cfg.meteo,
        all_in_turtle,
        sectors,
        pixel_size,
        cfg.area_ratio,
        scattering,
        cfg.scattering_max_iter,
        cfg.scattering_stop_ratio,
        cfg.scattering_coeff_par,
        cfg.scattering_coeff_nir,
        cfg.cache_radiation,
        raw,
    )
end

function _synthetic_quad_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float32}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    total_area_per_node = Dict{Int,Float64}()
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    java_item_id_per_node = Dict{Int,Int}()
    java_component_id_per_node = Dict{Int,Int}()

    xmins = Float64[]
    ymins = Float64[]
    xmaxs = Float64[]
    ymaxs = Float64[]

    for (i, spec) in enumerate(specs)
        p1 = ntuple(j -> Float64(spec.p1[j]), 3)
        p2 = ntuple(j -> Float64(spec.p2[j]), 3)
        p3 = ntuple(j -> Float64(spec.p3[j]), 3)
        p4 = ntuple(j -> Float64(spec.p4[j]), 3)
        item_id = get(spec, :item_id, i)
        component_id = get(spec, :component_id, 1)
        group = String(get(spec, :group, "plate"))
        type = String(get(spec, :type, "plate"))
        xs = (p1[1], p2[1], p3[1], p4[1])
        ys = (p1[2], p2[2], p3[2], p4[2])

        append!(xmins, xs)
        append!(ymins, ys)
        append!(xmaxs, xs)
        append!(ymaxs, ys)

        base = length(points)
        append!(
            points,
            GeometryBasics.Point{3,Float32}[
                GeometryBasics.Point{3,Float32}(Float32(p1[1]), Float32(p1[2]), Float32(p1[3])),
                GeometryBasics.Point{3,Float32}(Float32(p2[1]), Float32(p2[2]), Float32(p2[3])),
                GeometryBasics.Point{3,Float32}(Float32(p3[1]), Float32(p3[2]), Float32(p3[3])),
                GeometryBasics.Point{3,Float32}(Float32(p4[1]), Float32(p4[2]), Float32(p4[3])),
            ],
        )
        append!(faces, PlantGeom.Face3[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, [i, i])
        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        total_area_per_node[i] = area1 + area2
        barycenter_per_node[i] = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            (p1[3] + p2[3] + p3[3] + p4[3]) / 4,
        )
        node_group[i] = group
        node_type[i] = type
        java_item_id_per_node[i] = item_id
        java_component_id_per_node[i] = component_id
    end

    ArchimedLight.SceneGeometry(
        nothing,
        GeometryBasics.Mesh(points, faces),
        face2node,
        total_area_per_node,
        barycenter_per_node,
        node_group,
        node_type,
        java_item_id_per_node,
        java_component_id_per_node,
        "synthetic_horizontal_planes",
        (minimum(xmins), minimum(ymins), maximum(xmaxs), maximum(ymaxs)),
    )
end

function _synthetic_horizontal_scene(specs::AbstractVector{<:NamedTuple})
    quad_specs = map(specs) do spec
        (
            p1=(spec.x0, spec.y0, spec.z),
            p2=(spec.x1, spec.y0, spec.z),
            p3=(spec.x1, spec.y1, spec.z),
            p4=(spec.x0, spec.y1, spec.z),
            item_id=get(spec, :item_id, 1),
            component_id=get(spec, :component_id, 1),
            group=get(spec, :group, "plate"),
            type=get(spec, :type, "plate"),
        )
    end
    _synthetic_quad_scene(quad_specs)
end

function _synthetic_meteo_row(;
    date::Dates.Date=Dates.Date(2020, 6, 21),
    start_time::Dates.Time=Dates.Time(12),
    duration_seconds::Float64=1.0,
    latitude::Float64=0.0,
    relative_humidity::Float64=60.0,
    ri_par_f::Float64=100.0,
    ri_nir_f::Float64=0.0,
    direct_fraction::Float64=1.0,
    sun_azimut::Float64=180.0,
    sun_elevation::Float64=90.0,
    use::String="relativeHumidity RI_PAR_f",
)
    start_dt = Dates.DateTime(date, start_time)
    end_dt = start_dt + Dates.Millisecond(round(Int, duration_seconds * 1000))
    (
        date=date,
        hour_start=Dates.format(Dates.Time(start_dt), Dates.DateFormat("HH:MM:SS")),
        hour_end=Dates.format(Dates.Time(end_dt), Dates.DateFormat("HH:MM:SS")),
        step_duration=duration_seconds,
        latitude=latitude,
        relativeHumidity=relative_humidity,
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
        use=use,
    )
end

@testset "Single plate interception unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))

    scene = _synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate")])

    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
    total_area = sum(values(scene.total_area_per_node))
    total_projected_area = sum(values(first_order.projected_area_per_node))
    total_incident_power = sum(values(first_order.incident_par_power_per_node))

    @test isapprox(total_projected_area, total_area; atol=1e-12, rtol=1e-12)
    @test isapprox(total_incident_power, sky.ri_par_f * total_area; atol=1e-10, rtol=1e-10)

    budget = ArchimedLight.integrate_light(
        first_order,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    step = ArchimedLight.LightStepResult(sky, turtle, fluxes, first_order, nothing, budget, Dict{String,Float64}())
    rows = ArchimedLight.component_values_table(
        scene,
        step,
        cfg;
        step_duration_seconds=1.0,
        columns=["item_id", "component_id", "area", "Ri_PAR_0_f", "Ri_PAR_0_q"],
    ).rows
    @test length(rows) == length(scene.total_area_per_node)
    @test all(isapprox(Float64(r["Ri_PAR_0_f"]), sky.ri_par_f; atol=1e-10, rtol=1e-10) for r in rows)
    @test isapprox(sum(Float64(r["Ri_PAR_0_q"]) for r in rows), total_incident_power; atol=1e-10, rtol=1e-10)
end

@testset "Two stacked plates scattering unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))

    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate"),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate"),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    cfg_no_scat = _synthetic_plate_cfg(cfg_ref; sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01)
    turtle_no_scat = ArchimedLight.build_turtle(cfg_no_scat, sky)
    flux_no_scat = ArchimedLight.compute_directional_fluxes(sky, turtle_no_scat, cfg_no_scat)
    first_no_scat = ArchimedLight.compute_first_order(scene, turtle_no_scat, flux_no_scat, cfg_no_scat)

    @test isapprox(get(first_no_scat.projected_area_per_node, 1, 0.0), 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_no_scat.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_no_scat.incident_par_power_per_node, 1, 0.0), sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(get(first_no_scat.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    budget_no_scat = ArchimedLight.integrate_light(
        first_no_scat,
        nothing,
        cfg_no_scat;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget_no_scat.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(budget_no_scat.ri_par_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    cfg_scat = _synthetic_plate_cfg(cfg_ref; sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01)
    turtle_scat = ArchimedLight.build_turtle(cfg_scat, sky)
    flux_scat = ArchimedLight.compute_directional_fluxes(sky, turtle_scat, cfg_scat)
    first_scat = ArchimedLight.compute_first_order(scene, turtle_scat, flux_scat, cfg_scat)
    scat = ArchimedLight.compute_scattering(scene, turtle_scat, first_scat, cfg_scat)

    @test get(scat.added_par_power_per_node, 2, 0.0) > 0.0
    @test get(scat.added_par_power_per_node, 1, 0.0) > 0.0

    budget_scat = ArchimedLight.integrate_light(
        first_scat,
        scat,
        cfg_scat;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget_scat.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test get(budget_scat.ri_par_q_per_node, 2, 0.0) > 0.0
    @test get(budget_scat.ri_par_q_per_node, 2, 0.0) ≈ get(scat.added_par_power_per_node, 2, 0.0) atol=1e-10 rtol=1e-10
    @test get(budget_scat.ri_par_q_per_node, 1, 0.0) > get(budget_scat.ri_par_0_q_per_node, 1, 0.0)

    step_scat = ArchimedLight.LightStepResult(sky, turtle_scat, flux_scat, first_scat, scat, budget_scat, Dict{String,Float64}())
    rows = ArchimedLight.component_values_table(
        scene,
        step_scat,
        cfg_scat;
        step_duration_seconds=1.0,
        columns=["item_id", "component_id", "Ri_PAR_0_q", "Ri_PAR_q"],
    ).rows
    row_by_item = Dict(Int(r["item_id"]) => r for r in rows)
    @test isapprox(Float64(row_by_item[2]["Ri_PAR_0_q"]), 0.0; atol=1e-12, rtol=1e-12)
    @test Float64(row_by_item[2]["Ri_PAR_q"]) > 0.0
    @test Float64(row_by_item[1]["Ri_PAR_q"]) > Float64(row_by_item[1]["Ri_PAR_0_q"])
end

@testset "Partial overlap plates direct interception unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=0.5, y0=0.0, y1=1.0, z=1.0, group="upper_half", type="plate"),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower_full", type="plate"),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    @test isapprox(get(first_order.projected_area_per_node, 1, 0.0), 0.5; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.projected_area_per_node, 2, 0.0), 0.5; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.incident_par_power_per_node, 1, 0.0), 50.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(first_order.incident_par_power_per_node, 2, 0.0), 50.0; atol=1e-9, rtol=1e-9)

    budget = ArchimedLight.integrate_light(
        first_order,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), 50.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget.ri_par_0_q_per_node, 2, 0.0), 50.0; atol=1e-9, rtol=1e-9)

    step = ArchimedLight.LightStepResult(sky, turtle, fluxes, first_order, nothing, budget, Dict{String,Float64}())
    rows = ArchimedLight.component_values_table(
        scene,
        step,
        cfg;
        step_duration_seconds=1.0,
        columns=["item_id", "area", "Ri_PAR_0_q"],
    ).rows
    row_by_item = Dict(Int(r["item_id"]) => r for r in rows)
    @test isapprox(Float64(row_by_item[1]["area"]), 0.5; atol=1e-12, rtol=1e-12)
    @test isapprox(Float64(row_by_item[2]["area"]), 1.0; atol=1e-12, rtol=1e-12)
    @test isapprox(Float64(row_by_item[1]["Ri_PAR_0_q"]), 50.0; atol=1e-9, rtol=1e-9)
    @test isapprox(Float64(row_by_item[2]["Ri_PAR_0_q"]), 50.0; atol=1e-9, rtol=1e-9)
end

@testset "Tilted plate cosine projection unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    tilt_deg = 60.0
    c = cosd(tilt_deg)
    s = sind(tilt_deg)
    scene = _synthetic_quad_scene([
        (
            p1=(0.0, 0.0, 0.0),
            p2=(1.0, 0.0, 0.0),
            p3=(1.0, c, s),
            p4=(0.0, c, s),
            group="tilted",
            type="plate",
        ),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    @test isapprox(scene.total_area_per_node[1], 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.projected_area_per_node, 1, 0.0), c; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.incident_par_power_per_node, 1, 0.0), sky.ri_par_f * c; atol=1e-9, rtol=1e-9)

    budget = ArchimedLight.integrate_light(
        first_order,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget.ri_par_0_f_per_node, 1, 0.0), sky.ri_par_f * c; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), sky.ri_par_f * c; atol=1e-9, rtol=1e-9)
end

@testset "No-overlap stacked plates direct interception unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=1.5, x1=2.5, y0=0.0, y1=1.0, z=1.0, group="upper_offset", type="plate"),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower_full", type="plate"),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    @test isapprox(get(first_order.projected_area_per_node, 1, 0.0), 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.projected_area_per_node, 2, 0.0), 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.incident_par_power_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(first_order.incident_par_power_per_node, 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
end

@testset "Oblique sun lateral shadow unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=1.0, x1=2.0, y0=0.0, y1=1.0, z=1.0, group="upper_offset", type="plate", item_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0, group="lower_target", type="plate", item_id=2),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)

    sky_shadow = ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0)
    turtle_shadow = ArchimedLight.build_turtle(cfg, sky_shadow)
    flux_shadow = ArchimedLight.compute_directional_fluxes(sky_shadow, turtle_shadow, cfg)
    first_shadow = ArchimedLight.compute_first_order(scene, turtle_shadow, flux_shadow, cfg)

    expected_horizontal_power = sky_shadow.ri_par_f
    @test isapprox(get(first_shadow.incident_par_power_per_node, 1, 0.0), expected_horizontal_power; atol=1e-5, rtol=1e-7)
    @test isapprox(get(first_shadow.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)

    sky_clear = ArchimedLight.SkyState(0.0, 45.0, 100.0, 0.0, 1.0, 0.0)
    turtle_clear = ArchimedLight.build_turtle(cfg, sky_clear)
    flux_clear = ArchimedLight.compute_directional_fluxes(sky_clear, turtle_clear, cfg)
    first_clear = ArchimedLight.compute_first_order(scene, turtle_clear, flux_clear, cfg)

    @test isapprox(get(first_clear.incident_par_power_per_node, 1, 0.0), expected_horizontal_power; atol=1e-5, rtol=1e-7)
    @test isapprox(get(first_clear.incident_par_power_per_node, 2, 0.0), expected_horizontal_power; atol=1e-5, rtol=1e-7)
    @test get(first_shadow.incident_par_power_per_node, 2, 0.0) < 0.02 * get(first_clear.incident_par_power_per_node, 2, 0.0)
end

@testset "Vertical plate zenith unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_quad_scene([
        (
            p1=(0.0, 0.0, 0.0),
            p2=(0.0, 1.0, 0.0),
            p3=(0.0, 1.0, 1.0),
            p4=(0.0, 0.0, 1.0),
            group="vertical",
            type="plate",
            item_id=1,
        ),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    @test isapprox(scene.total_area_per_node[1], 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.projected_area_per_node, 1, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(first_order.incident_par_power_per_node, 1, 0.0), 0.0; atol=1e-10, rtol=1e-10)

    budget = ArchimedLight.integrate_light(
        first_order,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget.ri_par_0_f_per_node, 1, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), 0.0; atol=1e-10, rtol=1e-10)
end

@testset "Diffuse-only stacked plates unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", item_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", item_id=2),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=16, all_in_turtle=true, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 0.0, 1.0)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    upper_power = get(first_order.incident_par_power_per_node, 1, 0.0)
    lower_power = get(first_order.incident_par_power_per_node, 2, 0.0)
    @test upper_power > 90.0
    @test lower_power > 0.1
    @test lower_power < 1.0
    @test lower_power < upper_power

    budget = ArchimedLight.integrate_light(
        first_order,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), upper_power; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget.ri_par_0_q_per_node, 2, 0.0), lower_power; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget.ri_par_q_per_node, 1, 0.0), upper_power; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget.ri_par_q_per_node, 2, 0.0), lower_power; atol=1e-9, rtol=1e-9)
end

@testset "Virtual sensor transparency unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensors", type="pave", item_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", item_id=2),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    cfg_opaque = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    turtle_opaque = ArchimedLight.build_turtle(cfg_opaque, sky)
    flux_opaque = ArchimedLight.compute_directional_fluxes(sky, turtle_opaque, cfg_opaque)
    first_opaque = ArchimedLight.compute_first_order(scene, turtle_opaque, flux_opaque, cfg_opaque)
    budget_opaque = ArchimedLight.integrate_light(
        first_opaque,
        nothing,
        cfg_opaque;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )

    @test isapprox(get(first_opaque.incident_par_power_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(first_opaque.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(budget_opaque.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    mktempdir() do tmp
        model_sensor = joinpath(tmp, "model_sensor.yml")
        write(
            model_sensor,
            """
            ---
            Group: sensors
            Type:
              pave:
                  Interception:
                      model: VirtualSensor
            ...
            """,
        )

        cfg_sensor = _synthetic_plate_cfg(
            cfg_ref;
            sectors=1,
            all_in_turtle=false,
            scattering=false,
            pixel_size=0.01,
            models=[model_sensor],
        )
        turtle_sensor = ArchimedLight.build_turtle(cfg_sensor, sky)
        flux_sensor = ArchimedLight.compute_directional_fluxes(sky, turtle_sensor, cfg_sensor)
        first_sensor = ArchimedLight.compute_first_order(scene, turtle_sensor, flux_sensor, cfg_sensor)
        budget_sensor = ArchimedLight.integrate_light(
            first_sensor,
            nothing,
            cfg_sensor;
            step_duration_seconds=1.0,
            component_area_per_node=scene.total_area_per_node,
        )

        @test isapprox(get(first_sensor.incident_par_power_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(first_sensor.incident_par_power_per_node, 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(budget_sensor.ri_par_0_q_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(budget_sensor.ri_par_0_q_per_node, 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
        @test get(budget_sensor.ri_par_0_q_per_node, 2, 0.0) > get(budget_opaque.ri_par_0_q_per_node, 2, 0.0)

        step_sensor = ArchimedLight.LightStepResult(sky, turtle_sensor, flux_sensor, first_sensor, nothing, budget_sensor, Dict{String,Float64}())
        rows_sensor = ArchimedLight.component_values_table(
            scene,
            step_sensor,
            cfg_sensor;
            step_duration_seconds=1.0,
            columns=["item_id", "Ri_PAR_0_q"],
        ).rows
        row_by_item = Dict(Int(r["item_id"]) => r for r in rows_sensor)
        @test isapprox(Float64(row_by_item[1]["Ri_PAR_0_q"]), 100.0; atol=1e-9, rtol=1e-9)
        @test isapprox(Float64(row_by_item[2]["Ri_PAR_0_q"]), 100.0; atol=1e-9, rtol=1e-9)
    end
end

@testset "NIR-only stacked plates unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", item_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", item_id=2),
    ])

    cfg_direct = _synthetic_plate_cfg(cfg_ref; sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky_direct = ArchimedLight.SkyState(180.0, 90.0, 0.0, 100.0, 1.0, 0.0)
    turtle_direct = ArchimedLight.build_turtle(cfg_direct, sky_direct)
    flux_direct = ArchimedLight.compute_directional_fluxes(sky_direct, turtle_direct, cfg_direct)
    first_direct = ArchimedLight.compute_first_order(scene, turtle_direct, flux_direct, cfg_direct)
    budget_direct = ArchimedLight.integrate_light(
        first_direct,
        nothing,
        cfg_direct;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )

    @test isapprox(get(first_direct.incident_par_power_per_node, 1, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(first_direct.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(first_direct.incident_nir_power_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(first_direct.incident_nir_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(budget_direct.ri_par_0_q_per_node, 1, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(budget_direct.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(budget_direct.ri_nir_0_q_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget_direct.ri_nir_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    cfg_diffuse = _synthetic_plate_cfg(cfg_ref; sectors=16, all_in_turtle=true, scattering=false, pixel_size=0.01)
    sky_diffuse = ArchimedLight.SkyState(180.0, 90.0, 0.0, 100.0, 0.0, 1.0)
    turtle_diffuse = ArchimedLight.build_turtle(cfg_diffuse, sky_diffuse)
    flux_diffuse = ArchimedLight.compute_directional_fluxes(sky_diffuse, turtle_diffuse, cfg_diffuse)
    first_diffuse = ArchimedLight.compute_first_order(scene, turtle_diffuse, flux_diffuse, cfg_diffuse)
    budget_diffuse = ArchimedLight.integrate_light(
        first_diffuse,
        nothing,
        cfg_diffuse;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )

    upper_nir = get(first_diffuse.incident_nir_power_per_node, 1, 0.0)
    lower_nir = get(first_diffuse.incident_nir_power_per_node, 2, 0.0)
    @test isapprox(get(first_diffuse.incident_par_power_per_node, 1, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(first_diffuse.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test upper_nir > 90.0
    @test lower_nir > 0.1
    @test lower_nir < 1.0
    @test lower_nir < upper_nir
    @test isapprox(get(budget_diffuse.ri_nir_0_q_per_node, 1, 0.0), upper_nir; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget_diffuse.ri_nir_0_q_per_node, 2, 0.0), lower_nir; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget_diffuse.ri_nir_q_per_node, 1, 0.0), upper_nir; atol=1e-9, rtol=1e-9)
    @test isapprox(get(budget_diffuse.ri_nir_q_per_node, 2, 0.0), lower_nir; atol=1e-9, rtol=1e-9)
end

@testset "Mixed direct and diffuse decomposition unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", item_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", item_id=2),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01)
    sky = ArchimedLight.SkyState(180.0, 90.0, 60.0, 40.0, 0.75, 0.25)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_all = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)

    sun_ids = Set(findall(s -> s.source == :sun, turtle.sectors))
    sky_ids = Set(findall(s -> s.source == :sky, turtle.sectors))
    split_fluxes = function (keep_ids)
        ArchimedLight.DirectionalFluxes(
            fluxes.sector_ids,
            [i in keep_ids ? fluxes.par[i] : 0.0 for i in eachindex(fluxes.par)],
            [i in keep_ids ? fluxes.nir[i] : 0.0 for i in eachindex(fluxes.nir)],
        )
    end

    first_sun = ArchimedLight.compute_first_order(scene, turtle, split_fluxes(sun_ids), cfg)
    first_sky = ArchimedLight.compute_first_order(scene, turtle, split_fluxes(sky_ids), cfg)

    @test length(sun_ids) == 1
    @test !isempty(sky_ids)
    @test isapprox(get(first_sun.incident_par_power_per_node, 1, 0.0), 45.0; atol=1e-6, rtol=1e-9)
    @test isapprox(get(first_sun.incident_nir_power_per_node, 1, 0.0), 30.0; atol=1e-6, rtol=1e-9)
    @test isapprox(get(first_sun.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(first_sun.incident_nir_power_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test get(first_sky.incident_par_power_per_node, 2, 0.0) > 0.0
    @test get(first_sky.incident_nir_power_per_node, 2, 0.0) > 0.0

    for nid in keys(scene.total_area_per_node)
        @test isapprox(
            get(first_all.incident_par_power_per_node, nid, 0.0),
            get(first_sun.incident_par_power_per_node, nid, 0.0) + get(first_sky.incident_par_power_per_node, nid, 0.0);
            atol=1e-10,
            rtol=1e-10,
        )
        @test isapprox(
            get(first_all.incident_nir_power_per_node, nid, 0.0),
            get(first_sun.incident_nir_power_per_node, nid, 0.0) + get(first_sky.incident_nir_power_per_node, nid, 0.0);
            atol=1e-10,
            rtol=1e-10,
        )
    end

    budget_all = ArchimedLight.integrate_light(
        first_all,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    budget_sun = ArchimedLight.integrate_light(
        first_sun,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )
    budget_sky = ArchimedLight.integrate_light(
        first_sky,
        nothing,
        cfg;
        step_duration_seconds=1.0,
        component_area_per_node=scene.total_area_per_node,
    )

    for nid in keys(scene.total_area_per_node)
        @test isapprox(
            get(budget_all.ri_par_q_per_node, nid, 0.0),
            get(budget_sun.ri_par_q_per_node, nid, 0.0) + get(budget_sky.ri_par_q_per_node, nid, 0.0);
            atol=1e-10,
            rtol=1e-10,
        )
        @test isapprox(
            get(budget_all.ri_nir_q_per_node, nid, 0.0),
            get(budget_sun.ri_nir_q_per_node, nid, 0.0) + get(budget_sky.ri_nir_q_per_node, nid, 0.0);
            atol=1e-10,
            rtol=1e-10,
        )
    end
end

@testset "Model-driven absorptance unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="absorber", type="plate", item_id=1),
    ])
    meteo_row = _synthetic_meteo_row(; duration_seconds=2.0, ri_par_f=100.0, ri_nir_f=80.0, direct_fraction=1.0)

    mktempdir() do tmp
        model_path = joinpath(tmp, "model_absorber.yml")
        write(
            model_path,
            """
            ---
            Group: absorber
            Type:
              plate:
                  Interception:
                      model: Translucent
                      transparency: 0
                      optical_properties:
                          PAR: 0.2
                          NIR: 0.7
            ...
            """,
        )

        cfg = _synthetic_plate_cfg(
            cfg_ref;
            sectors=1,
            all_in_turtle=false,
            scattering=false,
            pixel_size=0.01,
            models=[model_path],
        )
        step = ArchimedLight.run_light_step(scene, meteo_row, cfg)

        @test isapprox(get(step.budget.ri_par_0_q_per_node, 1, 0.0), 200.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(step.budget.ri_nir_0_q_per_node, 1, 0.0), 160.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(step.budget.ra_par_0_q_per_node, 1, 0.0), 160.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(step.budget.ra_nir_0_q_per_node, 1, 0.0), 48.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(step.budget.ra_par_0_f_per_node, 1, 0.0), 80.0; atol=1e-9, rtol=1e-9)
        @test isapprox(get(step.budget.ra_nir_0_f_per_node, 1, 0.0), 24.0; atol=1e-9, rtol=1e-9)

        rows = ArchimedLight.component_values_table(
            scene,
            step,
            cfg;
            meteo_row=meteo_row,
            step_number=1,
            step_duration_seconds=2.0,
            columns=["item_id", "Ri_PAR_0_q", "Ra_PAR_0_q", "Ri_NIR_0_q", "Ra_NIR_0_q"],
        ).rows
        @test length(rows) == 1
        row = only(rows)
        @test isapprox(Float64(row["Ri_PAR_0_q"]), 200.0; atol=1e-9, rtol=1e-9)
        @test isapprox(Float64(row["Ra_PAR_0_q"]), 160.0; atol=1e-9, rtol=1e-9)
        @test isapprox(Float64(row["Ri_NIR_0_q"]), 160.0; atol=1e-9, rtol=1e-9)
        @test isapprox(Float64(row["Ra_NIR_0_q"]), 48.0; atol=1e-9, rtol=1e-9)
    end
end

@testset "Synthetic meteo series stability unit checks" begin
    cfg_ref = ArchimedLight.read_light_config(joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"))
    scene = _synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", item_id=1),
    ])
    cfg = _synthetic_plate_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    rows = [
        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; file="synthetic-meteo-series"))
    series = ArchimedLight.run_light_series(scene, meteo, cfg)

    @test length(series) == 3
    @test max_abs_float_dict_diff(series[1].budget.ri_par_f_per_node, series[2].budget.ri_par_f_per_node) == 0.0
    @test max_abs_float_dict_diff(series[1].budget.ri_nir_f_per_node, series[2].budget.ri_nir_f_per_node) == 0.0
    @test max_abs_float_dict_diff(series[1].budget.ri_par_f_per_node, series[3].budget.ri_par_f_per_node) == 0.0
    @test max_abs_float_dict_diff(series[1].budget.ri_nir_f_per_node, series[3].budget.ri_nir_f_per_node) == 0.0

    @test isapprox(get(series[1].budget.ri_par_q_per_node, 1, 0.0), 120.0 * 600.0; atol=1e-6, rtol=1e-12)
    @test isapprox(get(series[1].budget.ri_nir_q_per_node, 1, 0.0), 80.0 * 600.0; atol=1e-6, rtol=1e-12)
    @test isapprox(get(series[2].budget.ri_par_q_per_node, 1, 0.0), 3.0 * get(series[1].budget.ri_par_q_per_node, 1, 0.0); atol=1e-6, rtol=1e-12)
    @test isapprox(get(series[2].budget.ri_nir_q_per_node, 1, 0.0), 3.0 * get(series[1].budget.ri_nir_q_per_node, 1, 0.0); atol=1e-6, rtol=1e-12)
    @test isapprox(get(series[3].budget.ri_par_q_per_node, 1, 0.0), get(series[1].budget.ri_par_q_per_node, 1, 0.0); atol=1e-6, rtol=1e-12)
    @test isapprox(get(series[3].budget.ri_nir_q_per_node, 1, 0.0), get(series[1].budget.ri_nir_q_per_node, 1, 0.0); atol=1e-6, rtol=1e-12)

    step1 = ArchimedLight.run_light_step(scene, rows[1], cfg)
    step2 = ArchimedLight.run_light_step(scene, rows[2], cfg)
    @test max_abs_float_dict_diff(step1.budget.ri_par_q_per_node, series[1].budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step1.budget.ri_nir_q_per_node, series[1].budget.ri_nir_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step2.budget.ri_par_q_per_node, series[2].budget.ri_par_q_per_node) == 0.0
    @test max_abs_float_dict_diff(step2.budget.ri_nir_q_per_node, series[2].budget.ri_nir_q_per_node) == 0.0
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

end

if _RUN_SKY_MATRIX_TESTS
@testset "Sky matrix parity" begin
    fixtures = Dict(f.name => f for f in light_parity_fixtures())
    to_int(v) = v isa Number ? Int(round(v)) : parse(Int, string(v))

    function _row_float_value_local(row, col::String)
        sym = Symbol(col)
        sym in propertynames(row) || return nothing
        v = getproperty(row, sym)
        v === missing && return nothing
        v isa Number && return Float64(v)
        try
            return parse(Float64, string(v))
        catch
            return nothing
        end
    end

    function _filter_finite_rows(rows::Vector, key_cols::Vector{String}, value_col::String)
        keys = Set{Tuple}()
        filtered = Any[]
        for row in rows
            v = _row_float_value_local(row, value_col)
            if v !== nothing && isfinite(v)
                push!(filtered, row)
                push!(keys, Tuple(getproperty(row, Symbol(c)) for c in key_cols))
            end
        end
        return filtered, keys
    end

    function _assert_rows_match(
        expected_rows::Vector,
        observed_rows::Vector,
        key_cols::Vector{String},
        value_cols::Vector{String};
        atol::Float64=1e-6,
        rtol::Float64=1e-3,
        label::AbstractString="",
    )
        cmp = compare_rows_by_key(
            expected_rows,
            observed_rows;
            key_cols=key_cols,
            value_cols=value_cols,
            atol=atol,
            rtol=rtol,
            top_n=10,
        )
        if !cmp.ok
            @info "Sky matrix row mismatch" label missing_keys=cmp.missing_keys extra_keys=cmp.extra_keys mismatch_count=cmp.mismatches_total mismatches=cmp.mismatches
        end
        @test cmp.ok
    end

    function _format_meteo_date(v)
        v === missing && return ""
        if v isa Dates.Date
            return Dates.format(v, Dates.DateFormat("yyyy/mm/dd"))
        elseif v isa Dates.DateTime
            return Dates.format(Dates.Date(v), Dates.DateFormat("yyyy/mm/dd"))
        end
        return replace(string(v), '-' => '/')
    end

    function _observed_meteo_rows(series, meteo_rows)
        rows = Dict{String,Any}[]
        for i in eachindex(series)
            row = meteo_rows[i]
            step = series[i]
            push!(
                rows,
                Dict(
                    "step_number" => i - 1,
                    "date" => _format_meteo_date(getproperty(row, :date)),
                    "hour_start" => string(getproperty(row, :hour_start)),
                    "hour_end" => string(getproperty(row, :hour_end)),
                    "sun_elevation" => step.sky.sun_elevation_deg,
                    "sun_azimut" => step.sky.sun_azimuth_deg,
                ),
            )
        end
        rows
    end

    sky_matrix_fixtures = (
        "test-compare-sky06",
        "test-compare-sky16",
        "test-compare-sky46",
    )

    for fx_name in sky_matrix_fixtures
        fx = fixtures[fx_name]
        component_energy_atol = fx_name == "test-compare-sky46" ? 50.0 : 1e-3
        component_energy_rtol = fx_name == "test-compare-sky46" ? 3e-2 : 1e-3
        scat_atol = fx_name == "test-compare-sky46" ? 1.5e-5 : 1e-6
        scat_rtol = fx_name == "test-compare-sky46" ? 3e-2 : 1e-4
        expected_comp_path = _expected_component_values_path(fx)
        expected_scene_path = _expected_scene_values_path(fx)
        expected_sun_path = _expected_sun_log_path(fx)
        expected_summary_path = _expected_summary_path(fx)
        expected_scat_path = _expected_scattering_log_path(fx)
        expected_meteo_path = fx.expected_dir === nothing ? nothing : joinpath(fx.expected_dir, "meteo.csv")

        @test expected_comp_path !== nothing
        @test expected_scene_path !== nothing
        @test expected_sun_path !== nothing
        @test expected_summary_path !== nothing
        @test expected_scat_path !== nothing
        @test expected_meteo_path !== nothing
        @test isfile(expected_meteo_path)

        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        series = ArchimedLight.run_light_series(scene, selected, cfg)
        @test length(series) == length(selected.rows)

        mktempdir() do tmp
            out_paths = ArchimedLight.write_light_outputs(
                scene,
                series,
                cfg;
                meteo_rows=selected.rows,
                start_step_number=0,
                outdir=joinpath(tmp, "output"),
                write_component=true,
                write_scene=true,
                write_summary=true,
                write_sun_position_log=true,
                write_scattering_log=cfg.scattering,
                scattering_log_bands=["PAR"],
            )

            expected_comp = read_java_csv(expected_comp_path)
            observed_comp = read_java_csv(out_paths["component_values"])
            sort!(expected_comp; by=r -> (to_int(getproperty(r, :step_number)), to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id))))
            sort!(observed_comp; by=r -> (to_int(getproperty(r, :step_number)), to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id))))
            _assert_rows_match(expected_comp, observed_comp, ["step_number", "item_id", "component_id"], ["area", "barycentre_z"]; atol=1e-6, rtol=1e-4, label="component geometry $(fx_name)")
            expected_sky_fraction, finite_keys = _filter_finite_rows(expected_comp, ["step_number", "item_id", "component_id"], "sky_fraction")
            observed_sky_fraction = [r for r in observed_comp if (getproperty(r, :step_number), getproperty(r, :item_id), getproperty(r, :component_id)) in finite_keys]
            _assert_rows_match(expected_sky_fraction, observed_sky_fraction, ["step_number", "item_id", "component_id"], ["sky_fraction"]; atol=1e-4, rtol=1e-3, label="component sky fraction $(fx_name)")
            _assert_rows_match(
                expected_comp,
                observed_comp,
                ["step_number", "item_id", "component_id"],
                ["Ri_PAR_0_q", "Ri_PAR_q", "Ra_PAR_0_q", "Ra_PAR_q", "Ri_NIR_q", "Ra_NIR_q"];
                atol=component_energy_atol,
                rtol=component_energy_rtol,
                label="component energy $(fx_name)",
            )

            expected_scene = read_java_csv(expected_scene_path)
            observed_scene = read_java_csv(out_paths["scene_values"])
            sort!(expected_scene; by=r -> to_int(getproperty(r, :step_number)))
            sort!(observed_scene; by=r -> to_int(getproperty(r, :step_number)))
            scene_exact_cols = [c for c in ("date", "hour_start", "hour_end") if c in String.(propertynames(first(expected_scene))) && c in String.(propertynames(first(observed_scene)))]
            !isempty(scene_exact_cols) && _assert_rows_match(expected_scene, observed_scene, ["step_number"], scene_exact_cols; atol=0.0, rtol=0.0, label="scene exact $(fx_name)")
            _assert_rows_match(expected_scene, observed_scene, ["step_number"], ["RI_SW_f"]; atol=1e-4, rtol=1e-6, label="scene irradiance $(fx_name)")
            _assert_rows_match(expected_scene, observed_scene, ["step_number"], ["plot_area"]; atol=1e-5, rtol=0.0, label="scene plot area $(fx_name)")

            expected_summary = read_java_csv(expected_summary_path)
            observed_summary = read_java_csv(out_paths["summary"])
            sort!(expected_summary; by=r -> (to_int(getproperty(r, :step_number)), string(getproperty(r, :group)), string(getproperty(r, :type)), to_int(getproperty(r, :item_id))))
            sort!(observed_summary; by=r -> (to_int(getproperty(r, :step_number)), string(getproperty(r, :group)), string(getproperty(r, :type)), to_int(getproperty(r, :item_id))))
            summary_cols = [c for c in ("Ri_q", "Ra_q") if c in String.(propertynames(first(expected_summary))) && c in String.(propertynames(first(observed_summary)))]
            @test !isempty(summary_cols)
            _assert_rows_match(expected_summary, observed_summary, ["step_number", "group", "type", "item_id"], summary_cols; atol=1e-3, rtol=5e-4, label="summary $(fx_name)")

            expected_sun = read_java_csv(expected_sun_path)
            observed_sun = read_java_csv(out_paths["log_sun_position"])
            sort!(expected_sun; by=r -> to_int(getproperty(r, :stepNumber)))
            sort!(observed_sun; by=r -> to_int(getproperty(r, :stepNumber)))
            _assert_rows_match(expected_sun, observed_sun, ["stepNumber"], ["azimuthWeighted", "elevationWeighted"]; atol=1e-4, rtol=0.0, label="sun log $(fx_name)")

            expected_meteo = read_java_csv(expected_meteo_path)
            observed_meteo = _observed_meteo_rows(series, selected.rows)
            sort!(expected_meteo; by=r -> to_int(getproperty(r, :step_number)))
            sort!(observed_meteo; by=r -> Int(r["step_number"]))
            _assert_rows_match(expected_meteo, observed_meteo, ["step_number"], ["date", "hour_start", "hour_end"]; atol=0.0, rtol=0.0, label="meteo exact $(fx_name)")
            _assert_rows_match(expected_meteo, observed_meteo, ["step_number"], ["sun_elevation", "sun_azimut"]; atol=1e-4, rtol=0.0, label="meteo sun $(fx_name)")

            expected_scat = read_java_csv(expected_scat_path)
            observed_scat = read_java_csv(out_paths["log_iteration_scat_par"])
            sort!(expected_scat; by=r -> (to_int(getproperty(r, :step)), to_int(getproperty(r, :plantid)), to_int(getproperty(r, :nodeid)), to_int(getproperty(r, :iter))))
            sort!(observed_scat; by=r -> (to_int(getproperty(r, :step)), to_int(getproperty(r, :plantid)), to_int(getproperty(r, :nodeid)), to_int(getproperty(r, :iter))))
            scat_value_cols = [c for c in ("scat", "scat_W") if c in String.(propertynames(first(expected_scat))) && c in String.(propertynames(first(observed_scat)))]
            @test !isempty(scat_value_cols)
            scat_keys = ["step", "plantid", "nodeid", "iter"]
            observed_scat_keys = Set(Tuple(getproperty(r, Symbol(c)) for c in scat_keys) for r in observed_scat)
            expected_scat_keys = Set(Tuple(getproperty(r, Symbol(c)) for c in scat_keys) for r in expected_scat)
            @test isempty(setdiff(observed_scat_keys, expected_scat_keys))
            expected_scat_common = [r for r in expected_scat if Tuple(getproperty(r, Symbol(c)) for c in scat_keys) in observed_scat_keys]
            _assert_rows_match(
                expected_scat_common,
                observed_scat,
                scat_keys,
                scat_value_cols;
                atol=scat_atol,
                rtol=scat_rtol,
                label="scattering log $(fx_name)",
            )
        end
    end
end
end

if _RUN_PARITY_TESTS || _RUN_GUARD_TESTS
@testset "Parity reports" begin
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
    function _assert_rel_metric(actual, expected, metric::Symbol, label::AbstractString)
        err, lim = parity_rel_with_limit(actual, expected, metric)
        if !(err <= lim)
            @info "Parity relative mismatch" label metric actual expected err limit=lim
        end
        @test err <= lim
    end
    function _assert_abs_metric(actual, expected, metric::Symbol, label::AbstractString)
        err, lim = parity_abs_with_limit(actual, expected, metric)
        if !(err <= lim)
            @info "Parity absolute mismatch" label metric actual expected err limit=lim
        end
        @test err <= lim
    end
    function _assert_component_metric_limit(
        rows::Vector,
        field::Symbol,
        metric::Symbol,
        label::AbstractString;
        top_n::Int=8,
    )
        values = Float64[getproperty(r, field) for r in rows]
        @test !isempty(values)
        lim = parity_limit(metric)
        vmax = maximum(values)
        if !(vmax <= lim)
            top = sort(rows; by=r -> getproperty(r, field), rev=true)[1:min(top_n, length(rows))]
            @info "Parity component metric mismatch" label metric max_value=vmax limit=lim top
        end
        @test vmax <= lim
    end
    _PARITY_DEBUG_ROWS = lowercase(strip(get(ENV, "ARCHIMEDLIGHT_PARITY_DEBUG_ROWS", "0"))) in ("1", "true", "yes", "on")
    _PARITY_DEBUG_LABEL_FILTER = strip(get(ENV, "ARCHIMEDLIGHT_PARITY_DEBUG_LABEL_FILTER", ""))
    _PARITY_DEBUG_MAX_ROWS = let raw = strip(get(ENV, "ARCHIMEDLIGHT_PARITY_DEBUG_MAX_ROWS", "8"))
        parsed = try
            parse(Int, raw)
        catch
            8
        end
        max(parsed, 1)
    end

    function _parity_debug_enabled(label::AbstractString)
        _PARITY_DEBUG_ROWS || return false
        isempty(_PARITY_DEBUG_LABEL_FILTER) && return true
        occursin(_PARITY_DEBUG_LABEL_FILTER, label)
    end

    _rows_by_key(rows::Vector, key_cols::Vector{String}) = Dict(_row_key_from_dict(r, key_cols) => r for r in rows)
    _row_subset(row, cols::Vector{String}) = Dict(c => _row_get_by_col(row, c) for c in cols)

    function _log_row_mismatch_details(
        expected_rows::Vector,
        observed_rows::Vector,
        key_cols::Vector{String},
        cmp,
        label::AbstractString,
    )
        exp_map = _rows_by_key(expected_rows, key_cols)
        obs_map = _rows_by_key(observed_rows, key_cols)

        miss_n = min(length(cmp.missing_keys), _PARITY_DEBUG_MAX_ROWS)
        for i in 1:miss_n
            key = cmp.missing_keys[i]
            erow = get(exp_map, key, nothing)
            @info "Parity missing key detail" label rank=i key expected_row=(erow === nothing ? nothing : _row_subset(erow, key_cols))
        end

        extra_n = min(length(cmp.extra_keys), _PARITY_DEBUG_MAX_ROWS)
        for i in 1:extra_n
            key = cmp.extra_keys[i]
            orow = get(obs_map, key, nothing)
            @info "Parity extra key detail" label rank=i key observed_row=(orow === nothing ? nothing : _row_subset(orow, key_cols))
        end

        mm_n = min(length(cmp.mismatches), _PARITY_DEBUG_MAX_ROWS)
        for i in 1:mm_n
            mm = cmp.mismatches[i]
            key = mm.key
            col = string(mm.column)
            erow = get(exp_map, key, nothing)
            orow = get(obs_map, key, nothing)
            row_cols = unique(vcat(key_cols, [col]))
            @info "Parity mismatch detail" label rank=i key column=col expected=mm.expected observed=mm.observed abs_err=mm.abs_err rel_err=mm.rel_err expected_row=(erow === nothing ? nothing : _row_subset(erow, row_cols)) observed_row=(orow === nothing ? nothing : _row_subset(orow, row_cols))
        end
    end

    function _assert_rows_match(
        expected_rows::Vector,
        observed_rows::Vector,
        key_cols::Vector{String},
        value_cols::Vector{String};
        atol::Float64=1e-6,
        rtol::Float64=1e-3,
        top_n::Int=10,
        label::AbstractString="",
    )
        cmp = compare_rows_by_key(
            expected_rows,
            observed_rows;
            key_cols=key_cols,
            value_cols=value_cols,
            atol=atol,
            rtol=rtol,
            top_n=top_n,
        )
        if !cmp.ok
            @info "Parity row mismatch" label missing_keys=cmp.missing_keys extra_keys=cmp.extra_keys mismatch_count=cmp.mismatches_total mismatches=cmp.mismatches
            _parity_debug_enabled(label) && _log_row_mismatch_details(expected_rows, observed_rows, key_cols, cmp, label)
        end
        @test cmp.ok
    end

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
    function _component_snapshot_rows(scene, step, cfg, row)
        cols = ["step_number", "item_id", "component_id", "area", "Ri_PAR_0_f", "Ri_NIR_0_f", "barycentre_z"]
        raw_rows = ArchimedLight.component_values_table(
            scene,
            step,
            cfg;
            meteo_row=row,
            columns=cols,
        ).rows
        out = Vector{Dict{String,Any}}(undef, length(raw_rows))
        for i in eachindex(raw_rows)
            r = raw_rows[i]
            out[i] = Dict(
                "step_number" => Int(r["step_number"]),
                "item_id" => Int(r["item_id"]),
                "component_id" => Int(r["component_id"]),
                "area" => Float64(r["area"]),
                "irradiance_withoutScattering_PAR_NIR" => Float64(r["Ri_PAR_0_f"]) + Float64(r["Ri_NIR_0_f"]),
                "barycentre_z" => Float64(r["barycentre_z"]),
            )
        end
        sort!(out; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
        out
    end
    function _component_series_rows(scene, series, cfg, meteo_rows)
        cols = ["step_number", "item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"]
        out = Dict{Tuple{Int,Int,Int},Dict{String,Any}}()
        for i in eachindex(series)
            step = series[i]
            row = meteo_rows[i]
            rows = ArchimedLight.component_values_table(
                scene,
                step,
                cfg;
                meteo_row=row,
                step_number=i - 1,
                columns=cols,
            ).rows
            for r in rows
                key = (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"]))
                out[key] = Dict(
                    "step_number" => key[1],
                    "item_id" => key[2],
                    "component_id" => key[3],
                    "Ri_PAR_0_f" => Float64(r["Ri_PAR_0_f"]),
                    "Ri_NIR_0_f" => Float64(r["Ri_NIR_0_f"]),
                    "Ri_PAR_f" => Float64(r["Ri_PAR_f"]),
                    "Ri_NIR_f" => Float64(r["Ri_NIR_f"]),
                )
            end
        end
        rows = collect(values(out))
        sort!(rows; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
        rows
    end
    function _scene_snapshot_rows(scene, series, cfg, meteo_rows)
        cols = ["step_number", "date", "hour_start", "hour_end", "RI_SW_f"]
        rows = ArchimedLight.scene_values_table(
            scene,
            series,
            cfg;
            meteo_rows=meteo_rows,
            columns=cols,
        ).rows
        sort!(rows; by=r -> Int(r["step_number"]))
        rows
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
    test_root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests")
    function _stage_start(name::AbstractString)
        t0 = time()
        @info "Parity stage start" stage=name
        t0
    end
    function _stage_done(name::AbstractString, t0::Real)
        @info "Parity stage done" stage=name elapsed_s=round(time() - t0; digits=3)
        nothing
    end

    @testset "Hitcount parity" begin
        _stage_t0 = _stage_start("Hitcount parity")
    f_hit = fixtures["test-hitcount2"]
    r_false = fixture_parity_report(f_hit; all_in_turtle=false)
    r_true = fixture_parity_report(f_hit; all_in_turtle=true)
    @test r_false.expected_hitcount_total !== nothing
    @test r_true.expected_hitcount_total !== nothing
    @test r_false.snapshot.total_hits > 0
    @test r_true.snapshot.total_hits > 0
    _assert_rel_metric(r_false.snapshot.total_hits, r_false.expected_hitcount_total, :hitcount_total_rel, "test-hitcount2 total hits raycast")
    _assert_rel_metric(r_true.snapshot.total_hits, r_true.expected_hitcount_total, :hitcount_total_rel, "test-hitcount2 total hits turtle")

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
            _assert_rel_metric(got.hits, exp.hits, :hitcount_component_rel, "test-hitcount2 per-component hits key=$(key)")
        end

        # Item component IDs are now mapped explicitly; plant component irradiance is in tight parity.
        plant_key = (1, 2)
        @test haskey(expected, plant_key)
        @test haskey(observed, plant_key)
        expp = expected[plant_key]
        gotp = observed[plant_key]
        got_irr = gotp.power0 / max(expp.area, eps(Float64))
        _assert_rel_metric(got_irr, expp.irr0, :irr_component_rel_strict, "test-hitcount2 plant irradiance")
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

        irr_dense_hits_threshold = 400
        hit_rel_hi = NamedTuple{(:key,:value),Tuple{Tuple{Int,Int},Float64}}[]
        irr_rel_dense = NamedTuple{(:key,:value,:hits),Tuple{Tuple{Int,Int},Float64,Int}}[]
        irr_abs_sparse = NamedTuple{(:key,:value,:hits),Tuple{Tuple{Int,Int},Float64,Int}}[]
        for key in keys(expected)
            exp = expected[key]
            got = observed[key]
            rel_hit = abs(got.hits - exp.hits) / max(exp.hits, 1)
            if exp.hits >= 1000
                push!(hit_rel_hi, (key=key, value=rel_hit))
            end

            got_irr = got.power0 / max(exp.area, eps(Float64))
            abs_err_irr = abs(got_irr - exp.irr0)
            if abs(exp.irr0) >= 1.0
                if exp.hits >= irr_dense_hits_threshold
                    push!(irr_rel_dense, (key=key, value=abs_err_irr / max(abs(exp.irr0), eps(Float64)), hits=exp.hits))
                else
                    push!(irr_abs_sparse, (key=key, value=abs_err_irr, hits=exp.hits))
                end
            end
        end

        if all_in_turtle
            _assert_component_metric_limit(hit_rel_hi, :value, :hitcount_component_hi_rel_turtle, "test-hitcount3 high-hit rel (turtle)")
            _assert_component_metric_limit(irr_rel_dense, :value, :irr_component_rel_dense, "test-hitcount3 dense irradiance rel (turtle)")
            _assert_component_metric_limit(irr_abs_sparse, :value, :irr_component_abs_sparse, "test-hitcount3 sparse irradiance abs (turtle)")
        else
            _assert_component_metric_limit(hit_rel_hi, :value, :hitcount_component_hi_rel_raycast, "test-hitcount3 high-hit rel (raycast)")
            _assert_component_metric_limit(irr_rel_dense, :value, :irr_component_rel_dense, "test-hitcount3 dense irradiance rel (raycast)")
            _assert_component_metric_limit(irr_abs_sparse, :value, :irr_component_abs_sparse, "test-hitcount3 sparse irradiance abs (raycast)")
        end
    end

    for (all_in_turtle, abs_metric, rel_metric) in (
        (false, :hitcount_hist_abs_raycast, :hitcount_hist_rel_raycast),
        (true, :hitcount_hist_abs_turtle, :hitcount_hist_rel_turtle),
    )
        step = run_fixture_once(f_hit; all_in_turtle=all_in_turtle)
        expected_path = _expected_component_values_path(f_hit; all_in_turtle=all_in_turtle)
        @test expected_path !== nothing
        expected_hits = expected_component_hits(expected_path)
        got_hits = sort(collect(values(step.first_order.hits_per_node)))
        @test length(got_hits) == length(expected_hits)
        abs_err = abs.(got_hits .- expected_hits)
        rel_errs = abs_err ./ max.(abs.(Float64.(expected_hits)), 1.0)
        @test maximum(abs_err) <= parity_limit(abs_metric)
        @test maximum(rel_errs) <= parity_limit(rel_metric)
    end
        _stage_done("Hitcount parity", _stage_t0)
    end

    @testset "Sun and snapshot parity" begin
        _stage_t0 = _stage_start("Sun and snapshot parity")
    f_wsun = fixtures["test-weighted-sun"]
    r_wsun = fixture_parity_report(f_wsun)
    @test r_wsun.expected_sun_azimuth !== nothing
    @test r_wsun.expected_sun_elevation !== nothing
    @test isfinite(r_wsun.snapshot.sun_azimuth)
    @test isfinite(r_wsun.snapshot.sun_elevation)
    _assert_abs_metric(r_wsun.snapshot.sun_azimuth, r_wsun.expected_sun_azimuth, :sun_deg_snapshot, "test-weighted-sun snapshot azimuth")
    _assert_abs_metric(r_wsun.snapshot.sun_elevation, r_wsun.expected_sun_elevation, :sun_deg_snapshot, "test-weighted-sun snapshot elevation")

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
    @test max_az_err < parity_limit(:sun_deg_az_series)
    @test max_el_err < parity_limit(:sun_deg_el_series)

    if _RUN_PARITY_TESTS
    for fx_name in ("test-weighted-sun", "test-save_on_disk1", "test-save_on_disk6")
        max_abs_err = parity_limit(:scene_sw_series_max_abs; fixture=fx_name)
        max_mean_rel_err = parity_limit(:scene_sw_series_mean_rel; fixture=fx_name)
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

    # Frozen Java scene_values snapshots (full keyed rows).
    for fx_name in ("test-weighted-sun", "test-save_on_disk1", "test-save_on_disk6")
        sw_atol = parity_limit(:scene_snapshot_sw_atol; fixture=fx_name)
        sw_rtol = parity_limit(:scene_snapshot_sw_rtol; fixture=fx_name)
        fx = fixtures[fx_name]
        expected_path = _expected_scene_values_path(fx)
        @test expected_path !== nothing
        expected_rows = read_java_csv(expected_path)
        sort!(expected_rows; by=r -> Int(getproperty(r, :step_number)))

        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        observed_rows = _scene_snapshot_rows(scene, series, cfg, meteo.rows)

        _assert_rows_match(
            expected_rows,
            observed_rows,
            ["step_number"],
            ["date", "hour_start", "hour_end"];
            atol=0.0,
            rtol=0.0,
            label="scene snapshot fields fixture=$(fx_name)",
        )
        _assert_rows_match(
            expected_rows,
            observed_rows,
            ["step_number"],
            ["RI_SW_f"];
            atol=sw_atol,
            rtol=sw_rtol,
            label="scene snapshot RI_SW_f fixture=$(fx_name)",
        )
    end

    # Frozen Java component_values snapshots for single-step save fixtures.
    for fx_name in ("test-save_on_disk1", "test-save_on_disk6")
        fx = fixtures[fx_name]
        expected_path = _expected_component_values_path(fx)
        @test expected_path !== nothing
        expected_rows = read_java_csv(expected_path)
        sort!(expected_rows; by=r -> (Int(getproperty(r, :step_number)), Int(getproperty(r, :item_id)), Int(getproperty(r, :component_id))))

        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        row = first(meteo.rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        observed_rows = _component_snapshot_rows(scene, step, cfg, row)

        _assert_rows_match(
            expected_rows,
            observed_rows,
            ["step_number", "item_id", "component_id"],
            ["area", "barycentre_z"];
            atol=1e-6,
            rtol=1e-4,
            label="component snapshot geometry fixture=$(fx_name)",
        )
    end

    for fx_name in ("test-cafeier", "test-cafeier2", "test-cafeier_sensor", "test-cafeier_sensor2", "test-cafeier_sensor3")
        max_rel_err = parity_limit(:cafeier_total_rel; fixture=fx_name)
        fx = fixtures[fx_name]
        expected_path = _expected_component_values_path(fx)
        @test expected_path !== nothing
        expected_rows = read_java_csv(expected_path)
        exp_total = 0.0
        exp_keys = Set{Tuple{Int,Int}}()
        for r in expected_rows
            names = propertynames(r)
            (:area in names && :irradiance_withoutScattering_PAR_NIR in names && :item_id in names && :component_id in names) || continue
            push!(exp_keys, (to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id))))
            exp_total += Float64(getproperty(r, :area)) * Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR))
        end
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        key_by_node = ArchimedLight._interception_java_keys(scene, cfg)
        got_total = 0.0
        for (nid, key) in key_by_node
            key in exp_keys || continue
            got_total += get(step.first_order.incident_par_power_per_node, nid, 0.0)
            got_total += get(step.first_order.incident_nir_power_per_node, nid, 0.0)
        end
        @test relerr(got_total, exp_total) < max_rel_err

        summary_path = _expected_summary_path(fx)
        if summary_path !== nothing
            expected_summary = read_java_csv(summary_path)
            got_summary = ArchimedLight.summary_values_table(scene, [step], cfg; meteo_rows=[row], start_step_number=0)
            exp_by_item = Dict{Int,Float64}()
            for r in expected_summary
                item = to_int(getproperty(r, :item_id))
                exp_by_item[item] = get(exp_by_item, item, 0.0) + Float64(getproperty(r, :Ri_q))
            end
            got_by_item = Dict{Int,Float64}()
            for r in got_summary.rows
                item = Int(r["item_id"])
                got_by_item[item] = get(got_by_item, item, 0.0) + Float64(r["Ri_q"])
            end
            exp_riq = sum(values(exp_by_item))
            got_riq = sum(get(got_by_item, k, 0.0) for k in keys(exp_by_item))
            _assert_rel_metric(got_riq, exp_riq, :scene_riq_rel, "summary Ri_q by item fixture=$(fx_name)")
        end
    end
    end

        _stage_done("Sun and snapshot parity", _stage_t0)
    end

    @testset "Scattering links outputs and lightsource parity" begin
        _stage_t0 = _stage_start("Scattering links outputs and lightsource parity")
    f_scat = fixtures["test-scattering-one-plate"]
    r_scat = fixture_parity_report(f_scat)
    @test r_scat.expected_scattering_total !== nothing
    _assert_rel_metric(r_scat.julia_scattering_total, r_scat.expected_scattering_total, :scattering_total_rel_loose, "test-scattering-one-plate total scattering")

    if _RUN_PARITY_TESTS
    f_scat2 = fixtures["test-scattering-two-plates"]
    r_scat2 = fixture_parity_report(f_scat2)
    @test r_scat2.expected_scattering_total !== nothing
    _assert_rel_metric(r_scat2.julia_scattering_total, r_scat2.expected_scattering_total, :scattering_total_rel_strict, "test-scattering-two-plates total scattering")

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
        @test max_abs_float_dict_diff(bud_pix.ri_par_q_per_node, bud_lnk.ri_par_q_per_node) < parity_limit(:backend_identity_abs)
        @test max_abs_float_dict_diff(bud_pix.ri_nir_q_per_node, bud_lnk.ri_nir_q_per_node) < parity_limit(:backend_identity_abs)
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
        # Java uses Float32 turtle/pixel internals; with Java-aligned directions this
        # ratio stays effectively null but not bit-zero in Julia Float64 accumulation.
        @test got_std / max(abs(got_mean), eps(Float64)) < parity_limit(:area_ratio_uniformity_julia_relstd)
        @test exp_std / max(abs(exp_mean), eps(Float64)) < parity_limit(:area_ratio_uniformity_java_relstd)
        # Java sky conversion and units are matched closely.
        @test relerr(got_mean, exp_mean) < parity_limit(:area_ratio_mean_rel)

        area_ratio_metrics[area_name] = (got_mean, got_std, exp_mean, exp_std)
    end

    # Java runs this fixture family under four cache/all_in_turtle combinations.
    # Means should stay stable across those runtime switches.
    got_means = [v[1] for v in values(area_ratio_metrics)]
    mean_ref = max(abs(mean(got_means)), 1.0)
    @test (maximum(got_means) - minimum(got_means)) / mean_ref < parity_limit(:area_ratio_means_spread_rel)

    # Java test-area_ratio5 parity: area-ratio correction should keep leaf irradiance
    # stdev relatively stable across pixel sizes; without correction, stdev decreases
    # strongly with increased pixel resolution.
    function with_pixel_size_area_ratio(cfg::ArchimedLight.LightConfig, pixel_size_cm::Float64, area_ratio::Bool)
        raw = copy(cfg.raw)
        ArchimedLight.LightConfig(
            cfg.scene,
            cfg.meteo,
            cfg.all_in_turtle,
            cfg.turtle_sectors,
            pixel_size_cm / 100.0,
            area_ratio,
            cfg.scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            cfg.cache_radiation,
            raw,
        )
    end

    f_area5 = fixtures["test-area_ratio5"]
    cfg_area5 = ArchimedLight.read_light_config(f_area5.config_path)
    scene_area5 = ArchimedLight.read_scene(cfg_area5.scene)
    row_area5 = first(ArchimedLight.read_meteo(cfg_area5.meteo).rows)

    pixel_sweep = [
        (50000, 5.542562694341231),
        (100000, 3.919183666320266),
        (200000, 2.7712813471706156),
        (500000, 1.7527122188397935),
        (1000000, 1.2393546954101381),
    ]

    stdev_by_npix_ratio = Dict{Int,Float64}()
    stdev_by_npix_noratio = Dict{Int,Float64}()

    for (npix, px_cm) in pixel_sweep
        cfg_ratio = with_pixel_size_area_ratio(cfg_area5, px_cm, true)
        step_ratio = ArchimedLight.run_light_step(scene_area5, row_area5, cfg_ratio)
        rows_ratio = ArchimedLight.component_values_table(
            scene_area5,
            step_ratio,
            cfg_ratio;
            meteo_row=row_area5,
            columns=["type", "Ri_PAR_0_f", "Ri_NIR_0_f"],
        ).rows
        irr_ratio = Float64[
            Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_ratio if string(get(r, "type", "")) == "Leaf"
        ]
        @test !isempty(irr_ratio)
        stdev_by_npix_ratio[npix] = std(irr_ratio)

        cfg_noratio = with_pixel_size_area_ratio(cfg_area5, px_cm, false)
        step_noratio = ArchimedLight.run_light_step(scene_area5, row_area5, cfg_noratio)
        rows_noratio = ArchimedLight.component_values_table(
            scene_area5,
            step_noratio,
            cfg_noratio;
            meteo_row=row_area5,
            columns=["type", "Ri_PAR_0_f", "Ri_NIR_0_f"],
        ).rows
        irr_noratio = Float64[
            Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_noratio if string(get(r, "type", "")) == "Leaf"
        ]
        @test !isempty(irr_noratio)
        stdev_by_npix_noratio[npix] = std(irr_noratio)
    end

    expected_pix = read_java_csv(joinpath(fixture_path(f_area5), "expected", "pix.csv"))
    expected_stdev_ratio = Dict{Int,Float64}()
    expected_stdev_noratio = Dict{Int,Float64}()
    for r in expected_pix
        n = to_int(getproperty(r, :npix))
        expected_stdev_ratio[n] = Float64(getproperty(r, :irradiance_stdev_ratio))
        expected_stdev_noratio[n] = Float64(getproperty(r, :irradiance_stdev_noratio))
    end
    @test Set(keys(stdev_by_npix_ratio)) == Set(keys(expected_stdev_ratio))
    @test Set(keys(stdev_by_npix_noratio)) == Set(keys(expected_stdev_noratio))

    # Corrected branch is close to Java across the whole sweep.
    for n in keys(expected_stdev_ratio)
        @test abs(stdev_by_npix_ratio[n] - expected_stdev_ratio[n]) < parity_limit(:area_ratio5_ratio_abs)
    end
    # Uncorrected branch converges close to Java on high-resolution runs.
    for n in (500000, 1000000)
        @test abs(stdev_by_npix_noratio[n] - expected_stdev_noratio[n]) < parity_limit(:area_ratio5_noratio_abs_hi)
    end
    # Behavioral gate from Java README/run-test: correction stabilizes stdev,
    # while no-correction decreases stdev as resolution increases.
    ratio_vals = [stdev_by_npix_ratio[n] for (n, _) in pixel_sweep]
    noratio_vals = [stdev_by_npix_noratio[n] for (n, _) in pixel_sweep]
    @test maximum(ratio_vals) - minimum(ratio_vals) < parity_limit(:area_ratio5_ratio_spread)
    @test first(noratio_vals) > parity_limit(:area_ratio5_noratio_first_last_ratio) * last(noratio_vals)

    # Java test-ignore parity: Interception model=ignore drops ignored types from outputs,
    # and explicit/aliased ignore configurations are equivalent.
    ignore_root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-ignore")
    cfg_ignore0 = ArchimedLight.read_light_config(joinpath(ignore_root, "config.yml"))
    cfg_ignore1 = ArchimedLight.read_light_config(joinpath(ignore_root, "config_ignore1.yml"))
    cfg_ignore2 = ArchimedLight.read_light_config(joinpath(ignore_root, "config_ignore2.yml"))

    scene_ignore = ArchimedLight.read_scene(cfg_ignore0.scene)
    _, _, _, node_ids_ignore0, _, _ = ArchimedLight._scene_geometry_for_interception(scene_ignore, cfg_ignore0)
    _, _, _, node_ids_ignore1, _, _ = ArchimedLight._scene_geometry_for_interception(scene_ignore, cfg_ignore1)
    _, _, _, node_ids_ignore2, _, _ = ArchimedLight._scene_geometry_for_interception(scene_ignore, cfg_ignore2)

    @test any(get(scene_ignore.node_type, nid, "") == "Metamer" for nid in node_ids_ignore0)
    @test all(get(scene_ignore.node_type, nid, "") != "Metamer" for nid in node_ids_ignore1)
    @test length(node_ids_ignore1) < length(node_ids_ignore0)
    @test sort(node_ids_ignore1) == sort(node_ids_ignore2)

    keys_ignore1 = ArchimedLight._interception_java_keys(scene_ignore, cfg_ignore1)
    keys_ignore2 = ArchimedLight._interception_java_keys(scene_ignore, cfg_ignore2)
    @test keys_ignore1 == keys_ignore2

    # Java test-pavement parity: plot_paving creates synthetic pavement components.
    pavement_root = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-pavement")
    cfg_pav_yes = ArchimedLight.read_light_config(joinpath(pavement_root, "config.yml"))
    cfg_pav_no = ArchimedLight.read_light_config(joinpath(pavement_root, "config-nopavement.yml"))
    scene_pav = ArchimedLight.read_scene(cfg_pav_yes.scene)
    _, _, _, node_ids_yes, _, node_group_yes = ArchimedLight._scene_geometry_for_interception(scene_pav, cfg_pav_yes)
    _, _, _, node_ids_no, _, node_group_no = ArchimedLight._scene_geometry_for_interception(scene_pav, cfg_pav_no)
    pav_yes = Int[nid for nid in node_ids_yes if get(node_group_yes, nid, "") == "pavement"]
    pav_no = Int[nid for nid in node_ids_no if get(node_group_no, nid, "") == "pavement"]
    @test !isempty(pav_yes)
    @test isempty(pav_no)

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

    function _node_links_stats_tuples(rows)
        sort([(Int(getproperty(r, :plantid)), Int(getproperty(r, :nodeid)), Int(getproperty(r, :links)), Int(getproperty(r, :hits))) for r in rows])
    end
    function _node_links_dir_tuples(rows)
        sort([(Int(getproperty(r, :dir)), Int(getproperty(r, :plantid1)), Int(getproperty(r, :id1)), Int(getproperty(r, :plantid2)), Int(getproperty(r, :id2)), Int(getproperty(r, :n))) for r in rows])
    end
    function _node_links_stats_tuples_tbl(rows)
        sort([(Int(r["plantid"]), Int(r["nodeid"]), Int(r["links"]), Int(r["hits"])) for r in rows])
    end
    function _node_links_dir_tuples_tbl(rows)
        sort([(Int(r["dir"]), Int(r["plantid1"]), Int(r["id1"]), Int(r["plantid2"]), Int(r["id2"]), Int(r["n"])) for r in rows])
    end

    for links_name in ("test-links", "test-links2", "test-links3", "test-links4")
        fx = fixtures[links_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        sky = ArchimedLight.compute_sky(row, cfg)
        turtle = ArchimedLight.build_turtle(cfg, sky)

        got_stats = ArchimedLight.node_links_stats_alldirs_table(scene, turtle, cfg)
        exp_stats_path = joinpath(fixture_path(fx), "expected", "log-nodelinks-stats-alldirs.csv")
        @test isfile(exp_stats_path)
        exp_stats = read_java_csv(exp_stats_path)
        @test _node_links_stats_tuples_tbl(got_stats.rows) == _node_links_stats_tuples(exp_stats)

        got_dir0 = ArchimedLight.node_links_dir_table(scene, turtle, cfg; direction_index=0)
        exp_dir_path = joinpath(fixture_path(fx), "expected", "log-nodelinks-dir00.csv")
        @test isfile(exp_dir_path)
        exp_dir = read_java_csv(exp_dir_path)
        @test _node_links_dir_tuples_tbl(got_dir0.rows) == _node_links_dir_tuples(exp_dir)
    end

    function _component_par_by_item(scene, step, cfg, row)
        t = ArchimedLight.component_values_table(
            scene,
            step,
            cfg;
            meteo_row=row,
            columns=["item_id", "component_id", "group", "Ri_PAR_0_f", "Ri_PAR_f"],
        )
        return t.rows
    end

    function _sum_item(rows, item_id::Int, col::String)
        s = 0.0
        for r in rows
            Int(r["item_id"]) == item_id || continue
            s += Float64(r[col])
        end
        s
    end

    for links_name in ("test-links-sensor-plates", "test-links-sensor-plates2")
        fx = fixtures[links_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)

        step_pix = ArchimedLight.run_light_step(
            scene,
            row,
            cfg;
            scattering_backend=ArchimedLight.RaycastScatteringBackend(),
        )
        step_lnk = ArchimedLight.run_light_step(
            scene,
            row,
            cfg;
            scattering_backend=ArchimedLight.LinksScatteringBackend(),
        )

        rows_pix = _component_par_by_item(scene, step_pix, cfg, row)
        rows_lnk = _component_par_by_item(scene, step_lnk, cfg, row)
        @test length(rows_pix) == length(rows_lnk)

        map_pix_0 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"]) for r in rows_pix)
        map_lnk_0 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"]) for r in rows_lnk)
        map_pix_n = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_f"]) for r in rows_pix)
        map_lnk_n = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_f"]) for r in rows_lnk)
        @test maximum(abs(map_pix_0[k] - map_lnk_0[k]) for k in keys(map_pix_0); init=0.0) < parity_limit(:backend_identity_abs)
        @test maximum(abs(map_pix_n[k] - map_lnk_n[k]) for k in keys(map_pix_n); init=0.0) < parity_limit(:backend_identity_abs)

        if links_name == "test-links-sensor-plates"
            plate_scat = sum(
                Float64(r["Ri_PAR_f"] - r["Ri_PAR_0_f"]) for r in rows_pix if string(r["group"]) == "plates";
                init=0.0,
            )
            sensor_meas = sum(Float64(r["Ri_PAR_f"]) for r in rows_pix if string(r["group"]) == "sensors"; init=0.0)
            @test sensor_meas > 0.0
            @test abs(plate_scat - sensor_meas) / max(abs(plate_scat), eps(Float64)) < parity_limit(:links_sensor_plate_balance_rel)
        else
            coeffs = ArchimedLight._group_optical_coeffs(cfg)
            ro = get(get(coeffs, "plates", Dict{String,Float64}()), "PAR", cfg.scattering_coeff_par) / 2.0

            ip1 = _sum_item(rows_pix, 4, "Ri_PAR_f")
            ip2 = _sum_item(rows_pix, 1, "Ri_PAR_f")
            is2 = _sum_item(rows_pix, 2, "Ri_PAR_f")
            is1 = _sum_item(rows_pix, 3, "Ri_PAR_f")

            v_plate2 = ip1 * ro
            @test abs(v_plate2 - ip2) / max(abs(v_plate2), eps(Float64)) < parity_limit(:links_sensor_plates2_plate2_rel)

            v_sensor = ip1 * ro + ip2 * ro
            @test abs(is2 - v_sensor) / max(abs(v_sensor), eps(Float64)) < parity_limit(:links_sensor_plates2_is2_rel)
            @test abs(is1 - v_sensor) / max(abs(v_sensor), eps(Float64)) < parity_limit(:links_sensor_plates2_is1_rel)
        end
    end

    f_sensor4 = fixtures["test-cafeier_sensor4"]
    cfg_sensor4 = ArchimedLight.read_light_config(f_sensor4.config_path)
    scene_sensor4 = ArchimedLight.read_scene(cfg_sensor4.scene)
    row_sensor4 = first(ArchimedLight.read_meteo(cfg_sensor4.meteo).rows)
    step_sensor4 = ArchimedLight.run_light_step(scene_sensor4, row_sensor4, cfg_sensor4)

    summary4 = ArchimedLight.summary_values_table(
        scene_sensor4,
        [step_sensor4],
        cfg_sensor4;
        meteo_rows=[row_sensor4],
        start_step_number=0,
    )
    comp4 = ArchimedLight.component_values_table(
        scene_sensor4,
        step_sensor4,
        cfg_sensor4;
        meteo_row=row_sensor4,
        columns=["Ra_PAR_0_q", "Ra_NIR_0_q"],
    )

    ground_area = sum(
        Float64(r["area"]) for r in summary4.rows if Int(r["item_id"]) == -1;
        init=0.0,
    )
    @test ground_area > 0.0

    sky_irr = step_sensor4.sky.ri_sw_f
    stepdur = _step_duration_seconds(row_sensor4)
    tot_abs_energy = sum(Float64(r["Ra_PAR_0_q"] + r["Ra_NIR_0_q"]) for r in comp4.rows; init=0.0)
    abs_irr = tot_abs_energy / stepdur / ground_area
    @test sky_irr - abs_irr > 0.0

    vs_intercep = sum(Float64(r["Ri_q"]) for r in summary4.rows if string(r["group"]) == "sensors"; init=0.0)
    @test vs_intercep > 0.0
    vs_irr = vs_intercep / stepdur / ground_area
    @test abs(sky_irr - vs_irr) < parity_limit(:sensor4_vs_abs)

    scene_intercep = sum(Float64(r["Ri_q"]) for r in summary4.rows if string(r["group"]) != "sensors"; init=0.0)
    scene_irr = scene_intercep / stepdur / ground_area
    @test abs(sky_irr - scene_irr) < parity_limit(:sensor4_scene_abs)

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

    function _with_raw(cfg::ArchimedLight.LightConfig, raw::Dict{String,Any})
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

    function _counter_dirs(outdir::String)
        d = isdir(outdir) ? readdir(outdir) : String[]
        sort([x for x in d if occursin(r"^\d{6}$", x)])
    end

    function _model_radiance(path::String)
        txt = read(path, String)
        m = match(r"radiance\s*:\s*([-+0-9.eE]+)", txt)
        m === nothing && error("No radiance found in $path")
        parse(Float64, m.captures[1])
    end

    # Java run-test parity: output counter directories keep incrementing even when gaps exist.
    for cfg_path in (
        joinpath(test_root, "test-output-path", "config.yml"),
        joinpath(test_root, "test-output-path2", "project", "config.yml"),
    )
        cfg = ArchimedLight.read_light_config(cfg_path)
        mktempdir() do tmp
            out = joinpath(tmp, "output")
            raw = copy(cfg.raw)
            raw["output_directory"] = out
            haskey(raw, "simulation_directory") && delete!(raw, "simulation_directory")
            cfg2 = _with_raw(cfg, raw)

            @test basename(ArchimedLight.simulation_output_directory(cfg2)) == "000001"
            @test basename(ArchimedLight.simulation_output_directory(cfg2)) == "000002"
            @test length(_counter_dirs(out)) == 2

            rm(joinpath(out, "000002"); recursive=true, force=true)
            @test basename(ArchimedLight.simulation_output_directory(cfg2)) == "000002"
            @test basename(ArchimedLight.simulation_output_directory(cfg2)) == "000003"

            rm(joinpath(out, "000002"); recursive=true, force=true)
            @test basename(ArchimedLight.simulation_output_directory(cfg2)) == "000004"
            @test _counter_dirs(out) == ["000001", "000003", "000004"]
        end
    end

    # Java run-test parity: simulation_directory is reused and cleaned between runs.
    cfg_simdir = ArchimedLight.read_light_config(joinpath(test_root, "test-simulation-dir", "config.yml"))
    mktempdir() do tmp
        out = joinpath(tmp, "output")
        raw = copy(cfg_simdir.raw)
        raw["output_directory"] = out
        cfg2 = _with_raw(cfg_simdir, raw)

        sim1 = ArchimedLight.simulation_output_directory(cfg2)
        @test basename(sim1) == "simdir"
        dummy = joinpath(sim1, "dummy.txt")
        write(dummy, "x")
        @test isfile(dummy)

        sim2 = ArchimedLight.simulation_output_directory(cfg2)
        @test sim2 == sim1
        @test !isfile(dummy)
        @test count(==("simdir"), readdir(out)) == 1
    end

    function _run_fixture_series(cfg_path::String)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        return cfg, scene, meteo, series
    end

    function _xml_tag_value(path::String, tag::String)
        txt = read(path, String)
        m = match(Regex("<$(tag)>([^<]+)</$(tag)>"), txt)
        m === nothing && error("Missing tag <$tag> in $path")
        parse(Float64, strip(m.captures[1]))
    end

    # Java test-ops-output parity: export_ops selector and filenames.
    for (cfg_name, expected_count, expected_one_name) in (
        ("config1.yml", 1, "scene-step00004.ops"),
        ("config2.yml", 0, nothing),
        ("config3.yml", 4, nothing),
        ("config4.yml", 1, "scene-step00003.ops"),
    )
        cfg, scene, meteo, series = _run_fixture_series(joinpath(test_root, "test-ops-output", cfg_name))
        mktempdir() do out
            ArchimedLight.write_light_outputs(
                scene,
                series,
                cfg;
                meteo_rows=meteo.rows,
                outdir=out,
                write_component=false,
                write_scene=false,
                write_summary=false,
                write_sun_position_log=false,
                write_scattering_log=false,
            )
            ops_files = sort(filter(f -> endswith(f, ".ops"), readdir(out; join=true)))
            @test length(ops_files) == expected_count
            if expected_one_name !== nothing
                @test basename(only(ops_files)) == expected_one_name
            end
        end
    end

    for cfg_name in ("config5.yml", "config6.yml")
        cfg, scene, meteo, series = _run_fixture_series(joinpath(test_root, "test-ops-output", cfg_name))
        mktempdir() do out
            threw = false
            msg = ""
            try
                ArchimedLight.write_light_outputs(
                    scene,
                    series,
                    cfg;
                    meteo_rows=meteo.rows,
                    outdir=out,
                    write_component=false,
                    write_scene=false,
                    write_summary=false,
                    write_sun_position_log=false,
                    write_scattering_log=false,
                )
            catch e
                threw = e isa ErrorException
                msg = sprint(showerror, e)
            end
            @test threw
            @test occursin("invalid export_ops parameter", msg)
        end
    end

    # Java test-ops-output2/3/4 parity: exported GWA irradiance tags follow requested step(s).
    for (fx_name, expected_first_zero) in (("test-ops-output2", true), ("test-ops-output3", false))
        cfg, scene, meteo, series = _run_fixture_series(joinpath(test_root, fx_name, "config.yml"))
        mktempdir() do out
            ArchimedLight.write_light_outputs(
                scene,
                series,
                cfg;
                meteo_rows=meteo.rows,
                outdir=out,
                write_component=false,
                write_scene=false,
                write_summary=false,
                write_sun_position_log=false,
                write_scattering_log=false,
            )
            gwa_files = sort(filter(f -> endswith(lowercase(f), ".gwa"), readdir(out; join=true)))
            @test length(gwa_files) == 1
            ri_par = _xml_tag_value(only(gwa_files), "Ri_PAR_0_f")
            ri_nir = _xml_tag_value(only(gwa_files), "Ri_NIR_0_f")
            if expected_first_zero
                @test abs(ri_par) < parity_limit(:ops_export_zero_abs)
                @test abs(ri_nir) < parity_limit(:ops_export_zero_abs)
            else
                @test abs(ri_par) > parity_limit(:ops_export_zero_abs)
                @test abs(ri_nir) > parity_limit(:ops_export_zero_abs)
            end
        end
    end

    cfg_ops_all, scene_ops_all, meteo_ops_all, series_ops_all =
        _run_fixture_series(joinpath(test_root, "test-ops-output4", "config.yml"))
    mktempdir() do out
        ArchimedLight.write_light_outputs(
            scene_ops_all,
            series_ops_all,
            cfg_ops_all;
            meteo_rows=meteo_ops_all.rows,
            outdir=out,
            write_component=false,
            write_scene=false,
            write_summary=false,
            write_sun_position_log=false,
            write_scattering_log=false,
        )
        gwa_files = filter(f -> endswith(lowercase(f), ".gwa"), readdir(out; join=true))
        step_of(path) = parse(Int, only(match(r"step([0-9]+)", basename(path)).captures))
        ordered = sort(gwa_files; by=step_of)
        @test length(ordered) == length(meteo_ops_all.rows)
        @test [step_of(f) for f in ordered] == collect(1:length(ordered))

        for (i, f) in enumerate(ordered)
            ri_par = _xml_tag_value(f, "Ri_PAR_0_f")
            ri_nir = _xml_tag_value(f, "Ri_NIR_0_f")
            if i == 1
                @test abs(ri_par) < parity_limit(:ops_export_zero_abs)
                @test abs(ri_nir) < parity_limit(:ops_export_zero_abs)
            else
                @test abs(ri_par) > parity_limit(:ops_export_zero_abs)
                @test abs(ri_nir) > parity_limit(:ops_export_zero_abs)
            end
        end
    end

    # Java test-sortcsv parity.
    sortcsv_dir = joinpath(test_root, "test-sortcsv")
    sort_input = joinpath(sortcsv_dir, "input.csv")
    sort_expected = joinpath(sortcsv_dir, "expected.csv")
    _cell_equal(a, b) = begin
        sa = strip(string(a))
        sb = strip(string(b))
        pa = tryparse(Float64, sa)
        pb = tryparse(Float64, sb)
        if pa !== nothing && pb !== nothing
            return isapprox(pa, pb; atol=1e-9, rtol=1e-9)
        end
        sa == sb
    end
    mktempdir() do tmp
        out_sort = joinpath(tmp, "out.csv")
        ArchimedLight.sort_csv_file(sort_input, out_sort; sort_columns="stepNumber;plantId;nodeId")
        got_rows = read_java_csv(out_sort)
        exp_rows = read_java_csv(sort_expected)
        @test length(got_rows) == length(exp_rows)
        @test !isempty(got_rows)
        cols = propertynames(first(exp_rows))
        @test cols == propertynames(first(got_rows))
        for i in eachindex(exp_rows)
            for c in cols
                @test _cell_equal(getproperty(got_rows[i], c), getproperty(exp_rows[i], c))
            end
        end

        out_sort2 = joinpath(tmp, "out2.csv")
        ArchimedLight.sort_csv_file(sort_input, out_sort2; sort_columns="stepNumber,plantId,nodeId")
        got_rows2 = read_java_csv(out_sort2)
        @test length(got_rows2) == length(exp_rows)
        for i in eachindex(exp_rows)
            for c in cols
                @test _cell_equal(getproperty(got_rows2[i], c), getproperty(exp_rows[i], c))
            end
        end
    end

    function _run_fixture_step(cfg_path::String)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        return cfg, scene, row, step
    end

    # Java test-lightsource1/2 parity: sky-only run and emitter-only run match on receiver components.
    for fx_name in ("test-lightsource1", "test-lightsource2")
        tol = parity_limit(:lightsource12_max_abs; fixture=fx_name)
        cfg1, scene1, row1, step1 = _run_fixture_step(joinpath(test_root, fx_name, "config1.yml"))
        cfg2, scene2, row2, step2 = _run_fixture_step(joinpath(test_root, fx_name, "config2.yml"))

        rows1 = ArchimedLight.component_values_table(
            scene1,
            step1,
            cfg1;
            meteo_row=row1,
            columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f"],
        ).rows
        rows2 = ArchimedLight.component_values_table(
            scene2,
            step2,
            cfg2;
            meteo_row=row2,
            columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f"],
        ).rows

        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows1)
        d2 = Dict(
            (Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows2 if
            Int(r["item_id"]) != 2
        )
        ks = intersect(keys(d1), keys(d2))
        @test !isempty(ks)
        @test maximum(abs(d1[k] - d2[k]) for k in ks; init=0.0) < tol
    end

    # Java test-lightsource3 parity: linear additivity of sky and emitter contributions.
    cfg_l1, scene_l3, row_l1, step_l1 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config1.yml"))
    cfg_l2, _, row_l2, step_l2 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config2.yml"))
    cfg_l3, _, row_l3, step_l3 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config3.yml"))
    t1 = ArchimedLight.component_values_table(
        scene_l3,
        step_l1,
        cfg_l1;
        meteo_row=row_l1,
        columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"],
    ).rows
    t2 = ArchimedLight.component_values_table(
        scene_l3,
        step_l2,
        cfg_l2;
        meteo_row=row_l2,
        columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"],
    ).rows
    t3 = ArchimedLight.component_values_table(
        scene_l3,
        step_l3,
        cfg_l3;
        meteo_row=row_l3,
        columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"],
    ).rows
    for col in ("Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f")
        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t1)
        d2 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t2)
        d3 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t3)
        ks = intersect(intersect(keys(d1), keys(d2)), keys(d3))
        @test !isempty(ks)
        @test maximum(abs((d1[k] + d2[k]) - d3[k]) for k in ks; init=0.0) < parity_limit(:lightsource_additivity_abs)
    end

    # Java test-lightsource4/5/6 parity: emitted source power equals intercepted order-0 power.
    for (fx_name, model_name, strict_one_node) in (
        ("test-lightsource4", "model_pave2.yml", true),
        ("test-lightsource5", "model_lamp.yml", false),
        ("test-lightsource6", "model_lamp.yml", false),
    )
        cfg, scene, row, step = _run_fixture_step(joinpath(test_root, fx_name, "config.yml"))
        rows = ArchimedLight.component_values_table(
            scene,
            step,
            cfg;
            meteo_row=row,
            columns=["item_id", "component_id", "area", "Ri_PAR_0_f", "Ri_NIR_0_f"],
        ).rows
        source_power = _model_radiance(joinpath(test_root, fx_name, model_name))

        if strict_one_node
            vals = [Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows if Int(r["item_id"]) == 1]
            @test !isempty(vals)
            @test relerr(maximum(vals), source_power) < parity_limit(:lightsource_power_rel)
        else
            tot = sum(Float64(r["area"]) * Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows; init=0.0)
            @test relerr(tot, source_power) < parity_limit(:lightsource_power_rel)
        end
    end

    # Java test-lightsource7/8/9 parity: absorbed power from Ri_f and scat_factor matches source power.
    for (fx_name, model_name) in (
        ("test-lightsource7", "model_pave2.yml"),
        ("test-lightsource8", "model_lamp.yml"),
        ("test-lightsource9", "model_lamp.yml"),
    )
        cfg, scene, row, step = _run_fixture_step(joinpath(test_root, fx_name, "config.yml"))
        rows = ArchimedLight.component_values_table(
            scene,
            step,
            cfg;
            meteo_row=row,
            columns=["area", "Ri_PAR_f", "Ri_NIR_f", "scat_factor_PAR", "scat_factor_NIR"],
        ).rows
        source_power = _model_radiance(joinpath(test_root, fx_name, model_name))
        tot_abs = sum(
            (Float64(r["Ri_PAR_f"]) * (1.0 - Float64(r["scat_factor_PAR"])) +
            Float64(r["Ri_NIR_f"]) * (1.0 - Float64(r["scat_factor_NIR"]))) * Float64(r["area"]) for r in rows;
            init=0.0,
        )
        @test relerr(tot_abs, source_power) < parity_limit(:lightsource_power_rel)
    end

    # Java test-lightsource-box parity: closed box should intercept emitted lamp power.
    cfg_box, scene_box, row_box, step_box = _run_fixture_step(joinpath(test_root, "test-lightsource-box", "config.yml"))
    rows_box = ArchimedLight.component_values_table(
        scene_box,
        step_box,
        cfg_box;
        meteo_row=row_box,
        columns=["area", "Ri_PAR_0_f", "Ri_NIR_0_f"],
    ).rows
    source_box = _model_radiance(joinpath(test_root, "test-lightsource-box", "model_box.yml"))
    tot_box = sum(Float64(r["area"]) * (Float64(r["Ri_PAR_0_f"]) + Float64(r["Ri_NIR_0_f"])) for r in rows_box; init=0.0)
    @test relerr(tot_box, source_box) < parity_limit(:lightsource_power_rel)

    # Java test-lightsource-box2 parity: absorbed power with scattering matches emitted lamp power.
    cfg_box2, scene_box2, row_box2, step_box2 = _run_fixture_step(joinpath(test_root, "test-lightsource-box2", "config.yml"))
    rows_box2 = ArchimedLight.component_values_table(
        scene_box2,
        step_box2,
        cfg_box2;
        meteo_row=row_box2,
        columns=["area", "Ri_PAR_f", "Ri_NIR_f", "scat_factor_PAR", "scat_factor_NIR"],
    ).rows
    source_box2 = _model_radiance(joinpath(test_root, "test-lightsource-box2", "model_box.yml"))
    tot_box2 = sum(
        (
            Float64(r["Ri_PAR_f"]) * (1.0 - Float64(r["scat_factor_PAR"]) / 2.0) +
            Float64(r["Ri_NIR_f"]) * (1.0 - Float64(r["scat_factor_NIR"]) / 2.0)
        ) * Float64(r["area"]) for r in rows_box2;
        init=0.0,
    )
    @test relerr(tot_box2, source_box2) < parity_limit(:lightsource_box2_rel)

    function _check_summary_component_parity(cfg_path::String; use_scattering::Bool, max_scene_rel::Union{Nothing,Float64}=nothing)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        summary = ArchimedLight.summary_values_table(scene, series, cfg; meteo_rows=meteo.rows, start_step_number=0)

        max_rel = 0.0
        for i in eachindex(series)
            step = series[i]
            row = meteo.rows[i]
            dt = _step_duration_seconds(row)
            cols =
                if use_scattering
                    ["group", "type", "item_id", "area", "Ri_PAR_f", "Ri_NIR_f"]
                else
                    ["group", "type", "item_id", "area", "Ri_PAR_0_f", "Ri_NIR_0_f"]
                end
            crows = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, columns=cols).rows

            agg = Dict{Tuple{String,String,Int},Tuple{Float64,Float64}}()
            for r in crows
                key = (string(r["group"]), string(r["type"]), Int(r["item_id"]))
                area = Float64(r["area"])
                riq =
                    if use_scattering
                        (Float64(r["Ri_PAR_f"]) + Float64(r["Ri_NIR_f"])) * area * dt
                    else
                        (Float64(r["Ri_PAR_0_f"]) + Float64(r["Ri_NIR_0_f"])) * area * dt
                    end
                cur = get(agg, key, (0.0, 0.0))
                agg[key] = (cur[1] + area, cur[2] + riq)
            end

            srows = [r for r in summary.rows if Int(r["step_number"]) == i - 1]
            for r in srows
                key = (string(r["group"]), string(r["type"]), Int(r["item_id"]))
                area_exp, riq_exp = get(agg, key, (0.0, 0.0))
                area_got = Float64(r["area"])
                riq_got = Float64(r["Ri_q"])
                max_rel = max(max_rel, relerr(area_got, area_exp))
                max_rel = max(max_rel, relerr(riq_got, riq_exp))
            end
        end
        @test max_rel < parity_limit(:summary_component_rel)

        if max_scene_rel !== nothing
            ground = [
                Float64(r["area"]) for r in summary.rows if
                Int(r["step_number"]) == 0 && Int(r["item_id"]) == -1
            ]
            @test !isempty(ground)
            ground_area = first(ground)
            max_rel_scene = 0.0
            for i in eachindex(series)
                dt = _step_duration_seconds(meteo.rows[i])
                srows = [r for r in summary.rows if Int(r["step_number"]) == i - 1]
                ri_scene = sum(Float64(r["Ri_q"]) for r in srows) / dt / ground_area
                max_rel_scene = max(max_rel_scene, relerr(ri_scene, series[i].sky.ri_sw_f))
            end
            @test max_rel_scene < max_scene_rel
        end
    end

    # Java summary tests parity.
    _check_summary_component_parity(
        joinpath(test_root, "test-summary", "config.yml");
        use_scattering=false,
        max_scene_rel=parity_limit(:summary_scene_rel; fixture="test-summary"),
    )
    _check_summary_component_parity(
        joinpath(test_root, "test-summary2", "config.yml");
        use_scattering=false,
        max_scene_rel=parity_limit(:summary_scene_rel; fixture="test-summary2"),
    )
    _check_summary_component_parity(joinpath(test_root, "test-summary3", "config.yml"); use_scattering=true)
    _check_summary_component_parity(joinpath(test_root, "test-summary4", "config.yml"); use_scattering=true)

    # Java test-reducted-table parity: reduced and full pixel table paths must match Ri outputs.
    cfg_red1 = ArchimedLight.read_light_config(joinpath(test_root, "test-reducted-table", "config1.yml"))
    cfg_red2 = ArchimedLight.read_light_config(joinpath(test_root, "test-reducted-table", "config2.yml"))
    scene_red1 = ArchimedLight.read_scene(cfg_red1.scene)
    scene_red2 = ArchimedLight.read_scene(cfg_red2.scene)
    row_red1 = first(ArchimedLight.read_meteo(cfg_red1.meteo).rows)
    row_red2 = first(ArchimedLight.read_meteo(cfg_red2.meteo).rows)
    step_red1 = ArchimedLight.run_light_step(scene_red1, row_red1, cfg_red1)
    step_red2 = ArchimedLight.run_light_step(scene_red2, row_red2, cfg_red2)
    cols_red = ["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"]
    t_red1 = ArchimedLight.component_values_table(scene_red1, step_red1, cfg_red1; meteo_row=row_red1, columns=cols_red).rows
    t_red2 = ArchimedLight.component_values_table(scene_red2, step_red2, cfg_red2; meteo_row=row_red2, columns=cols_red).rows
    for col in cols_red[3:end]
        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t_red1)
        d2 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t_red2)
        @test maximum(abs(get(d1, k, 0.0) - get(d2, k, 0.0)) for k in union(keys(d1), keys(d2)); init=0.0) == 0.0
    end

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
    @test relerr(r_par_0, r_custom_0) < parity_limit(:customband_rel)
    @test relerr(r_par_n, r_custom_n) < parity_limit(:customband_rel)

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

    # Java test-cached-radiation4 parity (light-only subset): cached and uncached
    # runs must match on light outputs even when photo/energy sections are enabled.
    cached4_root = joinpath(test_root, "test-cached-radiation4")
    cfg4_1 = joinpath(cached4_root, "config.yml")
    cfg4_2 = joinpath(cached4_root, "config2.yml")
    @test isfile(cfg4_1)
    @test isfile(cfg4_2)
    s4_1 = run_fixture_series(cfg4_1)
    s4_2 = run_fixture_series(cfg4_2)
    @test length(s4_1) == length(s4_2)
    for i in eachindex(s4_1)
        @test max_abs_float_dict_diff(s4_1[i].first_order.projected_area_per_node, s4_2[i].first_order.projected_area_per_node) == 0.0
        @test max_abs_int_dict_diff(s4_1[i].first_order.hits_per_node, s4_2[i].first_order.hits_per_node) == 0
        @test max_abs_float_dict_diff(s4_1[i].budget.ri_par_q_per_node, s4_2[i].budget.ri_par_q_per_node) == 0.0
        @test max_abs_float_dict_diff(s4_1[i].budget.ri_nir_q_per_node, s4_2[i].budget.ri_nir_q_per_node) == 0.0
        @test max_abs_float_dict_diff(s4_1[i].budget.ra_par_q_per_node, s4_2[i].budget.ra_par_q_per_node) == 0.0
        @test max_abs_float_dict_diff(s4_1[i].budget.ra_nir_q_per_node, s4_2[i].budget.ra_nir_q_per_node) == 0.0
    end

    # Java test-parallel and test-parallel2 parity:
    # cache-radiation toggles must keep component light outputs stable.
    function _parallel_cfg_with_overrides(
        cfg::ArchimedLight.LightConfig;
        all_in_turtle::Union{Nothing,Bool}=nothing,
        cache_radiation::Union{Nothing,Bool}=nothing,
    )
        raw = copy(cfg.raw)
        a = all_in_turtle === nothing ? cfg.all_in_turtle : Bool(all_in_turtle)
        c = cache_radiation === nothing ? cfg.cache_radiation : Bool(cache_radiation)
        raw["all_in_turtle"] = a
        raw["cache_radiation"] = c
        ArchimedLight.LightConfig(
            cfg.scene,
            cfg.meteo,
            a,
            cfg.turtle_sectors,
            cfg.pixel_size,
            cfg.area_ratio,
            cfg.scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            c,
            raw,
        )
    end

    function _parallel_rows(fx_name::String; all_in_turtle::Union{Nothing,Bool}=nothing, cache_radiation::Bool=false)
        fx = fixtures[fx_name]
        cfg0 = ArchimedLight.read_light_config(fx.config_path)
        cfg = _parallel_cfg_with_overrides(cfg0; all_in_turtle=all_in_turtle, cache_radiation=cache_radiation)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        series = ArchimedLight.run_light_series(scene, selected, cfg)
        _component_series_rows(scene, series, cfg, selected.rows)
    end

    function _assert_parallel_fixture_cache_invariance(fx_name::String; all_in_turtle::Union{Nothing,Bool}=nothing)
        rows_uncached = _parallel_rows(fx_name; all_in_turtle=all_in_turtle, cache_radiation=false)
        rows_cached = _parallel_rows(fx_name; all_in_turtle=all_in_turtle, cache_radiation=true)
        _assert_rows_match(
            rows_uncached,
            rows_cached,
            ["step_number", "item_id", "component_id"],
            ["Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"];
            atol=parity_limit(:absorb_cached_abs),
            rtol=0.0,
            label="parallel cache invariance fixture=$(fx_name) all_in_turtle=$(all_in_turtle)",
        )
    end

    _assert_parallel_fixture_cache_invariance("test-parallel")
    for all_in_turtle in (true, false)
        _assert_parallel_fixture_cache_invariance("test-parallel2"; all_in_turtle=all_in_turtle)
    end

    # Frozen Java snapshots for compare fixtures (light-only scope).
    function _csv_columns(rows)
        isempty(rows) && return String[]
        String[string(n) for n in propertynames(first(rows))]
    end

    function _stack_component_rows_with_columns(scene, series, cfg, meteo_rows, cols::Vector{String})
        out = Dict{String,Any}[]
        for i in eachindex(series)
            rows = ArchimedLight.component_values_table(
                scene,
                series[i],
                cfg;
                meteo_row=meteo_rows[i],
                step_number=i - 1,
                columns=cols,
            ).rows
            append!(out, rows)
        end
        sort!(out; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
        out
    end

    compare_fx_names = (
        "test-compare-cafeier1",
        "test-compare-cafeier2",
        "test-compare-cafeier3",
        "test-compare-cafeier4",
        "test-compare-cafeier5",
        "test-compare-simpleplant",
        "test-compare-two_coffee",
        "test-compare-timestep",
    )
    for fx_name in compare_fx_names
        fx = fixtures[fx_name]
        expected_comp_path = _expected_component_values_path(fx)
        expected_scene_path = _expected_scene_values_path(fx)
        @test expected_comp_path !== nothing
        @test expected_scene_path !== nothing

        expected_comp = read_java_csv(expected_comp_path)
        expected_scene = read_java_csv(expected_scene_path)
        @test !isempty(expected_comp)
        @test !isempty(expected_scene)
        sort!(expected_comp; by=r -> (to_int(getproperty(r, :step_number)), to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id))))
        sort!(expected_scene; by=r -> to_int(getproperty(r, :step_number)))

        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        @test length(series) == length(selected.rows)

        comp_cols = _csv_columns(expected_comp)
        scene_cols = _csv_columns(expected_scene)
        observed_comp = _stack_component_rows_with_columns(scene, series, cfg, selected.rows, comp_cols)
        observed_scene = ArchimedLight.scene_values_table(scene, series, cfg; meteo_rows=selected.rows, columns=scene_cols).rows
        sort!(observed_scene; by=r -> Int(r["step_number"]))

        _assert_rows_match(
            expected_comp,
            observed_comp,
            ["step_number", "item_id", "component_id"],
            [c for c in comp_cols if !(c in ("step_number", "item_id", "component_id"))];
            atol=1e-3,
            rtol=3e-3,
            label="compare fixture component snapshot $(fx_name)",
        )
        _assert_rows_match(
            expected_scene,
            observed_scene,
            ["step_number"],
            [c for c in scene_cols if c != "step_number"];
            atol=1e-3,
            rtol=3e-3,
            label="compare fixture scene snapshot $(fx_name)",
        )

        expected_summary_path = _expected_summary_path(fx)
        if expected_summary_path !== nothing
            expected_summary = read_java_csv(expected_summary_path)
            !isempty(expected_summary) || continue
            sort!(
                expected_summary;
                by=r -> (
                    to_int(getproperty(r, :step_number)),
                    string(getproperty(r, :group)),
                    string(getproperty(r, :type)),
                    to_int(getproperty(r, :item_id)),
                ),
            )
            summary_cols = _csv_columns(expected_summary)
            summary_key_cols = [c for c in ("step_number", "group", "type", "item_id") if c in summary_cols]
            observed_summary = ArchimedLight.summary_values_table(
                scene,
                series,
                cfg;
                meteo_rows=selected.rows,
                start_step_number=0,
                columns=summary_cols,
            ).rows
            sort!(
                observed_summary;
                by=r -> (
                    Int(get(r, "step_number", 0)),
                    string(get(r, "group", "")),
                    string(get(r, "type", "")),
                    Int(get(r, "item_id", 0)),
                ),
            )
            _assert_rows_match(
                expected_summary,
                observed_summary,
                summary_key_cols,
                [c for c in summary_cols if !(c in summary_key_cols)];
                atol=1e-3,
                rtol=3e-3,
                label="compare fixture summary snapshot $(fx_name)",
            )
        end
    end
    end

        _stage_done("Scattering links outputs and lightsource parity", _stage_t0)
    end

    if _RUN_PARITY_TESTS
    @testset "Meteo and independent-step parity" begin
        _stage_t0 = _stage_start("Meteo and independent-step parity")
    function _light_cfg_with_meteo(cfg::ArchimedLight.LightConfig, meteo_path::String)
        raw = copy(cfg.raw)
        raw["meteo"] = meteo_path
        ArchimedLight.LightConfig(
            cfg.scene,
            meteo_path,
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

    function _light_cfg_with_meteo_controls(
        cfg::ArchimedLight.LightConfig;
        meteo_path::Union{Nothing,String}=nothing,
        meteo_range::Union{Nothing,String}=nothing,
    )
        raw = copy(cfg.raw)
        meteo_file = meteo_path === nothing ? cfg.meteo : meteo_path
        raw["meteo"] = meteo_file
        props0 = get(raw, "prop", Dict{String,Any}())
        props = Dict{String,Any}(string(k) => v for (k, v) in pairs(props0))
        delete!(props, "meteo_range")
        delete!(props, "meteoRange")
        delete!(props, "meteorange")
        if meteo_range !== nothing
            props["meteo_range"] = meteo_range
        end
        raw["prop"] = props
        ArchimedLight.LightConfig(
            cfg.scene,
            meteo_file,
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

    # Java test-meteo1/test-meteo2 parity: invalid meteo chronology must fail with explicit errors.
    for (fx_name, expected_msg) in (
        ("test-meteo1", "invalid overlapping meteo steps"),
        ("test-meteo2", "end is before start"),
    )
        cfg_bad = ArchimedLight.read_light_config(joinpath(test_root, fx_name, "config.yml"))
        scene_bad = ArchimedLight.read_scene(cfg_bad.scene)
        meteo_bad = ArchimedLight.read_meteo(cfg_bad.meteo)
        err = nothing
        try
            ArchimedLight.run_light_series(scene_bad, meteo_bad, cfg_bad)
        catch e
            err = e
        end
        @test err !== nothing
        @test occursin(expected_msg, sprint(showerror, err))
    end

    # Java test-meteoactive parity: active column filtering and failures.
    meteo_active_root = joinpath(test_root, "test-meteoactive")
    cfg_active_base = ArchimedLight.read_light_config(joinpath(meteo_active_root, "config.yml"))
    scene_active = ArchimedLight.read_scene(cfg_active_base.scene)

    function _active_steps(file_name::String; meteo_range::Union{Nothing,String}=nothing)
        cfg = _light_cfg_with_meteo_controls(
            cfg_active_base;
            meteo_path=joinpath(meteo_active_root, file_name),
            meteo_range=meteo_range,
        )
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        series = ArchimedLight.run_light_series(scene_active, meteo, cfg)
        @test length(series) == length(selected.rows)
        Int[round(Int, Float64(getproperty(r, :wind))) for r in selected.rows]
    end

    @test _active_steps("meteo1.csv") == [1, 2, 3, 4]
    @test _active_steps("meteo2.csv") == [1, 2, 3, 4, 5, 6]
    @test _active_steps("meteo3.csv") == [2, 5, 6]
    @test_throws ErrorException _active_steps("meteo4.csv")
    @test_throws ErrorException _active_steps("meteo5.csv")
    @test _active_steps("meteo2.csv"; meteo_range="2,5") == [2, 3, 4, 5]
    @test _active_steps("meteo2.csv"; meteo_range="2016/06/12-09:15:00,2016/06/12-10:45:00") == [3, 4, 5, 6]
    @test _active_steps("meteo3.csv"; meteo_range="2,5") == [2, 5]
    @test _active_steps("meteo3.csv"; meteo_range="2016/06/12-08:45:00,2016/06/12-10:45:00") == [2, 5, 6]

    # Java test-meteorange parity: index/date range selection and invalid range checks.
    meteo_range_root = joinpath(test_root, "test-meteorange")
    cfg_range_base = ArchimedLight.read_light_config(joinpath(meteo_range_root, "config.yml"))
    scene_range = ArchimedLight.read_scene(cfg_range_base.scene)
    meteo_range_all = ArchimedLight.read_meteo(cfg_range_base.meteo)

    function _steps_for_range(spec::String)
        cfg = _light_cfg_with_meteo_controls(cfg_range_base; meteo_range=spec)
        selected = ArchimedLight.prepare_meteo(meteo_range_all, cfg)
        series = ArchimedLight.run_light_series(scene_range, meteo_range_all, cfg)
        @test length(series) == length(selected.rows)
        Int[round(Int, Float64(getproperty(r, :wind))) for r in selected.rows]
    end

    @test _steps_for_range("2,5") == [2, 3, 4, 5]
    @test _steps_for_range("3,3") == [3]
    @test_throws ErrorException _steps_for_range("3,2")
    @test_throws ErrorException _steps_for_range("0,2")
    @test _steps_for_range("2016/06/12-08:40:00,2016/06/12-08:50:00") == [2]
    @test _steps_for_range("2016/06/12-08:30:00,2016/06/12-08:30:00") == [1, 2]
    @test _steps_for_range("2016/06/12-08:40:00,2016/06/12-09:15:00") == [2, 3]
    @test _steps_for_range("2016/06/12-08:45:00,2016/06/12-09:45:00") == [2, 3, 4]
    @test _steps_for_range("2016/06/12-07:45:00,2016/06/12-11:15:00") == [1, 2, 3, 4, 5, 6]

    # Java test-rhvpd parity: RH/VPD columns and use-tags consistency.
    rhvpd_root = joinpath(test_root, "test-rhvpd")
    cfg_rhvpd_base = ArchimedLight.read_light_config(joinpath(rhvpd_root, "config.yml"))
    scene_rhvpd = ArchimedLight.read_scene(cfg_rhvpd_base.scene)

    function _run_rhvpd(file_name::String)
        cfg = _light_cfg_with_meteo(cfg_rhvpd_base, joinpath(rhvpd_root, file_name))
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        ArchimedLight.run_light_step(scene_rhvpd, first(meteo.rows), cfg)
    end

    for file_name in ("meteo1.csv", "meteo2.csv", "meteo3.csv", "meteo4.csv", "meteo5.csv", "meteo6.csv")
        @test _run_rhvpd(file_name).sky.ri_sw_f > 0.0
    end
    for file_name in ("meteo7.csv", "meteo8.csv", "meteo9.csv", "meteo10.csv", "meteo11.csv", "meteo12.csv")
        @test_throws ErrorException _run_rhvpd(file_name)
    end

    # Java test-meteo-stepduration parity: hour_end and step_duration encodings are equivalent;
    # malformed step_duration must fail.
    meteo_step_root = joinpath(test_root, "test-meteo-stepduration")
    cfg_meteo_base = ArchimedLight.read_light_config(joinpath(meteo_step_root, "config.yml"))
    scene_meteo = ArchimedLight.read_scene(cfg_meteo_base.scene)

    function _run_meteo_step(file_name::String)
        cfg = _light_cfg_with_meteo(cfg_meteo_base, joinpath(meteo_step_root, file_name))
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        row = first(meteo.rows)
        step = ArchimedLight.run_light_step(scene_meteo, row, cfg)
        return row, step
    end

    row_m1, step_m1 = _run_meteo_step("meteo1.csv")
    row_m2, step_m2 = _run_meteo_step("meteo2.csv")
    @test isapprox(step_m1.sky.ri_sw_f, step_m2.sky.ri_sw_f; atol=1e-9, rtol=1e-9)
    @test isapprox(ArchimedLight._step_duration_seconds_local(row_m1), 1800.0; atol=1e-12, rtol=1e-12)
    @test isapprox(ArchimedLight._step_duration_seconds_local(row_m2), 1800.0; atol=1e-12, rtol=1e-12)

    cfg_m3 = _light_cfg_with_meteo(cfg_meteo_base, joinpath(meteo_step_root, "meteo3.csv"))
    row_m3 = first(ArchimedLight.read_meteo(cfg_m3.meteo).rows)
    @test_throws ErrorException ArchimedLight.run_light_step(scene_meteo, row_m3, cfg_m3)

    # Java test-independant-steps2 parity: after dropping the first non-common step and
    # renumbering, component outputs are identical.
    indep_root = joinpath(test_root, "test-independant-steps2")
    cfg_indep_base = ArchimedLight.read_light_config(joinpath(indep_root, "config.yml"))
    cfg_indep_1 = _light_cfg_with_meteo(cfg_indep_base, joinpath(indep_root, "meteo1.csv"))
    cfg_indep_2 = _light_cfg_with_meteo(cfg_indep_base, joinpath(indep_root, "meteo2.csv"))
    scene_indep = ArchimedLight.read_scene(cfg_indep_base.scene)
    meteo_indep_1 = ArchimedLight.read_meteo(cfg_indep_1.meteo)
    meteo_indep_2 = ArchimedLight.read_meteo(cfg_indep_2.meteo)
    series_indep_1 = ArchimedLight.run_light_series(scene_indep, meteo_indep_1, cfg_indep_1)
    series_indep_2 = ArchimedLight.run_light_series(scene_indep, meteo_indep_2, cfg_indep_2)

    component_cols = ArchimedLight.component_variable_names(cfg_indep_base)
    function _stack_component_rows(scene, series, cfg, meteo_rows, cols)
        out = Dict{String,Any}[]
        for i in eachindex(series)
            rows = ArchimedLight.component_values_table(
                scene,
                series[i],
                cfg;
                meteo_row=meteo_rows[i],
                step_number=i - 1,
                columns=cols,
            ).rows
            append!(out, rows)
        end
        sort!(out; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
        out
    end
    rows_indep_1 = _stack_component_rows(scene_indep, series_indep_1, cfg_indep_1, meteo_indep_1.rows, component_cols)
    rows_indep_2_raw = _stack_component_rows(scene_indep, series_indep_2, cfg_indep_2, meteo_indep_2.rows, component_cols)
    rows_indep_2 = Dict{String,Any}[]
    for r in rows_indep_2_raw
        step_no = Int(r["step_number"])
        step_no == 0 && continue
        rr = copy(r)
        rr["step_number"] = step_no - 1
        push!(rows_indep_2, rr)
    end
    sort!(rows_indep_2; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
    _assert_rows_match(
        rows_indep_1,
        rows_indep_2,
        ["step_number", "item_id", "component_id"],
        component_cols;
        atol=1e-9,
        rtol=1e-9,
        label="test-independant-steps2 component rows",
    )

    # Java test-independant-steps parity (light-only subset): after dropping the first
    # non-common step and renumbering, overlapping light outputs are identical.
    indep_root_light = joinpath(test_root, "test-independant-steps")
    cfg_indep_light_base = ArchimedLight.read_light_config(joinpath(indep_root_light, "config.yml"))
    cfg_indep_light_1 = _light_cfg_with_meteo(cfg_indep_light_base, joinpath(indep_root_light, "meteo1.csv"))
    cfg_indep_light_2 = _light_cfg_with_meteo(cfg_indep_light_base, joinpath(indep_root_light, "meteo2.csv"))
    scene_indep_light = ArchimedLight.read_scene(cfg_indep_light_base.scene)
    meteo_indep_light_1 = ArchimedLight.read_meteo(cfg_indep_light_1.meteo)
    meteo_indep_light_2 = ArchimedLight.read_meteo(cfg_indep_light_2.meteo)
    series_indep_light_1 = ArchimedLight.run_light_series(scene_indep_light, meteo_indep_light_1, cfg_indep_light_1)
    series_indep_light_2 = ArchimedLight.run_light_series(scene_indep_light, meteo_indep_light_2, cfg_indep_light_2)

    cols_indep_light = [
        "step_number",
        "item_id",
        "component_id",
        "sky_fraction",
        "Ri_PAR_0_f",
        "Ri_NIR_0_f",
        "Ri_PAR_0_q",
        "Ri_NIR_0_q",
        "Ri_PAR_f",
        "Ri_NIR_f",
        "Ri_PAR_q",
        "Ri_NIR_q",
        "Ra_PAR_0_f",
        "Ra_NIR_0_f",
        "Ra_PAR_0_q",
        "Ra_NIR_0_q",
        "Ra_PAR_f",
        "Ra_NIR_f",
        "Ra_PAR_q",
        "Ra_NIR_q",
    ]
    rows_indep_light_1 = _stack_component_rows(
        scene_indep_light,
        series_indep_light_1,
        cfg_indep_light_1,
        meteo_indep_light_1.rows,
        cols_indep_light,
    )
    rows_indep_light_2_raw = _stack_component_rows(
        scene_indep_light,
        series_indep_light_2,
        cfg_indep_light_2,
        meteo_indep_light_2.rows,
        cols_indep_light,
    )
    rows_indep_light_2 = Dict{String,Any}[]
    for r in rows_indep_light_2_raw
        step_no = Int(r["step_number"])
        step_no == 0 && continue
        rr = copy(r)
        rr["step_number"] = step_no - 1
        push!(rows_indep_light_2, rr)
    end
    sort!(rows_indep_light_2; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
    _assert_rows_match(
        rows_indep_light_1,
        rows_indep_light_2,
        ["step_number", "item_id", "component_id"],
        cols_indep_light;
        atol=1e-9,
        rtol=1e-9,
        label="test-independant-steps light component rows",
    )

        _stage_done("Meteo and independent-step parity", _stage_t0)
    end

    @testset "Absorption caching timestep parity" begin
        _stage_t0 = _stage_start("Absorption caching timestep parity")
    # Java absorb tests parity: within a type/domain, absorbed over incident ratios are constant.
    function _check_absorption_ratios(cfg_path::String; with_scattering::Bool)
        cfg = ArchimedLight.read_light_config(cfg_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        suffix = with_scattering ? "" : "_0"
        cols = [
            "type",
            "area",
            "step_duration",
            "Ri_PAR$(suffix)_f",
            "Ri_NIR$(suffix)_f",
            "Ra_PAR$(suffix)_f",
            "Ra_NIR$(suffix)_f",
            "Ra_PAR$(suffix)_q",
            "Ra_NIR$(suffix)_q",
        ]
        rows = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, columns=cols).rows
        for tname in ("Metamer", "Leaf"), band in ("PAR", "NIR")
            ri_f = "Ri_$(band)$(suffix)_f"
            ra_f = "Ra_$(band)$(suffix)_f"
            ra_q = "Ra_$(band)$(suffix)_q"

            vals_irr = Float64[]
            vals_energy = Float64[]
            for r in rows
                string(r["type"]) == tname || continue
                ri = Float64(r[ri_f])
                ri == 0.0 && continue
                area = Float64(r["area"])
                dt = Float64(r["step_duration"])
                push!(vals_irr, Float64(r[ra_f]) / ri)
                push!(vals_energy, Float64(r[ra_q]) / (ri * max(dt * area, eps(Float64))))
            end
            @test !isempty(vals_irr)
            ref_irr = first(vals_irr)
            ref_energy = first(vals_energy)
            @test maximum(abs(v - ref_irr) for v in vals_irr; init=0.0) < parity_limit(:absorb_cached_abs)
            @test maximum(abs(v - ref_energy) for v in vals_energy; init=0.0) < parity_limit(:absorb_cached_abs)
        end
    end
    _check_absorption_ratios(joinpath(test_root, "test-absorb", "config.yml"); with_scattering=false)
    _check_absorption_ratios(joinpath(test_root, "test-absorb2", "config.yml"); with_scattering=true)

    # Java test-cached-radiation parity: cached and uncached runs match across irradiance and absorbed outputs.
    cached_root = joinpath(test_root, "test-cached-radiation")
    cfg_cached_1 = ArchimedLight.read_light_config(joinpath(cached_root, "config.yml"))
    cfg_cached_2 = ArchimedLight.read_light_config(joinpath(cached_root, "config2.yml"))
    scene_cached_1 = ArchimedLight.read_scene(cfg_cached_1.scene)
    scene_cached_2 = ArchimedLight.read_scene(cfg_cached_2.scene)
    meteo_cached_1 = ArchimedLight.read_meteo(cfg_cached_1.meteo)
    meteo_cached_2 = ArchimedLight.read_meteo(cfg_cached_2.meteo)
    series_cached_1 = ArchimedLight.run_light_series(scene_cached_1, meteo_cached_1, cfg_cached_1)
    series_cached_2 = ArchimedLight.run_light_series(scene_cached_2, meteo_cached_2, cfg_cached_2)
    @test length(series_cached_1) == length(series_cached_2)

    function _cell_close(a, b; atol=1e-6, rtol=1e-3)
        sa = strip(string(a))
        sb = strip(string(b))
        pa = tryparse(Float64, sa)
        pb = tryparse(Float64, sb)
        if pa !== nothing && pb !== nothing
            return isapprox(pa, pb; atol=atol, rtol=rtol)
        end
        return sa == sb
    end

    cached_cols = [
        "item_id",
        "component_id",
        "Ri_PAR_0_f",
        "Ri_NIR_0_f",
        "Ri_PAR_f",
        "Ri_NIR_f",
        "Ra_PAR_0_f",
        "Ra_NIR_0_f",
        "Ra_PAR_0_q",
        "Ra_NIR_0_q",
        "Ra_PAR_f",
        "Ra_NIR_f",
        "Ra_PAR_q",
        "Ra_NIR_q",
    ]
    for i in eachindex(series_cached_1)
        rows1 = ArchimedLight.component_values_table(
            scene_cached_1,
            series_cached_1[i],
            cfg_cached_1;
            meteo_row=meteo_cached_1.rows[i],
            columns=cached_cols,
        ).rows
        rows2 = ArchimedLight.component_values_table(
            scene_cached_2,
            series_cached_2[i],
            cfg_cached_2;
            meteo_row=meteo_cached_2.rows[i],
            columns=cached_cols,
        ).rows
        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => r for r in rows1)
        d2 = Dict((Int(r["item_id"]), Int(r["component_id"])) => r for r in rows2)
        @test Set(keys(d1)) == Set(keys(d2))
        for k in keys(d1)
            for c in cached_cols[3:end]
                @test _cell_close(d1[k][c], d2[k][c])
            end
        end
    end

    # Java test-timestep{1,2,3} parity: one large meteo step and many small steps
    # must yield comparable average RI_SW_f.
    function _run_timestep_fixture(root_dir::String)
        cfg_one = ArchimedLight.read_light_config(joinpath(root_dir, "onestep", "config.yml"))
        cfg_many = ArchimedLight.read_light_config(joinpath(root_dir, "manysteps", "config.yml"))
        scene_one = ArchimedLight.read_scene(cfg_one.scene)
        scene_many = ArchimedLight.read_scene(cfg_many.scene)
        meteo_one = ArchimedLight.read_meteo(cfg_one.meteo)
        meteo_many = ArchimedLight.read_meteo(cfg_many.meteo)
        series_one = ArchimedLight.run_light_series(scene_one, meteo_one, cfg_one)
        series_many = ArchimedLight.run_light_series(scene_many, meteo_many, cfg_many)
        sv_one = ArchimedLight.scene_values_table(scene_one, series_one, cfg_one; meteo_rows=meteo_one.rows, columns=["RI_SW_f"]).rows
        sv_many = ArchimedLight.scene_values_table(scene_many, series_many, cfg_many; meteo_rows=meteo_many.rows, columns=["RI_SW_f"]).rows
        @test length(sv_one) == 1
        @test !isempty(sv_many)
        val_one = Float64(sv_one[1]["RI_SW_f"])
        val_many = mean(Float64(r["RI_SW_f"]) for r in sv_many)
        rel = relerr(val_many, val_one)
        return val_one, val_many, rel
    end

    _, _, rel_ts1 = _run_timestep_fixture(joinpath(test_root, "test-timestep1"))
    _, _, rel_ts2 = _run_timestep_fixture(joinpath(test_root, "test-timestep2"))
    _, _, rel_ts3 = _run_timestep_fixture(joinpath(test_root, "test-timestep3"))
    @test rel_ts1 < parity_limit(:timestep_rel_ts1)
    @test rel_ts2 < parity_limit(:timestep_rel_ts2)
    @test rel_ts3 < parity_limit(:timestep_rel_ts3)

    for disk_name in ("test-save_on_disk1", "test-save_on_disk2", "test-save_on_disk3", "test-save_on_disk4", "test-save_on_disk5")
        f_disk = fixtures[disk_name]
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
        if cfg_disk_on.scattering
            @test step_disk_2.scattering !== nothing
            @test step_mem.scattering !== nothing
            @test max_abs_float_dict_diff(step_disk_2.scattering.added_par_power_per_node, step_mem.scattering.added_par_power_per_node) == 0.0
            @test max_abs_float_dict_diff(step_disk_2.scattering.added_nir_power_per_node, step_mem.scattering.added_nir_power_per_node) == 0.0
        end
        rm(cache_dir; recursive=true, force=true)
    end

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
        _stage_done("Absorption caching timestep parity", _stage_t0)
    end
    end
end

end
