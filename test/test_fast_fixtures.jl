const _FAST_FIXTURE_CASE_FILTER = Set(
    filter(!isempty, strip.(lowercase.(split(get(ENV, "ARCHIMEDLIGHT_FAST_FIXTURE_CASE", ""), ",")))),
)

function _fast_case_enabled(name::String)
    isempty(_FAST_FIXTURE_CASE_FILTER) && return true
    return lowercase(name) in _FAST_FIXTURE_CASE_FILTER
end

function _render_ri_par_f_figure(scene, step, cfg; title::String)
    vertices, faces, face2node, _, _, _ = ArchimedLight._scene_geometry_for_interception(scene, cfg)
    metric = step.budget.ri_par_f_per_node

    v_sum = zeros(Float64, length(vertices))
    v_count = zeros(Int, length(vertices))
    for i in eachindex(faces)
        f = faces[i]
        v = get(metric, face2node[i], NaN)
        isfinite(v) || continue
        for vid in (Int(f[1]), Int(f[2]), Int(f[3]))
            v_sum[vid] += v
            v_count[vid] += 1
        end
    end
    vertex_values = Float64[v_count[i] > 0 ? (v_sum[i] / v_count[i]) : NaN for i in eachindex(vertices)]

    colorrange = (0.0, max(step.sky.ri_par_f, eps(Float64)))
    fig = Figure(size=(960, 720))
    ax = Axis3(
        fig[1, 1];
        title=title,
        aspect=:data,
        azimuth=1.45,
        elevation=0.35,
        perspectiveness=0.65,
    )
    p = mesh!(
        ax,
        vertices,
        faces;
        color=vertex_values,
        colormap=:viridis,
        colorrange=colorrange,
        nan_color=:lightgray,
        shading=false,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    Colorbar(fig[1, 2], p, label="Ri_PAR_f")
    return fig
end

@testset "Fast fixtures (manual, config-driven)" begin
    @testset "Sky only" begin
        if _fast_case_enabled("sky_06_direct")
            case_root = joinpath(@__DIR__, "fast_fixtures", "sky_06_direct", "input")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
            meteo = ArchimedLight.read_meteo(joinpath(case_root, "meteo.csv"))
            row = first(ArchimedLight.prepare_meteo(meteo, cfg).rows)
            sky = ArchimedLight.compute_sky(row, cfg)
            turtle = ArchimedLight.build_turtle(cfg, sky)
            flux = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)

            @test cfg.general.turtle_sectors == 6
            @test cfg.general.all_in_turtle == false
            @test length(turtle.sectors) == 7
            @test count(s -> s.source == :sun, turtle.sectors) == 1
            @test all(v -> v >= 0.0, flux.par)
            @test all(v -> v >= 0.0, flux.nir)
            @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
            @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
        end

        if _fast_case_enabled("sky_16_turtle")
            case_root = joinpath(@__DIR__, "fast_fixtures", "sky_16_turtle", "input")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
            meteo = ArchimedLight.read_meteo(joinpath(case_root, "meteo.csv"))
            row = first(ArchimedLight.prepare_meteo(meteo, cfg).rows)
            sky = ArchimedLight.compute_sky(row, cfg)
            turtle = ArchimedLight.build_turtle(cfg, sky)
            flux = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)

            @test cfg.general.turtle_sectors == 16
            @test cfg.general.all_in_turtle == true
            @test length(turtle.sectors) == 16
            @test count(s -> s.source == :sun, turtle.sectors) == 0
            @test all(v -> v >= 0.0, flux.par)
            @test all(v -> v >= 0.0, flux.nir)
            @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
            @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
        end

        if _fast_case_enabled("sky_46_direct")
            case_root = joinpath(@__DIR__, "fast_fixtures", "sky_46_direct", "input")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
            meteo = ArchimedLight.read_meteo(joinpath(case_root, "meteo.csv"))
            row = first(ArchimedLight.prepare_meteo(meteo, cfg).rows)
            sky = ArchimedLight.compute_sky(row, cfg)
            turtle = ArchimedLight.build_turtle(cfg, sky)
            flux = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)

            @test cfg.general.turtle_sectors == 46
            @test cfg.general.all_in_turtle == false
            @test length(turtle.sectors) == 47
            @test count(s -> s.source == :sun, turtle.sectors) == 1
            @test all(v -> v >= 0.0, flux.par)
            @test all(v -> v >= 0.0, flux.nir)
            @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
            @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
        end
    end

    @testset "Simple plant (numeric + visual refs)" begin
        if _fast_case_enabled("simpleplant_16_notoric")
            case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "input", "config.yml"))
            scene = ArchimedLight.read_scene(cfg.paths.scene)
            meteo = ArchimedLight.read_meteo(cfg.paths.meteo)
            selected = ArchimedLight.prepare_meteo(meteo, cfg)
            series = ArchimedLight.run_light_series(scene, meteo, cfg)
            step = first(series)
            meteo_row = first(selected.rows)

                tmpdir = mktempdir()
                observed_csv = joinpath(tmpdir, "component_values.csv")
                ArchimedLight.write_component_values_csv(
                    observed_csv,
                    scene,
                    step,
                    cfg;
                    meteo_row=meteo_row,
                    step_number=0,
                    columns=["step_number", "source_topology_id", "object_id", "area", "Ri_PAR_0_q"],
                    strict=false,
                )
                expected_csv = joinpath(case_root, "expected", "component_values.csv")
                expected_rows = collect(Tables.rowtable(CSV.File(expected_csv; delim=';', normalizenames=false)))
                observed_rows = collect(Tables.rowtable(CSV.File(observed_csv; delim=';', normalizenames=false)))

                @test length(observed_rows) == length(expected_rows)
                key_notoric(row) = (Int(row.step_number), Int(row.source_topology_id))
                expected_map = Dict(key_notoric(r) => r for r in expected_rows)
                observed_map = Dict(key_notoric(r) => r for r in observed_rows)
                @test Set(keys(observed_map)) == Set(keys(expected_map))
                for k in keys(expected_map)
                    er = expected_map[k]
                    or = observed_map[k]
                    @test isapprox(Float64(or.area), Float64(er.area); atol=1e-8, rtol=1e-6)
                    @test isapprox(Float64(or.Ri_PAR_0_q), Float64(er.Ri_PAR_0_q); atol=1e-4, rtol=1e-6)
                end

                fig = _render_ri_par_f_figure(scene, step, cfg; title="simpleplant_16_notoric | Ri_PAR_f")
                ref_png = joinpath(case_root, "expected", "ri_par_f_step0.png")
                @test_reference relpath(ref_png, @__DIR__) fig by = ReferenceTests.psnr_equality(35)
        end

        if _fast_case_enabled("simpleplant_16_toric")
            case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_toric")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "input", "config.yml"))
            scene = ArchimedLight.read_scene(cfg.paths.scene)
            meteo = ArchimedLight.read_meteo(cfg.paths.meteo)
            selected = ArchimedLight.prepare_meteo(meteo, cfg)
            series = ArchimedLight.run_light_series(scene, meteo, cfg)
            step = first(series)
            meteo_row = first(selected.rows)

                tmpdir = mktempdir()
                observed_csv = joinpath(tmpdir, "component_values.csv")
                ArchimedLight.write_component_values_csv(
                    observed_csv,
                    scene,
                    step,
                    cfg;
                    meteo_row=meteo_row,
                    step_number=0,
                    columns=["step_number", "source_topology_id", "object_id", "area", "Ri_PAR_0_q"],
                    strict=false,
                )
                expected_csv = joinpath(case_root, "expected", "component_values.csv")
                expected_rows = collect(Tables.rowtable(CSV.File(expected_csv; delim=';', normalizenames=false)))
                observed_rows = collect(Tables.rowtable(CSV.File(observed_csv; delim=';', normalizenames=false)))

                @test length(observed_rows) == length(expected_rows)
                key_toric(row) = (Int(row.step_number), Int(row.source_topology_id))
                expected_map = Dict(key_toric(r) => r for r in expected_rows)
                observed_map = Dict(key_toric(r) => r for r in observed_rows)
                @test Set(keys(observed_map)) == Set(keys(expected_map))
                for k in keys(expected_map)
                    er = expected_map[k]
                    or = observed_map[k]
                    @test isapprox(Float64(or.area), Float64(er.area); atol=1e-8, rtol=1e-6)
                    @test isapprox(Float64(or.Ri_PAR_0_q), Float64(er.Ri_PAR_0_q); atol=1e-4, rtol=1e-6)
                end

                fig = _render_ri_par_f_figure(scene, step, cfg; title="simpleplant_16_toric | Ri_PAR_f")
                ref_png = joinpath(case_root, "expected", "ri_par_f_step0.png")
                @test_reference relpath(ref_png, @__DIR__) fig by = ReferenceTests.psnr_equality(35)
        end
    end
end
