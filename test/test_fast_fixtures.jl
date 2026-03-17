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

            @test get(cfg.general, "sky_sectors", 46) == 6
            @test get(cfg.general, "all_in_turtle", false) == false
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

            @test get(cfg.general, "sky_sectors", 46) == 16
            @test get(cfg.general, "all_in_turtle", false) == true
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

            @test get(cfg.general, "sky_sectors", 46) == 46
            @test get(cfg.general, "all_in_turtle", false) == false
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
            scene = ArchimedLight.read_scene(cfg.source_files.scene)
            meteo = ArchimedLight.read_meteo(cfg.source_files.meteo)
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
                    columns=["step_number", "item_id", "component_id", "area", "Ri_PAR_0_q"],
                    strict=false,
                )
                observed_rows = collect(Tables.rowtable(CSV.File(observed_csv; delim=';', normalizenames=false)))

                @test !isempty(observed_rows)
                @test all(Float64(row.area) > 0.0 for row in observed_rows)
                @test any(Float64(row.Ri_PAR_0_q) > 0.0 for row in observed_rows)
                @test sum(Float64(row.Ri_PAR_0_q) for row in observed_rows) > 0.0
                @test _render_ri_par_f_figure(scene, step, cfg; title="simpleplant_16_notoric | Ri_PAR_f") isa Figure
        end

        if _fast_case_enabled("simpleplant_16_toric")
            case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_toric")
            cfg = ArchimedLight.read_light_config(joinpath(case_root, "input", "config.yml"))
            scene = ArchimedLight.read_scene(cfg.source_files.scene)
            meteo = ArchimedLight.read_meteo(cfg.source_files.meteo)
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
                    columns=["step_number", "item_id", "component_id", "area", "Ri_PAR_0_q"],
                    strict=false,
                )
                observed_rows = collect(Tables.rowtable(CSV.File(observed_csv; delim=';', normalizenames=false)))

                @test !isempty(observed_rows)
                @test all(Float64(row.area) > 0.0 for row in observed_rows)
                @test any(Float64(row.Ri_PAR_0_q) > 0.0 for row in observed_rows)
                @test sum(Float64(row.Ri_PAR_0_q) for row in observed_rows) > 0.0
                @test _render_ri_par_f_figure(scene, step, cfg; title="simpleplant_16_toric | Ri_PAR_f") isa Figure
        end
    end
end
