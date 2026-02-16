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
    scene_par_0 = sum(values(step_custom.budget.ri_par_0_q_per_node)) * step_seconds
    scene_par_n = sum(values(step_custom.budget.ri_par_q_per_node)) * step_seconds
    scene_custom_0 = sum(values(step_custom.budget.extra_0_q_per_band["CUSTOM"])) * step_seconds
    scene_custom_n = sum(values(step_custom.budget.extra_q_per_band["CUSTOM"])) * step_seconds

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
end
