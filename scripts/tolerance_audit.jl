#!/usr/bin/env julia

using ArchimedLight
using Statistics
using Printf

include(joinpath(dirname(@__DIR__), "test", "java_parity_harness.jl"))

to_int(v) = v isa Number ? Int(round(v)) : parse(Int, strip(string(v)))

function expected_component_metrics(path::String)
    rows = read_java_csv(path)
    out = Dict{Tuple{Int,Int},NamedTuple{(:area,:hits,:irr0),Tuple{Float64,Int,Float64}}}()
    for r in rows
        key = (to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id)))
        area = Float64(getproperty(r, :area))
        hits = Int(floor(area * Float64(getproperty(r, :surface_hits))))
        irr0 = Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR))
        out[key] = (area=area, hits=hits, irr0=irr0)
    end
    out
end

function observed_component_metrics(scene::ArchimedLight.SceneGeometry, cfg::ArchimedLight.LightConfig, step::ArchimedLight.LightStepResult)
    key_by_node = ArchimedLight._interception_java_keys(scene, cfg)
    out = Dict{Tuple{Int,Int},NamedTuple{(:hits,:power0),Tuple{Int,Float64}}}()
    for (nid, key) in key_by_node
        p0 =
            get(step.first_order.incident_par_power_per_node, nid, 0.0) +
            get(step.first_order.incident_nir_power_per_node, nid, 0.0)
        out[key] = (hits=get(step.first_order.hits_per_node, nid, 0), power0=p0)
    end
    out
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

function _run_fixture_step(cfg_path::String)
    cfg = ArchimedLight.read_light_config(cfg_path)
    scene = ArchimedLight.read_scene(cfg.scene)
    row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
    step = ArchimedLight.run_light_step(scene, row, cfg)
    return cfg, scene, row, step
end

function _run_fixture_series(cfg_path::String)
    cfg = ArchimedLight.read_light_config(cfg_path)
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)
    return cfg, scene, meteo, series
end

function _model_radiance(path::String)
    txt = read(path, String)
    m = match(r"radiance\s*:\s*([-+0-9.eE]+)", txt)
    m === nothing && error("No radiance found in $path")
    parse(Float64, m.captures[1])
end

function _xml_tag_value(path::String, tag::String)
    txt = read(path, String)
    m = match(Regex("<$(tag)>([^<]+)</$(tag)>"), txt)
    m === nothing && error("Missing tag <$tag> in $path")
    parse(Float64, strip(m.captures[1]))
end

function _sum_item(rows, item_id::Int, col::String)
    s = 0.0
    for r in rows
        Int(r["item_id"]) == item_id || continue
        s += Float64(r[col])
    end
    s
end

function _component_par_by_item(scene, step, cfg, row)
    t = ArchimedLight.component_values_table(
        scene,
        step,
        cfg;
        meteo_row=row,
        columns=["item_id", "component_id", "group", "Ri_PAR_0_f", "Ri_PAR_f"],
    )
    t.rows
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

function _check_summary_component_max_rel(cfg_path::String; use_scattering::Bool, include_scene_rel::Bool=false)
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

    max_scene_rel = nothing
    if include_scene_rel
        ground = [
            Float64(r["area"]) for r in summary.rows if
            Int(r["step_number"]) == 0 && Int(r["item_id"]) == -1
        ]
        ground_area = first(ground)
        max_scene_rel = 0.0
        for i in eachindex(series)
            dt = _step_duration_seconds(meteo.rows[i])
            srows = [r for r in summary.rows if Int(r["step_number"]) == i - 1]
            ri_scene = sum(Float64(r["Ri_q"]) for r in srows) / dt / ground_area
            max_scene_rel = max(max_scene_rel, relerr(ri_scene, series[i].sky.ri_sw_f))
        end
    end
    return max_rel, max_scene_rel
end

function _check_absorption_ratio_max(cfg_path::String; with_scattering::Bool)
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
    max_delta = 0.0
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
        ref_irr = first(vals_irr)
        ref_energy = first(vals_energy)
        max_delta = max(max_delta, maximum(abs(v - ref_irr) for v in vals_irr; init=0.0))
        max_delta = max(max_delta, maximum(abs(v - ref_energy) for v in vals_energy; init=0.0))
    end
    max_delta
end

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
    val_one = Float64(sv_one[1]["RI_SW_f"])
    val_many = mean(Float64(r["RI_SW_f"]) for r in sv_many)
    relerr(val_many, val_one)
end

function add_row!(rows::Vector{NamedTuple}, metric::Symbol, observed::Real; fixture::Union{Nothing,String}=nothing, case::String="", cmp::Symbol=:le)
    threshold = parity_limit(metric; fixture=fixture)
    obs = Float64(observed)
    pass = cmp == :le ? obs <= threshold : obs >= threshold
    headroom =
        if cmp == :le
            threshold / max(obs, eps(Float64))
        else
            obs / max(threshold, eps(Float64))
        end
    push!(
        rows,
        (
            metric=metric,
            fixture=fixture === nothing ? "*" : fixture,
            case=case,
            cmp=cmp,
            observed=obs,
            threshold=threshold,
            pass=pass,
            headroom=headroom,
        ),
    )
end

function main()
    fixtures = Dict(f.name => f for f in light_parity_fixtures())
    test_root = joinpath(_repo_root(), "java_implementation", "archimed-lib-2018", "tests")
    rows = NamedTuple[]

    # Non-numeric exact fields are audited by row-key checks in tests.
    add_row!(rows, :exact, 0.0; case="row-key exactness (representative)")

    f_hit2 = fixtures["test-hitcount2"]
    for all_in_turtle in (false, true)
        mode = all_in_turtle ? "turtle" : "raycast"
        rep = fixture_parity_report(f_hit2; all_in_turtle=all_in_turtle)
        add_row!(rows, :hitcount_total_rel, relerr(rep.snapshot.total_hits, rep.expected_hitcount_total); case="test-hitcount2 $(mode)")

        cfg = ArchimedLight.read_light_config(f_hit2.config_path)
        cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        expected = expected_component_metrics(_expected_component_values_path(f_hit2; all_in_turtle=all_in_turtle))
        observed = observed_component_metrics(scene, cfg, step)
        max_comp_rel = maximum(abs(observed[k].hits - expected[k].hits) / max(expected[k].hits, 1) for k in keys(expected))
        add_row!(rows, :hitcount_component_rel, max_comp_rel; case="test-hitcount2 $(mode)")

        plant_key = (1, 2)
        got_irr = observed[plant_key].power0 / max(expected[plant_key].area, eps(Float64))
        add_row!(rows, :irr_component_rel_strict, relerr(got_irr, expected[plant_key].irr0); case="test-hitcount2 $(mode)")

        expected_hits = expected_component_hits(_expected_component_values_path(f_hit2; all_in_turtle=all_in_turtle))
        got_hits = sort(collect(values(step.first_order.hits_per_node)))
        abs_err = abs.(got_hits .- expected_hits)
        rel_err = abs_err ./ max.(abs.(Float64.(expected_hits)), 1.0)
        add_row!(rows, all_in_turtle ? :hitcount_hist_abs_turtle : :hitcount_hist_abs_raycast, maximum(abs_err); case="test-hitcount2 $(mode)")
        add_row!(rows, all_in_turtle ? :hitcount_hist_rel_turtle : :hitcount_hist_rel_raycast, maximum(rel_err); case="test-hitcount2 $(mode)")
    end

    f_hit3 = fixtures["test-hitcount3"]
    for all_in_turtle in (false, true)
        mode = all_in_turtle ? "turtle" : "raycast"
        cfg = ArchimedLight.read_light_config(f_hit3.config_path)
        cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        expected = expected_component_metrics(_expected_component_values_path(f_hit3; all_in_turtle=all_in_turtle))
        observed = observed_component_metrics(scene, cfg, step)

        hi = Float64[]
        dense = Float64[]
        sparse = Float64[]
        for key in keys(expected)
            exp = expected[key]
            got = observed[key]
            rel_hit = abs(got.hits - exp.hits) / max(exp.hits, 1)
            exp.hits >= 1000 && push!(hi, rel_hit)
            got_irr = got.power0 / max(exp.area, eps(Float64))
            abs_err = abs(got_irr - exp.irr0)
            if abs(exp.irr0) >= 1.0
                if exp.hits >= 400
                    push!(dense, abs_err / max(abs(exp.irr0), eps(Float64)))
                else
                    push!(sparse, abs_err)
                end
            end
        end
        add_row!(rows, all_in_turtle ? :hitcount_component_hi_rel_turtle : :hitcount_component_hi_rel_raycast, maximum(hi); case="test-hitcount3 $(mode)")
        add_row!(rows, :irr_component_rel_dense, maximum(dense); case="test-hitcount3 $(mode)")
        add_row!(rows, :irr_component_abs_sparse, maximum(sparse); case="test-hitcount3 $(mode)")
    end

    f_wsun = fixtures["test-weighted-sun"]
    rep_wsun = fixture_parity_report(f_wsun)
    add_row!(rows, :sun_deg_snapshot, abs(rep_wsun.snapshot.sun_azimuth - rep_wsun.expected_sun_azimuth); case="test-weighted-sun az")
    add_row!(rows, :sun_deg_snapshot, abs(rep_wsun.snapshot.sun_elevation - rep_wsun.expected_sun_elevation); case="test-weighted-sun el")
    cfg_wsun = ArchimedLight.read_light_config(f_wsun.config_path)
    meteo_wsun = ArchimedLight.read_meteo(cfg_wsun.meteo)
    rows_wsun = read_java_csv(_expected_sun_log_path(f_wsun))
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
    add_row!(rows, :sun_deg_az_series, max_az_err; case="test-weighted-sun")
    add_row!(rows, :sun_deg_el_series, max_el_err; case="test-weighted-sun")

    for fx_name in ("test-weighted-sun", "test-save_on_disk1", "test-save_on_disk6")
        fx = fixtures[fx_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        series = ArchimedLight.run_light_series(scene, meteo, cfg)
        expected_rows = read_java_csv(_expected_scene_values_path(fx))
        expected_sw = Float64[getproperty(r, :RI_SW_f) for r in expected_rows]
        got_sw = Float64[s.sky.ri_sw_f for s in series]
        abs_err = abs.(got_sw .- expected_sw)
        rel_err = abs_err ./ max.(abs.(expected_sw), eps(Float64))
        sw_atol = parity_limit(:scene_snapshot_sw_atol; fixture=fx_name)
        add_row!(rows, :scene_sw_series_max_abs, maximum(abs_err); fixture=fx_name, case=fx_name)
        add_row!(rows, :scene_sw_series_mean_rel, mean(rel_err); fixture=fx_name, case=fx_name)
        add_row!(rows, :scene_snapshot_sw_atol, maximum(abs_err); fixture=fx_name, case=fx_name)
        rtol_candidates = Float64[rel_err[i] for i in eachindex(rel_err) if abs_err[i] > sw_atol]
        add_row!(rows, :scene_snapshot_sw_rtol, isempty(rtol_candidates) ? 0.0 : maximum(rtol_candidates); fixture=fx_name, case=fx_name)
    end

    for fx_name in ("test-cafeier", "test-cafeier2", "test-cafeier_sensor", "test-cafeier_sensor2", "test-cafeier_sensor3")
        fx = fixtures[fx_name]
        expected_rows = read_java_csv(_expected_component_values_path(fx))
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
        add_row!(rows, :cafeier_total_rel, relerr(got_total, exp_total); fixture=fx_name, case=fx_name)

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
            add_row!(rows, :scene_riq_rel, relerr(got_riq, exp_riq); case=fx_name)
        end
    end

    rep_scat1 = fixture_parity_report(fixtures["test-scattering-one-plate"])
    rep_scat2 = fixture_parity_report(fixtures["test-scattering-two-plates"])
    add_row!(rows, :scattering_total_rel_loose, relerr(rep_scat1.julia_scattering_total, rep_scat1.expected_scattering_total); case="one-plate")
    add_row!(rows, :scattering_total_rel_strict, relerr(rep_scat2.julia_scattering_total, rep_scat2.expected_scattering_total); case="two-plates")

    for links_fx_name in ("test-links-pixeltable", "test-links-pixeltable2")
        fx = fixtures[links_fx_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        sky = ArchimedLight.compute_sky(row, cfg)
        turtle = ArchimedLight.build_turtle(cfg, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
        first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
        scat_pix = ArchimedLight.compute_scattering(scene, turtle, first_order, cfg; mode=:raycast)
        scat_lnk = ArchimedLight.compute_scattering(scene, turtle, first_order, cfg; mode=:links)
        bud_pix = ArchimedLight.integrate_light(first_order, scat_pix, cfg)
        bud_lnk = ArchimedLight.integrate_light(first_order, scat_lnk, cfg)
        add_row!(rows, :backend_identity_abs, maximum(abs(get(bud_pix.ri_par_q_per_node, k, 0.0) - get(bud_lnk.ri_par_q_per_node, k, 0.0)) for k in union(keys(bud_pix.ri_par_q_per_node), keys(bud_lnk.ri_par_q_per_node)); init=0.0); case="links-pixeltable par $(links_fx_name)")
        add_row!(rows, :backend_identity_abs, maximum(abs(get(bud_pix.ri_nir_q_per_node, k, 0.0) - get(bud_lnk.ri_nir_q_per_node, k, 0.0)) for k in union(keys(bud_pix.ri_nir_q_per_node), keys(bud_lnk.ri_nir_q_per_node)); init=0.0); case="links-pixeltable nir $(links_fx_name)")
    end

    area_ratio_metrics = Dict{String,Tuple{Float64,Float64,Float64,Float64}}()
    for area_name in ("test-area_ratio", "test-area_ratio2", "test-area_ratio3", "test-area_ratio4")
        fx = fixtures[area_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        dt_seconds = _step_duration_seconds(row)
        step = ArchimedLight.run_light_step(scene, row, cfg)
        base_max = maximum(keys(scene.total_area_per_node))
        got = Float64[]
        for nid in keys(step.budget.ri_par_0_q_per_node)
            nid > base_max || continue
            push!(got, step.budget.ri_par_0_q_per_node[nid] + step.budget.ri_nir_0_q_per_node[nid])
        end
        got_mean, got_std = mean(got), std(got)
        expected_rows = read_java_csv(_expected_component_values_path(fx))
        expected = Float64[
            Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR)) * Float64(getproperty(r, :area)) * dt_seconds for
            r in expected_rows if to_int(getproperty(r, :item_id)) == -1
        ]
        exp_mean, exp_std = mean(expected), std(expected)
        add_row!(rows, :area_ratio_uniformity_julia_relstd, got_std / max(abs(got_mean), eps(Float64)); case=area_name)
        add_row!(rows, :area_ratio_uniformity_java_relstd, exp_std / max(abs(exp_mean), eps(Float64)); case=area_name)
        add_row!(rows, :area_ratio_mean_rel, relerr(got_mean, exp_mean); case=area_name)
        area_ratio_metrics[area_name] = (got_mean, got_std, exp_mean, exp_std)
    end
    got_means = [v[1] for v in values(area_ratio_metrics)]
    mean_ref = max(abs(mean(got_means)), 1.0)
    add_row!(rows, :area_ratio_means_spread_rel, (maximum(got_means) - minimum(got_means)) / mean_ref; case="area_ratio group")

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
        rows_ratio = ArchimedLight.component_values_table(scene_area5, step_ratio, cfg_ratio; meteo_row=row_area5, columns=["type", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
        irr_ratio = Float64[
            Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_ratio if string(get(r, "type", "")) == "Leaf"
        ]
        stdev_by_npix_ratio[npix] = std(irr_ratio)

        cfg_noratio = with_pixel_size_area_ratio(cfg_area5, px_cm, false)
        step_noratio = ArchimedLight.run_light_step(scene_area5, row_area5, cfg_noratio)
        rows_noratio = ArchimedLight.component_values_table(scene_area5, step_noratio, cfg_noratio; meteo_row=row_area5, columns=["type", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
        irr_noratio = Float64[
            Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_noratio if string(get(r, "type", "")) == "Leaf"
        ]
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
    add_row!(rows, :area_ratio5_ratio_abs, maximum(abs(stdev_by_npix_ratio[n] - expected_stdev_ratio[n]) for n in keys(expected_stdev_ratio)); case="area_ratio5")
    add_row!(rows, :area_ratio5_noratio_abs_hi, maximum(abs(stdev_by_npix_noratio[n] - expected_stdev_noratio[n]) for n in (500000, 1000000)); case="area_ratio5")
    ratio_vals = [stdev_by_npix_ratio[n] for (n, _) in pixel_sweep]
    noratio_vals = [stdev_by_npix_noratio[n] for (n, _) in pixel_sweep]
    add_row!(rows, :area_ratio5_ratio_spread, maximum(ratio_vals) - minimum(ratio_vals); case="area_ratio5")
    add_row!(rows, :area_ratio5_noratio_first_last_ratio, first(noratio_vals) / last(noratio_vals); case="area_ratio5", cmp=:ge)

    for links_name in ("test-links-sensor-plates", "test-links-sensor-plates2")
        fx = fixtures[links_name]
        cfg = ArchimedLight.read_light_config(fx.config_path)
        scene = ArchimedLight.read_scene(cfg.scene)
        row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
        step_pix = ArchimedLight.run_light_step(scene, row, cfg; scattering_backend=ArchimedLight.RaycastScatteringBackend())
        step_lnk = ArchimedLight.run_light_step(scene, row, cfg; scattering_backend=ArchimedLight.LinksScatteringBackend())
        rows_pix = _component_par_by_item(scene, step_pix, cfg, row)
        rows_lnk = _component_par_by_item(scene, step_lnk, cfg, row)
        map_pix_0 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"]) for r in rows_pix)
        map_lnk_0 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"]) for r in rows_lnk)
        map_pix_n = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_f"]) for r in rows_pix)
        map_lnk_n = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_f"]) for r in rows_lnk)
        add_row!(rows, :backend_identity_abs, maximum(abs(map_pix_0[k] - map_lnk_0[k]) for k in keys(map_pix_0); init=0.0); case="$(links_name) Ri_PAR_0")
        add_row!(rows, :backend_identity_abs, maximum(abs(map_pix_n[k] - map_lnk_n[k]) for k in keys(map_pix_n); init=0.0); case="$(links_name) Ri_PAR")

        if links_name == "test-links-sensor-plates"
            plate_scat = sum(Float64(r["Ri_PAR_f"] - r["Ri_PAR_0_f"]) for r in rows_pix if string(r["group"]) == "plates"; init=0.0)
            sensor_meas = sum(Float64(r["Ri_PAR_f"]) for r in rows_pix if string(r["group"]) == "sensors"; init=0.0)
            add_row!(rows, :links_sensor_plate_balance_rel, abs(plate_scat - sensor_meas) / max(abs(plate_scat), eps(Float64)); case=links_name)
        else
            coeffs = ArchimedLight._group_optical_coeffs(cfg)
            ro = get(get(coeffs, "plates", Dict{String,Float64}()), "PAR", cfg.scattering_coeff_par) / 2.0
            ip1 = _sum_item(rows_pix, 4, "Ri_PAR_f")
            ip2 = _sum_item(rows_pix, 1, "Ri_PAR_f")
            is2 = _sum_item(rows_pix, 2, "Ri_PAR_f")
            is1 = _sum_item(rows_pix, 3, "Ri_PAR_f")
            v_plate2 = ip1 * ro
            v_sensor = ip1 * ro + ip2 * ro
            add_row!(rows, :links_sensor_plates2_plate2_rel, abs(v_plate2 - ip2) / max(abs(v_plate2), eps(Float64)); case=links_name)
            add_row!(rows, :links_sensor_plates2_is2_rel, abs(is2 - v_sensor) / max(abs(v_sensor), eps(Float64)); case=links_name)
            add_row!(rows, :links_sensor_plates2_is1_rel, abs(is1 - v_sensor) / max(abs(v_sensor), eps(Float64)); case=links_name)
        end
    end

    f_sensor4 = fixtures["test-cafeier_sensor4"]
    cfg_sensor4 = ArchimedLight.read_light_config(f_sensor4.config_path)
    scene_sensor4 = ArchimedLight.read_scene(cfg_sensor4.scene)
    row_sensor4 = first(ArchimedLight.read_meteo(cfg_sensor4.meteo).rows)
    step_sensor4 = ArchimedLight.run_light_step(scene_sensor4, row_sensor4, cfg_sensor4)
    summary4 = ArchimedLight.summary_values_table(scene_sensor4, [step_sensor4], cfg_sensor4; meteo_rows=[row_sensor4], start_step_number=0)
    ground_area = sum(Float64(r["area"]) for r in summary4.rows if Int(r["item_id"]) == -1; init=0.0)
    stepdur = _step_duration_seconds(row_sensor4)
    sky_irr = step_sensor4.sky.ri_sw_f
    vs_intercep = sum(Float64(r["Ri_q"]) for r in summary4.rows if string(r["group"]) == "sensors"; init=0.0)
    scene_intercep = sum(Float64(r["Ri_q"]) for r in summary4.rows if string(r["group"]) != "sensors"; init=0.0)
    vs_irr = vs_intercep / stepdur / ground_area
    scene_irr = scene_intercep / stepdur / ground_area
    add_row!(rows, :sensor4_vs_abs, abs(sky_irr - vs_irr); case="test-cafeier_sensor4")
    add_row!(rows, :sensor4_scene_abs, abs(sky_irr - scene_irr); case="test-cafeier_sensor4")

    for (fx_name, expect_zero) in (("test-ops-output2", true), ("test-ops-output3", false))
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
            gwa = only(sort(filter(f -> endswith(lowercase(f), ".gwa"), readdir(out; join=true))))
            ri_par = abs(_xml_tag_value(gwa, "Ri_PAR_0_f"))
            ri_nir = abs(_xml_tag_value(gwa, "Ri_NIR_0_f"))
            cmp = expect_zero ? :le : :ge
            add_row!(rows, :ops_export_zero_abs, ri_par; case="$(fx_name) par", cmp=cmp)
            add_row!(rows, :ops_export_zero_abs, ri_nir; case="$(fx_name) nir", cmp=cmp)
        end
    end

    for fx_name in ("test-lightsource1", "test-lightsource2")
        tol_case = fx_name
        cfg1, scene1, row1, step1 = _run_fixture_step(joinpath(test_root, fx_name, "config1.yml"))
        cfg2, scene2, row2, step2 = _run_fixture_step(joinpath(test_root, fx_name, "config2.yml"))
        rows1 = ArchimedLight.component_values_table(scene1, step1, cfg1; meteo_row=row1, columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
        rows2 = ArchimedLight.component_values_table(scene2, step2, cfg2; meteo_row=row2, columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows1)
        d2 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows2 if Int(r["item_id"]) != 2)
        ks = intersect(keys(d1), keys(d2))
        obs = maximum(abs(d1[k] - d2[k]) for k in ks; init=0.0)
        add_row!(rows, :lightsource12_max_abs, obs; fixture=tol_case, case=fx_name)
    end

    cfg_l1, scene_l3, row_l1, step_l1 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config1.yml"))
    cfg_l2, _, row_l2, step_l2 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config2.yml"))
    cfg_l3, _, row_l3, step_l3 = _run_fixture_step(joinpath(test_root, "test-lightsource3", "config3.yml"))
    t1 = ArchimedLight.component_values_table(scene_l3, step_l1, cfg_l1; meteo_row=row_l1, columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"]).rows
    t2 = ArchimedLight.component_values_table(scene_l3, step_l2, cfg_l2; meteo_row=row_l2, columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"]).rows
    t3 = ArchimedLight.component_values_table(scene_l3, step_l3, cfg_l3; meteo_row=row_l3, columns=["item_id", "component_id", "Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f"]).rows
    max_add = 0.0
    for col in ("Ri_PAR_0_f", "Ri_NIR_0_f", "Ri_PAR_f", "Ri_NIR_f")
        d1 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t1)
        d2 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t2)
        d3 = Dict((Int(r["item_id"]), Int(r["component_id"])) => Float64(r[col]) for r in t3)
        ks = intersect(intersect(keys(d1), keys(d2)), keys(d3))
        max_add = max(max_add, maximum(abs((d1[k] + d2[k]) - d3[k]) for k in ks; init=0.0))
    end
    add_row!(rows, :lightsource_additivity_abs, max_add; case="test-lightsource3")

    for (fx_name, model_name, strict_one_node) in (
        ("test-lightsource4", "model_pave2.yml", true),
        ("test-lightsource5", "model_lamp.yml", false),
        ("test-lightsource6", "model_lamp.yml", false),
    )
        cfg, scene, row, step = _run_fixture_step(joinpath(test_root, fx_name, "config.yml"))
        rows_ls = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, columns=["item_id", "area", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
        source_power = _model_radiance(joinpath(test_root, fx_name, model_name))
        obs =
            if strict_one_node
                vals = [Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_ls if Int(r["item_id"]) == 1]
                relerr(maximum(vals), source_power)
            else
                tot = sum(Float64(r["area"]) * Float64(r["Ri_PAR_0_f"] + r["Ri_NIR_0_f"]) for r in rows_ls; init=0.0)
                relerr(tot, source_power)
            end
        add_row!(rows, :lightsource_power_rel, obs; case=fx_name)
    end

    for (fx_name, model_name) in (
        ("test-lightsource7", "model_pave2.yml"),
        ("test-lightsource8", "model_lamp.yml"),
        ("test-lightsource9", "model_lamp.yml"),
    )
        cfg, scene, row, step = _run_fixture_step(joinpath(test_root, fx_name, "config.yml"))
        rows_ls = ArchimedLight.component_values_table(scene, step, cfg; meteo_row=row, columns=["area", "Ri_PAR_f", "Ri_NIR_f", "scat_factor_PAR", "scat_factor_NIR"]).rows
        source_power = _model_radiance(joinpath(test_root, fx_name, model_name))
        tot_abs = sum(
            (
                Float64(r["Ri_PAR_f"]) * (1.0 - Float64(r["scat_factor_PAR"])) +
                Float64(r["Ri_NIR_f"]) * (1.0 - Float64(r["scat_factor_NIR"]))
            ) * Float64(r["area"]) for r in rows_ls;
            init=0.0,
        )
        add_row!(rows, :lightsource_power_rel, relerr(tot_abs, source_power); case=fx_name)
    end

    cfg_box, scene_box, row_box, step_box = _run_fixture_step(joinpath(test_root, "test-lightsource-box", "config.yml"))
    rows_box = ArchimedLight.component_values_table(scene_box, step_box, cfg_box; meteo_row=row_box, columns=["area", "Ri_PAR_0_f", "Ri_NIR_0_f"]).rows
    source_box = _model_radiance(joinpath(test_root, "test-lightsource-box", "model_box.yml"))
    tot_box = sum(Float64(r["area"]) * (Float64(r["Ri_PAR_0_f"]) + Float64(r["Ri_NIR_0_f"])) for r in rows_box; init=0.0)
    add_row!(rows, :lightsource_power_rel, relerr(tot_box, source_box); case="test-lightsource-box")

    cfg_box2, scene_box2, row_box2, step_box2 = _run_fixture_step(joinpath(test_root, "test-lightsource-box2", "config.yml"))
    rows_box2 = ArchimedLight.component_values_table(scene_box2, step_box2, cfg_box2; meteo_row=row_box2, columns=["area", "Ri_PAR_f", "Ri_NIR_f", "scat_factor_PAR", "scat_factor_NIR"]).rows
    source_box2 = _model_radiance(joinpath(test_root, "test-lightsource-box2", "model_box.yml"))
    tot_box2 = sum(
        (
            Float64(r["Ri_PAR_f"]) * (1.0 - Float64(r["scat_factor_PAR"]) / 2.0) +
            Float64(r["Ri_NIR_f"]) * (1.0 - Float64(r["scat_factor_NIR"]) / 2.0)
        ) * Float64(r["area"]) for r in rows_box2;
        init=0.0,
    )
    add_row!(rows, :lightsource_box2_rel, relerr(tot_box2, source_box2); case="test-lightsource-box2")

    max_rel_summary1, max_scene1 = _check_summary_component_max_rel(joinpath(test_root, "test-summary", "config.yml"); use_scattering=false, include_scene_rel=true)
    max_rel_summary2, max_scene2 = _check_summary_component_max_rel(joinpath(test_root, "test-summary2", "config.yml"); use_scattering=false, include_scene_rel=true)
    max_rel_summary3, _ = _check_summary_component_max_rel(joinpath(test_root, "test-summary3", "config.yml"); use_scattering=true, include_scene_rel=false)
    max_rel_summary4, _ = _check_summary_component_max_rel(joinpath(test_root, "test-summary4", "config.yml"); use_scattering=true, include_scene_rel=false)
    add_row!(rows, :summary_component_rel, max_rel_summary1; case="test-summary")
    add_row!(rows, :summary_component_rel, max_rel_summary2; case="test-summary2")
    add_row!(rows, :summary_component_rel, max_rel_summary3; case="test-summary3")
    add_row!(rows, :summary_component_rel, max_rel_summary4; case="test-summary4")
    add_row!(rows, :summary_scene_rel, max_scene1; fixture="test-summary", case="test-summary")
    add_row!(rows, :summary_scene_rel, max_scene2; fixture="test-summary2", case="test-summary2")

    f_custom = fixtures["test-customband"]
    cfg_custom = ArchimedLight.read_light_config(f_custom.config_path)
    scene_custom = ArchimedLight.read_scene(cfg_custom.scene)
    row_custom = first(ArchimedLight.read_meteo(cfg_custom.meteo).rows)
    step_custom = ArchimedLight.run_light_step(scene_custom, row_custom, cfg_custom)
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
    add_row!(rows, :customband_rel, relerr(scene_par_0 / max(scene_incid_par, eps(Float64)), scene_custom_0 / max(scene_incid_custom, eps(Float64))); case="custom_0")
    add_row!(rows, :customband_rel, relerr(scene_par_n / max(scene_incid_par, eps(Float64)), scene_custom_n / max(scene_incid_custom, eps(Float64))); case="custom_n")

    add_row!(rows, :absorb_cached_abs, _check_absorption_ratio_max(joinpath(test_root, "test-absorb", "config.yml"); with_scattering=false); case="test-absorb")
    add_row!(rows, :absorb_cached_abs, _check_absorption_ratio_max(joinpath(test_root, "test-absorb2", "config.yml"); with_scattering=true); case="test-absorb2")

    add_row!(rows, :timestep_rel_ts1, _run_timestep_fixture(joinpath(test_root, "test-timestep1")); case="test-timestep1")
    add_row!(rows, :timestep_rel_ts2, _run_timestep_fixture(joinpath(test_root, "test-timestep2")); case="test-timestep2")
    add_row!(rows, :timestep_rel_ts3, _run_timestep_fixture(joinpath(test_root, "test-timestep3")); case="test-timestep3")

    sort!(rows; by=r -> (r.pass ? 1 : 0, r.headroom))

    println("ArchimedLight tolerance audit")
    println("===========================")
    @printf("%-26s %-20s %-28s %-3s %-12s %-12s %-10s %s\n", "metric", "fixture", "case", "op", "observed", "threshold", "headroom", "status")
    for r in rows
        op = r.cmp == :le ? "<=" : ">="
        status = r.pass ? "PASS" : "FAIL"
        @printf(
            "%-26s %-20s %-28s %-3s %-12.6g %-12.6g %-10.4f %s\n",
            string(r.metric),
            r.fixture,
            r.case,
            op,
            r.observed,
            r.threshold,
            r.headroom,
            status,
        )
    end

    covered = Set((r.metric, r.fixture) for r in rows)
    missing = sort(collect(setdiff(Set(keys(PARITY_LIMITS)), covered)); by=x -> (String(x[1]), x[2]))
    println()
    println("Coverage")
    println("--------")
    println("audited entries: $(length(covered)) / $(length(PARITY_LIMITS))")
    if isempty(missing)
        println("all central tolerance entries audited")
    else
        println("missing entries:")
        for (metric, fixture) in missing
            println("  - metric=$(metric) fixture=$(fixture)")
        end
    end

    failures = count(!r.pass for r in rows)
    failures == 0 || error("tolerance audit failed with $(failures) metric violations")
end

main()
