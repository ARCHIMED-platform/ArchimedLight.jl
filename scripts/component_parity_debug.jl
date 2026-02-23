#!/usr/bin/env julia

using ArchimedLight
using Statistics

include(joinpath(dirname(@__DIR__), "test", "java_parity_harness.jl"))

to_int(v) = v isa Number ? Int(round(v)) : parse(Int, strip(string(v)))

function expected_component_metric(component_csv::AbstractString, key::Tuple{Int,Int})
    rows = read_java_csv(component_csv)
    for r in rows
        k = (to_int(getproperty(r, :item_id)), to_int(getproperty(r, :component_id)))
        k == key || continue
        area = Float64(getproperty(r, :area))
        hits = Int(floor(area * Float64(getproperty(r, :surface_hits))))
        irr0 = Float64(getproperty(r, :irradiance_withoutScattering_PAR_NIR))
        return (area=area, hits=hits, irr0=irr0)
    end
    return nothing
end

function component_sector_report(
    fx::ParityFixture,
    key::Tuple{Int,Int};
    all_in_turtle::Bool,
)
    cfg = ArchimedLight.read_light_config(fx.config_path)
    cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
    scene = ArchimedLight.read_scene(cfg.scene)
    row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
    sky = ArchimedLight.compute_sky(row, cfg)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)

    vertices, faces, face2node, _, plotbox, node_group = ArchimedLight._scene_geometry_for_interception(scene, cfg)
    key_by_node = ArchimedLight._interception_java_keys(scene, cfg)
    node_ids = Int[nid for (nid, k) in key_by_node if k == key]
    isempty(node_ids) && error("component $(key) not found in scene")
    length(node_ids) == 1 || error("component $(key) maps to $(length(node_ids)) nodes, expected one")
    target_nid = node_ids[1]

    virtual_nodes = ArchimedLight._virtual_sensor_node_ids(node_group, cfg)
    upper_hit = ArchimedLight._use_upper_hit_pixel_table(cfg)
    cache_ctx = ArchimedLight._projection_cache_context(vertices, faces, face2node, plotbox, cfg)

    rows = NamedTuple[]
    total_hits = 0
    total_power = 0.0
    for i in eachindex(turtle.sectors)
        _, _, projected_mesh_area, projected_pixels_area = ArchimedLight._direction_projection_cached(
            vertices,
            faces,
            face2node,
            turtle.sectors[i].direction,
            cfg,
            plotbox,
            cache_ctx;
            upper_hit=upper_hit,
        )
        vis, node_hits = ArchimedLight._rasterize_direction_java(
            vertices,
            faces,
            face2node,
            turtle.sectors[i].direction,
            cfg,
            plotbox;
            cache_ctx=cache_ctx,
            virtual_nodes=virtual_nodes,
            upper_hit=upper_hit,
        )
        hits = get(node_hits, target_nid, 0)
        area = get(vis, target_nid, 0.0)
        flux_sw = fluxes.par[i] + fluxes.nir[i]
        power = area * flux_sw
        pma = get(projected_mesh_area, target_nid, 0.0)
        ppa = get(projected_pixels_area, target_nid, 0.0)
        ratio = cfg.area_ratio && ppa > 0.0 ? pma / ppa : 1.0
        total_hits += hits
        total_power += power
        push!(
            rows,
            (
                dir=i,
                source=String(turtle.sectors[i].source),
                hits=hits,
                visible_area=area,
                flux_sw=flux_sw,
                power_sw=power,
                area_ratio=ratio,
            ),
        )
    end

    step = ArchimedLight.run_light_step(scene, row, cfg)
    got_power_step =
        get(step.first_order.incident_par_power_per_node, target_nid, 0.0) +
        get(step.first_order.incident_nir_power_per_node, target_nid, 0.0)
    got_hits_step = get(step.first_order.hits_per_node, target_nid, 0)

    expected_path = _expected_component_values_path(fx; all_in_turtle=all_in_turtle)
    expected = expected_path === nothing ? nothing : expected_component_metric(expected_path, key)

    (
        key=key,
        node_id=target_nid,
        all_in_turtle=all_in_turtle,
        rows=rows,
        sum_hits=total_hits,
        sum_power=total_power,
        step_hits=got_hits_step,
        step_power=got_power_step,
        expected=expected,
    )
end

function _resolve_java_debug_dir(path::AbstractString)
    isdir(path) || return nothing
    has_direct_logs = !isempty(readdir(path; join=true))
    if has_direct_logs && !isempty(filter(p -> startswith(basename(p), "log-arearatio-"), readdir(path; join=true)))
        return path
    end
    step_dirs = filter(p -> isdir(p), readdir(path; join=true))
    for d in sort(step_dirs)
        if !isempty(filter(p -> startswith(basename(p), "log-arearatio-"), readdir(d; join=true)))
            return d
        end
    end
    return nothing
end

function java_debug_metrics(debug_dir::AbstractString, key::Tuple{Int,Int})
    data_dir = _resolve_java_debug_dir(debug_dir)
    data_dir === nothing && return nothing

    ratio = Dict{Int,Float64}()
    vis = Dict{Int,Float64}()

    for p in sort(readdir(data_dir; join=true))
        b = basename(p)
        if startswith(b, "log-arearatio-")
            rows = read_java_csv(p)
            for r in rows
                item = to_int(getproperty(r, :plantid))
                comp = to_int(getproperty(r, :nodeid))
                (item, comp) == key || continue
                dir = to_int(getproperty(r, :dir))
                ratio[dir + 1] = Float64(getproperty(r, :arearatio))
            end
        elseif startswith(b, "log-visiblearea-")
            rows = read_java_csv(p)
            for r in rows
                item = to_int(getproperty(r, :plantid))
                comp = to_int(getproperty(r, :nodeid))
                (item, comp) == key || continue
                dir = to_int(getproperty(r, :dir))
                vis[dir + 1] = Float64(getproperty(r, :visiblearea))
            end
        end
    end

    return (dir_area_ratio=ratio, dir_visible_area=vis, source_dir=data_dir)
end

function java_debug_metrics_all(debug_dir::AbstractString)
    data_dir = _resolve_java_debug_dir(debug_dir)
    data_dir === nothing && return nothing

    ratio = Dict{Tuple{Int,Int,Int},Float64}()
    vis = Dict{Tuple{Int,Int,Int},Float64}()

    for p in sort(readdir(data_dir; join=true))
        b = basename(p)
        if startswith(b, "log-arearatio-")
            rows = read_java_csv(p)
            for r in rows
                dir = to_int(getproperty(r, :dir)) + 1
                item = to_int(getproperty(r, :plantid))
                comp = to_int(getproperty(r, :nodeid))
                ratio[(dir, item, comp)] = Float64(getproperty(r, :arearatio))
            end
        elseif startswith(b, "log-visiblearea-")
            rows = read_java_csv(p)
            for r in rows
                dir = to_int(getproperty(r, :dir)) + 1
                item = to_int(getproperty(r, :plantid))
                comp = to_int(getproperty(r, :nodeid))
                vis[(dir, item, comp)] = Float64(getproperty(r, :visiblearea))
            end
        end
    end

    return (dir_area_ratio=ratio, dir_visible_area=vis, source_dir=data_dir)
end

function scan_java_debug_mismatches(
    fx::ParityFixture,
    debug_dir::AbstractString;
    all_in_turtle::Bool=false,
    topn::Int=40,
)
    jdbg = java_debug_metrics_all(debug_dir)
    jdbg === nothing && error("No Java debug logs found under $(repr(debug_dir))")

    cfg = ArchimedLight.read_light_config(fx.config_path)
    cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle)
    scene = ArchimedLight.read_scene(cfg.scene)
    row = first(ArchimedLight.read_meteo(cfg.meteo).rows)
    sky = ArchimedLight.compute_sky(row, cfg)
    turtle = ArchimedLight.build_turtle(cfg, sky)

    vertices, faces, face2node, _, plotbox, node_group = ArchimedLight._scene_geometry_for_interception(scene, cfg)
    key_by_node = ArchimedLight._interception_java_keys(scene, cfg)
    node_by_key = Dict{Tuple{Int,Int},Int}(k => nid for (nid, k) in key_by_node)

    virtual_nodes = ArchimedLight._virtual_sensor_node_ids(node_group, cfg)
    upper_hit = ArchimedLight._use_upper_hit_pixel_table(cfg)
    cache_ctx = ArchimedLight._projection_cache_context(vertices, faces, face2node, plotbox, cfg)

    rows = NamedTuple[]
    for d in eachindex(turtle.sectors)
        _, _, projected_mesh_area, projected_pixels_area = ArchimedLight._direction_projection_cached(
            vertices,
            faces,
            face2node,
            turtle.sectors[d].direction,
            cfg,
            plotbox,
            cache_ctx;
            upper_hit=upper_hit,
        )
        vis, node_hits = ArchimedLight._rasterize_direction_java(
            vertices,
            faces,
            face2node,
            turtle.sectors[d].direction,
            cfg,
            plotbox;
            cache_ctx=cache_ctx,
            virtual_nodes=virtual_nodes,
            upper_hit=upper_hit,
        )

        for ((dir, item, comp), jratio) in jdbg.dir_area_ratio
            dir == d || continue
            nid = get(node_by_key, (item, comp), 0)
            nid == 0 && continue

            ppa = get(projected_pixels_area, nid, 0.0)
            pma = get(projected_mesh_area, nid, 0.0)
            gratio = cfg.area_ratio && ppa > 0.0 ? pma / ppa : 1.0

            jvis = get(jdbg.dir_visible_area, (dir, item, comp), NaN)
            gvis = get(vis, nid, 0.0)
            hits = get(node_hits, nid, 0)

            ratio_abs = abs(gratio - jratio)
            ratio_rel = ratio_abs / max(abs(jratio), eps(Float64))
            vis_abs = if isfinite(jvis)
                abs(gvis - jvis)
            else
                NaN
            end
            vis_pix = if isfinite(vis_abs)
                vis_abs / max(plotbox.pixel_area, eps(Float64))
            else
                NaN
            end

            push!(
                rows,
                (
                    dir=dir,
                    item=item,
                    comp=comp,
                    nid=nid,
                    source=String(turtle.sectors[d].source),
                    hits=hits,
                    jratio=jratio,
                    gratio=gratio,
                    ratio_abs=ratio_abs,
                    ratio_rel=ratio_rel,
                    jvis=jvis,
                    gvis=gvis,
                    vis_abs=vis_abs,
                    vis_pix=vis_pix,
                ),
            )
        end
    end

    sort!(rows; by=r -> (
        isfinite(r.ratio_rel) ? r.ratio_rel : -Inf,
        isfinite(r.vis_pix) ? r.vis_pix : -Inf,
        r.vis_abs,
    ), rev=true)

    n = min(topn, length(rows))
    println()
    println("scan fixture=$(fx.name) all_in_turtle=$(all_in_turtle) java_debug=$(jdbg.source_dir)")
    println("rows compared: $(length(rows))")
    if !isempty(rows)
        rrels = [r.ratio_rel for r in rows if isfinite(r.ratio_rel)]
        vpix = [r.vis_pix for r in rows if isfinite(r.vis_pix)]
        !isempty(rrels) && println("ratio_rel mean=$(mean(rrels)) p99=$(quantile(rrels, 0.99)) max=$(maximum(rrels))")
        !isempty(vpix) && println("vis_pixel_abs mean=$(mean(vpix)) p99=$(quantile(vpix, 0.99)) max=$(maximum(vpix))")
    end

    bydir = Dict{Int,Int}()
    for r in rows
        if (r.ratio_rel > 0.02) || (isfinite(r.vis_pix) && r.vis_pix > 0.1)
            bydir[r.dir] = get(bydir, r.dir, 0) + 1
        end
    end
    println("dirs with significant mismatches:")
    for d in sort(collect(keys(bydir)))
        println("  dir=$(d) count=$(bydir[d])")
    end

    println()
    println("top $(n) mismatches:")
    for r in rows[1:n]
        println(
            "  dir=$(lpad(r.dir, 3)) key=($(r.item),$(r.comp)) hits=$(lpad(r.hits, 4)) " *
            "ratio julia=$(round(r.gratio; digits=9)) java=$(round(r.jratio; digits=9)) rel=$(round(r.ratio_rel; digits=6)) " *
            "vis julia=$(round(r.gvis; digits=10)) java=$(round(r.jvis; digits=10)) pixΔ=$(round(r.vis_pix; digits=4))",
        )
    end
end

function print_report(rep; topn::Int=12, java_debug=nothing)
    println()
    println("mode all_in_turtle=$(rep.all_in_turtle) key=$(rep.key) node=$(rep.node_id)")
    println("----------------------------------------------------------------")
    println("hits (sector sum): $(rep.sum_hits)")
    println("hits (run_light_step): $(rep.step_hits)")
    println("power SW (sector sum): $(round(rep.sum_power; digits=12))")
    println("power SW (run_light_step): $(round(rep.step_power; digits=12))")
    if rep.expected !== nothing
        exp = rep.expected
        exp_power = exp.irr0 * exp.area
        println("expected hits: $(exp.hits)")
        println("expected power SW: $(round(exp_power; digits=12))")
        println("expected irradiance SW: $(round(exp.irr0; digits=9))")
        println("julia irradiance SW: $(round(rep.step_power / max(exp.area, eps(Float64)); digits=9))")
    end

    rows = sort(rep.rows; by=r -> (r.power_sw, r.hits), rev=true)
    n = min(topn, length(rows))
    println()
    println("top $(n) contributing directions:")
    for r in rows[1:n]
        println(
            "  dir=$(lpad(r.dir, 3)) source=$(r.source) hits=$(lpad(r.hits, 4)) " *
            "area=$(round(r.visible_area; digits=10)) flux=$(round(r.flux_sw; digits=6)) " *
            "power=$(round(r.power_sw; digits=10))",
        )
    end

    if java_debug !== nothing
        println()
        println("java debug dir: $(java_debug.source_dir)")
        println("dir-by-dir area ratio / visible area:")
        for r in sort(rep.rows; by=x -> x.dir)
            j_ratio = get(java_debug.dir_area_ratio, r.dir, NaN)
            j_vis = get(java_debug.dir_visible_area, r.dir, NaN)
            println(
                "  dir=$(lpad(r.dir, 3)) " *
                "ratio julia=$(round(r.area_ratio; digits=9)) java=$(round(j_ratio; digits=9)) " *
                "vis julia=$(round(r.visible_area; digits=10)) java=$(round(j_vis; digits=10))",
            )
        end
    end
end

function parse_mode(arg::String)
    s = lowercase(strip(arg))
    if s in ("both", "all")
        return (false, true)
    elseif s in ("raycast", "false", "0")
        return (false,)
    elseif s in ("turtle", "true", "1")
        return (true,)
    end
    error("invalid mode $(repr(arg)); use one of: both, raycast, turtle")
end

function main()
    if !isempty(ARGS) && lowercase(strip(ARGS[1])) == "scan"
        length(ARGS) >= 3 || error(
            "usage: julia --project=. scripts/component_parity_debug.jl scan <fixture> <java_debug_dir> [mode=raycast|turtle] [topn=40]",
        )
        fixture_name = String(ARGS[2])
        java_debug_dir = String(ARGS[3])
        mode = length(ARGS) >= 4 ? parse_mode(ARGS[4]) : (false,)
        topn = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 40

        fixtures = Dict(f.name => f for f in light_parity_fixtures())
        haskey(fixtures, fixture_name) || error("unknown fixture $(repr(fixture_name))")
        fx = fixtures[fixture_name]

        for all_in_turtle in mode
            scan_java_debug_mismatches(fx, java_debug_dir; all_in_turtle=all_in_turtle, topn=topn)
        end
        return
    end

    length(ARGS) >= 3 || error(
        "usage: julia --project=. scripts/component_parity_debug.jl <fixture> <item_id> <component_id> [mode=both] [topn=12] [java_debug_dir]",
    )
    fixture_name = String(ARGS[1])
    item_id = parse(Int, ARGS[2])
    component_id = parse(Int, ARGS[3])
    mode = length(ARGS) >= 4 ? parse_mode(ARGS[4]) : (false, true)
    topn = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 12
    java_debug_dir = length(ARGS) >= 6 ? String(ARGS[6]) : nothing

    fixtures = Dict(f.name => f for f in light_parity_fixtures())
    haskey(fixtures, fixture_name) || error("unknown fixture $(repr(fixture_name))")
    fx = fixtures[fixture_name]
    key = (item_id, component_id)

    for all_in_turtle in mode
        rep = component_sector_report(fx, key; all_in_turtle=all_in_turtle)
        jdbg = java_debug_dir === nothing ? nothing : java_debug_metrics(java_debug_dir, key)
        print_report(rep; topn=topn, java_debug=jdbg)
    end
end

main()
