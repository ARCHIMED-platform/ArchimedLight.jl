function _dict_zero(node_ids)
    Dict{Int,Float64}(nid => 0.0 for nid in node_ids)
end

function _sum_dict_values(d::Dict{Int,Float64})
    s = 0.0
    for v in values(d)
        s += v
    end
    s
end

function _pair_counts_for_scattering(scene::SceneGeometry, turtle::TurtleGrid, cfg::LightConfig)
    vertices, faces, face2node, node_ids, plotbox, node_group = _scene_geometry_for_interception(scene, cfg)
    pair_counts = Dict{Tuple{Int,Int},Int}()
    sun_hits = Dict{Int,Int}()

    for sector in turtle.sectors
        pixel_hits, node_hits, _, _ = _direction_projection(vertices, faces, face2node, sector.direction, cfg, plotbox)

        if sector.source == :sun
            for (nid, h) in node_hits
                sun_hits[nid] = get(sun_hits, nid, 0) + h
            end
            continue
        end

        for stack in values(pixel_hits)
            length(stack) <= 1 && continue
            sort!(stack, by=x -> x[1], rev=true)
            for h in 1:(length(stack) - 1)
                above = stack[h][2]
                below = stack[h + 1][2]
                pair_counts[(above, below)] = get(pair_counts, (above, below), 0) + 1
                pair_counts[(below, above)] = get(pair_counts, (below, above), 0) + 1
            end
        end
    end

    return pair_counts, sun_hits, node_ids, node_group
end

function _all_dir_hits_for_scattering(first::FirstOrderResult, sun_hits::Dict{Int,Int}, cfg::LightConfig, node_ids)
    all_hits = Dict{Int,Int}(nid => get(first.hits_per_node, nid, 0) for nid in node_ids)
    if !cfg.all_in_turtle
        for (nid, hsun) in sun_hits
            all_hits[nid] = max(0, get(all_hits, nid, 0) - hsun)
        end
    end
    all_hits
end

function _group_optical_coeffs(cfg::LightConfig)
    models = get(cfg.raw, "models", nothing)
    models isa AbstractVector || return Dict{String,Dict{String,Float64}}()

    base = get(cfg.raw, "__base_dir", dirname(cfg.scene))
    coeffs = Dict{String,Dict{String,Float64}}()

    for m in models
        mp = String(m)
        path = isabspath(mp) ? mp : normpath(joinpath(base, mp))
        isfile(path) || continue

        d = try
            YAML.load_file(path)
        catch
            nothing
        end
        d isa AbstractDict || continue
        d = _to_string_dict(d)

        group = haskey(d, "Group") ? string(d["Group"]) : ""
        isempty(group) && continue

        types = get(d, "Type", nothing)
        types isa AbstractDict || continue
        isempty(types) && continue

        tkey = first(keys(types))
        tconf = _to_string_dict(types[tkey])
        tconf isa AbstractDict || continue

        inter = get(tconf, "Interception", nothing)
        inter isa AbstractDict || continue
        inter = _to_string_dict(inter)
        iuse = get(inter, "use", nothing)
        iconf =
            if iuse !== nothing && haskey(inter, string(iuse))
                inter[string(iuse)]
            else
                inter
            end
        iconf isa AbstractDict || continue
        iconf = _to_string_dict(iconf)
        op = get(iconf, "optical_properties", nothing)
        op isa AbstractDict || continue
        op = _to_string_dict(op)

        c = Dict{String,Float64}()
        for (k, v) in op
            try
                c[uppercase(string(k))] = Float64(v)
            catch
            end
        end
        coeffs[group] = c
    end

    coeffs
end

function _scattering_context(scene::SceneGeometry, turtle::TurtleGrid, first::FirstOrderResult, cfg::LightConfig)
    pair_counts, sun_hits, geom_node_ids, node_group = _pair_counts_for_scattering(scene, turtle, cfg)
    group_coeffs = _group_optical_coeffs(cfg)

    node_set = Set{Int}()
    for nid in geom_node_ids
        push!(node_set, nid)
    end
    for nid in keys(first.incident_par_power_per_node)
        push!(node_set, nid)
    end
    for nid in keys(first.incident_nir_power_per_node)
        push!(node_set, nid)
    end
    node_ids = collect(node_set)
    all_hits = _all_dir_hits_for_scattering(first, sun_hits, cfg, node_ids)
    return pair_counts, all_hits, node_ids, node_group, group_coeffs
end

function _default_band_coeff(cfg::LightConfig, band_key::String)
    bk = uppercase(band_key)
    bk == "NIR" && return cfg.scattering_coeff_nir
    return cfg.scattering_coeff_par
end

function _scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    pair_counts,
    all_hits::Dict{Int,Int},
    node_ids,
    node_group::Dict{Int,String},
    group_coeffs::Dict{String,Dict{String,Float64}},
    cfg::LightConfig,
    band_key::String,
    default_coeff::Float64,
)
    coeff_by_node = Dict{Int,Float64}()
    for nid in node_ids
        g = get(node_group, nid, "")
        c = get(group_coeffs, g, Dict{String,Float64}())
        coeff_by_node[nid] = get(c, uppercase(band_key), default_coeff)
    end

    current = Dict{Int,Float64}(nid => get(initial_power_per_node, nid, 0.0) for nid in node_ids)
    added = _dict_zero(node_ids)
    ref = _sum_dict_values(current)
    thr = cfg.scattering_stop_ratio * max(ref, eps(Float64))
    iterations = 0

    converged = false
    for it in 1:cfg.scattering_max_iter
        iterations = it

        hit_energy = _dict_zero(node_ids)
        for nid in node_ids
            nh = get(all_hits, nid, 0)
            if nh > 0
                hit_energy[nid] = get(current, nid, 0.0) * get(coeff_by_node, nid, default_coeff) / nh / 2.0
            end
        end

        next = _dict_zero(node_ids)
        for ((to, from), cnt) in pair_counts
            next[to] = get(next, to, 0.0) + cnt * get(hit_energy, from, 0.0)
        end

        total_next = _sum_dict_values(next)

        for nid in node_ids
            added[nid] = get(added, nid, 0.0) + get(next, nid, 0.0)
        end

        current = next

        if total_next <= thr
            converged = true
            break
        end
    end
    return added, iterations, converged
end

function compute_scattering_band(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode")
    pair_counts, all_hits, node_ids, node_group, group_coeffs = _scattering_context(scene, turtle, first, cfg)
    initial =
        if initial_power_per_node === nothing
            Dict{Int,Float64}(nid => get(first.incident_par_power_per_node, nid, 0.0) for nid in node_ids)
        else
            Dict{Int,Float64}(nid => get(initial_power_per_node, nid, 0.0) for nid in node_ids)
        end
    dflt = isnothing(default_coeff) ? _default_band_coeff(cfg, String(band)) : default_coeff
    added, iterations, converged = _scattering_one_band(
        initial,
        pair_counts,
        all_hits,
        node_ids,
        node_group,
        group_coeffs,
        cfg,
        String(band),
        dflt,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged)
end

function compute_scattering(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode")
    pair_counts, all_hits, node_ids, node_group, group_coeffs = _scattering_context(scene, turtle, first, cfg)

    initial_par = Dict{Int,Float64}(nid => get(first.incident_par_power_per_node, nid, 0.0) for nid in node_ids)
    initial_nir = Dict{Int,Float64}(nid => get(first.incident_nir_power_per_node, nid, 0.0) for nid in node_ids)

    added_par, it_par, conv_par = _scattering_one_band(
        initial_par,
        pair_counts,
        all_hits,
        node_ids,
        node_group,
        group_coeffs,
        cfg,
        "PAR",
        cfg.scattering_coeff_par,
    )
    added_nir, it_nir, conv_nir = _scattering_one_band(
        initial_nir,
        pair_counts,
        all_hits,
        node_ids,
        node_group,
        group_coeffs,
        cfg,
        "NIR",
        cfg.scattering_coeff_nir,
    )

    ScatteringResult(added_par, added_nir, max(it_par, it_nir), conv_par && conv_nir)
end
