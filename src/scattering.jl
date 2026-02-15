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

function _next_scattering_iteration(current::Dict{Int,Float64}, all_hits::Dict{Int,Int}, pair_counts, coeff::Float64, node_ids)
    hit_energy = _dict_zero(node_ids)
    for nid in node_ids
        nh = get(all_hits, nid, 0)
        hit_energy[nid] = nh > 0 ? get(current, nid, 0.0) * coeff / nh / 2.0 : 0.0
    end

    next = _dict_zero(node_ids)
    for ((to, from), cnt) in pair_counts
        next[to] = get(next, to, 0.0) + cnt * get(hit_energy, from, 0.0)
    end
    next
end

function _group_optical_coeffs(cfg::LightConfig)
    models = get(cfg.raw, "models", nothing)
    models isa AbstractVector || return Dict{String,Tuple{Float64,Float64}}()

    base = get(cfg.raw, "__base_dir", dirname(cfg.scene))
    coeffs = Dict{String,Tuple{Float64,Float64}}()

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

        par = haskey(op, "PAR") ? Float64(op["PAR"]) : cfg.scattering_coeff_par
        nir = haskey(op, "NIR") ? Float64(op["NIR"]) : cfg.scattering_coeff_nir
        coeffs[group] = (par, nir)
    end

    coeffs
end

function compute_scattering(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode")

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

    current_par = Dict{Int,Float64}(nid => get(first.incident_par_power_per_node, nid, 0.0) for nid in node_ids)
    current_nir = Dict{Int,Float64}(nid => get(first.incident_nir_power_per_node, nid, 0.0) for nid in node_ids)
    added_par = _dict_zero(node_ids)
    added_nir = _dict_zero(node_ids)

    ref_par = _sum_dict_values(current_par)
    ref_nir = _sum_dict_values(current_nir)
    thr_par = cfg.scattering_stop_ratio * max(ref_par, eps(Float64))
    thr_nir = cfg.scattering_stop_ratio * max(ref_nir, eps(Float64))

    converged = false
    iterations = 0
    for it in 1:cfg.scattering_max_iter
        iterations = it
        coeff_par = Dict{Int,Float64}()
        coeff_nir = Dict{Int,Float64}()
        for nid in node_ids
            g = get(node_group, nid, "")
            cp, cn = get(group_coeffs, g, (cfg.scattering_coeff_par, cfg.scattering_coeff_nir))
            coeff_par[nid] = cp
            coeff_nir[nid] = cn
        end

        # Build per-node hit energies with optical properties.
        hit_energy_par = _dict_zero(node_ids)
        hit_energy_nir = _dict_zero(node_ids)
        for nid in node_ids
            nh = get(all_hits, nid, 0)
            if nh > 0
                hit_energy_par[nid] = get(current_par, nid, 0.0) * get(coeff_par, nid, cfg.scattering_coeff_par) / nh / 2.0
                hit_energy_nir[nid] = get(current_nir, nid, 0.0) * get(coeff_nir, nid, cfg.scattering_coeff_nir) / nh / 2.0
            end
        end

        next_par = _dict_zero(node_ids)
        next_nir = _dict_zero(node_ids)
        for ((to, from), cnt) in pair_counts
            next_par[to] = get(next_par, to, 0.0) + cnt * get(hit_energy_par, from, 0.0)
            next_nir[to] = get(next_nir, to, 0.0) + cnt * get(hit_energy_nir, from, 0.0)
        end

        total_next_par = _sum_dict_values(next_par)
        total_next_nir = _sum_dict_values(next_nir)

        for nid in node_ids
            added_par[nid] = get(added_par, nid, 0.0) + get(next_par, nid, 0.0)
            added_nir[nid] = get(added_nir, nid, 0.0) + get(next_nir, nid, 0.0)
        end

        current_par = next_par
        current_nir = next_nir

        if total_next_par <= thr_par && total_next_nir <= thr_nir
            converged = true
            break
        end
    end

    ScatteringResult(added_par, added_nir, iterations, converged)
end
