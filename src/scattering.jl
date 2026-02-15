function _dict_copy_float(d::Dict{Int,Float64})
    Dict{Int,Float64}(k => v for (k, v) in d)
end

function _normalized_capture(projected_area_per_node::Dict{Int,Float64})
    s = sum(values(projected_area_per_node))
    s <= 0.0 && return Dict{Int,Float64}(k => 0.0 for k in keys(projected_area_per_node))
    Dict{Int,Float64}(k => v / s for (k, v) in projected_area_per_node)
end

function compute_scattering(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast
)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode")

    capture = _normalized_capture(first.projected_area_per_node)

    added_par = Dict{Int,Float64}(k => 0.0 for k in keys(capture))
    added_nir = Dict{Int,Float64}(k => 0.0 for k in keys(capture))

    current_par = _dict_copy_float(first.incident_par_power_per_node)
    current_nir = _dict_copy_float(first.incident_nir_power_per_node)

    initial_total = sum(values(current_par)) + sum(values(current_nir))
    threshold = cfg.scattering_stop_ratio * max(initial_total, eps(Float64))

    converged = false
    iter = 0
    for k in 1:cfg.scattering_max_iter
        iter = k
        emitted_par_total = cfg.scattering_coeff_par * sum(values(current_par))
        emitted_nir_total = cfg.scattering_coeff_nir * sum(values(current_nir))
        emitted_total = emitted_par_total + emitted_nir_total

        if emitted_total <= threshold
            converged = true
            break
        end

        next_par = Dict{Int,Float64}(nid => 0.0 for nid in keys(capture))
        next_nir = Dict{Int,Float64}(nid => 0.0 for nid in keys(capture))

        for (nid, w) in capture
            p = emitted_par_total * w
            n = emitted_nir_total * w
            next_par[nid] = p
            next_nir[nid] = n
            added_par[nid] = get(added_par, nid, 0.0) + p
            added_nir[nid] = get(added_nir, nid, 0.0) + n
        end

        current_par = next_par
        current_nir = next_nir
    end

    ScatteringResult(added_par, added_nir, iter, converged)
end

