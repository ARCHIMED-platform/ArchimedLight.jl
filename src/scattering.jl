function _dict_zero(node_ids)
    Dict{Int,Float64}(nid => 0.0 for nid in node_ids)
end

function _copy_node_values(source::Dict{Int,Float64}, node_ids)
    Dict{Int,Float64}(nid => get(source, nid, 0.0) for nid in node_ids)
end

function _sum_dict_values(d::Dict{Int,Float64})
    s = 0.0
    for v in values(d)
        s += v
    end
    s
end

function _add_pair_count!(pair_counts::Dict{Tuple{Int,Int},Int}, to::Int, from::Int)
    pair_counts[(to, from)] = get(pair_counts, (to, from), 0) + 1
    return nothing
end

function _accept_scattering_link(node_group::Dict{Int,String}, to::Int, from::Int)
    !(get(node_group, to, "") == "pavement" && get(node_group, from, "") == "pavement")
end

function _stack_transfer_pairs!(
    pair_counts::Dict{Tuple{Int,Int},Int},
    node_ids_stack::Vector{Int},
    scatt_up::Vector{Int},
    scatt_down::Vector{Int},
    node_group::Dict{Int,String},
)
    n_hits = length(node_ids_stack)
    @inbounds for h in (n_hits - 1):-1:1
        to_above = node_ids_stack[h]
        from_below = scatt_up[h + 1]
        if from_below != 0 && _accept_scattering_link(node_group, to_above, from_below)
            _add_pair_count!(pair_counts, to_above, from_below)
        end

        to_below = node_ids_stack[h + 1]
        from_above = scatt_down[h]
        if from_above != 0 && _accept_scattering_link(node_group, to_below, from_above)
            _add_pair_count!(pair_counts, to_below, from_above)
        end
    end
    return nothing
end

function _pair_counts_for_scattering(scene::SceneGeometry, turtle::TurtleGrid, cfg::LightConfig)
    vertices, faces, face2node, node_ids, plotbox, node_group = _scene_geometry_for_interception(scene, cfg)
    virtual_nodes = _virtual_sensor_node_ids(node_group, cfg)
    cache_ctx = _projection_cache_context(vertices, faces, face2node, plotbox, cfg)
    pair_counts = Dict{Tuple{Int,Int},Int}()
    sun_hits = Dict{Int,Int}()

    for sector in turtle.sectors
        pixel_hits, node_hits, _, _ =
            _direction_projection_cached(vertices, faces, face2node, sector.direction, cfg, plotbox, cache_ctx)

        if sector.source == :sun
            for (nid, h) in node_hits
                sun_hits[nid] = get(sun_hits, nid, 0) + h
            end
            continue
        end

        for stack in values(pixel_hits)
            length(stack) <= 1 && continue
            # Java PixelTable uses a stable height sort (Collections.sort on comparator).
            sort!(stack, by=x -> x[1], rev=true, alg=Base.Sort.MergeSort)

            n_hits = length(stack)
            node_ids_stack = Vector{Int}(undef, n_hits)
            @inbounds for i in 1:n_hits
                node_ids_stack[i] = stack[i][2]
            end

            # Mirror EnergyTransferTask: carry nearest non-virtual diffuser upward/downward
            # across virtual sensors, then transfer between adjacent stack positions.
            scatt_up = Vector{Int}(undef, n_hits)
            up = 0
            @inbounds for h in n_hits:-1:1
                nid = node_ids_stack[h]
                if !(nid in virtual_nodes)
                    up = nid
                end
                scatt_up[h] = up
            end

            scatt_down = Vector{Int}(undef, n_hits)
            down = 0
            @inbounds for h in 1:n_hits
                nid = node_ids_stack[h]
                if !(nid in virtual_nodes)
                    down = nid
                end
                scatt_down[h] = down
            end

            _stack_transfer_pairs!(pair_counts, node_ids_stack, scatt_up, scatt_down, node_group)
        end
    end

    return pair_counts, sun_hits, node_ids, node_group
end

function _all_dir_hits_for_scattering(first::FirstOrderResult, sun_hits::Dict{Int,Int}, cfg::LightConfig, node_ids)
    all_hits = Dict{Int,Int}(nid => get(first.hits_per_node, nid, 0) for nid in node_ids)
    if !cfg.general.all_in_turtle
        for (nid, hsun) in sun_hits
            all_hits[nid] = max(0, get(all_hits, nid, 0) - hsun)
        end
    end
    all_hits
end

function _group_optical_coeffs(cfg::LightConfig)
    coeffs = Dict{Tuple{String,String},Dict{String,Float64}}()

    for group_model in values(cfg.models)
        group = group_model.group
        isempty(group) && continue

        isempty(group_model.types) && continue

        for (type_name, type_model) in group_model.types
            interception = type_model.interception
            interception === nothing && continue

            coeff =
                if lowercase(strip(interception.model)) == "virtualsensor"
                    Dict{String,Float64}("PAR" => 0.0, "NIR" => 0.0)
                else
                    props = interception.optical_properties
                    props === nothing && continue
                    Dict{String,Float64}("PAR" => props.par, "NIR" => props.nir)
                end

            coeffs[(group, type_name)] = coeff
            # Group-level fallback for nodes without explicit type metadata.
            if !haskey(coeffs, (group, "*"))
                coeffs[(group, "*")] = coeff
            end
        end
    end

    coeffs
end

struct ScatteringSceneContext
    pair_counts::Dict{Tuple{Int,Int},Int}
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
end

function _scattering_context(scene::SceneGeometry, turtle::TurtleGrid, first::FirstOrderResult, cfg::LightConfig)
    pair_counts, sun_hits, geom_node_ids, node_group = _pair_counts_for_scattering(scene, turtle, cfg)
    group_type_coeffs = _group_optical_coeffs(cfg)

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
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(node_group, nid, get(scene.node_group, nid, ""))
        default_type = g == "pavement" ? "Cobblestone" : ""
        node_type[nid] = get(scene.node_type, nid, default_type)
    end
    return ScatteringSceneContext(pair_counts, all_hits, node_ids, node_group, node_type, group_type_coeffs)
end

function _scattering_backend_from_mode(mode::Symbol)
    mode == :raycast && return RaycastScatteringBackend()
    mode == :links && return LinksScatteringBackend()
    error("Unsupported scattering mode: $mode (supported: :raycast, :links)")
end

function _resolve_scattering_backend(mode::Symbol, backend::Nothing)
    return _scattering_backend_from_mode(mode)
end

function _resolve_scattering_backend(mode::Symbol, backend::ScatteringBackend)
    mode in (:raycast, :links) || error("Unsupported scattering mode: $mode (supported: :raycast, :links)")
    return backend
end

function _resolve_scattering_backend(mode::Symbol, backend)
    error(
        "Unsupported scattering backend selector type: $(typeof(backend)). " *
        "Use `nothing`, `RaycastScatteringBackend()`, or `LinksScatteringBackend()`.",
    )
end

"""
    build_scattering_transfer_graph(scene, turtle, first, cfg; mode=:raycast, backend=nothing)::ScatteringTransferGraph

Build the scattering transfer graph (pair links and per-node hit normalization data)
independently from iterative scattering propagation.
"""
function build_scattering_transfer_graph(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return build_scattering_transfer_graph(scene, turtle, first, cfg, b)
end

function build_scattering_transfer_graph(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::RaycastScatteringBackend,
)
    context = _scattering_context(scene, turtle, first, cfg)
    return ScatteringTransferGraph(
        context.pair_counts,
        context.all_hits,
        context.node_ids,
        context.node_group,
        context.node_type,
        context.group_type_coeffs,
    )
end

function build_scattering_transfer_graph(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::LinksScatteringBackend,
)
    # CPU reference currently uses the same transfer-graph construction for both modes.
    return build_scattering_transfer_graph(scene, turtle, first, cfg, RaycastScatteringBackend())
end

function _default_band_coeff(cfg::LightConfig, band_key::String)
    bk = uppercase(band_key)
    bk == "NIR" && return cfg.general.scattering_coeff_nir
    return cfg.general.scattering_coeff_par
end

function _coeff_by_node(
    graph::ScatteringTransferGraph,
    band_key::String,
    default_coeff::Float64,
)
    coeff_by_node = Dict{Int,Float64}()
    band = uppercase(band_key)
    for nid in graph.node_ids
        g = get(graph.node_group, nid, "")
        t = get(graph.node_type, nid, "")
        c = get(
            graph.group_type_coeffs,
            (g, t),
            get(graph.group_type_coeffs, (g, "*"), Dict{String,Float64}()),
        )
        coeff_by_node[nid] = get(c, band, default_coeff)
    end
    coeff_by_node
end

function _initial_scattering_power(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}},
    band::AbstractString,
)
    if initial_power_per_node !== nothing
        return _copy_node_values(initial_power_per_node, graph.node_ids)
    end
    source = uppercase(String(band)) == "NIR" ? first.incident_nir_power_per_node : first.incident_par_power_per_node
    return _copy_node_values(source, graph.node_ids)
end

function _propagate_scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    cfg::LightConfig,
    default_coeff::Float64,
)
    node_ids = graph.node_ids
    pair_counts = graph.pair_counts
    all_hits = graph.all_hits
    current = _copy_node_values(initial_power_per_node, node_ids)
    added = _dict_zero(node_ids)
    ref = _sum_dict_values(current)
    thr = cfg.general.scattering_stop_ratio * max(ref, eps(Float64))
    iterations = 0

    converged = false
    for it in 1:cfg.general.scattering_max_iter
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

        if total_next < thr
            converged = true
            break
        end
    end
    return added, iterations, converged
end

function _scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    cfg::LightConfig,
    band_key::String,
    default_coeff::Float64,
)
    coeffs = _coeff_by_node(graph, band_key, default_coeff)
    return _propagate_scattering_one_band(initial_power_per_node, graph, coeffs, cfg, default_coeff)
end

"""
    compute_scattering_band(graph, first, cfg; mode=:raycast, backend=nothing, band="PAR", initial_power_per_node=nothing, default_coeff=nothing)

Run one-band iterative scattering from a pre-built transfer graph.
"""
function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return compute_scattering_band(
        graph,
        first,
        cfg,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::RaycastScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    initial = _initial_scattering_power(graph, first, initial_power_per_node, band)
    dflt = isnothing(default_coeff) ? _default_band_coeff(cfg, String(band)) : default_coeff
    added, iterations, converged = _scattering_one_band(
        initial,
        graph,
        cfg,
        String(band),
        dflt,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged)
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::LinksScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering_band(
        graph,
        first,
        cfg,
        RaycastScatteringBackend();
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, turtle, first, cfg, b)
    return compute_scattering_band(
        graph,
        first,
        cfg,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

"""
    compute_scattering(scene, turtle, first, cfg; mode=:raycast, backend=nothing)::ScatteringResult
    compute_scattering(graph, first, cfg; mode=:raycast, backend=nothing)::ScatteringResult

Compute iterative multiple scattering for PAR and NIR from first-order incident power.
When a `ScatteringTransferGraph` is provided, transfer-link construction is skipped and only
iterative propagation is run.
"""
function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return compute_scattering(graph, first, cfg, b)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::RaycastScatteringBackend,
)
    initial_par = _initial_scattering_power(graph, first, nothing, "PAR")
    initial_nir = _initial_scattering_power(graph, first, nothing, "NIR")

    added_par, it_par, conv_par = _scattering_one_band(
        initial_par,
        graph,
        cfg,
        "PAR",
        cfg.general.scattering_coeff_par,
    )
    added_nir, it_nir, conv_nir = _scattering_one_band(
        initial_nir,
        graph,
        cfg,
        "NIR",
        cfg.general.scattering_coeff_nir,
    )

    ScatteringResult(added_par, added_nir, max(it_par, it_nir), conv_par && conv_nir)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    cfg::LightConfig,
    ::LinksScatteringBackend,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering(graph, first, cfg, RaycastScatteringBackend())
end

function compute_scattering(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, turtle, first, cfg, b)
    return compute_scattering(graph, first, cfg, b)
end
