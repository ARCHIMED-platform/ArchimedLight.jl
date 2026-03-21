function _dict_zero(node_ids)
    Dict{Int,Float64}(nid => 0.0 for nid in node_ids)
end

struct ScatteringTopologyCache
    pair_counts::ScatteringPairCounts
    sun_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
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

"""
    _pack_scattering_edge(to, from) -> UInt64

Pack one scattering edge `(to, from)` into a single 64-bit key.

We reserve 32 bits for each node id:

    [to ........ 32 bits][from ...... 32 bits]

The explicit `UInt32` conversion constrains each id to one half of the packed key, and the
`UInt64` widening makes the shift/OR operations safe. This is cheaper to hash in the topology
builder than a `Tuple{Int,Int}` key.
"""
@inline function _pack_scattering_edge(to::Int, from::Int)
    return (UInt64(UInt32(to)) << 32) | UInt64(UInt32(from))
end

# Recover the upper and lower 32-bit halves of a packed `(to, from)` edge key.
@inline _unpack_scattering_to(edge::UInt64) = Int(UInt32(edge >> 32))
@inline _unpack_scattering_from(edge::UInt64) = Int(UInt32(edge & 0xffffffff))

function _edge_counts_from_packed(edge_counts::Dict{UInt64,Int})
    isempty(edge_counts) && return ScatteringPairCounts(Int[], Int[], Int[])

    # Sort once when materializing the final graph so iteration order is deterministic.
    packed_edges = sort!(collect(keys(edge_counts)))
    to_nodes = Int[]
    from_nodes = Int[]
    counts = Int[]
    sizehint!(to_nodes, length(packed_edges))
    sizehint!(from_nodes, length(packed_edges))
    sizehint!(counts, length(packed_edges))

    for edge in packed_edges
        count = edge_counts[edge]
        if count > 0
            push!(to_nodes, _unpack_scattering_to(edge))
            push!(from_nodes, _unpack_scattering_from(edge))
            push!(counts, count)
        end
    end
    return ScatteringPairCounts(to_nodes, from_nodes, counts)
end

@inline function _add_packed_edge_count!(edge_counts::Dict{UInt64,Int}, to::Int, from::Int)
    edge = _pack_scattering_edge(to, from)
    edge_counts[edge] = get(edge_counts, edge, 0) + 1
    return nothing
end

function _accept_scattering_link(node_group::Dict{Int,String}, to::Int, from::Int)
    !(get(node_group, to, "") == "pavement" && get(node_group, from, "") == "pavement")
end

mutable struct ScatteringStackScratch
    nearest_nonvirtual_below::Vector{Int}
end

ScatteringStackScratch() = ScatteringStackScratch(Int[])

function _stack_transfer_pairs!(
    edge_counts::Dict{UInt64,Int},
    stack,
    virtual_nodes::Set{Int},
    scratch::ScatteringStackScratch,
    node_group::Dict{Int,String},
)
    n_hits = length(stack)
    below = scratch.nearest_nonvirtual_below
    resize!(below, n_hits)

    nearest_below = 0
    @inbounds for h in n_hits:-1:1
        nid = _stack_hit_node(stack, h)
        if !(nid in virtual_nodes)
            nearest_below = nid
        end
        below[h] = nearest_below
    end

    nearest_above = 0
    @inbounds for h in 1:(n_hits - 1)
        to_above = _stack_hit_node(stack, h)
        if !(to_above in virtual_nodes)
            nearest_above = to_above
        end

        from_below = below[h + 1]
        if from_below != 0 && _accept_scattering_link(node_group, to_above, from_below)
            _add_packed_edge_count!(edge_counts, to_above, from_below)
        end

        to_below = _stack_hit_node(stack, h + 1)
        from_above = nearest_above
        if from_above != 0 && _accept_scattering_link(node_group, to_below, from_above)
            _add_packed_edge_count!(edge_counts, to_below, from_above)
        end
    end
    return nothing
end

function _accumulate_scattering_counts!(
    edge_counts::Dict{UInt64,Int},
    sun_hits::Dict{Int,Int},
    sector::TurtleSector,
    projection::DirectionProjectionResult,
    virtual_nodes::Set{Int},
    node_group::Dict{Int,String},
    scratch::ScatteringStackScratch;
    stacks_sorted::Bool=false,
)
    if sector.source == :sun
        for (nid, h) in projection.node_hits
            sun_hits[nid] = get(sun_hits, nid, 0) + h
        end
        return nothing
    end

    for stack in values(projection.pixel_hits)
        length(stack) <= 1 && continue
        stacks_sorted || _sort_hit_stack!(stack)
        _stack_transfer_pairs!(edge_counts, stack, virtual_nodes, scratch, node_group)
    end
    return nothing
end

function _pair_counts_for_scattering(scene::SceneGeometry, models::LightModels, turtle::TurtleGrid, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    virtual_nodes = _virtual_sensor_node_ids(geometry.node_group, geometry.node_type, models)
    cache_ctx = _projection_cache_context(geometry.vertices, geometry.faces, geometry.face2node, geometry.plotbox, options)
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for sector in turtle.sectors
        projection =
            _direction_projection_cached(geometry.vertices, geometry.faces, geometry.face2node, sector.direction, options, geometry.plotbox, cache_ctx)
        _accumulate_scattering_counts!(
            edge_counts,
            sun_hits,
            sector,
            projection,
            virtual_nodes,
            geometry.node_group,
            scratch,
        )
    end

    return _edge_counts_from_packed(edge_counts), sun_hits, geometry.node_ids, geometry.node_group
end

function _pair_counts_from_projections(
    turtle::TurtleGrid,
    projections::AbstractVector{DirectionProjectionResult},
    virtual_nodes::Set{Int},
    node_group::Dict{Int,String},
    stacks_sorted::Bool=false,
)
    edge_counts = Dict{UInt64,Int}()
    sun_hits = Dict{Int,Int}()
    scratch = ScatteringStackScratch()

    for i in eachindex(turtle.sectors)
        _accumulate_scattering_counts!(
            edge_counts,
            sun_hits,
            turtle.sectors[i],
            projections[i],
            virtual_nodes,
            node_group,
            scratch;
            stacks_sorted=stacks_sorted,
        )
    end

    return _edge_counts_from_packed(edge_counts), sun_hits
end

function _all_dir_hits_for_scattering(first::FirstOrderResult, sun_hits::Dict{Int,Int}, options::LightOptions, node_ids)
    all_hits = Dict{Int,Int}(nid => get(first.hits_per_node, nid, 0) for nid in node_ids)
    if !options.all_in_turtle
        for (nid, hsun) in sun_hits
            all_hits[nid] = max(0, get(all_hits, nid, 0) - hsun)
        end
    end
    all_hits
end

function _group_optical_coeffs(models::LightModels)
    coeffs = Dict{Tuple{String,String},Dict{String,Float64}}()

    for group_model in values(models)
        group = group_model.group
        isempty(group) && continue

        isempty(group_model.types) && continue

        for (type_name, type_model) in group_model.types
            interception = type_model.interception
            interception === nothing && continue

            coeff =
                if _is_sensor_interception(interception)
                    Dict{String,Float64}("PAR" => 0.0, "NIR" => 0.0)
                else
                    props = interception.optical_properties
                    props === nothing && continue
                    band_coeffs = Dict{String,Float64}()
                    has_par = get(props.extras, "__has_par", true)
                    has_nir = get(props.extras, "__has_nir", true)
                    has_par && (band_coeffs["PAR"] = props.par)
                    has_nir && (band_coeffs["NIR"] = props.nir)
                    band_coeffs
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
    pair_counts::ScatteringPairCounts
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
end

function _build_scattering_topology_cache(
    scene::SceneGeometry,
    models::LightModels,
    prepared::PreparedInterceptionData,
    turtle::TurtleGrid,
    projections::AbstractVector{DirectionProjectionResult},
    stacks_sorted::Bool=false,
)
    pair_counts, sun_hits = _pair_counts_from_projections(
        turtle,
        projections,
        prepared.virtual_nodes,
        prepared.geometry.node_group,
        stacks_sorted,
    )
    group_type_coeffs = _group_optical_coeffs(models)
    node_ids = copy(prepared.geometry.node_ids)
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(prepared.geometry.node_group, nid, _scene_group(scene, nid, ""))
        default_type = g == "pavement" ? "Cobblestone" : ""
        node_type[nid] = _scene_type(scene, nid, default_type)
    end
    return ScatteringTopologyCache(
        pair_counts,
        sun_hits,
        node_ids,
        prepared.geometry.node_group,
        node_type,
        group_type_coeffs,
    )
end

function _node_ids_for_scattering(topology::ScatteringTopologyCache, first::FirstOrderResult)
    node_set = Set{Int}()
    for nid in topology.node_ids
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.par)
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.nir)
        push!(node_set, nid)
    end
    return collect(node_set)
end

function _transfer_graph_from_topology(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
)
    node_ids = _node_ids_for_scattering(topology, first)
    all_hits = _all_dir_hits_for_scattering(first, topology.sun_hits, options, node_ids)
    return ScatteringTransferGraph(
        topology.pair_counts,
        all_hits,
        node_ids,
        topology.node_group,
        topology.node_type,
        topology.group_type_coeffs,
    )
end

function _scattering_context(scene::SceneGeometry, models::LightModels, turtle::TurtleGrid, first::FirstOrderResult, options::LightOptions)
    pair_counts, sun_hits, geom_node_ids, node_group = _pair_counts_for_scattering(scene, models, turtle, options)
    group_type_coeffs = _group_optical_coeffs(models)

    node_set = Set{Int}()
    for nid in geom_node_ids
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.par)
        push!(node_set, nid)
    end
    for nid in keys(first.incident_power.nir)
        push!(node_set, nid)
    end
    node_ids = collect(node_set)
    all_hits = _all_dir_hits_for_scattering(first, sun_hits, options, node_ids)
    node_type = Dict{Int,String}()
    for nid in node_ids
        g = get(node_group, nid, _scene_group(scene, nid, ""))
        default_type = g == "pavement" ? "Cobblestone" : ""
        node_type[nid] = _scene_type(scene, nid, default_type)
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
    build_scattering_transfer_graph(scene, models, turtle, first, options; mode=:raycast, backend=nothing)::ScatteringTransferGraph

Build the scattering transfer graph (pair links and per-node hit normalization data)
independently from iterative scattering propagation.
"""
function build_scattering_transfer_graph(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return build_scattering_transfer_graph(scene, models, turtle, first, options, b)
end

function build_scattering_transfer_graph(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    prepared = _prepare_interception_data(scene, models, options)
    projections = _build_direction_projections(prepared, turtle, options)
    topology = _build_scattering_topology_cache(scene, models, prepared, turtle, projections)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    # CPU reference currently uses the same transfer-graph construction for both modes.
    return build_scattering_transfer_graph(scene, models, turtle, first, options, RaycastScatteringBackend())
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return build_scattering_transfer_graph(topology, first, options, b)
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    return _transfer_graph_from_topology(topology, first, options)
end

function build_scattering_transfer_graph(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    return build_scattering_transfer_graph(topology, first, options, RaycastScatteringBackend())
end

function _default_band_coeff(options::LightOptions, band_key::String)
    bk = uppercase(band_key)
    bk == "NIR" && return options.scattering_coeff_nir
    return options.scattering_coeff_par
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
    source = uppercase(String(band)) == "NIR" ? first.incident_power.nir : first.incident_power.par
    return _copy_node_values(source, graph.node_ids)
end

function _propagate_scattering_one_band(
    initial_power_per_node::Dict{Int,Float64},
    graph::ScatteringTransferGraph,
    coeff_by_node::Dict{Int,Float64},
    options::LightOptions,
    default_coeff::Float64,
)
    node_ids = graph.node_ids
    pair_counts = graph.pair_counts
    all_hits = graph.all_hits
    current = _copy_node_values(initial_power_per_node, node_ids)
    added = _dict_zero(node_ids)
    ref = _sum_dict_values(current)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    iterations = 0

    converged = false
    for it in 1:options.scattering_max_iter
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
    options::LightOptions,
    band_key::String,
    default_coeff::Float64,
)
    coeffs = _coeff_by_node(graph, band_key, default_coeff)
    return _propagate_scattering_one_band(initial_power_per_node, graph, coeffs, options, default_coeff)
end

"""
    compute_scattering_band(graph, first, options; mode=:raycast, backend=nothing, band="PAR", initial_power_per_node=nothing, default_coeff=nothing)

Run one-band iterative scattering from a pre-built transfer graph.
"""
function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions;
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
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    initial = _initial_scattering_power(graph, first, initial_power_per_node, band)
    dflt = isnothing(default_coeff) ? _default_band_coeff(options, String(band)) : default_coeff
    added, iterations, converged = _scattering_one_band(
        initial,
        graph,
        options,
        String(band),
        dflt,
    )
    return (added_power_per_node=added, iterations=iterations, converged=converged)
end

function compute_scattering_band(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend;
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering_band(
        graph,
        first,
        options,
        RaycastScatteringBackend();
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, models, turtle, first, options, b)
    return compute_scattering_band(
        graph,
        first,
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

function compute_scattering_band(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
    band::AbstractString="PAR",
    initial_power_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    default_coeff::Union{Nothing,Float64}=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(topology, first, options, b)
    return compute_scattering_band(
        graph,
        first,
        options,
        b;
        band=band,
        initial_power_per_node=initial_power_per_node,
        default_coeff=default_coeff,
    )
end

"""
    compute_scattering(scene, models, turtle, first, options; mode=:raycast, backend=nothing)::ScatteringResult
    compute_scattering(graph, first, options; mode=:raycast, backend=nothing)::ScatteringResult

Compute iterative multiple scattering for PAR and NIR from first-order incident power.
When a `ScatteringTransferGraph` is provided, transfer-link construction is skipped and only
iterative propagation is run.
"""
function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    return compute_scattering(graph, first, options, b)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::RaycastScatteringBackend,
)
    initial_par = _initial_scattering_power(graph, first, nothing, "PAR")
    initial_nir = _initial_scattering_power(graph, first, nothing, "NIR")

    added_par, it_par, conv_par = _scattering_one_band(
        initial_par,
        graph,
        options,
        "PAR",
        options.scattering_coeff_par,
    )
    added_nir, it_nir, conv_nir = _scattering_one_band(
        initial_nir,
        graph,
        options,
        "NIR",
        options.scattering_coeff_nir,
    )

    ScatteringResult(SpectralNodeValues(added_par, added_nir), max(it_par, it_nir), conv_par && conv_nir)
end

function compute_scattering(
    graph::ScatteringTransferGraph,
    first::FirstOrderResult,
    options::LightOptions,
    ::LinksScatteringBackend,
)
    # CPU reference currently shares the same iterative propagation path.
    return compute_scattering(graph, first, options, RaycastScatteringBackend())
end

function compute_scattering(
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(scene, models, turtle, first, options, b)
    return compute_scattering(graph, first, options, b)
end

function compute_scattering(
    topology::ScatteringTopologyCache,
    first::FirstOrderResult,
    options::LightOptions;
    mode=:raycast,
    backend=nothing,
)
    b = _resolve_scattering_backend(mode, backend)
    graph = build_scattering_transfer_graph(topology, first, options, b)
    return compute_scattering(graph, first, options, b)
end
