using ArchimedLight
using KernelAbstractions
using Printf

const STACK_DIAG_PAVING = parse(Int, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_COFFEE_PAVING", "2500"))
const STACK_DIAG_STEP = parse(Int, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_STEP", "1"))
const STACK_DIAG_SECTOR = lowercase(get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_SECTOR", "first_scattering"))
const STACK_DIAG_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_MAX_HITS", "64"))
const STACK_DIAG_HIT_EPSILON = parse(Float32, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_HIT_EPSILON", "-1.0"))
const STACK_DIAG_EXAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_EXAMPLES", "8"))
const STACK_DIAG_STACK_PREVIEW = parse(Int, get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_STACK_PREVIEW", "12"))
const STACK_DIAG_SCAN_ALL =
    lowercase(get(ENV, "ARCHIMEDLIGHT_STACK_DIAG_SCAN_ALL", "false")) in ("1", "true", "yes", "on") ||
    STACK_DIAG_SECTOR in ("all", "summary")

function _coffee_fixture()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config; plot_paving_override=STACK_DIAG_PAVING)
    options = ArchimedLight.LightOptions(options; cache_radiation=true, meteo_range=nothing, scattering=true)
    return scene, models, meteo, options
end

function _sector_index(turtle)
    if STACK_DIAG_SECTOR in ("first_scattering", "all", "summary")
        idx = findfirst(sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0, turtle.sectors)
        idx === nothing && error("No non-sun upward sector found in turtle.")
        return idx
    end
    return parse(Int, STACK_DIAG_SECTOR)
end

function _raster_stack(pixel_hits, pixel_idx)
    stack = get(pixel_hits, pixel_idx, nothing)
    stack === nothing && return Int[], Float64[]
    ArchimedLight._sort_hit_stack!(stack)
    n = length(stack)
    nodes = Vector{Int}(undef, n)
    heights = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        nodes[i] = ArchimedLight._stack_hit_node(stack, i)
        heights[i] = ArchimedLight._stack_hit_height(stack, i)
    end
    return nodes, heights
end

function _raycore_stack(traced, pixel_idx)
    n = Int(traced.counts[pixel_idx])
    n == 0 && return Int[], Float64[]
    nodes = Vector{Int}(undef, n)
    heights = Vector{Float64}(undef, n)
    offset = (pixel_idx - 1) * traced.max_hits
    @inbounds for i in 1:n
        nodes[i] = Int(traced.nodes[offset + i])
        heights[i] = Float64(traced.heights[offset + i])
    end
    return nodes, heights
end

function _preview(values)
    limit = min(length(values), STACK_DIAG_STACK_PREVIEW)
    suffix = length(values) > limit ? ", ..." : ""
    return "[" * join(values[1:limit], ", ") * suffix * "]"
end

function _node_label(geometry, node_idx::Int)
    if 1 <= node_idx <= length(geometry.node_ids)
        node_id = geometry.node_ids[node_idx]
        group = get(geometry.node_group, node_id, "")
        type_name = get(geometry.node_type, node_id, "")
        return "$(node_id) ($(group)/$(type_name))"
    end
    return "invalid-index-$node_idx"
end

function _stack_length(pixel_hits, pixel_idx)
    stack = get(pixel_hits, pixel_idx, nothing)
    return stack === nothing ? 0 : length(stack)
end

function _sector_count_summary(prepared, raycore_data, sector, options)
    geometry = prepared.geometry
    n_pixels = geometry.plotbox.nx * geometry.plotbox.ny
    raster_projection = ArchimedLight._prepared_direction_projection(prepared, sector.direction, options)
    traced = ArchimedLight._raycore_trace_direction_stacks(raycore_data, sector.direction, options)
    raster_total = 0
    raster_occupied = 0
    equal_counts = 0
    raycore_less = 0
    raycore_more = 0
    @inbounds for pixel_idx in 1:n_pixels
        nr = _stack_length(raster_projection.pixel_hits, pixel_idx)
        nk = Int(traced.counts[pixel_idx])
        raster_total += nr
        nr > 0 && (raster_occupied += 1)
        if nr == nk
            equal_counts += 1
        elseif nk < nr
            raycore_less += 1
        else
            raycore_more += 1
        end
    end
    raycore_total = sum(x -> Int(x), traced.counts; init=0)
    raycore_occupied = count(!=(Int32(0)), traced.counts)
    return (
        raster_total=raster_total,
        raycore_total=raycore_total,
        raster_occupied=raster_occupied,
        raycore_occupied=raycore_occupied,
        equal_counts=equal_counts,
        raycore_less=raycore_less,
        raycore_more=raycore_more,
        overflow=count(traced.overflow),
    )
end

function _print_top_node_deltas(label, deltas, geometry; positive::Bool)
    order = sortperm(deltas; rev=positive)
    printed = 0
    println(label)
    for idx in order
        delta = deltas[idx]
        positive ? (delta > 0 || break) : (delta < 0 || break)
        @printf("  %+8d  node_index=%-6d node=%s\n", delta, idx, _node_label(geometry, idx))
        printed += 1
        printed >= STACK_DIAG_EXAMPLES && break
    end
    printed == 0 && println("  none")
end

function main()
    scene, models, meteo, options = _coffee_fixture()
    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_budget_maps=true)
    geometry = prepared.geometry

    rows = ArchimedLight.prepare_meteo(meteo, options).rows
    1 <= STACK_DIAG_STEP <= length(rows) ||
        error("ARCHIMEDLIGHT_STACK_DIAG_STEP=$(STACK_DIAG_STEP) is outside 1:$(length(rows)).")
    row = rows[STACK_DIAG_STEP]
    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    sector_idx = _sector_index(turtle)
    1 <= sector_idx <= length(turtle.sectors) ||
        error("ARCHIMEDLIGHT_STACK_DIAG_SECTOR=$(STACK_DIAG_SECTOR) resolved outside 1:$(length(turtle.sectors)).")
    sector = turtle.sectors[sector_idx]

    ib = ArchimedLight.RaycoreInterceptionBackend(
        backend=KernelAbstractions.CPU(),
        max_hits_per_pixel=STACK_DIAG_MAX_HITS,
        hit_epsilon=STACK_DIAG_HIT_EPSILON,
    )
    raycore_data, was_chunked =
        ArchimedLight._raycore_initial_scene_data(prepared, ib.config, options; toricity=options.toricity)

    raster_projection = ArchimedLight._prepared_direction_projection(prepared, sector.direction, options)
    traced = ArchimedLight._raycore_trace_direction_stacks(raycore_data, sector.direction, options)

    n_pixels = geometry.plotbox.nx * geometry.plotbox.ny
    raster_total = 0
    raycore_total = 0
    raster_occupied = 0
    raycore_occupied = 0
    equal_counts = 0
    raycore_less = 0
    raycore_more = 0
    same_nodes = 0
    mismatched_nodes = 0
    examples = NamedTuple[]
    raster_node_counts = zeros(Int, length(geometry.node_ids))
    raycore_node_counts = zeros(Int, length(geometry.node_ids))

    @inbounds for pixel_idx in 1:n_pixels
        raster_nodes, raster_heights = _raster_stack(raster_projection.pixel_hits, pixel_idx)
        raycore_nodes, raycore_heights = _raycore_stack(traced, pixel_idx)
        nr = length(raster_nodes)
        nk = length(raycore_nodes)
        raster_total += nr
        raycore_total += nk
        nr > 0 && (raster_occupied += 1)
        nk > 0 && (raycore_occupied += 1)
        if nr == nk
            equal_counts += 1
        elseif nk < nr
            raycore_less += 1
        else
            raycore_more += 1
        end
        if raster_nodes == raycore_nodes
            same_nodes += 1
        else
            mismatched_nodes += 1
            if length(examples) < STACK_DIAG_EXAMPLES
                push!(
                    examples,
                    (
                        pixel=pixel_idx,
                        raster_nodes=copy(raster_nodes),
                        raycore_nodes=copy(raycore_nodes),
                        raster_heights=copy(raster_heights),
                        raycore_heights=copy(raycore_heights),
                    ),
                )
            end
        end
        for idx in raster_nodes
            1 <= idx <= length(raster_node_counts) && (raster_node_counts[idx] += 1)
        end
        for idx in raycore_nodes
            1 <= idx <= length(raycore_node_counts) && (raycore_node_counts[idx] += 1)
        end
    end

    overflow_count = count(traced.overflow)

    println("Raycore stack parity diagnostic")
    println("workload: coffee")
    println("paving: ", STACK_DIAG_PAVING)
    println("step: ", STACK_DIAG_STEP, " / ", length(rows))
    println("sector: ", sector_idx, " id=", sector.id, " source=", sector.source, " direction=", sector.direction)
    println("raycore max_hits_per_pixel: ", STACK_DIAG_MAX_HITS)
    println("raycore hit_epsilon: ", STACK_DIAG_HIT_EPSILON)
    println("raycore scene chunked: ", was_chunked)
    println("production projection uses raster compat: ", ArchimedLight._raycore_use_raster_compat_projection(raycore_data, sector.direction, options))
    println()
    @printf("%-28s %12s %12s\n", "metric", "raster", "raycore")
    @printf("%-28s %12d %12d\n", "occupied pixels", raster_occupied, raycore_occupied)
    @printf("%-28s %12d %12d\n", "total stack hits", raster_total, raycore_total)
    @printf("%-28s %12d %12d\n", "max stack capacity", 0, n_pixels * STACK_DIAG_MAX_HITS)
    println("raycore overflow pixels: ", overflow_count)
    println("equal count pixels: ", equal_counts)
    println("raycore count lower pixels: ", raycore_less)
    println("raycore count higher pixels: ", raycore_more)
    println("identical node-stack pixels: ", same_nodes)
    println("mismatched node-stack pixels: ", mismatched_nodes)
    println()

    for ex in examples
        println("pixel ", ex.pixel)
        println("  raster count=", length(ex.raster_nodes), " nodes=", _preview(ex.raster_nodes))
        println("  raycore count=", length(ex.raycore_nodes), " nodes=", _preview(ex.raycore_nodes))
        println("  raster heights=", _preview(round.(ex.raster_heights; digits=4)))
        println("  raycore heights=", _preview(round.(ex.raycore_heights; digits=4)))
    end
    println()

    deltas = raster_node_counts .- raycore_node_counts
    _print_top_node_deltas("Top nodes with raster > raycore hit counts", deltas, geometry; positive=true)
    _print_top_node_deltas("Top nodes with raycore > raster hit counts", deltas, geometry; positive=false)

    if STACK_DIAG_SCAN_ALL
        println()
        println("Per-sector stack-count summary")
        @printf(
            "%6s %8s %12s %12s %10s %10s %10s %10s\n",
            "sector",
            "source",
            "raster",
            "raycore",
            "ratio",
            "lower",
            "higher",
            "overflow",
        )
        worst_idx = 0
        worst_ratio = Inf
        for (idx, candidate) in pairs(turtle.sectors)
            candidate.source == :sun && continue
            Float32(candidate.direction[3]) > 0.0f0 || continue
            summary = _sector_count_summary(prepared, raycore_data, candidate, options)
            ratio =
                summary.raster_total == 0 ?
                (summary.raycore_total == 0 ? 1.0 : Inf) :
                summary.raycore_total / summary.raster_total
            if ratio < worst_ratio
                worst_ratio = ratio
                worst_idx = idx
            end
            @printf(
                "%6d %8s %12d %12d %9.3f%% %10d %10d %10d\n",
                idx,
                string(candidate.source),
                summary.raster_total,
                summary.raycore_total,
                100 * ratio,
                summary.raycore_less,
                summary.raycore_more,
                summary.overflow,
            )
        end
        worst_idx == 0 || println("lowest raycore/raster stack-hit ratio sector: ", worst_idx, " (", round(100 * worst_ratio; digits=3), "%)")
    end
end

main()
