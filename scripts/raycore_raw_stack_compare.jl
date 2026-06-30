if lowercase(get(ENV, "ARCHIMEDLIGHT_RAW_STACK_BACKEND", "metal")) == "metal" &&
   Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using GeometryBasics
using KernelAbstractions
using MultiScaleTreeGraph
using PlantGeom
using Printf

const RAW_WORKLOAD = lowercase(get(ENV, "ARCHIMEDLIGHT_RAW_STACK_WORKLOAD", "coffee"))
const RAW_BACKEND = lowercase(get(ENV, "ARCHIMEDLIGHT_RAW_STACK_BACKEND", "metal"))
const RAW_MAX_HITS = [
    parse(Int, strip(value)) for value in split(get(ENV, "ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP", "32,64,128,256"), ',')
    if !isempty(strip(value))
]
const RAW_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_WORKGROUPSIZE", "256"))
const RAW_HIT_EPSILON = parse(Float32, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_HIT_EPSILON", "1.0e-4"))
const RAW_DISTANCE_ATOL = parse(Float32, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_DISTANCE_ATOL", "1.0e-4"))
const RAW_SECTOR = lowercase(get(ENV, "ARCHIMEDLIGHT_RAW_STACK_SECTOR", "first_scattering"))
const RAW_STEP = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_STEP", "1"))
const RAW_PIXEL_SIZE = get(ENV, "ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE", "")
const RAW_COFFEE_PAVING = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_COFFEE_PAVING", "2500"))
const RAW_SYNTHETIC_STACK_LAYERS = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_LAYERS", "384"))
const RAW_SYNTHETIC_STACK_SPACING = parse(Float64, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_SPACING", "0.002"))
const RAW_FACE_CHUNK_LIMITS_RAW = get(ENV, "ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS", "auto")
const RAW_OUTPUT = get(ENV, "ARCHIMEDLIGHT_RAW_STACK_OUTPUT", "")
const RAYCORE_INVALID_NODE = ArchimedLight.RAYCORE_INVALID_NODE
const RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY = ArchimedLight._raycore_software_traversal_stack_capacity()

isempty(RAW_MAX_HITS) && error("ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP selected no hit-buffer sizes.")
RAW_SYNTHETIC_STACK_LAYERS > 0 || error("ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_LAYERS must be positive.")
RAW_SYNTHETIC_STACK_SPACING > 0 || error("ARCHIMEDLIGHT_RAW_STACK_SYNTHETIC_SPACING must be positive.")

function _raw_face_chunk_limits()
    limits = Union{Symbol,Int}[]
    for raw in split(RAW_FACE_CHUNK_LIMITS_RAW, ',')
        value = lowercase(strip(raw))
        isempty(value) && continue
        if value in ("auto", "default")
            push!(limits, :auto)
        elseif value in ("unchunked", "none", "merged", "merged_mesh")
            push!(limits, :unchunked)
        else
            limit = parse(Int, value)
            limit > 0 || error("ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS values must be positive integers, auto, or unchunked.")
            push!(limits, limit)
        end
    end
    isempty(limits) && error("ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS selected no chunk modes.")
    return limits
end

_raw_chunk_mode_label(mode::Symbol) = string(mode)
_raw_chunk_mode_label(mode::Int) = string(mode)

function _load_optional_package(name::Symbol)
    Base.find_package(String(name)) === nothing &&
        error("Package $name is not available in the active environment.")
    return Base.require(Main, name)
end

function _backend_from_env()
    RAW_BACKEND == "cpu" && return KernelAbstractions.CPU()
    if RAW_BACKEND == "metal"
        mod = _load_optional_package(:Metal)
        array_type = getproperty(mod, :MtlArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif RAW_BACKEND == "cuda"
        mod = _load_optional_package(:CUDA)
        functional = getproperty(mod, :functional)
        Base.invokelatest(functional) || error("CUDA is available but CUDA.functional() is false.")
        array_type = getproperty(mod, :CuArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif RAW_BACKEND == "oneapi"
        mod = _load_optional_package(:oneAPI)
        array_type = getproperty(mod, :oneArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif RAW_BACKEND == "amdgpu"
        mod = _load_optional_package(:AMDGPU)
        array_type = getproperty(mod, :ROCArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    end
    error("Unsupported ARCHIMEDLIGHT_RAW_STACK_BACKEND=$RAW_BACKEND. Use cpu, metal, cuda, oneapi, or amdgpu.")
end

function _with_common_options(options)
    kwargs = (; cache_radiation=true, scattering=true, meteo_range=nothing)
    if !isempty(strip(RAW_PIXEL_SIZE))
        kwargs = merge(kwargs, (; pixel_size=parse(Float64, RAW_PIXEL_SIZE)))
    end
    return ArchimedLight.LightOptions(options; kwargs...)
end

function _panel_mesh(width, length, height, inclination_deg)
    angle = deg2rad(inclination_deg)
    yhalf = 0.5f0 * Float32(length * cos(angle))
    zhalf = 0.5f0 * Float32(length * sin(angle))
    w = Float32(width)
    h = Float32(height)
    points = GeometryBasics.Point3f[
        GeometryBasics.Point3f(0, -yhalf, h - zhalf),
        GeometryBasics.Point3f(w, -yhalf, h - zhalf),
        GeometryBasics.Point3f(w, yhalf, h + zhalf),
        GeometryBasics.Point3f(0, yhalf, h + zhalf),
    ]
    return GeometryBasics.Mesh(points, GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)])
end

function _load_simple_workload()
    config = joinpath(@__DIR__, "..", "test", "fast_fixtures", "simpleplant_16_toric", "input", "config.yml")
    options, scene, _meteo, models = ArchimedLight.read_config(config)
    options = _with_common_options(options)
    sky = ArchimedLight.SkyState(180.0, 55.0, 420.0, 220.0, 0.70, 0.30)
    return (; name="simple", scene, models, options, sky)
end

function _load_coffee_workload()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config; plot_paving_override=RAW_COFFEE_PAVING)
    options = _with_common_options(options)
    rows = ArchimedLight.prepare_meteo(meteo, options).rows
    1 <= RAW_STEP <= length(rows) ||
        error("ARCHIMEDLIGHT_RAW_STACK_STEP=$RAW_STEP is outside 1:$(length(rows)).")
    sky = ArchimedLight.compute_sky(rows[RAW_STEP], options)
    return (; name="coffee", scene, models, options, sky)
end

function _load_agripv_workload()
    root = get(
        ENV,
        "ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT",
        "/Users/rvezy/Documents/cirad/Articles_Rapports_Communications/Conférences/2026_fspm/coutellier_agripv",
    )
    opf_path = joinpath(root, "0_simulations", "archicrop", "wheat", "plant_1995-06-24.opf")
    isfile(opf_path) || error("Agrivoltaics wheat OPF not found at $opf_path")

    n_rows = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_ROWS", "2"))
    plants_per_row = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_PLANTS_PER_ROW", "8"))
    interrow = parse(Float64, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_INTERROW", "0.25"))
    intrarow = parse(Float64, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_INTRAROW", "0.30"))
    ground_nx = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX", "18"))
    ground_ny = parse(Int, get(ENV, "ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY", "24"))
    ground_nx > 0 || error("ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NX must be positive.")
    ground_ny > 0 || error("ARCHIMEDLIGHT_RAW_STACK_AGRIPV_GROUND_NY must be positive.")
    width = max(interrow * n_rows, interrow)
    scene_length = max(intrarow * plants_per_row, intrarow)

    wheat = PlantGeom.read_opf(opf_path, mtg_type=MultiScaleTreeGraph.NodeMTG)
    panel = _panel_mesh(width, 3.0, 1.6, 25.0)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, width, scene_length)) do s
        PlantGeom.add_object!(s, panel; group="panel", type="Panel", id=1, at=(0.0, scene_length / 2, 0.0))
        id = 2
        for row in 0:(n_rows - 1), col in 0:(plants_per_row - 1)
            PlantGeom.add_plant!(
                s,
                wheat;
                group="wheat",
                id=id,
                at=((row + 0.5) * interrow, (col + 0.5) * intrarow, 0.0),
                rotate=(z=5.0 * sin(id),),
                deg=true,
            )
            id += 1
        end
        PlantGeom.add_ground!(s; nx=ground_nx, ny=ground_ny, group="pavement", type="Cobblestone")
    end
    models = ArchimedLight.models_for(
        "wheat" => ("*" => ArchimedLight.translucent(par=0.15, nir=0.90),),
        "panel" => ("Panel" => ArchimedLight.translucent(par=0.0, nir=0.0),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
    options = _with_common_options(ArchimedLight.LightOptions(scattering=true, cache_radiation=true, all_in_turtle=true, toricity=true))
    sky = ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30)
    return (; name="agripv", scene, models, options, sky)
end

function _stacked_plate_mesh(layers::Int, spacing::Float64)
    points = GeometryBasics.Point3f[]
    faces = GeometryBasics.TriangleFace{Int}[]
    for layer in 0:(layers - 1)
        z = Float32(1.0 + layer * spacing)
        base = length(points)
        append!(
            points,
            GeometryBasics.Point3f[
                GeometryBasics.Point3f(0, 0, z),
                GeometryBasics.Point3f(1, 0, z),
                GeometryBasics.Point3f(1, 1, z),
                GeometryBasics.Point3f(0, 1, z),
            ],
        )
        append!(faces, GeometryBasics.TriangleFace{Int}[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
    end
    return GeometryBasics.Mesh(points, faces)
end

function _load_synthetic_stack_workload()
    mesh = _stacked_plate_mesh(RAW_SYNTHETIC_STACK_LAYERS, RAW_SYNTHETIC_STACK_SPACING)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, 1.0, 1.0), source_path="raycore_raw_synthetic_stack") do s
        PlantGeom.add_object!(s, mesh; group="synthetic", type="Stack", id=1)
    end
    models = ArchimedLight.models_for(
        "synthetic" => ("Stack" => ArchimedLight.translucent(par=0.15, nir=0.90),),
    )
    options = _with_common_options(
        ArchimedLight.LightOptions(
            scattering=true,
            cache_radiation=true,
            all_in_turtle=true,
            toricity=false,
            pixel_size=0.05,
            turtle_sectors=6,
        ),
    )
    sky = ArchimedLight.SkyState(180.0, 90.0, 420.0, 220.0, 0.70, 0.30)
    return (; name="synthetic_stack", scene, models, options, sky)
end

function _load_workload()
    RAW_WORKLOAD == "simple" && return _load_simple_workload()
    RAW_WORKLOAD == "coffee" && return _load_coffee_workload()
    RAW_WORKLOAD == "agripv" && return _load_agripv_workload()
    RAW_WORKLOAD in ("synthetic_stack", "synthetic_deep_blas") && return _load_synthetic_stack_workload()
    error("Unsupported ARCHIMEDLIGHT_RAW_STACK_WORKLOAD=$RAW_WORKLOAD. Use simple, coffee, agripv, or synthetic_stack.")
end

function _selected_sectors(turtle)
    if RAW_SECTOR in ("all", "summary")
        return [idx for (idx, sector) in pairs(turtle.sectors) if sector.source != :sun && Float32(sector.direction[3]) > 0.0f0]
    elseif RAW_SECTOR == "first_scattering"
        idx = findfirst(sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0, turtle.sectors)
        idx === nothing && error("No non-sun upward sector found in turtle.")
        return [idx]
    end
    idx = parse(Int, RAW_SECTOR)
    1 <= idx <= length(turtle.sectors) || error("ARCHIMEDLIGHT_RAW_STACK_SECTOR=$idx is outside 1:$(length(turtle.sectors)).")
    return [idx]
end

function _mismatch_summary(data, direction, options, mismatch)
    mismatch === nothing && return "none"
    parts = ["kind=$(mismatch.kind)"]
    hasproperty(mismatch, :pixel) && push!(parts, "pixel=$(mismatch.pixel)")
    hasproperty(mismatch, :slot) && push!(parts, "slot=$(mismatch.slot)")
    if hasproperty(mismatch, :pixel)
        far = ArchimedLight._raycore_projection_far(data, direction, options)
        ray = ArchimedLight._raycore_sky_to_ground_ray(data.prepared.geometry.plotbox, mismatch.pixel, direction, far)
        push!(parts, "origin=$(Tuple(ray.o))")
        push!(parts, "direction=$(Tuple(ray.d))")
    end
    for name in (:reference_count, :candidate_count, :reference_overflow, :candidate_overflow,
        :metadata, :reference_instance, :candidate_instance, :reference_distance, :candidate_distance,
        :reference, :candidate, :diagnosis, :unmatched_reference_count, :unmatched_candidate_count,
        :first_unmatched_reference, :first_unmatched_candidate)
        hasproperty(mismatch, name) && push!(parts, "$(name)=$(getproperty(mismatch, name))")
    end
    return join(parts, ";")
end

function _raw_config(backend, max_hits)
    return ArchimedLight.RaycoreBackendConfig(
        ;
        backend=backend,
        max_hits_per_pixel=max_hits,
        workgroupsize=RAW_WORKGROUPSIZE,
        hit_epsilon=RAW_HIT_EPSILON,
        max_prechunk_instances=0,
        allow_fallback=false,
    )
end

function _scene_data_for(prepared, config, options, chunk_mode::Union{Symbol,Int})
    chunk_mode === :auto &&
        return ArchimedLight._raycore_initial_scene_data(prepared, config, options; toricity=options.toricity)
    chunk_mode === :unchunked &&
        return (
            ArchimedLight._raycore_scene_data(prepared, config; toricity=options.toricity),
            false,
        )
    return (
        ArchimedLight._raycore_scene_data(
            prepared,
            config;
            toricity=options.toricity,
            face_chunk_limit=chunk_mode,
        ),
        true,
    )
end

function _same_raycore_topology(left, right)
    return left.geometry_mode == right.geometry_mode &&
           left.chunked_tlas == right.chunked_tlas &&
           left.tlas_instances == right.tlas_instances &&
           left.tlas_geometries == right.tlas_geometries &&
           left.tlas_depth == right.tlas_depth &&
           left.blas_depth == right.blas_depth &&
           left.blas_node_count == right.blas_node_count &&
           left.expanded_face_count == right.expanded_face_count &&
           left.expanded_face_instance_upper_bound == right.expanded_face_instance_upper_bound &&
           left.reference_compact_face_count == right.reference_compact_face_count
end

function _reference_result(data, candidate_topology, chunk_label)
    reference_topology = _raycore_topology_summary(data)
    return (data, _raw_chunk_mode_label(chunk_label), _same_raycore_topology(reference_topology, candidate_topology))
end

function _reference_data_matching_candidate(prepared, config, options, candidate_topology)
    if !candidate_topology.chunked_tlas
        data = ArchimedLight._raycore_scene_data(prepared, config; toricity=options.toricity)
        return _reference_result(data, candidate_topology, :unchunked)
    end

    for limit in ArchimedLight._raycore_retry_chunk_limits()
        data = ArchimedLight._raycore_scene_data(
            prepared,
            config;
            toricity=options.toricity,
            face_chunk_limit=limit,
        )
        reference_topology = _raycore_topology_summary(data)
        _same_raycore_topology(reference_topology, candidate_topology) &&
            return (data, string(limit), true)
    end

    data = ArchimedLight._raycore_scene_data(
        prepared,
        config;
        toricity=options.toricity,
        face_chunk_limit=ArchimedLight._raycore_face_chunk_limit(),
    )
    return _reference_result(data, candidate_topology, ArchimedLight._raycore_face_chunk_limit())
end

function _reference_data_for(prepared, candidate_data, config, options, chunk_mode::Union{Symbol,Int})
    candidate_topology = _raycore_topology_summary(candidate_data)
    chunk_mode isa Int &&
        return _reference_result(
            ArchimedLight._raycore_scene_data(
                prepared,
                config;
                toricity=options.toricity,
                face_chunk_limit=chunk_mode,
            ),
            candidate_topology,
            chunk_mode,
        )
    chunk_mode === :unchunked &&
        return _reference_result(
            ArchimedLight._raycore_scene_data(prepared, config; toricity=options.toricity),
            candidate_topology,
            :unchunked,
        )
    return _reference_data_matching_candidate(prepared, config, options, candidate_topology)
end

function _bvh_depth(nodes)
    isempty(nodes) && return 0
    function visit(node_idx::Int, depth::Int)
        1 <= node_idx <= length(nodes) || return depth
        node = nodes[node_idx]
        node.child0 == RAYCORE_INVALID_NODE && return depth
        return max(visit(Int(node.child0), depth + 1), visit(Int(node.child1), depth + 1))
    end
    return visit(1, 1)
end

function _max_blas_depth_and_nodes(data)
    static = data.tlas.static_tlas
    nodes = Array(static.all_blas_nodes)
    descriptors = Array(static.blas_descriptors)
    isempty(nodes) && return (depth=0, node_count=0)
    isempty(descriptors) && return (depth=_bvh_depth(nodes), node_count=length(nodes))

    max_depth = 0
    for idx in eachindex(descriptors)
        first_node = Int(descriptors[idx].nodes_offset) + 1
        last_node = idx == length(descriptors) ? length(nodes) : Int(descriptors[idx+1].nodes_offset)
        first_node <= last_node || continue
        max_depth = max(max_depth, _bvh_depth(@view nodes[first_node:last_node]))
    end
    return (depth=max_depth, node_count=length(nodes))
end

function _raycore_topology_summary(data)
    static = data.tlas.static_tlas
    blas = _max_blas_depth_and_nodes(data)
    shape = ArchimedLight._raycore_scene_shape_summary(data)
    return (
        geometry_mode=shape.geometry_mode,
        chunked_tlas=shape.chunked_tlas,
        tlas_instances=shape.tlas_instances,
        tlas_geometries=shape.tlas_geometries,
        tlas_depth=_bvh_depth(Array(static.nodes)),
        blas_depth=blas.depth,
        blas_node_count=blas.node_count,
        expanded_face_count=shape.expanded_face_count,
        expanded_face_instance_upper_bound=shape.expanded_face_instance_upper_bound,
        reference_compact_face_count=shape.reference_compact_face_count,
    )
end

function _run_row(prepared, options, direction, backend, max_hits, chunk_mode::Union{Symbol,Int})
    candidate_config = _raw_config(backend, max_hits)
    reference_config = _raw_config(KernelAbstractions.CPU(), max_hits)
    candidate_data, candidate_was_chunked =
        _scene_data_for(prepared, candidate_config, options, chunk_mode)
    reference_data, reference_face_chunk_limit, reference_matches_candidate_topology =
        _reference_data_for(prepared, candidate_data, reference_config, options, chunk_mode)
    comparison = ArchimedLight._raycore_raw_stack_comparison(
        reference_data,
        candidate_data,
        direction,
        options;
        distance_atol=RAW_DISTANCE_ATOL,
    )
    reference_topology = _raycore_topology_summary(reference_data)
    candidate_topology = _raycore_topology_summary(candidate_data)
    return (
        comparison...,
        requested_face_chunk_limit=_raw_chunk_mode_label(chunk_mode),
        reference_face_chunk_limit=reference_face_chunk_limit,
        reference_matches_candidate_topology=reference_matches_candidate_topology,
        candidate_chunked=candidate_was_chunked,
        reference_geometry_mode=reference_topology.geometry_mode,
        candidate_geometry_mode=candidate_topology.geometry_mode,
        reference_chunked=reference_topology.chunked_tlas,
        candidate_chunked_tlas=candidate_topology.chunked_tlas,
        reference_instances=reference_topology.tlas_instances,
        candidate_instances=candidate_topology.tlas_instances,
        reference_geometries=reference_topology.tlas_geometries,
        candidate_geometries=candidate_topology.tlas_geometries,
        reference_tlas_depth=reference_topology.tlas_depth,
        candidate_tlas_depth=candidate_topology.tlas_depth,
        reference_blas_depth=reference_topology.blas_depth,
        candidate_blas_depth=candidate_topology.blas_depth,
        raycore_traversal_stack_capacity=RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
        reference_blas_depth_over_stack_capacity=reference_topology.blas_depth > RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
        candidate_blas_depth_over_stack_capacity=candidate_topology.blas_depth > RAYCORE_SOFTWARE_TRAVERSAL_STACK_CAPACITY,
        reference_blas_node_count=reference_topology.blas_node_count,
        candidate_blas_node_count=candidate_topology.blas_node_count,
        reference_expanded_face_count=reference_topology.expanded_face_count,
        candidate_expanded_face_count=candidate_topology.expanded_face_count,
        reference_expanded_face_instance_upper_bound=reference_topology.expanded_face_instance_upper_bound,
        candidate_expanded_face_instance_upper_bound=candidate_topology.expanded_face_instance_upper_bound,
        reference_compact_face_count=reference_topology.reference_compact_face_count,
        candidate_compact_face_count=candidate_topology.reference_compact_face_count,
        mismatch_summary=_mismatch_summary(reference_data, direction, options, comparison.mismatch),
    )
end

function _write_rows(rows)
    isempty(RAW_OUTPUT) && return nothing
    open(RAW_OUTPUT, "w") do io
        fields = propertynames(rows[1])
        println(io, join(string.(fields), "\t"))
        for row in rows
            println(io, join((string(getproperty(row, field)) for field in fields), "\t"))
        end
    end
    @info "Wrote raw stack comparison rows" path=RAW_OUTPUT rows=length(rows)
    return nothing
end

function main()
    backend = _backend_from_env()
    workload = _load_workload()
    options = workload.options
    turtle = ArchimedLight.build_turtle(options, workload.sky)
    sector_indices = _selected_sectors(turtle)
    chunk_modes = _raw_face_chunk_limits()
    prepared = ArchimedLight._prepare_interception_data(
        workload.scene,
        workload.models,
        options;
        include_budget_maps=true,
        include_raycore_instancing=true,
    )

    rows = NamedTuple[]
    println("Raycore raw stack comparison")
    println("workload: ", workload.name)
    println("candidate backend: ", typeof(backend))
    println("hit_epsilon: ", RAW_HIT_EPSILON)
    println("distance_atol: ", RAW_DISTANCE_ATOL)
    println("sectors: ", join(sector_indices, ", "))
    println("face_chunk_limits: ", join((_raw_chunk_mode_label(mode) for mode in chunk_modes), ", "))
    @printf(
        "%8s %8s %8s %9s %10s %12s %12s %9s %10s %10s %9s %8s %8s %8s %8s %-18s %-20s %-16s %-16s %-16s %8s %8s\n",
        "maxhits",
        "sector",
        "chunk",
        "refchunk",
        "topomatch",
        "ref_hits",
        "got_hits",
        "hit_pct",
        "ref_occ",
        "got_occ",
        "occ_pct",
        "ref_ovf",
        "got_ovf",
        "ref_max",
        "got_max",
        "ok",
        "diagnosis",
        "mismatch",
        "ref_mode",
        "got_mode",
        "ref_blas",
        "got_blas",
    )
    for max_hits in RAW_MAX_HITS, chunk_mode in chunk_modes, sector_idx in sector_indices
        sector = turtle.sectors[sector_idx]
        row = _run_row(prepared, options, sector.direction, backend, max_hits, chunk_mode)
        full_row = merge(
            (
                workload=workload.name,
                backend=RAW_BACKEND,
                max_hits=max_hits,
                sector=sector_idx,
                sector_id=sector.id,
                sector_source=sector.source,
            ),
            row,
        )
        push!(rows, full_row)
        @printf(
            "%8d %8d %8s %9s %10s %12d %12d %8.3f%% %10d %10d %8.3f%% %8d %8d %8d %8d %-18s %-20s %-16s %-16s %-16s %8d %8d\n",
            max_hits,
            sector_idx,
            row.requested_face_chunk_limit,
            row.reference_face_chunk_limit,
            string(row.reference_matches_candidate_topology),
            row.reference_hits,
            row.candidate_hits,
            100 * row.hit_ratio,
            row.reference_occupied,
            row.candidate_occupied,
            100 * row.occupied_ratio,
            row.reference_overflow,
            row.candidate_overflow,
            row.reference_max_seen,
            row.candidate_max_seen,
            string(row.ok),
            string(row.diagnosis),
            row.mismatch === nothing ? "none" : string(row.mismatch.kind),
            string(row.reference_geometry_mode),
            string(row.candidate_geometry_mode),
            row.reference_blas_depth,
            row.candidate_blas_depth,
        )
        row.mismatch === nothing || println("  first mismatch: ", row.mismatch_summary)
    end
    _write_rows(rows)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
