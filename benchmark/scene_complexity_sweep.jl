const _BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", ""))
const _BENCH_WANTS_METAL = _BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using GeometryBasics
using KernelAbstractions
using MultiScaleTreeGraph
using PlantGeom
using Printf
using Statistics

const BENCH_SCENES_RAW = get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SCENES", "simple,coffee")
const BENCH_SCENES = [strip(lowercase(value)) for value in split(BENCH_SCENES_RAW, ',') if !isempty(strip(value))]
const BENCH_SECTORS = [parse(Int, strip(value)) for value in split(get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SECTORS", "6,16,46"), ',') if !isempty(strip(value))]
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_SAMPLES", "1"))
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_OUTPUT", "/tmp/archimed_scene_complexity_sweep.tsv")
const BENCH_BACKENDS = [strip(lowercase(value)) for value in split(get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_BACKENDS", "normal_cpu,raycore_ka_cpu,raycore_metal_gpu"), ',') if !isempty(strip(value))]
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "32"))
const BENCH_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", "256"))
const BENCH_RASTERGPU_TILE_SIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_SIZE", "1"))
const BENCH_RASTERGPU_TILE_FACE_CAPACITY =
    parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_FACE_CAPACITY", "32"))
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto"))
const BENCH_MAX_PRECHUNK_INSTANCES = haskey(ENV, "ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES") ?
                                     parse(Int, ENV["ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES"]) :
                                     nothing
const BENCH_VALIDATE = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_VALIDATE", "1")) in ("1", "true", "yes", "on")

function _split_float_env(name::AbstractString, default::Vector{Float64})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    return [parse(Float64, strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

function _flag_required(name::AbstractString)
    return lowercase(get(ENV, name, "")) in ("required", "force", "error")
end

function _optional_metal_backend()
    _BENCH_WANTS_METAL || return nothing
    if !Sys.isapple() || Sys.ARCH != :aarch64
        msg = "Metal benchmark requires Apple Silicon macOS."
        _flag_required("ARCHIMEDLIGHT_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end
    if Base.find_package("Metal") === nothing
        msg = "Metal.jl is not available in this environment."
        _flag_required("ARCHIMEDLIGHT_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _raycore_interception_backend(backend)
    kwargs = (
        backend=backend,
        max_hits_per_pixel=BENCH_MAX_HITS,
        workgroupsize=BENCH_WORKGROUPSIZE,
        edge_accumulation=BENCH_EDGE_ACCUMULATION,
        validate=BENCH_VALIDATE,
    )
    BENCH_MAX_PRECHUNK_INSTANCES !== nothing &&
        (kwargs = merge(kwargs, (; max_prechunk_instances=BENCH_MAX_PRECHUNK_INSTANCES)))
    return ArchimedLight.RaycoreInterceptionBackend(; kwargs...)
end

function _rastergpu_interception_backend(backend)
    return ArchimedLight.RasterGPUBackend(
        backend=backend,
        max_hits_per_pixel=BENCH_MAX_HITS,
        workgroupsize=BENCH_WORKGROUPSIZE,
        tile_size=BENCH_RASTERGPU_TILE_SIZE,
        tile_face_capacity=BENCH_RASTERGPU_TILE_FACE_CAPACITY,
        edge_accumulation=BENCH_EDGE_ACCUMULATION,
    )
end

function _backend_cases()
    cases = NamedTuple[
        (name="normal_cpu", interception_backend=:raster_cpu, scattering_backend=nothing),
    ]
    if "raycore_ka_cpu" in BENCH_BACKENDS || "all" in BENCH_BACKENDS
        ib = _raycore_interception_backend(KernelAbstractions.CPU())
        push!(cases, (name="raycore_ka_cpu", interception_backend=ib, scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib)))
    end
    if "rastergpu_ka_cpu" in BENCH_BACKENDS || "all" in BENCH_BACKENDS
        ib = _rastergpu_interception_backend(KernelAbstractions.CPU())
        push!(cases, (name="rastergpu_ka_cpu", interception_backend=ib, scattering_backend=ArchimedLight.RasterGPUScatteringBackend(ib)))
    end
    if "raycore_metal_gpu" in BENCH_BACKENDS || "all" in BENCH_BACKENDS
        metal = _optional_metal_backend()
        if metal !== nothing
            ib = _raycore_interception_backend(metal)
            push!(cases, (name="raycore_metal_gpu", interception_backend=ib, scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib)))
        end
    end
    if "rastergpu_metal_gpu" in BENCH_BACKENDS || "all" in BENCH_BACKENDS
        metal = _optional_metal_backend()
        if metal !== nothing
            ib = _rastergpu_interception_backend(metal)
            push!(cases, (name="rastergpu_metal_gpu", interception_backend=ib, scattering_backend=ArchimedLight.RasterGPUScatteringBackend(ib)))
        end
    end
    selected = Set(BENCH_BACKENDS)
    "all" in selected && return cases
    return [case for case in cases if lowercase(case.name) in selected]
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
    pixels = _split_float_env("ARCHIMEDLIGHT_COMPLEXITY_SIMPLE_PIXELS", [0.40, 0.20])
    return (
        name="simple_plant",
        scene=scene,
        models=models,
        base_options=options,
        pixels=pixels,
        sky=ArchimedLight.SkyState(180.0, 55.0, 420.0, 220.0, 0.70, 0.30),
        step_seconds=900.0,
    )
end

function _load_coffee_workload()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    paving = parse(Int, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PAVING", "2500"))
    options, scene, _meteo, models = ArchimedLight.read_config(config; plot_paving_override=paving)
    pixels = _split_float_env("ARCHIMEDLIGHT_COMPLEXITY_COFFEE_PIXELS", [0.006, 0.003])
    return (
        name="coffee_docs",
        scene=scene,
        models=models,
        base_options=options,
        pixels=pixels,
        sky=ArchimedLight.SkyState(180.0, 52.0, 360.0, 260.0, 0.70, 0.30),
        step_seconds=300.0,
    )
end

function _load_agripv_workload()
    root = get(ENV, "ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT", "")
    !isempty(root) || error("Set ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT to enable the agrivoltaics workload.")
    opf_path = joinpath(root, "0_simulations", "archicrop", "wheat", "plant_1995-06-24.opf")
    isfile(opf_path) || error("Agrivoltaics wheat OPF not found at $opf_path")

    n_rows = parse(Int, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_ROWS", "2"))
    plants_per_row = parse(Int, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_PLANTS_PER_ROW", "8"))
    interrow = parse(Float64, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_INTERROW", "0.25"))
    intrarow = parse(Float64, get(ENV, "ARCHIMEDLIGHT_COMPLEXITY_AGRIPV_INTRAROW", "0.30"))
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
        PlantGeom.add_ground!(s; nx=18, ny=24, group="pavement", type="Cobblestone")
    end
    models = ArchimedLight.models_for(
        "wheat" => ("*" => ArchimedLight.translucent(par=0.15, nir=0.90),),
        "panel" => ("Panel" => ArchimedLight.translucent(par=0.0, nir=0.0),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
    options = ArchimedLight.LightOptions(scattering=true, cache_radiation=true, all_in_turtle=true, toricity=true)
    pixels = _split_float_env("ARCHIMEDLIGHT_COMPLEXITY_COMPLEX_PIXELS", [0.10, 0.05])
    return (
        name="agripv_wheat_panel",
        scene=scene,
        models=models,
        base_options=options,
        pixels=pixels,
        sky=ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
        step_seconds=1800.0,
    )
end

function _load_generated_complex_workload()
    leaf = _panel_mesh(0.12, 0.35, 0.55, 35.0)
    scene = PlantGeom.make_scene(domain=(-0.5, -0.5, 3.5, 3.5)) do s
        id = 1
        for i in 0:7, j in 0:7
            PlantGeom.add_object!(
                s,
                leaf;
                group="synthetic_canopy",
                type="Leaf",
                id=id,
                at=(0.35 + 0.35 * i, 0.35 + 0.35 * j, 0.05 * sin(i + j)),
                rotate=(z=13.0 * id, x=20.0 + 5.0 * sin(id)),
                deg=true,
            )
            id += 1
        end
        PlantGeom.add_ground!(s; nx=24, ny=24, group="pavement", type="Cobblestone")
    end
    models = ArchimedLight.models_for(
        "synthetic_canopy" => ("Leaf" => ArchimedLight.translucent(par=0.15, nir=0.90),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
    options = ArchimedLight.LightOptions(scattering=true, cache_radiation=true, all_in_turtle=true, toricity=true)
    pixels = _split_float_env("ARCHIMEDLIGHT_COMPLEXITY_COMPLEX_PIXELS", [0.10, 0.05])
    return (
        name="generated_canopy",
        scene=scene,
        models=models,
        base_options=options,
        pixels=pixels,
        sky=ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
        step_seconds=1800.0,
    )
end

function _load_complex_workload()
    try
        return _load_agripv_workload()
    catch err
        @warn "Falling back to generated complex scene." reason=sprint(showerror, err)
        return _load_generated_complex_workload()
    end
end

function _selected_workloads()
    selected = Set(BENCH_SCENES)
    "all" in selected && (selected = Set(["simple", "coffee", "complex"]))
    workloads = NamedTuple[]
    "simple" in selected && push!(workloads, _load_simple_workload())
    "coffee" in selected && push!(workloads, _load_coffee_workload())
    "complex" in selected && push!(workloads, _load_complex_workload())
    isempty(workloads) && error("No scene selected. Use simple, coffee, complex, or all.")
    return workloads
end

function _options_for(workload, pixel_size::Float64, sectors::Int)
    return ArchimedLight.LightOptions(
        workload.base_options;
        pixel_size=pixel_size,
        turtle_sectors=sectors,
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
    )
end

function _scene_metrics(scene, models, options)
    summary = ArchimedLight.summarize_scene(scene; models=models)
    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_budget_maps=true)
    return (
        nodes=summary.node_count,
        faces=summary.face_count,
        objects=summary.object_count,
        pixel_cells=prepared.geometry.plotbox.nx * prepared.geometry.plotbox.ny,
        plot_nx=prepared.geometry.plotbox.nx,
        plot_ny=prepared.geometry.plotbox.ny,
    )
end

function _new_simulation(workload, options, case)
    if case.name == "normal_cpu"
        return ArchimedLight.LightSimulation(
            workload.scene,
            workload.models;
            options=options,
            interception_backend=:raster_cpu,
            scattering_mode=:raycast,
        )
    end
    return ArchimedLight.LightSimulation(
        workload.scene,
        workload.models;
        options=options,
        interception_backend=case.interception_backend,
        scattering_backend=case.scattering_backend,
    )
end

function _step_totals(step)
    totals = ArchimedLight._light_step_energy_totals(step)
    return (incident_par=totals.incident_par_total, incident_nir=totals.incident_nir_total)
end

function _mib(bytes)
    return Float64(bytes) / 1024.0^2
end

function _reduction_label(caps)
    caps === nothing && return ""
    return string(
        "atomics=", Int(caps.supports_atomics),
        ",edge=", Int(caps.dense_edge_accumulation),
        ",count=", Int(caps.node_count_reduction),
        ",area=", Int(caps.sector_area_reduction),
        ",fused=", Int(caps.fused_count_area_reduction),
    )
end

function _cache_shape_and_reductions(cache)
    if cache === nothing
        shape = ArchimedLight._raycore_scene_shape_summary(nothing)
        return shape, nothing
    end

    if cache.raycore_data !== nothing
        shape = ArchimedLight._raycore_scene_shape_summary(cache.raycore_data)
        caps = ArchimedLight._raycore_device_reduction_capabilities(cache.raycore_data)
        return shape, caps
    end

    rastergpu_data = cache.rastergpu_data
    if rastergpu_data !== nothing
        face_count = length(rastergpu_data.face_i_dev)
        shape = (
            geometry_mode=:raster_gpu,
            chunked_tlas=false,
            tlas_instances=0,
            tlas_geometries=0,
            node_count=length(rastergpu_data.node_counts_host),
            expanded_face_count=face_count,
            expanded_face_instance_upper_bound=face_count,
            reference_compact_face_count=face_count,
            dense_edge_pairs=length(rastergpu_data.dense_edge_counts_host),
            edge_key_capacity=length(rastergpu_data.edge_keys_host),
        )
        caps = (
            supports_atomics=KernelAbstractions.supports_atomics(rastergpu_data.backend),
            dense_edge_accumulation=ArchimedLight._rastergpu_dense_edge_counts_fits(rastergpu_data),
            node_count_reduction=true,
            sector_area_reduction=true,
            fused_count_area_reduction=true,
        )
        return shape, caps
    end

    shape = ArchimedLight._raycore_scene_shape_summary(nothing)
    return shape, nothing
end

function _time_once(workload, options, case)
    sim = _new_simulation(workload, options, case)
    prepare = @timed ArchimedLight._ensure_light_cache!(sim)
    populate = @timed ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=workload.step_seconds)
    reuse = @timed ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=workload.step_seconds)
    cache = sim.cache
    shape, caps = _cache_shape_and_reductions(cache)
    totals = _step_totals(reuse.value)
    return (
        status="ok",
        error_type="",
        error_reason=:none,
        error_stage=:none,
        error_message="",
        validation_hit_ratio=missing,
        validation_occupied_ratio=missing,
        validation_reference_hits=missing,
        validation_raycore_hits=missing,
        validation_reference_occupied=missing,
        validation_raycore_occupied=missing,
        validation_reference_overflow=missing,
        validation_raycore_overflow=missing,
        validation_directions_tested=missing,
        validation_direction_count=missing,
        prepare_seconds=prepare.time,
        prepare_mib=_mib(prepare.bytes),
        populate_seconds=populate.time,
        populate_mib=_mib(populate.bytes),
        reuse_seconds=reuse.time,
        reuse_mib=_mib(reuse.bytes),
        resolved_backend=cache === nothing ? "unprepared" : string(nameof(typeof(cache.resolved_interception_backend))),
        geometry_mode=shape.geometry_mode,
        raycore_chunked=shape.chunked_tlas,
        raycore_tlas_instances=shape.tlas_instances,
        raycore_tlas_geometries=shape.tlas_geometries,
        raycore_node_count=shape.node_count,
        raycore_expanded_face_count=shape.expanded_face_count,
        raycore_expanded_face_instance_upper_bound=shape.expanded_face_instance_upper_bound,
        raycore_reference_compact_face_count=shape.reference_compact_face_count,
        raycore_dense_edge_pairs=shape.dense_edge_pairs,
        raycore_edge_key_capacity=shape.edge_key_capacity,
        reduction_capabilities=_reduction_label(caps),
        incident_par=totals.incident_par,
        incident_nir=totals.incident_nir,
    )
end

_validation_field(validation, field::Symbol, default=missing) =
    validation !== nothing && hasproperty(validation, field) ? getproperty(validation, field) : default

function _error_case_result(err, case)
    is_validation = err isa ArchimedLight.RaycoreValidationError
    validation = is_validation ? err.validation : nothing
    return (
        status="error",
        error_type=string(nameof(typeof(err))),
        error_reason=is_validation ? err.reason : :exception,
        error_stage=is_validation ? err.stage : :benchmark,
        error_message=sprint(showerror, err),
        validation_hit_ratio=_validation_field(validation, :hit_ratio),
        validation_occupied_ratio=_validation_field(validation, :occupied_ratio),
        validation_reference_hits=_validation_field(validation, :reference_hits),
        validation_raycore_hits=_validation_field(validation, :raycore_hits),
        validation_reference_occupied=_validation_field(validation, :reference_occupied),
        validation_raycore_occupied=_validation_field(validation, :raycore_occupied),
        validation_reference_overflow=_validation_field(validation, :reference_overflow),
        validation_raycore_overflow=_validation_field(validation, :raycore_overflow),
        validation_directions_tested=_validation_field(validation, :directions_tested),
        validation_direction_count=_validation_field(validation, :direction_count),
        prepare_seconds=NaN,
        prepare_mib=NaN,
        populate_seconds=NaN,
        populate_mib=NaN,
        reuse_seconds=NaN,
        reuse_mib=NaN,
        resolved_backend=case.name == "normal_cpu" ? "RasterCPUBackend" :
                         startswith(case.name, "rastergpu_") ? "RasterGPUBackend" :
                         "RaycoreInterceptionBackend",
        geometry_mode=:error,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_node_count=0,
        raycore_expanded_face_count=0,
        raycore_expanded_face_instance_upper_bound=0,
        raycore_reference_compact_face_count=0,
        raycore_dense_edge_pairs=0,
        raycore_edge_key_capacity=0,
        reduction_capabilities="",
        incident_par=NaN,
        incident_nir=NaN,
    )
end

function _median(values)
    return median(collect(values))
end

function _benchmark_case(workload, options, case)
    # One unreported warm-up keeps method and kernel compilation out of the row.
    _time_once(workload, options, case)
    samples = [_time_once(workload, options, case) for _ in 1:BENCH_SAMPLES]
    last = samples[end]
    return merge(
        (
            prepare_seconds=_median(s.prepare_seconds for s in samples),
            prepare_mib=_median(s.prepare_mib for s in samples),
            populate_seconds=_median(s.populate_seconds for s in samples),
            populate_mib=_median(s.populate_mib for s in samples),
            reuse_seconds=_median(s.reuse_seconds for s in samples),
            reuse_mib=_median(s.reuse_mib for s in samples),
        ),
        NamedTuple{(
            :status,
            :error_type,
            :error_reason,
            :error_stage,
            :error_message,
            :validation_hit_ratio,
            :validation_occupied_ratio,
            :validation_reference_hits,
            :validation_raycore_hits,
            :validation_reference_occupied,
            :validation_raycore_occupied,
            :validation_reference_overflow,
            :validation_raycore_overflow,
            :validation_directions_tested,
            :validation_direction_count,
            :resolved_backend,
            :geometry_mode,
            :raycore_chunked,
            :raycore_tlas_instances,
            :raycore_tlas_geometries,
            :raycore_node_count,
            :raycore_expanded_face_count,
            :raycore_expanded_face_instance_upper_bound,
            :raycore_reference_compact_face_count,
            :raycore_dense_edge_pairs,
            :raycore_edge_key_capacity,
            :reduction_capabilities,
            :incident_par,
            :incident_nir,
        )}((
            last.status,
            last.error_type,
            last.error_reason,
            last.error_stage,
            last.error_message,
            last.validation_hit_ratio,
            last.validation_occupied_ratio,
            last.validation_reference_hits,
            last.validation_raycore_hits,
            last.validation_reference_occupied,
            last.validation_raycore_occupied,
            last.validation_reference_overflow,
            last.validation_raycore_overflow,
            last.validation_directions_tested,
            last.validation_direction_count,
            last.resolved_backend,
            last.geometry_mode,
            last.raycore_chunked,
            last.raycore_tlas_instances,
            last.raycore_tlas_geometries,
            last.raycore_node_count,
            last.raycore_expanded_face_count,
            last.raycore_expanded_face_instance_upper_bound,
            last.raycore_reference_compact_face_count,
            last.raycore_dense_edge_pairs,
            last.raycore_edge_key_capacity,
            last.reduction_capabilities,
            last.incident_par,
            last.incident_nir,
        )),
    )
end

function _relative_delta(value, reference)
    reference == 0.0 && return value == 0.0 ? 0.0 : Inf
    return (value - reference) / abs(reference)
end

function _tsv_cell(value)
    return replace(string(value), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _write_results(results)
    isempty(BENCH_OUTPUT) && return results
    columns = Symbol[]
    for row in results, key in keys(row)
        key in columns || push!(columns, key)
    end
    dir = dirname(BENCH_OUTPUT)
    isempty(dir) || dir == "." || mkpath(dir)
    open(BENCH_OUTPUT, "w") do io
        println(io, join(string.(columns), '\t'))
        for row in results
            println(io, join((_tsv_cell(key in keys(row) ? getproperty(row, key) : missing) for key in columns), '\t'))
        end
    end
    @info "Wrote scene complexity benchmark results" path=BENCH_OUTPUT rows=length(results)
    return results
end

function main()
    workloads = _selected_workloads()
    cases = _backend_cases()
    results = NamedTuple[]
    baselines = Dict{Tuple{String,Float64,Int},NamedTuple}()

    println("ArchimedLight scene complexity sweep")
    println("scenes: ", join((w.name for w in workloads), ", "))
    println("sectors: ", join(BENCH_SECTORS, ", "))
    println("samples per row after warmup: ", BENCH_SAMPLES)
    println("backends: ", join((c.name for c in cases), ", "))
    println("rastergpu tile size: ", BENCH_RASTERGPU_TILE_SIZE)
    println("rastergpu tile face capacity: ", BENCH_RASTERGPU_TILE_FACE_CAPACITY)
    println("output: ", BENCH_OUTPUT)
    println()
    @printf("%-20s %-16s %-7s %8s %7s %8s %8s %9s %11s %-26s %-24s %10s %10s %10s %10s %10s %10s %10s %9s %9s\n",
        "scene", "backend", "status", "pixel", "dirs", "nodes", "faces", "pixels", "resolved", "mode", "reductions",
        "prepare", "populate", "reuse", "prep MiB", "pop MiB", "reuse MiB", "PAR d%", "tlas", "scratch")

    for workload in workloads
        for pixel_size in workload.pixels, sectors in BENCH_SECTORS
            options = _options_for(workload, pixel_size, sectors)
            metrics = _scene_metrics(workload.scene, workload.models, options)
            combo_rows = NamedTuple[]
            for case in cases
                timing = try
                    _benchmark_case(workload, options, case)
                catch err
                    _error_case_result(err, case)
                end
                row = merge(
                    (
                        scene=workload.name,
                        backend=case.name,
                        pixel_size=pixel_size,
                        sectors=sectors,
                        nodes=metrics.nodes,
                        faces=metrics.faces,
                        objects=metrics.objects,
                        pixel_cells=metrics.pixel_cells,
                        plot_nx=metrics.plot_nx,
                        plot_ny=metrics.plot_ny,
                    ),
                    timing,
                )
                push!(combo_rows, row)
                case.name == "normal_cpu" && row.status == "ok" && (baselines[(workload.name, pixel_size, sectors)] = row)
            end
            base = get(baselines, (workload.name, pixel_size, sectors), nothing)
            for row in combo_rows
                par_delta = base === nothing || row.status != "ok" ? NaN : 100 * _relative_delta(row.incident_par, base.incident_par)
                nir_delta = base === nothing || row.status != "ok" ? NaN : 100 * _relative_delta(row.incident_nir, base.incident_nir)
                full_row = merge(row, (par_delta_percent=par_delta, nir_delta_percent=nir_delta))
                push!(results, full_row)
                @printf("%-20s %-16s %-7s %8.4f %7d %8d %8d %9d %-11s %-26s %-24s %10.3f %10.3f %10.3f %10.1f %10.1f %10.1f %10.3f %9s %9s\n",
                    full_row.scene,
                    full_row.backend,
                    full_row.status,
                    full_row.pixel_size,
                    full_row.sectors,
                    full_row.nodes,
                    full_row.faces,
                    full_row.pixel_cells,
                    full_row.resolved_backend,
                    string(full_row.geometry_mode),
                    full_row.reduction_capabilities,
                    full_row.prepare_seconds,
                    full_row.populate_seconds,
                    full_row.reuse_seconds,
                    full_row.prepare_mib,
                    full_row.populate_mib,
                    full_row.reuse_mib,
                    full_row.par_delta_percent,
                    string(full_row.raycore_tlas_instances, "/", full_row.raycore_tlas_geometries),
                    string(full_row.raycore_dense_edge_pairs, "/", full_row.raycore_edge_key_capacity),
                )
                full_row.status == "error" && @info "Scene complexity row failed" scene=full_row.scene backend=full_row.backend reason=full_row.error_reason stage=full_row.error_stage message=full_row.error_message
            end
        end
    end
    return _write_results(results)
end

main()
