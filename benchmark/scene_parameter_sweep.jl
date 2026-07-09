const _BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_PARAM_METAL", get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", "")))
const _BENCH_WANTS_METAL = _BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using GeometryBasics
using KernelAbstractions
using PlantGeom
using Printf
using Statistics

const BENCH_SCENES = [
    strip(lowercase(value))
    for value in split(get(ENV, "ARCHIMEDLIGHT_PARAM_SCENES", "single_plate,stacked_plates,canopy_grid"), ',')
    if !isempty(strip(value))
]
const BENCH_BACKENDS = [
    strip(lowercase(value))
    for value in split(get(ENV, "ARCHIMEDLIGHT_PARAM_BACKENDS", "normal_cpu,raycore_ka_cpu"), ',')
    if !isempty(strip(value))
]
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_PARAM_SAMPLES", "1"))
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_PARAM_OUTPUT", "/tmp/archimed_scene_parameter_sweep.tsv")
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_PARAM_MAX_HITS", get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "32")))
const BENCH_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_PARAM_WORKGROUPSIZE", get(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", "256")))
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_PARAM_EDGE_ACCUMULATION", get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto")))
const BENCH_VALIDATE = lowercase(get(ENV, "ARCHIMEDLIGHT_PARAM_VALIDATE", "1")) in ("1", "true", "yes", "on")
const BENCH_SCATTERING = lowercase(get(ENV, "ARCHIMEDLIGHT_PARAM_SCATTERING", "1")) in ("1", "true", "yes", "on")
const BENCH_CACHE_RADIATION = lowercase(get(ENV, "ARCHIMEDLIGHT_PARAM_CACHE_RADIATION", "1")) in ("1", "true", "yes", "on")

function _split_int_env(name::AbstractString, default::Vector{Int})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    return [parse(Int, strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

function _split_float_env(name::AbstractString, default::Vector{Float64})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    return [parse(Float64, strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

function _scene_env_prefix(name::AbstractString)
    clean = replace(uppercase(name), r"[^A-Z0-9]+" => "_")
    return "ARCHIMEDLIGHT_PARAM_$(clean)"
end

function _scene_pixels(name::AbstractString, default::Vector{Float64})
    global_pixels = _split_float_env("ARCHIMEDLIGHT_PARAM_PIXELS", default)
    return _split_float_env("$(_scene_env_prefix(name))_PIXELS", global_pixels)
end

function _scene_sectors(name::AbstractString, default::Vector{Int})
    global_sectors = _split_int_env("ARCHIMEDLIGHT_PARAM_SECTORS", default)
    return _split_int_env("$(_scene_env_prefix(name))_SECTORS", global_sectors)
end

function _flag_required(name::AbstractString)
    return lowercase(get(ENV, name, "")) in ("required", "force", "error")
end

function _optional_metal_backend()
    _BENCH_WANTS_METAL || return nothing
    if !Sys.isapple() || Sys.ARCH != :aarch64
        msg = "Metal benchmark requires Apple Silicon macOS."
        _flag_required("ARCHIMEDLIGHT_PARAM_METAL") ? error(msg) : (@warn msg; return nothing)
    end
    if Base.find_package("Metal") === nothing
        msg = "Metal.jl is not available in this environment."
        _flag_required("ARCHIMEDLIGHT_PARAM_METAL") ? error(msg) : (@warn msg; return nothing)
    end
    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _raycore_interception_backend(backend)
    return ArchimedLight.RaycoreInterceptionBackend(
        backend=backend,
        max_hits_per_pixel=BENCH_MAX_HITS,
        workgroupsize=BENCH_WORKGROUPSIZE,
        edge_accumulation=BENCH_EDGE_ACCUMULATION,
        validate=BENCH_VALIDATE,
    )
end

function _backend_cases()
    selected = Set(BENCH_BACKENDS)
    cases = NamedTuple[]
    if "all" in selected || "normal_cpu" in selected
        push!(cases, (name="normal_cpu", interception_backend=:raster_cpu, scattering_backend=nothing))
    end
    if "all" in selected || "raycore_ka_cpu" in selected
        ib = _raycore_interception_backend(KernelAbstractions.CPU())
        push!(cases, (name="raycore_ka_cpu", interception_backend=ib, scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib)))
    end
    if "all" in selected || "raycore_metal_gpu" in selected
        metal = _optional_metal_backend()
        if metal !== nothing
            ib = _raycore_interception_backend(metal)
            push!(cases, (name="raycore_metal_gpu", interception_backend=ib, scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib)))
        end
    end
    isempty(cases) && error("No backend selected. Use normal_cpu, raycore_ka_cpu, raycore_metal_gpu, or all.")
    return cases
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

function _default_models()
    return ArchimedLight.models_for(
        "leaf" => ("Leaf" => ArchimedLight.translucent(par=0.15, nir=0.85),),
        "panel" => ("Panel" => ArchimedLight.translucent(par=0.0, nir=0.0),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
end

function _single_plate_workload()
    leaf = _panel_mesh(1.0, 1.0, 0.8, 0.0)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, 1.2, 1.2)) do s
        PlantGeom.add_object!(s, leaf; group="leaf", type="Leaf", id=1, at=(0.1, 0.6, 0.0))
        PlantGeom.add_ground!(s; nx=8, ny=8, group="pavement", type="Cobblestone")
    end
    name = "single_plate"
    return (
        name=name,
        scene=scene,
        models=_default_models(),
        pixels=_scene_pixels(name, [0.20, 0.10]),
        sectors=_scene_sectors(name, [1, 6, 16]),
        sky=ArchimedLight.SkyState(180.0, 70.0, 300.0, 100.0, 0.75, 0.25),
        step_seconds=900.0,
        toricity=false,
    )
end

function _stacked_plates_workload()
    leaf = _panel_mesh(1.0, 1.0, 0.0, 0.0)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, 1.2, 1.2)) do s
        for (id, z) in enumerate((0.35, 0.65, 0.95))
            PlantGeom.add_object!(
                s,
                leaf;
                group="leaf",
                type="Leaf",
                id=id,
                at=(0.1 + 0.03 * id, 0.6, z),
                rotate=(z=12.0 * id,),
                deg=true,
            )
        end
        PlantGeom.add_ground!(s; nx=10, ny=10, group="pavement", type="Cobblestone")
    end
    name = "stacked_plates"
    return (
        name=name,
        scene=scene,
        models=_default_models(),
        pixels=_scene_pixels(name, [0.15, 0.075]),
        sectors=_scene_sectors(name, [1, 6, 16]),
        sky=ArchimedLight.SkyState(170.0, 55.0, 360.0, 180.0, 0.70, 0.30),
        step_seconds=900.0,
        toricity=false,
    )
end

function _canopy_grid_workload()
    leaf = _panel_mesh(0.16, 0.36, 0.75, 35.0)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, 3.0, 3.0)) do s
        id = 1
        for i in 0:5, j in 0:5
            PlantGeom.add_object!(
                s,
                leaf;
                group="leaf",
                type="Leaf",
                id=id,
                at=(0.28 + 0.42 * i, 0.30 + 0.42 * j, 0.04 * sin(i + j)),
                rotate=(z=17.0 * id, x=8.0 * sin(id)),
                deg=true,
            )
            id += 1
        end
        PlantGeom.add_ground!(s; nx=20, ny=20, group="pavement", type="Cobblestone")
    end
    name = "canopy_grid"
    return (
        name=name,
        scene=scene,
        models=_default_models(),
        pixels=_scene_pixels(name, [0.12, 0.06]),
        sectors=_scene_sectors(name, [6, 16, 46]),
        sky=ArchimedLight.SkyState(210.0, 52.0, 420.0, 220.0, 0.70, 0.30),
        step_seconds=1800.0,
        toricity=true,
    )
end

function _panel_canopy_workload()
    leaf = _panel_mesh(0.12, 0.30, 0.55, 25.0)
    panel = _panel_mesh(2.8, 0.9, 1.8, 22.0)
    scene = PlantGeom.make_scene(domain=(0.0, 0.0, 3.0, 3.0)) do s
        PlantGeom.add_object!(s, panel; group="panel", type="Panel", id=1, at=(0.1, 1.5, 0.0))
        id = 2
        for i in 0:4, j in 0:4
            PlantGeom.add_object!(
                s,
                leaf;
                group="leaf",
                type="Leaf",
                id=id,
                at=(0.40 + 0.45 * i, 0.35 + 0.45 * j, 0.02 * sin(id)),
                rotate=(z=11.0 * id, x=15.0),
                deg=true,
            )
            id += 1
        end
        PlantGeom.add_ground!(s; nx=20, ny=20, group="pavement", type="Cobblestone")
    end
    name = "panel_canopy"
    return (
        name=name,
        scene=scene,
        models=_default_models(),
        pixels=_scene_pixels(name, [0.12, 0.06]),
        sectors=_scene_sectors(name, [6, 16, 46]),
        sky=ArchimedLight.SkyState(210.0, 48.0, 420.0, 260.0, 0.68, 0.32),
        step_seconds=1800.0,
        toricity=true,
    )
end

function _coffee_workload()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    paving = parse(Int, get(ENV, "ARCHIMEDLIGHT_PARAM_COFFEE_PAVING", "2500"))
    options, scene, _meteo, models = ArchimedLight.read_config(config; plot_paving_override=paving)
    name = "coffee_docs"
    return (
        name=name,
        scene=scene,
        models=models,
        pixels=_scene_pixels(name, [0.04, 0.02]),
        sectors=_scene_sectors(name, [6, 16]),
        sky=ArchimedLight.SkyState(180.0, 52.0, 360.0, 260.0, 0.70, 0.30),
        step_seconds=300.0,
        toricity=options.toricity,
    )
end

function _workload_by_name(name::AbstractString)
    name in ("single", "single_plate") && return _single_plate_workload()
    name in ("stacked", "stacked_plates") && return _stacked_plates_workload()
    name in ("canopy", "canopy_grid") && return _canopy_grid_workload()
    name in ("panel", "panel_canopy") && return _panel_canopy_workload()
    name in ("coffee", "coffee_docs") && return _coffee_workload()
    error("Unknown scene '$name'. Use single_plate, stacked_plates, canopy_grid, panel_canopy, coffee_docs, or all.")
end

function _selected_workloads()
    selected = "all" in BENCH_SCENES ?
               ["single_plate", "stacked_plates", "canopy_grid", "panel_canopy", "coffee_docs"] :
               BENCH_SCENES
    return [_workload_by_name(name) for name in selected]
end

function _options_for(workload, pixel_size::Float64, sectors::Int)
    return ArchimedLight.LightOptions(
        pixel_size=pixel_size,
        turtle_sectors=sectors,
        scattering=BENCH_SCATTERING,
        cache_radiation=BENCH_CACHE_RADIATION,
        all_in_turtle=true,
        toricity=workload.toricity,
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

_mib(bytes) = Float64(bytes) / 1024.0^2

_median(values) = median(collect(values))

function _time_once(workload, options, case)
    sim = _new_simulation(workload, options, case)
    prepare = @timed ArchimedLight._ensure_light_cache!(sim)
    populate = @timed ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=workload.step_seconds)
    reuse = @timed ArchimedLight.run_light(sim, workload.sky; step_duration_seconds=workload.step_seconds)
    cache = sim.cache
    raycore_data = cache === nothing ? nothing : cache.raycore_data
    shape = ArchimedLight._raycore_scene_shape_summary(raycore_data)
    totals = _step_totals(reuse.value)
    return (
        status="ok",
        error_type="",
        error_reason=:none,
        error_stage=:none,
        error_message="",
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
        raycore_expanded_face_count=shape.expanded_face_count,
        raycore_dense_edge_pairs=shape.dense_edge_pairs,
        incident_par=totals.incident_par,
        incident_nir=totals.incident_nir,
    )
end

function _error_case_result(err, case)
    is_validation = err isa ArchimedLight.RaycoreValidationError
    return (
        status="error",
        error_type=string(nameof(typeof(err))),
        error_reason=is_validation ? err.reason : :exception,
        error_stage=is_validation ? err.stage : :benchmark,
        error_message=sprint(showerror, err),
        prepare_seconds=NaN,
        prepare_mib=NaN,
        populate_seconds=NaN,
        populate_mib=NaN,
        reuse_seconds=NaN,
        reuse_mib=NaN,
        resolved_backend=case.name == "normal_cpu" ? "RasterCPUBackend" : "RaycoreInterceptionBackend",
        geometry_mode=:error,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_expanded_face_count=0,
        raycore_dense_edge_pairs=0,
        incident_par=NaN,
        incident_nir=NaN,
    )
end

function _benchmark_case(workload, options, case)
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
        (
            status=last.status,
            error_type=last.error_type,
            error_reason=last.error_reason,
            error_stage=last.error_stage,
            error_message=last.error_message,
            resolved_backend=last.resolved_backend,
            geometry_mode=last.geometry_mode,
            raycore_chunked=last.raycore_chunked,
            raycore_tlas_instances=last.raycore_tlas_instances,
            raycore_tlas_geometries=last.raycore_tlas_geometries,
            raycore_expanded_face_count=last.raycore_expanded_face_count,
            raycore_dense_edge_pairs=last.raycore_dense_edge_pairs,
            incident_par=last.incident_par,
            incident_nir=last.incident_nir,
        ),
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
    @info "Wrote scene parameter sweep results" path=BENCH_OUTPUT rows=length(results)
    return results
end

function main()
    workloads = _selected_workloads()
    cases = _backend_cases()
    results = NamedTuple[]
    baselines = Dict{Tuple{String,Float64,Int},NamedTuple}()

    println("ArchimedLight scene parameter sweep")
    println("scenes: ", join((w.name for w in workloads), ", "))
    println("backends: ", join((c.name for c in cases), ", "))
    println("samples per row after warmup: ", BENCH_SAMPLES)
    println("output: ", BENCH_OUTPUT)
    println()
    @printf("%-16s %-16s %-7s %8s %7s %8s %8s %9s %-11s %-22s %10s %10s %10s %10s\n",
        "scene", "backend", "status", "pixel", "dirs", "nodes", "faces", "pixels", "resolved", "mode",
        "prepare", "populate", "reuse", "PAR d%")

    for workload in workloads
        for pixel_size in workload.pixels, sectors in workload.sectors
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
                @printf("%-16s %-16s %-7s %8.4f %7d %8d %8d %9d %-11s %-22s %10.3f %10.3f %10.3f %10.3f\n",
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
                    full_row.prepare_seconds,
                    full_row.populate_seconds,
                    full_row.reuse_seconds,
                    full_row.par_delta_percent,
                )
                @info "Scene parameter benchmark row complete" scene=full_row.scene backend=full_row.backend status=full_row.status pixel_size=full_row.pixel_size sectors=full_row.sectors prepare_seconds=full_row.prepare_seconds populate_seconds=full_row.populate_seconds reuse_seconds=full_row.reuse_seconds
                full_row.status == "error" && @info "Scene parameter row failed" scene=full_row.scene backend=full_row.backend reason=full_row.error_reason stage=full_row.error_stage message=full_row.error_message
            end
        end
    end
    return _write_results(results)
end

main()
