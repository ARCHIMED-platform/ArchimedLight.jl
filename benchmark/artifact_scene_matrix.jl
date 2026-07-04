using ArchimedLight
using Artifacts: artifact_hash, artifact_path
using KernelAbstractions
using Pkg.Artifacts: ensure_artifact_installed
using Printf
using Statistics

const REPO_ROOT = dirname(@__DIR__)
const ARTIFACT_NAME = "archimedlight-benchmark-scenes"
const ARTIFACTS_TOML = joinpath(REPO_ROOT, "Artifacts.toml")
const DEFAULT_OUTPUT = "/tmp/archimedlight_artifact_scene_matrix.tsv"

const DEFAULT_SCENES = [
    "simple_plant",
    "coffee",
    "wheat",
    "agrivoltaics_wheat",
    "elaeis",
    "elaeis_two_plants",
]
const DEFAULT_BACKENDS = ["normal_cpu", "rastergpu_metal", "raycore_metal"]
const DEFAULT_DIRECTIONS = [16, 48, 136]
const DEFAULT_PIXELS_CM = [0.1, 0.5, 1.0, 10.0]
const DEFAULT_BOOL = [false, true]
const DEFAULT_STEPS = [1, 12, 24, 48]
const DEFAULT_RAYCORE_TORIC_TRAVERSAL = [:replicated, :periodic]

function _split_string_env(name::AbstractString, default::Vector{String})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    values = [strip(value) for value in split(raw, ',') if !isempty(strip(value))]
    "all" in lowercase.(values) && return default
    return values
end

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

function _env_bool(name::AbstractString, default::Bool)
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    value = lowercase(strip(raw))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("Unsupported boolean value for $name=$(repr(raw)). Use true or false.")
end

function _split_bool_env(name::AbstractString, default::Vector{Bool})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    return [_parse_bool_token(name, strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

function _parse_bool_token(name::AbstractString, raw::AbstractString)
    value = lowercase(raw)
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    error("Unsupported boolean value in $name: $(repr(raw)). Use true or false.")
end

function _split_symbol_env(name::AbstractString, default::Vector{Symbol})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    return [Symbol(strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

const BENCH_SCENES = _split_string_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_SCENES", DEFAULT_SCENES)
const BENCH_BACKENDS = _split_string_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_BACKENDS", DEFAULT_BACKENDS)
const BENCH_DIRECTIONS = _split_int_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_DIRECTIONS", DEFAULT_DIRECTIONS)
const BENCH_PIXELS_CM = _split_float_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_PIXELS_CM", DEFAULT_PIXELS_CM)
const BENCH_SCATTERING = _split_bool_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_SCATTERING", DEFAULT_BOOL)
const BENCH_CACHE_RADIATION = _split_bool_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_CACHE_RADIATION", DEFAULT_BOOL)
const BENCH_STEPS = _split_int_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_STEPS", DEFAULT_STEPS)
const BENCH_RAYCORE_TORIC_TRAVERSAL =
    _split_symbol_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_RAYCORE_TORIC_TRAVERSAL", DEFAULT_RAYCORE_TORIC_TRAVERSAL)
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_SAMPLES", "1"))
const BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_WARMUPS", "0"))
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_OUTPUT", DEFAULT_OUTPUT)
const BENCH_DRY_RUN = _env_bool("ARCHIMEDLIGHT_ARTIFACT_BENCH_DRY_RUN", false)
const BENCH_MAX_PIXEL_CELLS = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_MAX_PIXEL_CELLS", "25000000"))
const BENCH_STEP_SECONDS = parse(Float64, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_STEP_SECONDS", "1800.0"))
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "32"))
const BENCH_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", "256"))
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto"))
const BENCH_RASTERGPU_TILE_SIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_SIZE", "1"))
const BENCH_RASTERGPU_TILE_FACE_CAPACITY =
    parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_FACE_CAPACITY", "32"))
const BENCH_VALIDATE = _env_bool("ARCHIMEDLIGHT_BENCH_VALIDATE", true)

function _validate_raycore_toric_traversal_values(values)
    allowed = Set(DEFAULT_RAYCORE_TORIC_TRAVERSAL)
    unknown = setdiff(Set(values), allowed)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_ARTIFACT_BENCH_RAYCORE_TORIC_TRAVERSAL value(s): " *
        "$(join(string.(sort!(collect(unknown))), ", ")). Use replicated, periodic, or both.",
    )
    return values
end

_validate_raycore_toric_traversal_values(BENCH_RAYCORE_TORIC_TRAVERSAL)

function _benchmark_scenes_root()
    override = get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_SCENES_ROOT", "")
    if !isempty(strip(override))
        root = normpath(override)
        isdir(root) || error("ARCHIMEDLIGHT_ARTIFACT_BENCH_SCENES_ROOT does not exist: $root")
        return _normalize_scene_root(root)
    end

    try
        hash = artifact_hash(ARTIFACT_NAME, ARTIFACTS_TOML)
        hash === nothing && error("artifact hash not found")
        installed_path = artifact_path(hash)
        isdir(installed_path) && return _normalize_scene_root(installed_path)

        local_root = joinpath(@__DIR__, "benchmark_scenes")
        isdir(local_root) && return _normalize_scene_root(local_root)

        ensure_artifact_installed(ARTIFACT_NAME, ARTIFACTS_TOML)
        return _normalize_scene_root(artifact_path(hash))
    catch err
        local_root = joinpath(@__DIR__, "benchmark_scenes")
        isdir(local_root) && return _normalize_scene_root(local_root)
        error(
            "Could not resolve artifact $(repr(ARTIFACT_NAME)) from $ARTIFACTS_TOML. " *
            "Set ARCHIMEDLIGHT_ARTIFACT_BENCH_SCENES_ROOT to a benchmark_scenes directory. " *
            "Original error: $(sprint(showerror, err))",
        )
    end
end

function _normalize_scene_root(root::AbstractString)
    if any(endswith(name, ".ops") for name in readdir(root))
        return root
    end
    nested = joinpath(root, "benchmark_scenes")
    if isdir(nested) && any(endswith(name, ".ops") for name in readdir(nested))
        return nested
    end
    error("Scene root does not contain OPS files: $root")
end

function _scene_specs(root::AbstractString)
    ops = Dict(splitext(basename(path))[1] => path for path in filter(path -> endswith(path, ".ops"), readdir(root; join=true)))
    missing = setdiff(BENCH_SCENES, collect(keys(ops)))
    isempty(missing) || error(
        "Requested benchmark scene(s) not found in $root: $(join(sort!(missing), ", ")). " *
        "Available scenes: $(join(sort!(collect(keys(ops))), ", ")).",
    )
    return [(name=name, path=ops[name]) for name in BENCH_SCENES]
end

function _benchmark_models()
    vegetation = ArchimedLight.translucent(par=0.15, nir=0.90)
    panel = ArchimedLight.translucent(par=0.0, nir=0.0)
    return ArchimedLight.models_for(
        "panel" => (
            "*" => panel,
            "Panel" => panel,
        ),
        "*" => (
            "*" => vegetation,
            "Panel" => panel,
        ),
    )
end

function _sky_series(n::Int)
    n > 0 || error("Step count must be positive, got $n.")
    return [
        ArchimedLight.SkyState(
            90.0 + 180.0 * (i - 1) / max(n - 1, 1),
            25.0 + 35.0 * sin(pi * (i - 0.5) / n),
            320.0 + 80.0 * sin(2pi * (i - 1) / max(n, 1)),
            180.0 + 50.0 * cos(2pi * (i - 1) / max(n, 1)),
            0.70,
            0.30,
        )
        for i in 1:n
    ]
end

function _metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        return nothing, "Metal benchmarks require Apple Silicon macOS."
    end
    if Base.find_package("Metal") === nothing
        return nothing, "Metal.jl is not available in this Julia environment."
    end
    try
        metal_mod = Base.require(Main, :Metal)
        array_type = getproperty(metal_mod, :MtlArray)
        dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
        backend = Base.invokelatest(KernelAbstractions.get_backend, dev_array)
        return backend, ""
    catch err
        return nothing, "Could not initialize Metal backend: $(sprint(showerror, err))"
    end
end

function _backend_cases()
    selected = lowercase.(BENCH_BACKENDS)
    unknown = setdiff(selected, DEFAULT_BACKENDS)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_ARTIFACT_BENCH_BACKENDS value(s): $(join(sort!(unknown), ", ")). " *
        "Use normal_cpu, rastergpu_metal, raycore_metal, or all.",
    )

    if BENCH_DRY_RUN
        cases = NamedTuple[]
        for name in selected
            if name == "raycore_metal"
                for toric_traversal in BENCH_RAYCORE_TORIC_TRAVERSAL
                    push!(
                        cases,
                        (
                            name=name,
                            available=true,
                            unavailable_reason="",
                            interception_backend=nothing,
                            scattering_backend=nothing,
                            toric_traversal=toric_traversal,
                        ),
                    )
                end
            else
                push!(
                    cases,
                    (
                        name=name,
                        available=true,
                        unavailable_reason="",
                        interception_backend=nothing,
                        scattering_backend=nothing,
                        toric_traversal=:not_applicable,
                    ),
                )
            end
        end
        return cases
    end

    metal, metal_reason = _metal_backend()
    cases = NamedTuple[]
    for name in selected
        if name == "normal_cpu"
            push!(
                cases,
                (
                    name=name,
                    available=true,
                    unavailable_reason="",
                    interception_backend=:raster_cpu,
                    scattering_backend=nothing,
                    toric_traversal=:not_applicable,
                ),
            )
        elseif name == "rastergpu_metal"
            if metal === nothing
                push!(
                    cases,
                    (
                        name=name,
                        available=false,
                        unavailable_reason=metal_reason,
                        interception_backend=nothing,
                        scattering_backend=nothing,
                        toric_traversal=:not_applicable,
                    ),
                )
            else
                ib = ArchimedLight.RasterGPUBackend(
                    backend=metal,
                    max_hits_per_pixel=BENCH_MAX_HITS,
                    workgroupsize=BENCH_WORKGROUPSIZE,
                    tile_size=BENCH_RASTERGPU_TILE_SIZE,
                    tile_face_capacity=BENCH_RASTERGPU_TILE_FACE_CAPACITY,
                    edge_accumulation=BENCH_EDGE_ACCUMULATION,
                )
                push!(
                    cases,
                    (
                        name=name,
                        available=true,
                        unavailable_reason="",
                        interception_backend=ib,
                        scattering_backend=ArchimedLight.RasterGPUScatteringBackend(
                            ib;
                            edge_accumulation=BENCH_EDGE_ACCUMULATION,
                        ),
                        toric_traversal=:not_applicable,
                    ),
                )
            end
        elseif name == "raycore_metal"
            if metal === nothing
                for toric_traversal in BENCH_RAYCORE_TORIC_TRAVERSAL
                    push!(
                        cases,
                        (
                            name=name,
                            available=false,
                            unavailable_reason=metal_reason,
                            interception_backend=nothing,
                            scattering_backend=nothing,
                            toric_traversal=toric_traversal,
                        ),
                    )
                end
            else
                for toric_traversal in BENCH_RAYCORE_TORIC_TRAVERSAL
                    ib = ArchimedLight.RaycoreInterceptionBackend(
                        backend=metal,
                        max_hits_per_pixel=BENCH_MAX_HITS,
                        workgroupsize=BENCH_WORKGROUPSIZE,
                        edge_accumulation=BENCH_EDGE_ACCUMULATION,
                        toric_traversal=toric_traversal,
                        validate=BENCH_VALIDATE,
                    )
                    push!(
                        cases,
                        (
                            name=name,
                            available=true,
                            unavailable_reason="",
                            interception_backend=ib,
                            scattering_backend=ArchimedLight.RaycoreScatteringBackend(
                                ib;
                                edge_accumulation=BENCH_EDGE_ACCUMULATION,
                            ),
                            toric_traversal=toric_traversal,
                        ),
                    )
                end
            end
        end
    end
    return cases
end

function _options(pixel_size_m::Float64, directions::Int, scattering::Bool, cache_radiation::Bool)
    return ArchimedLight.LightOptions(
        pixel_size=pixel_size_m,
        turtle_sectors=directions,
        scattering=scattering,
        cache_radiation=cache_radiation,
        all_in_turtle=true,
        toricity=true,
    )
end

function _scene_metrics(scene, models, pixel_size_m::Float64)
    summary = ArchimedLight.summarize_scene(scene; models=models)
    domain = summary.domain
    if domain === nothing
        return (
            nodes=summary.node_count,
            faces=summary.face_count,
            objects=summary.object_count,
            domain_x=NaN,
            domain_y=NaN,
            pixel_cells_estimate=typemax(Int),
            plot_nx_estimate=typemax(Int),
            plot_ny_estimate=typemax(Int),
        )
    end
    xmin, ymin, xmax, ymax = domain
    domain_x = xmax - xmin
    domain_y = ymax - ymin
    nx = max(1, ceil(Int, domain_x / pixel_size_m))
    ny = max(1, ceil(Int, domain_y / pixel_size_m))
    return (
        nodes=summary.node_count,
        faces=summary.face_count,
        objects=summary.object_count,
        domain_x=domain_x,
        domain_y=domain_y,
        pixel_cells_estimate=nx * ny,
        plot_nx_estimate=nx,
        plot_ny_estimate=ny,
    )
end

function _ops_domain_estimate(path::AbstractString)
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        parts = split(stripped)
        if !isempty(parts) && parts[1] == "T" && length(parts) >= 6
            vals = parse.(Float64, parts[2:6])
            domain_x = abs(vals[4] - vals[1])
            domain_y = abs(vals[5] - vals[2])
            domain_x == 0.0 && (domain_x = abs(vals[4]))
            domain_y == 0.0 && (domain_y = abs(vals[5]))
            return (domain_x=domain_x, domain_y=domain_y)
        end
    end
    return (domain_x=NaN, domain_y=NaN)
end

function _ops_object_count_estimate(path::AbstractString)
    count = 0
    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        parts = split(stripped)
        isempty(parts) && continue
        parts[1] == "T" && continue
        parts[1] == "-1" && continue
        count += 1
    end
    return count
end

function _dry_run_scene_metrics(path::AbstractString, pixel_size_m::Float64)
    domain = _ops_domain_estimate(path)
    if !isfinite(domain.domain_x) || !isfinite(domain.domain_y)
        pixel_cells = typemax(Int)
        nx = typemax(Int)
        ny = typemax(Int)
    else
        nx = max(1, ceil(Int, domain.domain_x / pixel_size_m))
        ny = max(1, ceil(Int, domain.domain_y / pixel_size_m))
        pixel_cells = nx * ny
    end
    return (
        nodes=-1,
        faces=-1,
        objects=_ops_object_count_estimate(path),
        domain_x=domain.domain_x,
        domain_y=domain.domain_y,
        pixel_cells_estimate=pixel_cells,
        plot_nx_estimate=nx,
        plot_ny_estimate=ny,
    )
end

function _new_simulation(scene, models, options, case)
    if case.name == "normal_cpu"
        return ArchimedLight.LightSimulation(
            scene,
            models;
            options=options,
            interception_backend=:raster_cpu,
            scattering_mode=:raycast,
        )
    end
    return ArchimedLight.LightSimulation(
        scene,
        models;
        options=options,
        interception_backend=case.interception_backend,
        scattering_backend=case.scattering_backend,
    )
end

function _step_totals(step)
    totals = ArchimedLight._light_step_energy_totals(step)
    return (
        incident_par=Float64(totals.incident_par_total),
        incident_nir=Float64(totals.incident_nir_total),
        absorbed_par=Float64(totals.absorbed_par_total),
        absorbed_nir=Float64(totals.absorbed_nir_total),
    )
end

function _series_totals(steps)
    totals = (incident_par=0.0, incident_nir=0.0, absorbed_par=0.0, absorbed_nir=0.0)
    for step in steps
        current = _step_totals(step)
        totals = (
            incident_par=totals.incident_par + current.incident_par,
            incident_nir=totals.incident_nir + current.incident_nir,
            absorbed_par=totals.absorbed_par + current.absorbed_par,
            absorbed_nir=totals.absorbed_nir + current.absorbed_nir,
        )
    end
    return totals
end

function _run_series(sim, skies)
    return [ArchimedLight.run_light(sim, sky; step_duration_seconds=BENCH_STEP_SECONDS) for sky in skies]
end

function _cache_info(sim)
    cache = sim.cache
    cache === nothing && return (
        resolved_backend="",
        geometry_mode=:unprepared,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_expanded_face_count=0,
    )
    if cache.raycore_data !== nothing
        shape = ArchimedLight._raycore_scene_shape_summary(cache.raycore_data)
        return (
            resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
            geometry_mode=shape.geometry_mode,
            raycore_chunked=shape.chunked_tlas,
            raycore_tlas_instances=shape.tlas_instances,
            raycore_tlas_geometries=shape.tlas_geometries,
            raycore_expanded_face_count=shape.expanded_face_count,
        )
    end
    if hasproperty(cache, :rastergpu_data) && cache.rastergpu_data !== nothing
        return (
            resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
            geometry_mode=:raster_gpu,
            raycore_chunked=false,
            raycore_tlas_instances=0,
            raycore_tlas_geometries=0,
            raycore_expanded_face_count=0,
        )
    end
    return (
        resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
        geometry_mode=:raster_cpu,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_expanded_face_count=0,
    )
end

function _measure(scene, models, options, case, skies)
    for _ in 1:BENCH_WARMUPS
        sim = _new_simulation(scene, models, options, case)
        _run_series(sim, skies)
    end

    times_ms = Float64[]
    alloc_mb = Float64[]
    totals = (incident_par=NaN, incident_nir=NaN, absorbed_par=NaN, absorbed_nir=NaN)
    info = (
        resolved_backend="",
        geometry_mode=:unknown,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_expanded_face_count=0,
    )
    for _ in 1:BENCH_SAMPLES
        GC.gc()
        sim = _new_simulation(scene, models, options, case)
        series = Ref{Any}()
        t0 = time_ns()
        bytes = @allocated begin
            series[] = _run_series(sim, skies)
        end
        push!(times_ms, (time_ns() - t0) / 1e6)
        push!(alloc_mb, Float64(bytes) / 1024.0^2)
        totals = _series_totals(series[])
        info = _cache_info(sim)
    end
    return merge(
        (
            status="ok",
            error_type="",
            error_message="",
            median_ms=median(times_ms),
            min_ms=minimum(times_ms),
            max_ms=maximum(times_ms),
            median_alloc_mb=median(alloc_mb),
        ),
        info,
        totals,
    )
end

function _empty_measurement(status::AbstractString; error_type::AbstractString="", error_message::AbstractString="", resolved_backend::AbstractString="", geometry_mode=:none)
    return (
        status=String(status),
        error_type=String(error_type),
        error_message=String(error_message),
        median_ms=NaN,
        min_ms=NaN,
        max_ms=NaN,
        median_alloc_mb=NaN,
        resolved_backend=String(resolved_backend),
        geometry_mode=geometry_mode,
        raycore_chunked=false,
        raycore_tlas_instances=0,
        raycore_tlas_geometries=0,
        raycore_expanded_face_count=0,
        incident_par=NaN,
        incident_nir=NaN,
        absorbed_par=NaN,
        absorbed_nir=NaN,
    )
end

function _relative_delta(value, reference)
    isfinite(value) || return NaN
    isfinite(reference) || return NaN
    reference == 0.0 && return value == 0.0 ? 0.0 : Inf
    return 100.0 * (value - reference) / abs(reference)
end

function _tsv_cell(value)
    return replace(string(value), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _write_results(rows)
    isempty(BENCH_OUTPUT) && return rows
    columns = Symbol[]
    for row in rows, key in keys(row)
        key in columns || push!(columns, key)
    end
    dir = dirname(BENCH_OUTPUT)
    isempty(dir) || dir == "." || mkpath(dir)
    open(BENCH_OUTPUT, "w") do io
        println(io, join(string.(columns), '\t'))
        for row in rows
            println(io, join((_tsv_cell(key in keys(row) ? getproperty(row, key) : missing) for key in columns), '\t'))
        end
    end
    @info "Wrote artifact scene benchmark results" path=BENCH_OUTPUT rows=length(rows)
    return rows
end

function _print_header(scene_root, cases)
    println("ArchimedLight artifact scene benchmark")
    println("scene root: ", scene_root)
    println("scenes: ", join(BENCH_SCENES, ", "))
    println("backends: ", join((case.name for case in cases), ", "))
    println("directions: ", join(BENCH_DIRECTIONS, ", "))
    println("pixels cm: ", join(BENCH_PIXELS_CM, ", "))
    println("scattering: ", join(BENCH_SCATTERING, ", "))
    println("cache_radiation: ", join(BENCH_CACHE_RADIATION, ", "))
    println("steps: ", join(BENCH_STEPS, ", "))
    println("raycore toric traversal: ", join(BENCH_RAYCORE_TORIC_TRAVERSAL, ", "))
    println("samples: ", BENCH_SAMPLES)
    println("dry run: ", BENCH_DRY_RUN)
    println("max pixel cells: ", BENCH_MAX_PIXEL_CELLS == 0 ? "disabled" : string(BENCH_MAX_PIXEL_CELLS))
    println("output: ", BENCH_OUTPUT)
    println()
    @printf("%-20s %-16s %-10s %-7s %-8s %-5s %-5s %5s %12s %-24s %10s %10s %10s\n",
        "scene", "backend", "toric", "dirs", "pixelcm", "scat", "cache", "steps", "pixels", "status", "median", "alloc", "PAR d%")
end

function main()
    BENCH_SAMPLES > 0 || error("ARCHIMEDLIGHT_ARTIFACT_BENCH_SAMPLES must be positive.")
    scene_root = _benchmark_scenes_root()
    scene_specs = _scene_specs(scene_root)
    cases = _backend_cases()
    models = _benchmark_models()
    rows = NamedTuple[]
    baselines = Dict{Tuple{String,Int,Float64,Bool,Bool,Int},NamedTuple}()

    _print_header(scene_root, cases)

    for spec in scene_specs
        scene = BENCH_DRY_RUN ? nothing : ArchimedLight.read_scene(spec.path)
        scene_name = spec.name
        for directions in BENCH_DIRECTIONS, pixel_cm in BENCH_PIXELS_CM, scattering in BENCH_SCATTERING,
            cache_radiation in BENCH_CACHE_RADIATION, step_count in BENCH_STEPS

            pixel_m = pixel_cm / 100.0
            options = _options(pixel_m, directions, scattering, cache_radiation)
            metrics = BENCH_DRY_RUN ? _dry_run_scene_metrics(spec.path, pixel_m) : _scene_metrics(scene, models, pixel_m)
            skies = BENCH_DRY_RUN ? ArchimedLight.SkyState[] : _sky_series(step_count)
            combo_key = (scene_name, directions, pixel_cm, scattering, cache_radiation, step_count)
            combo_rows = NamedTuple[]

            for case in cases
                measurement =
                    if BENCH_DRY_RUN
                        _empty_measurement("dry_run")
                    elseif BENCH_MAX_PIXEL_CELLS > 0 && metrics.pixel_cells_estimate > BENCH_MAX_PIXEL_CELLS
                        _empty_measurement(
                            "skipped_resource_limit";
                            error_message="estimated pixel cells $(metrics.pixel_cells_estimate) exceed limit $BENCH_MAX_PIXEL_CELLS",
                        )
                    elseif !case.available
                        _empty_measurement("skipped_backend_unavailable"; error_message=case.unavailable_reason)
                    else
                        try
                            _measure(scene, models, options, case, skies)
                        catch err
                            _empty_measurement(
                                "error";
                                error_type=string(nameof(typeof(err))),
                                error_message=sprint(showerror, err),
                            )
                        end
                    end

                row = merge(
                    (
                        scene=scene_name,
                        scene_path=spec.path,
                        backend=case.name,
                        toric_traversal=case.toric_traversal,
                        directions=directions,
                        pixel_size_cm=pixel_cm,
                        pixel_size_m=pixel_m,
                        scattering=scattering,
                        cache_radiation=cache_radiation,
                        steps=step_count,
                        samples=BENCH_SAMPLES,
                    ),
                    metrics,
                    measurement,
                )
                case.name == "normal_cpu" && row.status == "ok" && (baselines[combo_key] = row)
                push!(combo_rows, row)
            end

            base = get(baselines, combo_key, nothing)
            for row in combo_rows
                par_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.incident_par, base.incident_par) : NaN
                nir_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.incident_nir, base.incident_nir) : NaN
                absorbed_par_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.absorbed_par, base.absorbed_par) : NaN
                absorbed_nir_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.absorbed_nir, base.absorbed_nir) : NaN
                full_row = merge(
                    row,
                    (
                        par_delta_percent=par_delta,
                        nir_delta_percent=nir_delta,
                        absorbed_par_delta_percent=absorbed_par_delta,
                        absorbed_nir_delta_percent=absorbed_nir_delta,
                    ),
                )
                push!(rows, full_row)
                @printf("%-20s %-16s %-10s %7d %8.3f %-5s %-5s %5d %12d %-24s %10.3f %10.3f %10.3f\n",
                    full_row.scene,
                    full_row.backend,
                    string(full_row.toric_traversal),
                    full_row.directions,
                    full_row.pixel_size_cm,
                    string(full_row.scattering),
                    string(full_row.cache_radiation),
                    full_row.steps,
                    full_row.pixel_cells_estimate,
                    full_row.status,
                    full_row.median_ms,
                    full_row.median_alloc_mb,
                    full_row.par_delta_percent,
                )
                full_row.status == "error" && @info "Artifact scene benchmark row failed" scene=full_row.scene backend=full_row.backend message=full_row.error_message
            end
        end
    end

    return _write_results(rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
