using ArchimedLight
using Artifacts: artifact_hash, artifact_path
using Dates
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
const DEFAULT_BACKENDS = ["normal_cpu", "rastergpu_metal"]
const DEFAULT_DIRECTIONS = [16, 46, 136]
const DEFAULT_PIXELS_CM = [0.1, 0.5, 1.0, 10.0]
const DEFAULT_BOOL = [false, true]
const DEFAULT_STEPS = [1, 12, 24, 48]

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

function _int_or_auto_env(name::AbstractString, default)
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    value = lowercase(strip(raw))
    value == "auto" && return :auto
    return parse(Int, value)
end

_timestamp() = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")

function _log_progress(message::AbstractString; kwargs...)
    print("[", _timestamp(), "] ", message)
    for (key, value) in kwargs
        print(" ", key, "=", value)
    end
    println()
    flush(stdout)
    return nothing
end

function _case_label(case)
    return string(case.name, case.toric_traversal == :not_applicable ? "" : "[$(case.toric_traversal)]")
end

const BENCH_SCENES = _split_string_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_SCENES", DEFAULT_SCENES)
const BENCH_BACKENDS = _split_string_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_BACKENDS", DEFAULT_BACKENDS)
const BENCH_DIRECTIONS = _split_int_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_DIRECTIONS", DEFAULT_DIRECTIONS)
const BENCH_PIXELS_CM = _split_float_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_PIXELS_CM", DEFAULT_PIXELS_CM)
const BENCH_SCATTERING = _split_bool_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_SCATTERING", DEFAULT_BOOL)
const BENCH_CACHE_RADIATION = _split_bool_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_CACHE_RADIATION", DEFAULT_BOOL)
const BENCH_STEPS = _split_int_env("ARCHIMEDLIGHT_ARTIFACT_BENCH_STEPS", DEFAULT_STEPS)
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_SAMPLES", "1"))
const BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_WARMUPS", "0"))
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_OUTPUT", DEFAULT_OUTPUT)
const BENCH_DRY_RUN = _env_bool("ARCHIMEDLIGHT_ARTIFACT_BENCH_DRY_RUN", false)
const BENCH_MAX_PIXEL_CELLS = parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_MAX_PIXEL_CELLS", "25000000"))
const BENCH_STEP_SECONDS = parse(Float64, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_STEP_SECONDS", "1800.0"))
const BENCH_MAX_HITS = _int_or_auto_env("ARCHIMEDLIGHT_BENCH_MAX_HITS", :auto)
const BENCH_WORKGROUPSIZE = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", "256"))
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto"))
const BENCH_RASTERGPU_TILE_SIZE = _int_or_auto_env("ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_SIZE", :auto)
const BENCH_RASTERGPU_TILE_FACE_CAPACITY =
    _int_or_auto_env("ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_FACE_CAPACITY", :auto)
const BENCH_RASTERGPU_AUTO_RETRIES = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_RASTERGPU_AUTO_RETRIES", "4"))
const BENCH_VALIDATE = _env_bool("ARCHIMEDLIGHT_BENCH_VALIDATE", true)
const BENCH_SKIP_HUGE_CASES = _env_bool("ARCHIMEDLIGHT_ARTIFACT_BENCH_SKIP_HUGE_CASES", true)
const BENCH_SKIP_HUGE_MIN_DIRECTIONS =
    parse(Int, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_SKIP_HUGE_MIN_DIRECTIONS", "136"))
const BENCH_SKIP_HUGE_MIN_HOURS =
    parse(Float64, get(ENV, "ARCHIMEDLIGHT_ARTIFACT_BENCH_SKIP_HUGE_MIN_HOURS", "2.0"))
const RESULT_COLUMNS = Symbol[
    :scene,
    :scene_path,
    :backend,
    :toric_traversal,
    :directions,
    :pixel_size_cm,
    :pixel_size_m,
    :scattering,
    :cache_radiation,
    :steps,
    :samples,
    :nodes,
    :faces,
    :objects,
    :domain_x,
    :domain_y,
    :pixel_cells_estimate,
    :plot_nx_estimate,
    :plot_ny_estimate,
    :status,
    :error_type,
    :error_message,
    :median_ms,
    :min_ms,
    :max_ms,
    :median_alloc_mb,
    :resolved_backend,
    :geometry_mode,
    :rastergpu_fused_dense,
    :rastergpu_max_hits_per_pixel,
    :rastergpu_tile_size,
    :rastergpu_tile_face_capacity,
    :rastergpu_edge_accumulation,
    :rastergpu_retry_count,
    :incident_par,
    :incident_nir,
    :absorbed_par,
    :absorbed_nir,
    :par_delta_percent,
    :nir_delta_percent,
    :absorbed_par_delta_percent,
    :absorbed_nir_delta_percent,
]

function _preload_metal_if_requested()
    selected = lowercase.(BENCH_BACKENDS)
    "rastergpu_metal" in selected || return nothing
    Sys.isapple() && Sys.ARCH == :aarch64 || return nothing
    Base.find_package("Metal") === nothing && return nothing
    return Base.require(Main, :Metal)
end

const BENCH_METAL_MODULE = _preload_metal_if_requested()

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
        metal_mod = BENCH_METAL_MODULE === nothing ? Base.require(Main, :Metal) : BENCH_METAL_MODULE
        array_type = getproperty(metal_mod, :MtlArray)
        dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
        backend = Base.invokelatest(KernelAbstractions.get_backend, dev_array)
        return backend, ""
    catch err
        return nothing, "Could not initialize Metal backend: $(sprint(showerror, err))"
    end
end

_setting_is_auto(value) = value === :auto
_setting_value(value, fallback::Int) = _setting_is_auto(value) ? fallback : Int(value)

function _rastergpu_auto_requested()
    return _setting_is_auto(BENCH_MAX_HITS) ||
           _setting_is_auto(BENCH_RASTERGPU_TILE_SIZE) ||
           _setting_is_auto(BENCH_RASTERGPU_TILE_FACE_CAPACITY)
end

function _next_bucket(value::Integer)
    buckets = (32, 64, 96, 128, 160, 192, 256, 384, 512, 768, 1024, 1536, 2048, 4096)
    target = max(1, ceil(Int, 1.15 * value))
    for bucket in buckets
        bucket >= target && return bucket
    end
    return target
end

function _observed_capacity(message::AbstractString, field::AbstractString)
    m = match(Regex(field * "=(\\d+)"), message)
    m === nothing && return nothing
    return parse(Int, m.captures[1])
end

function _rastergpu_initial_config(scene_name::AbstractString, metrics)
    dense_canopy = occursin("elaeis", lowercase(scene_name))
    default_max_hits = dense_canopy ? 192 : 128
    default_tile_face_capacity = dense_canopy ? 256 : 128
    return (
        max_hits_per_pixel=_setting_value(BENCH_MAX_HITS, default_max_hits),
        tile_size=_setting_value(BENCH_RASTERGPU_TILE_SIZE, 1),
        tile_face_capacity=_setting_value(BENCH_RASTERGPU_TILE_FACE_CAPACITY, default_tile_face_capacity),
        edge_accumulation=BENCH_EDGE_ACCUMULATION,
    )
end

function _rastergpu_config_with_retry(config, message::AbstractString)
    observed_tile_capacity = _observed_capacity(message, "observed_max_tile_face_capacity")
    observed_hits = _observed_capacity(message, "observed_max_hits_per_pixel")
    next_tile_capacity = config.tile_face_capacity
    next_max_hits = config.max_hits_per_pixel

    if observed_tile_capacity !== nothing && _setting_is_auto(BENCH_RASTERGPU_TILE_FACE_CAPACITY)
        next_tile_capacity = max(config.tile_face_capacity, _next_bucket(observed_tile_capacity))
    end
    if observed_hits !== nothing && _setting_is_auto(BENCH_MAX_HITS)
        next_max_hits = max(config.max_hits_per_pixel, _next_bucket(observed_hits))
    end

    if next_tile_capacity == config.tile_face_capacity && next_max_hits == config.max_hits_per_pixel
        return nothing
    end
    return (
        max_hits_per_pixel=next_max_hits,
        tile_size=config.tile_size,
        tile_face_capacity=next_tile_capacity,
        edge_accumulation=config.edge_accumulation,
    )
end

function _make_rastergpu_backends(metal, config)
    ib = ArchimedLight.RasterGPUBackend(
        backend=metal,
        max_hits_per_pixel=config.max_hits_per_pixel,
        workgroupsize=BENCH_WORKGROUPSIZE,
        tile_size=config.tile_size,
        tile_face_capacity=config.tile_face_capacity,
        edge_accumulation=config.edge_accumulation,
    )
    return ib, ArchimedLight.RasterGPUScatteringBackend(
        ib;
        edge_accumulation=config.edge_accumulation,
    )
end

function _with_rastergpu_config(case, config)
    ib, sb = _make_rastergpu_backends(case.metal_backend, config)
    return merge(
        case,
        (
            interception_backend=ib,
            scattering_backend=sb,
            rastergpu_max_hits_per_pixel=config.max_hits_per_pixel,
            rastergpu_tile_size=config.tile_size,
            rastergpu_tile_face_capacity=config.tile_face_capacity,
            rastergpu_edge_accumulation=config.edge_accumulation,
        ),
    )
end

function _rastergpu_result_config(config, retry_count::Int; fused_dense=missing)
    return (
        rastergpu_max_hits_per_pixel=config.max_hits_per_pixel,
        rastergpu_tile_size=config.tile_size,
        rastergpu_tile_face_capacity=config.tile_face_capacity,
        rastergpu_edge_accumulation=config.edge_accumulation,
        rastergpu_fused_dense=fused_dense,
        rastergpu_retry_count=retry_count,
    )
end

function _backend_cases()
    selected = lowercase.(BENCH_BACKENDS)
    unknown = setdiff(selected, DEFAULT_BACKENDS)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_ARTIFACT_BENCH_BACKENDS value(s): $(join(sort!(unknown), ", ")). " *
        "Use normal_cpu, rastergpu_metal, or all.",
    )

    if BENCH_DRY_RUN
        cases = NamedTuple[]
        for name in selected
            push!(
                cases,
                (
                    name=name,
                    available=true,
                    unavailable_reason="",
                    interception_backend=nothing,
                    scattering_backend=nothing,
                    toric_traversal=:not_applicable,
                    metal_backend=nothing,
                ),
            )
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
                    metal_backend=nothing,
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
                        metal_backend=nothing,
                    ),
                )
            else
                config = _rastergpu_initial_config("", nothing)
                ib, sb = _make_rastergpu_backends(metal, config)
                push!(
                    cases,
                    (
                        name=name,
                        available=true,
                        unavailable_reason="",
                        interception_backend=ib,
                        scattering_backend=sb,
                        toric_traversal=:not_applicable,
                        metal_backend=metal,
                        rastergpu_max_hits_per_pixel=config.max_hits_per_pixel,
                        rastergpu_tile_size=config.tile_size,
                        rastergpu_tile_face_capacity=config.tile_face_capacity,
                        rastergpu_edge_accumulation=config.edge_accumulation,
                    ),
                )
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
        rastergpu_fused_dense=missing,
    )
    if hasproperty(cache, :rastergpu_data) && cache.rastergpu_data !== nothing
        fused_dense = try
            ArchimedLight._rastergpu_use_fused_dense_edges(cache.rastergpu_data)
        catch
            missing
        end
        return (
            resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
            geometry_mode=:raster_gpu,
            rastergpu_fused_dense=fused_dense,
        )
    end
    return (
        resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
        geometry_mode=:raster_cpu,
        rastergpu_fused_dense=missing,
    )
end

function _measure(scene, models, options, case, skies; row_context=nothing)
    for warmup in 1:BENCH_WARMUPS
        row_context !== nothing && _log_progress(
            "warmup start";
            row_context...,
            warmup=warmup,
            warmups=BENCH_WARMUPS,
        )
        sim = _new_simulation(scene, models, options, case)
        _run_series(sim, skies)
        row_context !== nothing && _log_progress(
            "warmup done";
            row_context...,
            warmup=warmup,
            warmups=BENCH_WARMUPS,
        )
    end

    times_ms = Float64[]
    alloc_mb = Float64[]
    totals = (incident_par=NaN, incident_nir=NaN, absorbed_par=NaN, absorbed_nir=NaN)
    info = (
        resolved_backend="",
        geometry_mode=:unknown,
        rastergpu_fused_dense=missing,
    )
    for sample in 1:BENCH_SAMPLES
        GC.gc()
        sim = _new_simulation(scene, models, options, case)
        series = Ref{Any}()
        row_context !== nothing && _log_progress(
            "sample start";
            row_context...,
            sample=sample,
            samples=BENCH_SAMPLES,
        )
        t0 = time_ns()
        bytes = @allocated begin
            series[] = _run_series(sim, skies)
        end
        elapsed_ms = (time_ns() - t0) / 1e6
        push!(times_ms, elapsed_ms)
        push!(alloc_mb, Float64(bytes) / 1024.0^2)
        totals = _series_totals(series[])
        info = _cache_info(sim)
        row_context !== nothing && _log_progress(
            "sample done";
            row_context...,
            sample=sample,
            samples=BENCH_SAMPLES,
            elapsed_ms=round(elapsed_ms; digits=3),
            alloc_mb=round(Float64(bytes) / 1024.0^2; digits=3),
        )
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

function _measure_rastergpu_adaptive(
    scene,
    models,
    options,
    case,
    skies,
    scene_name,
    metrics;
    initial_config=nothing,
    row_context=nothing,
)
    config = initial_config === nothing ? _rastergpu_initial_config(scene_name, metrics) : initial_config
    last_error = nothing
    last_message = ""

    for retry_count in 0:BENCH_RASTERGPU_AUTO_RETRIES
        configured_case = _with_rastergpu_config(case, config)
        row_context !== nothing && _log_progress(
            "rastergpu config attempt";
            row_context...,
            retry=retry_count,
            max_hits=config.max_hits_per_pixel,
            tile_size=config.tile_size,
            tile_face_capacity=config.tile_face_capacity,
            edge_accumulation=config.edge_accumulation,
        )
        try
            measurement = _measure(scene, models, options, configured_case, skies; row_context=row_context)
            return merge(
                measurement,
                _rastergpu_result_config(
                    config,
                    retry_count;
                    fused_dense=measurement.rastergpu_fused_dense,
                ),
            )
        catch err
            last_error = err
            last_message = sprint(showerror, err)
            next_config = _rastergpu_config_with_retry(config, last_message)
            if next_config === nothing || retry_count == BENCH_RASTERGPU_AUTO_RETRIES
                row_context !== nothing && _log_progress(
                    "rastergpu config failed";
                    row_context...,
                    retry=retry_count,
                    max_hits=config.max_hits_per_pixel,
                    tile_size=config.tile_size,
                    tile_face_capacity=config.tile_face_capacity,
                    error=last_message,
                )
                return merge(
                    _empty_measurement(
                        "error";
                        error_type=string(nameof(typeof(err))),
                        error_message=last_message,
                    ),
                    _rastergpu_result_config(config, retry_count),
                )
            end
            row_context !== nothing && _log_progress(
                "rastergpu retry";
                row_context...,
                retry=retry_count + 1,
                previous_max_hits=config.max_hits_per_pixel,
                previous_tile_face_capacity=config.tile_face_capacity,
                next_max_hits=next_config.max_hits_per_pixel,
                next_tile_face_capacity=next_config.tile_face_capacity,
                error=last_message,
            )
            config = next_config
        end
    end

    return merge(
        _empty_measurement(
            "error";
            error_type=last_error === nothing ? "" : string(nameof(typeof(last_error))),
            error_message=last_message,
        ),
        _rastergpu_result_config(config, BENCH_RASTERGPU_AUTO_RETRIES),
    )
end

function _measure_case(
    scene,
    models,
    options,
    case,
    skies,
    scene_name,
    metrics;
    rastergpu_config_cache=nothing,
    rastergpu_config_key=nothing,
    row_context=nothing,
)
    if case.name == "rastergpu_metal" && _rastergpu_auto_requested()
        initial_config =
            rastergpu_config_cache !== nothing && rastergpu_config_key !== nothing ?
            get(rastergpu_config_cache, rastergpu_config_key, nothing) :
            nothing
        measurement = _measure_rastergpu_adaptive(
            scene,
            models,
            options,
            case,
            skies,
            scene_name,
            metrics;
            initial_config=initial_config,
            row_context=row_context,
        )
        if measurement.status == "ok" &&
           rastergpu_config_cache !== nothing &&
           rastergpu_config_key !== nothing
            rastergpu_config_cache[rastergpu_config_key] = (
                max_hits_per_pixel=measurement.rastergpu_max_hits_per_pixel,
                tile_size=measurement.rastergpu_tile_size,
                tile_face_capacity=measurement.rastergpu_tile_face_capacity,
                edge_accumulation=measurement.rastergpu_edge_accumulation,
            )
        end
        return measurement
    end

    measurement = _measure(scene, models, options, case, skies; row_context=row_context)
    if case.name == "rastergpu_metal"
        config = (
            max_hits_per_pixel=case.rastergpu_max_hits_per_pixel,
            tile_size=case.rastergpu_tile_size,
            tile_face_capacity=case.rastergpu_tile_face_capacity,
            edge_accumulation=case.rastergpu_edge_accumulation,
        )
        return merge(
            measurement,
            _rastergpu_result_config(config, 0; fused_dense=measurement.rastergpu_fused_dense),
        )
    end
    return measurement
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
        rastergpu_fused_dense=missing,
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

function _huge_case_skip_reason(directions::Int, scattering::Bool, cache_radiation::Bool, step_count::Int)
    BENCH_SKIP_HUGE_CASES || return ""
    directions >= BENCH_SKIP_HUGE_MIN_DIRECTIONS || return ""
    simulated_hours = step_count * BENCH_STEP_SECONDS / 3600.0
    simulated_hours >= BENCH_SKIP_HUGE_MIN_HOURS || return ""
    (scattering && cache_radiation) || !cache_radiation || return ""
    return "skipped by ARCHIMEDLIGHT_ARTIFACT_BENCH_SKIP_HUGE_CASES: " *
           "scattering=$scattering, cache_radiation=$cache_radiation, directions=$directions, " *
           "simulated_hours=$(round(simulated_hours; digits=3))"
end

function _with_deltas(row, base)
    par_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.incident_par, base.incident_par) : NaN
    nir_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.incident_nir, base.incident_nir) : NaN
    absorbed_par_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.absorbed_par, base.absorbed_par) : NaN
    absorbed_nir_delta = row.status == "ok" && base !== nothing ? _relative_delta(row.absorbed_nir, base.absorbed_nir) : NaN
    return merge(
        row,
        (
            par_delta_percent=par_delta,
            nir_delta_percent=nir_delta,
            absorbed_par_delta_percent=absorbed_par_delta,
            absorbed_nir_delta_percent=absorbed_nir_delta,
        ),
    )
end

function _print_result_row(row)
    @printf("%-20s %-16s %-10s %7d %8.3f %-5s %-5s %5d %12d %-24s %10.3f %10.3f %10.3f\n",
        row.scene,
        row.backend,
        string(row.toric_traversal),
        row.directions,
        row.pixel_size_cm,
        string(row.scattering),
        string(row.cache_radiation),
        row.steps,
        row.pixel_cells_estimate,
        row.status,
        row.median_ms,
        row.median_alloc_mb,
        row.par_delta_percent,
    )
    flush(stdout)
    return nothing
end

function _tsv_cell(value)
    return replace(string(value), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function _write_result_row(io, row)
    println(
        io,
        join((_tsv_cell(key in keys(row) ? getproperty(row, key) : missing) for key in RESULT_COLUMNS), '\t'),
    )
    return nothing
end

function _prepare_results_file()
    isempty(BENCH_OUTPUT) && return nothing
    dir = dirname(BENCH_OUTPUT)
    isempty(dir) || dir == "." || mkpath(dir)
    open(BENCH_OUTPUT, "w") do io
        println(io, join(string.(RESULT_COLUMNS), '\t'))
    end
    @info "Initialized artifact scene benchmark results" path=BENCH_OUTPUT
    return nothing
end

function _append_result(row)
    isempty(BENCH_OUTPUT) && return row
    open(BENCH_OUTPUT, "a") do io
        _write_result_row(io, row)
        flush(io)
    end
    return row
end

function _write_results(rows)
    isempty(BENCH_OUTPUT) && return rows
    dir = dirname(BENCH_OUTPUT)
    isempty(dir) || dir == "." || mkpath(dir)
    open(BENCH_OUTPUT, "w") do io
        println(io, join(string.(RESULT_COLUMNS), '\t'))
        for row in rows
            _write_result_row(io, row)
        end
    end
    @info "Wrote artifact scene benchmark results" path=BENCH_OUTPUT rows=length(rows)
    return rows
end

function _print_header(scene_root, cases)
    _log_progress("benchmark start"; output=BENCH_OUTPUT)
    println("ArchimedLight artifact scene benchmark")
    println("scene root: ", scene_root)
    println("scenes: ", join(BENCH_SCENES, ", "))
    println("backends: ", join((case.name for case in cases), ", "))
    println("directions: ", join(BENCH_DIRECTIONS, ", "))
    println("pixels cm: ", join(BENCH_PIXELS_CM, ", "))
    println("scattering: ", join(BENCH_SCATTERING, ", "))
    println("cache_radiation: ", join(BENCH_CACHE_RADIATION, ", "))
    println("steps: ", join(BENCH_STEPS, ", "))
    println("samples: ", BENCH_SAMPLES)
    println("dry run: ", BENCH_DRY_RUN)
    println("max pixel cells: ", BENCH_MAX_PIXEL_CELLS == 0 ? "disabled" : string(BENCH_MAX_PIXEL_CELLS))
    println("output: ", BENCH_OUTPUT)
    println()
    @printf("%-20s %-16s %-10s %-7s %-8s %-5s %-5s %5s %12s %-24s %10s %10s %10s\n",
        "scene", "backend", "toric", "dirs", "pixelcm", "scat", "cache", "steps", "pixels", "status", "median", "alloc", "PAR d%")
end

function main()
    _log_progress("initializing benchmark")
    BENCH_SAMPLES > 0 || error("ARCHIMEDLIGHT_ARTIFACT_BENCH_SAMPLES must be positive.")
    scene_root = _benchmark_scenes_root()
    scene_specs = _scene_specs(scene_root)
    cases = _backend_cases()
    models = _benchmark_models()
    rows = NamedTuple[]
    baselines = Dict{Tuple{String,Int,Float64,Bool,Bool,Int},NamedTuple}()
    rastergpu_config_cache = Dict{Tuple{String,Int,Float64,Bool,Bool},Any}()

    _print_header(scene_root, cases)
    _prepare_results_file()

    for spec in scene_specs
        _log_progress("scene load start"; scene=spec.name, path=spec.path)
        scene = BENCH_DRY_RUN ? nothing : ArchimedLight.read_scene(spec.path)
        _log_progress("scene load done"; scene=spec.name)
        scene_name = spec.name
        for directions in BENCH_DIRECTIONS, pixel_cm in BENCH_PIXELS_CM, scattering in BENCH_SCATTERING,
            cache_radiation in BENCH_CACHE_RADIATION, step_count in BENCH_STEPS

            pixel_m = pixel_cm / 100.0
            options = _options(pixel_m, directions, scattering, cache_radiation)
            metrics = BENCH_DRY_RUN ? _dry_run_scene_metrics(spec.path, pixel_m) : _scene_metrics(scene, models, pixel_m)
            skies = BENCH_DRY_RUN ? ArchimedLight.SkyState[] : _sky_series(step_count)
            combo_key = (scene_name, directions, pixel_cm, scattering, cache_radiation, step_count)
            rastergpu_config_key = (scene_name, directions, pixel_cm, scattering, cache_radiation)
            skip_reason = _huge_case_skip_reason(directions, scattering, cache_radiation, step_count)
            combo_row_count = 0
            _log_progress(
                "combo start";
                scene=scene_name,
                directions=directions,
                pixel_cm=pixel_cm,
                scattering=scattering,
                cache_radiation=cache_radiation,
                steps=step_count,
                pixel_cells=metrics.pixel_cells_estimate,
            )

            for case in cases
                row_context = (
                    scene=scene_name,
                    backend=_case_label(case),
                    directions=directions,
                    pixel_cm=pixel_cm,
                    scattering=scattering,
                    cache_radiation=cache_radiation,
                    steps=step_count,
                    pixel_cells=metrics.pixel_cells_estimate,
                )
                _log_progress("row start"; row_context...)
                measurement =
                    if BENCH_DRY_RUN
                        _empty_measurement("dry_run")
                    elseif !isempty(skip_reason)
                        _empty_measurement("skipped_huge_case"; error_message=skip_reason)
                    elseif BENCH_MAX_PIXEL_CELLS > 0 && metrics.pixel_cells_estimate > BENCH_MAX_PIXEL_CELLS
                        _empty_measurement(
                            "skipped_resource_limit";
                            error_message="estimated pixel cells $(metrics.pixel_cells_estimate) exceed limit $BENCH_MAX_PIXEL_CELLS",
                        )
                    elseif !case.available
                        _empty_measurement("skipped_backend_unavailable"; error_message=case.unavailable_reason)
                    else
                        try
                            _measure_case(
                                scene,
                                models,
                                options,
                                case,
                                skies,
                                scene_name,
                                metrics;
                                rastergpu_config_cache=rastergpu_config_cache,
                                rastergpu_config_key=rastergpu_config_key,
                                row_context=row_context,
                            )
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
                base = get(baselines, combo_key, nothing)
                full_row = _with_deltas(row, base)
                push!(rows, full_row)
                combo_row_count += 1
                _print_result_row(full_row)
                _append_result(full_row)
                full_row.status == "error" && @info "Artifact scene benchmark row failed" scene=full_row.scene backend=full_row.backend message=full_row.error_message
                _log_progress(
                    "row done";
                    row_context...,
                    status=row.status,
                    median_ms=round(row.median_ms; digits=3),
                )
            end
            _log_progress(
                "combo done";
                scene=scene_name,
                directions=directions,
                pixel_cm=pixel_cm,
                scattering=scattering,
                cache_radiation=cache_radiation,
                steps=step_count,
                rows=combo_row_count,
            )
        end
    end

    result = _write_results(rows)
    _log_progress("benchmark done"; rows=length(rows), output=BENCH_OUTPUT)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
