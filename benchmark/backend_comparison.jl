const _BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", ""))
const _BENCH_WANTS_METAL = _BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")
const _BENCH_CUDA_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_CUDA", ""))
const _BENCH_WANTS_CUDA = _BENCH_CUDA_FLAG in ("1", "true", "yes", "on", "required", "force")
const _BENCH_ONEAPI_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_ONEAPI", ""))
const _BENCH_WANTS_ONEAPI = _BENCH_ONEAPI_FLAG in ("1", "true", "yes", "on", "required", "force")
const _BENCH_AMDGPU_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_AMDGPU", ""))
const _BENCH_WANTS_AMDGPU = _BENCH_AMDGPU_FLAG in ("1", "true", "yes", "on", "required", "force")

if _BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end
if _BENCH_WANTS_CUDA && Base.find_package("CUDA") !== nothing
    Base.require(Main, :CUDA)
end
if _BENCH_WANTS_ONEAPI && Base.find_package("oneAPI") !== nothing
    Base.require(Main, :oneAPI)
end
if _BENCH_WANTS_AMDGPU && Base.find_package("AMDGPU") !== nothing
    Base.require(Main, :AMDGPU)
end

include(joinpath(@__DIR__, "..", "scripts", "raycore_device_smoke.jl"))

using Dates
using Printf
using Statistics

function _split_env_values(name::AbstractString, parse_value)
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return nothing
    return [parse_value(strip(value)) for value in split(raw, ',') if !isempty(strip(value))]
end

const BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WARMUPS", "2"))
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_SAMPLES", "5"))
const BENCH_SECTORS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_SECTORS", "46"))
const BENCH_PIXEL_SIZE = parse(Float64, get(ENV, "ARCHIMEDLIGHT_BENCH_PIXEL_SIZE", "0.01"))
const BENCH_WORKLOAD = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_WORKLOAD", "smoke"))
const BENCH_PLOT_PAVING = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_PLOT_PAVING", "2500"))
const BENCH_BACKENDS = split(lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_BACKENDS", "all")), ',')
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "32"))
const BENCH_MAX_PRECHUNK_INSTANCES = haskey(ENV, "ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES") ?
                                     parse(Int, ENV["ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES"]) :
                                     nothing
const BENCH_WORKGROUPSIZE = haskey(ENV, "ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE") ?
                            parse(Int, ENV["ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE"]) :
                            nothing
const BENCH_EDGE_ACCUMULATION = Symbol(get(ENV, "ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION", "auto"))
const BENCH_EXECUTION = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_EXECUTION", "oneshot"))
const BENCH_COMPONENTS = split(lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_COMPONENTS", "first_order,scattering,topology,propagation")), ',')

const BENCH_SECTOR_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_SECTOR_SWEEP", x -> parse(Int, x)), [BENCH_SECTORS])
const BENCH_PIXEL_SIZE_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_PIXEL_SIZE_SWEEP", x -> parse(Float64, x)), [BENCH_PIXEL_SIZE])
const BENCH_PLOT_PAVING_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_PLOT_PAVING_SWEEP", x -> parse(Int, x)), [BENCH_PLOT_PAVING])
const BENCH_MAX_HIT_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_MAX_HITS_SWEEP", x -> parse(Int, x)), [BENCH_MAX_HITS])
const BENCH_MAX_PRECHUNK_INSTANCE_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES_SWEEP", x -> parse(Int, x)), [BENCH_MAX_PRECHUNK_INSTANCES])
const BENCH_WORKGROUPSIZE_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE_SWEEP", x -> parse(Int, x)), [BENCH_WORKGROUPSIZE])
const BENCH_EDGE_ACCUMULATION_VALUES =
    something(_split_env_values("ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION_SWEEP", Symbol), [BENCH_EDGE_ACCUMULATION])
const BENCH_HAS_SECTOR_OVERRIDE = haskey(ENV, "ARCHIMEDLIGHT_BENCH_SECTORS") || haskey(ENV, "ARCHIMEDLIGHT_BENCH_SECTOR_SWEEP")
const BENCH_HAS_PIXEL_SIZE_OVERRIDE = haskey(ENV, "ARCHIMEDLIGHT_BENCH_PIXEL_SIZE") || haskey(ENV, "ARCHIMEDLIGHT_BENCH_PIXEL_SIZE_SWEEP")

_workgroupsize_label(workgroupsize) = workgroupsize === nothing ? "default" : string(workgroupsize)
_workgroupsize_values_label(values) = join((_workgroupsize_label(value) for value in values), ',')
_max_prechunk_instances_label(value) = value === nothing ? "default" : string(value <= 0 ? "disabled" : value)
_max_prechunk_instances_values_label(values) = join((_max_prechunk_instances_label(value) for value in values), ',')

function _validate_edge_accumulation_values(values)
    allowed = Set([:auto, :sparse_host_reduce, :dense_atomic])
    unknown = setdiff(Set(values), allowed)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_BENCH_EDGE_ACCUMULATION value(s): $(join(string.(sort!(collect(unknown))), ',')). " *
        "Use auto, sparse_host_reduce, or dense_atomic.",
    )
    return values
end

_validate_edge_accumulation_values(BENCH_EDGE_ACCUMULATION_VALUES)

function _bench_components()
    allowed = Set(["first_order", "scattering", "topology", "propagation", "integration", "stack_trace"])
    components = Set(component for component in strip.(BENCH_COMPONENTS) if !isempty(component))
    "all" in components && return allowed
    unknown = setdiff(components, allowed)
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_BENCH_COMPONENTS=$(join(sort!(collect(unknown)), ',')). " *
        "Use all, first_order, scattering, topology, propagation, integration, or stack_trace.",
    )
    isempty(components) && error(
        "ARCHIMEDLIGHT_BENCH_COMPONENTS selected no components. " *
        "Use all, first_order, scattering, topology, propagation, integration, or stack_trace.",
    )
    return components
end

function _benchmark_meteo_row(;
    date::Dates.Date=Dates.Date(2020, 6, 21),
    start_time::Dates.Time=Dates.Time(12),
    duration_seconds::Float64=1.0,
    ri_par_f::Float64=100.0,
    ri_nir_f::Float64=50.0,
    direct_fraction::Float64=1.0,
    sun_azimut::Float64=180.0,
    sun_elevation::Float64=90.0,
)
    start_dt = Dates.DateTime(date, start_time)
    end_dt = start_dt + Dates.Millisecond(round(Int, duration_seconds * 1000))
    return (
        date=date,
        hour_start=Dates.format(Dates.Time(start_dt), Dates.DateFormat("HH:MM:SS")),
        hour_end=Dates.format(Dates.Time(end_dt), Dates.DateFormat("HH:MM:SS")),
        step_duration=duration_seconds,
        latitude=0.0,
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
    )
end

function _flag_required(name::AbstractString)
    return lowercase(get(ENV, name, "")) in ("required", "force", "error")
end

function _optional_metal_backend()
    if !_BENCH_WANTS_METAL
        return nothing
    end
    if !Sys.isapple() || Sys.ARCH != :aarch64
        msg = "Metal benchmark requires Apple Silicon macOS."
        _flag_required("ARCHIMEDLIGHT_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end
    if Base.find_package("Metal") === nothing
        msg = "Metal.jl is not available. Run this with `--project=test/gpu` after instantiating that environment."
        _flag_required("ARCHIMEDLIGHT_BENCH_METAL") ? error(msg) : (@warn msg; return nothing)
    end

    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _optional_array_backend(package::Symbol, array_type_name::Symbol, env_name::AbstractString, wants_backend::Bool)
    wants_backend || return nothing
    if Base.find_package(string(package)) === nothing
        msg = "$(package).jl is not available. Run this with an environment that includes $(package).jl."
        _flag_required(env_name) ? error(msg) : (@warn msg; return nothing)
    end

    try
        mod = Base.require(Main, package)
        array_type = getproperty(mod, array_type_name)
        dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
        return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
    catch err
        msg = "Could not initialize $(package).jl benchmark backend: $(sprint(showerror, err))"
        _flag_required(env_name) ? error(msg) : (@warn msg; return nothing)
    end
end

_optional_cuda_backend() = _optional_array_backend(:CUDA, :CuArray, "ARCHIMEDLIGHT_BENCH_CUDA", _BENCH_WANTS_CUDA)
_optional_oneapi_backend() = _optional_array_backend(:oneAPI, :oneArray, "ARCHIMEDLIGHT_BENCH_ONEAPI", _BENCH_WANTS_ONEAPI)
_optional_amdgpu_backend() = _optional_array_backend(:AMDGPU, :ROCArray, "ARCHIMEDLIGHT_BENCH_AMDGPU", _BENCH_WANTS_AMDGPU)

function _raycore_interception_backend(backend; max_hits::Int, workgroupsize, edge_accumulation::Symbol, max_prechunk_instances)
    kwargs = (
        backend=backend,
        max_hits_per_pixel=max_hits,
        edge_accumulation=edge_accumulation,
    )
    max_prechunk_instances !== nothing && (kwargs = merge(kwargs, (; max_prechunk_instances=max_prechunk_instances)))
    if workgroupsize === nothing
        return ArchimedLight.RaycoreInterceptionBackend(; kwargs...)
    end
    kwargs = merge(kwargs, (; workgroupsize=workgroupsize))
    return ArchimedLight.RaycoreInterceptionBackend(; kwargs...)
end

function _backend_cases(;
    max_hits::Int=BENCH_MAX_HITS,
    workgroupsize=BENCH_WORKGROUPSIZE,
    edge_accumulation::Symbol=BENCH_EDGE_ACCUMULATION,
    max_prechunk_instances=BENCH_MAX_PRECHUNK_INSTANCES,
)
    cases = [
        (
            name="normal_cpu",
            backend=nothing,
            interception_backend=nothing,
            scattering_backend=nothing,
        ),
        (
            name="raycore_ka_cpu",
            backend=KernelAbstractions.CPU(),
            interception_backend=_raycore_interception_backend(
                KernelAbstractions.CPU();
                max_hits=max_hits,
                workgroupsize=workgroupsize,
                edge_accumulation=edge_accumulation,
                max_prechunk_instances=max_prechunk_instances,
            ),
            scattering_backend=nothing,
        ),
    ]
    cases[2] = merge(cases[2], (; scattering_backend=ArchimedLight.RaycoreScatteringBackend(cases[2].interception_backend; edge_accumulation=edge_accumulation)))

    metal_backend = _optional_metal_backend()
    if metal_backend !== nothing
        ib = _raycore_interception_backend(
            metal_backend;
            max_hits=max_hits,
            workgroupsize=workgroupsize,
            edge_accumulation=edge_accumulation,
            max_prechunk_instances=max_prechunk_instances,
        )
        push!(
            cases,
            (
                name="raycore_metal_gpu",
                backend=metal_backend,
                interception_backend=ib,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=edge_accumulation),
            ),
        )
    end
    cuda_backend = _optional_cuda_backend()
    if cuda_backend !== nothing
        ib = _raycore_interception_backend(
            cuda_backend;
            max_hits=max_hits,
            workgroupsize=workgroupsize,
            edge_accumulation=edge_accumulation,
            max_prechunk_instances=max_prechunk_instances,
        )
        push!(
            cases,
            (
                name="raycore_cuda_gpu",
                backend=cuda_backend,
                interception_backend=ib,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=edge_accumulation),
            ),
        )
    end
    oneapi_backend = _optional_oneapi_backend()
    if oneapi_backend !== nothing
        ib = _raycore_interception_backend(
            oneapi_backend;
            max_hits=max_hits,
            workgroupsize=workgroupsize,
            edge_accumulation=edge_accumulation,
            max_prechunk_instances=max_prechunk_instances,
        )
        push!(
            cases,
            (
                name="raycore_oneapi_gpu",
                backend=oneapi_backend,
                interception_backend=ib,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=edge_accumulation),
            ),
        )
    end
    amdgpu_backend = _optional_amdgpu_backend()
    if amdgpu_backend !== nothing
        ib = _raycore_interception_backend(
            amdgpu_backend;
            max_hits=max_hits,
            workgroupsize=workgroupsize,
            edge_accumulation=edge_accumulation,
            max_prechunk_instances=max_prechunk_instances,
        )
        push!(
            cases,
            (
                name="raycore_amdgpu_gpu",
                backend=amdgpu_backend,
                interception_backend=ib,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=edge_accumulation),
            ),
        )
    end
    selected = Set(backend for backend in strip.(BENCH_BACKENDS) if !isempty(backend))
    if "all" in selected
        return cases
    end
    filtered = [case for case in cases if lowercase(case.name) in selected]
    isempty(filtered) && error(
        "ARCHIMEDLIGHT_BENCH_BACKENDS selected no backends. " *
        "Use all, normal_cpu, raycore_ka_cpu, raycore_metal_gpu, raycore_cuda_gpu, raycore_oneapi_gpu, or raycore_amdgpu_gpu. " *
        "GPU backends also require the matching ARCHIMEDLIGHT_BENCH_METAL, ARCHIMEDLIGHT_BENCH_CUDA, " *
        "ARCHIMEDLIGHT_BENCH_ONEAPI, or ARCHIMEDLIGHT_BENCH_AMDGPU flag.",
    )
    return filtered
end

function _with_scattering(options::ArchimedLight.LightOptions, scattering::Bool)
    return ArchimedLight.LightOptions(options; scattering=scattering)
end

function _benchmark_fixture(; sectors::Int=BENCH_SECTORS, pixel_size::Float64=BENCH_PIXEL_SIZE, plot_paving::Int=BENCH_PLOT_PAVING)
    if BENCH_WORKLOAD == "smoke"
        scene = _smoke_scene()
        models = _smoke_models()
        meteo_row = _benchmark_meteo_row(; ri_par_f=100.0, ri_nir_f=50.0)
        first_order_options = ArchimedLight.LightOptions(
            turtle_sectors=sectors,
            all_in_turtle=false,
            scattering=false,
            pixel_size=pixel_size,
            toricity=false,
        )
        scattering_options = ArchimedLight.LightOptions(
            turtle_sectors=sectors,
            all_in_turtle=false,
            scattering=true,
            pixel_size=pixel_size,
            toricity=false,
        )
        return (
            name="smoke_stacked_plates",
            scene=scene,
            models=models,
            meteo_row=meteo_row,
            first_order_options=first_order_options,
            scattering_options=scattering_options,
        )
    elseif BENCH_WORKLOAD in ("coffee", "coffee_home", "home")
        repo_root = dirname(@__DIR__)
        config_path = joinpath(repo_root, "example_2", "config.yml")
        options, scene, meteo, models = ArchimedLight.read_config(config_path; plot_paving_override=plot_paving)
        if BENCH_HAS_SECTOR_OVERRIDE || BENCH_HAS_PIXEL_SIZE_OVERRIDE
            options = ArchimedLight.LightOptions(
                options;
                turtle_sectors=BENCH_HAS_SECTOR_OVERRIDE ? sectors : options.turtle_sectors,
                pixel_size=BENCH_HAS_PIXEL_SIZE_OVERRIDE ? pixel_size : options.pixel_size,
            )
        end
        meteo_row = first(ArchimedLight.prepare_meteo(meteo, options).rows)
        return (
            name="coffee_home_figure",
            scene=scene,
            models=models,
            meteo_row=meteo_row,
            first_order_options=_with_scattering(options, false),
            scattering_options=_with_scattering(options, true),
        )
    end

    error("Unsupported ARCHIMEDLIGHT_BENCH_WORKLOAD=$BENCH_WORKLOAD. Use smoke or coffee_home.")
end

function _measure(label::AbstractString, f; warmups::Int=BENCH_WARMUPS, samples::Int=BENCH_SAMPLES)
    for _ in 1:warmups
        f()
    end

    times_ms = Float64[]
    allocations = Int[]
    for _ in 1:samples
        GC.gc()
        t0 = time_ns()
        bytes = @allocated f()
        push!(times_ms, (time_ns() - t0) / 1e6)
        push!(allocations, bytes)
    end

    return (
        label=label,
        median_ms=median(times_ms),
        min_ms=minimum(times_ms),
        max_ms=maximum(times_ms),
        median_alloc_mb=median(allocations) / 1024^2,
    )
end

function _first_order_call(scene, models, meteo_row, options, case)
    if case.name == "normal_cpu"
        return ArchimedLight.run_light_step(scene, models, meteo_row, options)
    else
        return ArchimedLight.run_light_step(
            scene,
            models,
            meteo_row,
            options;
            interception_backend=case.interception_backend,
        )
    end
end

function _scattering_call(scene, models, meteo_row, options, case)
    if case.name == "normal_cpu"
        return ArchimedLight.run_light_step(scene, models, meteo_row, options; scattering_mode=:raycast)
    else
        return ArchimedLight.run_light_step(
            scene,
            models,
            meteo_row,
            options;
            interception_backend=case.interception_backend,
            scattering_backend=case.scattering_backend,
        )
    end
end

function _cached_call(scene, models, meteo_row, options, case; scattering_mode::Symbol=:raycast)
    sim =
        if case.name == "normal_cpu"
            ArchimedLight.LightSimulation(
                scene,
                models;
                options=options,
                interception_backend=:raster_cpu,
                scattering_mode=scattering_mode,
            )
        else
            ArchimedLight.LightSimulation(
                scene,
                models;
                options=options,
                interception_backend=case.interception_backend,
                scattering_backend=case.scattering_backend,
            )
        end
    cache = ArchimedLight._ensure_light_cache!(sim)
    backend_info = (
        resolved_backend=string(nameof(typeof(cache.resolved_interception_backend))),
        fallback_reason=cache.fallback_reason,
    )
    return () -> ArchimedLight.run_light(sim, meteo_row), backend_info
end

function _backend_result_label(case, backend_info)
    backend_info === nothing && return case.name
    resolved_backend = backend_info.resolved_backend
    expected = case.name == "normal_cpu" ? "RasterCPUBackend" : "RaycoreInterceptionBackend"
    resolved_backend == expected && return case.name
    return "$(case.name)->$(resolved_backend)"
end

_fallback_result_label(backend_info) =
    backend_info === nothing ? "" : string(backend_info.fallback_reason)

function _first_order_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return () -> _first_order_call(scene, models, meteo_row, options, case), nothing
    elseif BENCH_EXECUTION == "cached"
        return _cached_call(scene, models, meteo_row, options, case)
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _scattering_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return () -> _scattering_call(scene, models, meteo_row, options, case), nothing
    elseif BENCH_EXECUTION == "cached"
        return _cached_call(scene, models, meteo_row, options, case; scattering_mode=:raycast)
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _sky_turtle_fluxes(meteo_row, options)
    sky = ArchimedLight.compute_sky(meteo_row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(meteo_row, sky, turtle, options)
    return sky, turtle, fluxes
end

function _component_cache(scene, models, options, case)
    if case.name == "normal_cpu"
        return ArchimedLight.prepare_light_cache(
            scene,
            models,
            options;
            interception_backend=:raster_cpu,
        )
    end
    return ArchimedLight.prepare_light_cache(
        scene,
        models,
        options;
        interception_backend=case.interception_backend,
        scattering_backend=case.scattering_backend,
    )
end

function _prepared_topology_state(scene, models, meteo_row, options, case)
    _sky, turtle, fluxes = _sky_turtle_fluxes(meteo_row, options)
    cache = _component_cache(scene, models, options, case)
    backend =
        cache.raycore_data !== nothing && cache.scattering_backend isa ArchimedLight.RaycoreScatteringBackend ?
        cache.scattering_backend :
        nothing
    return (
        turtle=turtle,
        fluxes=fluxes,
        prepared=cache.prepared,
        raycore_data=cache.raycore_data,
        backend=backend,
        resolved_interception_backend=cache.resolved_interception_backend,
        fallback_reason=cache.fallback_reason,
    )
end

function _build_first_topology(scene, models, options, state, case)
    if state.raycore_data === nothing
        return ArchimedLight._stream_first_order_with_scattering_topology(
            state.prepared,
            scene,
            models,
            state.turtle,
            state.fluxes,
            options,
        )
    end
    return ArchimedLight._stream_first_order_with_raycore_scattering_topology(
        state.raycore_data,
        scene,
        models,
        state.turtle,
        state.fluxes,
        options,
    )
end

function _topology_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return function ()
            state = _prepared_topology_state(scene, models, meteo_row, options, case)
            return _build_first_topology(scene, models, options, state, case)
        end, nothing
    elseif BENCH_EXECUTION == "cached"
        state = _prepared_topology_state(scene, models, meteo_row, options, case)
        backend_info = (
            resolved_backend=string(nameof(typeof(state.resolved_interception_backend))),
            fallback_reason=state.fallback_reason,
        )
        return () -> _build_first_topology(scene, models, options, state, case), backend_info
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _propagation_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return function ()
            state = _prepared_topology_state(scene, models, meteo_row, options, case)
            first, topology = _build_first_topology(scene, models, options, state, case)
            backend = case.name == "normal_cpu" ? nothing : state.backend
            return ArchimedLight.compute_scattering(topology, first, options; backend=backend)
        end, nothing
    elseif BENCH_EXECUTION == "cached"
        state = _prepared_topology_state(scene, models, meteo_row, options, case)
        first, topology = _build_first_topology(scene, models, options, state, case)
        backend = case.name == "normal_cpu" ? nothing : state.backend
        backend_info = (
            resolved_backend=string(nameof(typeof(state.resolved_interception_backend))),
            fallback_reason=state.fallback_reason,
        )
        return () -> ArchimedLight.compute_scattering(topology, first, options; backend=backend), backend_info
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _integration_benchmark_call(scene, models, meteo_row, options, case)
    integrate_from_state(first, scat, state) = ArchimedLight.integrate_light(
        scene,
        models,
        first,
        scat,
        options;
        meteo_row=meteo_row,
        component_area_per_node=state.prepared.component_area_per_node,
        absorption_par_per_node=state.prepared.absorption_par_per_node,
        absorption_nir_per_node=state.prepared.absorption_nir_per_node,
    )
    if BENCH_EXECUTION == "oneshot"
        return function ()
            state = _prepared_topology_state(scene, models, meteo_row, options, case)
            first, topology = _build_first_topology(scene, models, options, state, case)
            backend = case.name == "normal_cpu" ? nothing : state.backend
            scat = ArchimedLight.compute_scattering(topology, first, options; backend=backend)
            return integrate_from_state(first, scat, state)
        end, nothing
    elseif BENCH_EXECUTION == "cached"
        state = _prepared_topology_state(scene, models, meteo_row, options, case)
        first, topology = _build_first_topology(scene, models, options, state, case)
        backend = case.name == "normal_cpu" ? nothing : state.backend
        scat = ArchimedLight.compute_scattering(topology, first, options; backend=backend)
        backend_info = (
            resolved_backend=string(nameof(typeof(state.resolved_interception_backend))),
            fallback_reason=state.fallback_reason,
        )
        return () -> integrate_from_state(first, scat, state), backend_info
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _stack_trace_state(scene, models, meteo_row, options, case)
    case.name == "normal_cpu" && return nothing
    _sky, turtle, _fluxes = _sky_turtle_fluxes(meteo_row, options)
    cache = _component_cache(scene, models, options, case)
    directions = [sector.direction for sector in turtle.sectors if Float32(sector.direction[3]) > 0.0f0]
    return (
        raycore_data=cache.raycore_data,
        directions=directions,
        resolved_interception_backend=cache.resolved_interception_backend,
        fallback_reason=cache.fallback_reason,
    )
end

function _stack_trace_summary(state, options)
    state === nothing && return ""
    state.raycore_data === nothing && return ",resolved=$(nameof(typeof(state.resolved_interception_backend))),fallback=$(state.fallback_reason),stack_trace=skipped"
    total_pixels = 0
    occupied_pixels = 0
    total_hits = 0
    max_seen = 0
    overflow = false
    max_hits = state.raycore_data.max_hits_per_pixel
    for direction in state.directions
        traced = ArchimedLight._raycore_trace_direction_stack_nodes(
            state.raycore_data,
            direction,
            options,
        )
        total_pixels += length(traced.counts)
        overflow |= any(traced.overflow)
        @inbounds for count in traced.counts
            c = Int(count)
            occupied_pixels += c > 0
            total_hits += c
            max_seen = max(max_seen, c)
        end
    end
    capacity = total_pixels * max_hits
    hit_util = capacity == 0 ? 0.0 : 100 * total_hits / capacity
    occupied = total_pixels == 0 ? 0.0 : 100 * occupied_pixels / total_pixels
    return @sprintf(
        ",hit_util=%.2f%%,max_seen=%d,occupied=%.2f%%,overflow=%s",
        hit_util,
        max_seen,
        occupied,
        string(overflow),
    )
end

function _stack_trace_benchmark_call(scene, models, meteo_row, options, case)
    case.name == "normal_cpu" && return nothing, ""
    if BENCH_EXECUTION == "oneshot"
        state = _stack_trace_state(scene, models, meteo_row, options, case)
        summary = _stack_trace_summary(state, options)
        state.raycore_data === nothing && return nothing, summary
        return function ()
            local_state = _stack_trace_state(scene, models, meteo_row, options, case)
            local_state.raycore_data === nothing && return nothing
            for direction in local_state.directions
                ArchimedLight._raycore_trace_direction_stack_nodes_device(
                    local_state.raycore_data,
                    direction,
                    options,
                )
            end
            return nothing
        end, summary
    elseif BENCH_EXECUTION == "cached"
        state = _stack_trace_state(scene, models, meteo_row, options, case)
        summary = _stack_trace_summary(state, options)
        state.raycore_data === nothing && return nothing, summary
        return function ()
            for direction in state.directions
                ArchimedLight._raycore_trace_direction_stack_nodes_device(
                    state.raycore_data,
                    direction,
                    options,
                )
            end
            return nothing
        end, summary
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _print_markdown_table(results)
    println("| workload | backend | fallback | median ms | min ms | max ms | median alloc MiB |")
    println("|---|---:|---:|---:|---:|---:|---:|")
    for row in results
        fallback = hasproperty(row, :fallback_reason) ? string(row.fallback_reason) : ""
        @printf(
            "| %s | %s | %s | %.3f | %.3f | %.3f | %.3f |\n",
            row.workload,
            row.backend,
            fallback,
            row.median_ms,
            row.min_ms,
            row.max_ms,
            row.median_alloc_mb,
        )
    end
end

function _variant_label(
    component::AbstractString,
    fixture,
    first_order_options,
    plot_paving::Int,
    max_hits::Int,
    max_prechunk_instances,
    workgroupsize,
    edge_accumulation::Symbol;
    extra::AbstractString="",
)
    return string(
        component,
        "[scene=", fixture.name,
        ",sectors=", first_order_options.turtle_sectors,
        ",pixel=", first_order_options.pixel_size,
        ",paving=", BENCH_WORKLOAD == "smoke" ? "n/a" : string(plot_paving),
        ",max_hits=", max_hits,
        ",max_prechunk=", _max_prechunk_instances_label(max_prechunk_instances),
        ",wg=", _workgroupsize_label(workgroupsize),
        ",edge=", edge_accumulation,
        extra,
        "]",
    )
end

function _configured_or_values(values, has_override::Bool)
    return (BENCH_WORKLOAD == "smoke" || has_override) ? join(values, ',') : "config"
end

function _append_variant_results!(
    results,
    components,
    fixture,
    cases,
    plot_paving::Int,
    max_hits::Int,
    max_prechunk_instances,
    workgroupsize,
    edge_accumulation::Symbol,
)
    scene = fixture.scene
    models = fixture.models
    meteo_row = fixture.meteo_row
    first_order_options = fixture.first_order_options
    scattering_options = fixture.scattering_options
    for case in cases
        if "first_order" in components
            call, backend_info =
                _first_order_benchmark_call(scene, models, meteo_row, first_order_options, case)
            push!(
                results,
                merge(
                    _measure("first_order", call),
                    (;
                        workload=_variant_label("first_order", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation),
                        backend=_backend_result_label(case, backend_info),
                        fallback_reason=_fallback_result_label(backend_info),
                    ),
                ),
            )
        end
        if "scattering" in components
            call, backend_info =
                _scattering_benchmark_call(scene, models, meteo_row, scattering_options, case)
            push!(
                results,
                merge(
                    _measure("scattering", call),
                    (;
                        workload=_variant_label("scattering", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation),
                        backend=_backend_result_label(case, backend_info),
                        fallback_reason=_fallback_result_label(backend_info),
                    ),
                ),
            )
        end
        if "topology" in components
            call, backend_info =
                _topology_benchmark_call(scene, models, meteo_row, scattering_options, case)
            push!(
                results,
                merge(
                    _measure("topology", call),
                    (;
                        workload=_variant_label("topology", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation),
                        backend=_backend_result_label(case, backend_info),
                        fallback_reason=_fallback_result_label(backend_info),
                    ),
                ),
            )
        end
        if "propagation" in components
            call, backend_info =
                _propagation_benchmark_call(scene, models, meteo_row, scattering_options, case)
            push!(
                results,
                merge(
                    _measure("propagation", call),
                    (;
                        workload=_variant_label("propagation", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation),
                        backend=_backend_result_label(case, backend_info),
                        fallback_reason=_fallback_result_label(backend_info),
                    ),
                ),
            )
        end
        if "integration" in components
            call, backend_info =
                _integration_benchmark_call(scene, models, meteo_row, scattering_options, case)
            push!(
                results,
                merge(
                    _measure("integration", call),
                    (;
                        workload=_variant_label("integration", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation),
                        backend=_backend_result_label(case, backend_info),
                        fallback_reason=_fallback_result_label(backend_info),
                    ),
                ),
            )
        end
        if "stack_trace" in components
            trace_call, trace_summary = _stack_trace_benchmark_call(scene, models, meteo_row, scattering_options, case)
            if trace_call !== nothing
                push!(
                    results,
                    merge(
                        _measure("stack_trace", trace_call),
                        (;
                            workload=_variant_label("stack_trace", fixture, first_order_options, plot_paving, max_hits, max_prechunk_instances, workgroupsize, edge_accumulation; extra=trace_summary),
                            backend=case.name,
                            fallback_reason="",
                        ),
                    ),
                )
            end
        end
    end
    return results
end

function _run_raycore_smokes(cases)
    for case in cases
        if case.name != "normal_cpu"
            run_raycore_device_smoke(;
                backend=case.backend,
                first_order_area_atol=case.backend isa KernelAbstractions.CPU ? 1e-10 : 1e-5,
                first_order_area_rtol=case.backend isa KernelAbstractions.CPU ? 1e-10 : 1e-5,
                first_order_power_atol=case.backend isa KernelAbstractions.CPU ? 1e-8 : 1e-4,
                first_order_power_rtol=case.backend isa KernelAbstractions.CPU ? 1e-8 : 1e-5,
                scattering_atol=case.backend isa KernelAbstractions.CPU ? 1e-10 : 1e-4,
                scattering_rtol=case.backend isa KernelAbstractions.CPU ? 1e-10 : 1e-4,
            )
        end
    end
end

function main()
    components = _bench_components()
    results = NamedTuple[]

    for plot_paving in BENCH_PLOT_PAVING_VALUES,
        sectors in BENCH_SECTOR_VALUES,
        pixel_size in BENCH_PIXEL_SIZE_VALUES,
        max_hits in BENCH_MAX_HIT_VALUES,
        max_prechunk_instances in BENCH_MAX_PRECHUNK_INSTANCE_VALUES,
        workgroupsize in BENCH_WORKGROUPSIZE_VALUES,
        edge_accumulation in BENCH_EDGE_ACCUMULATION_VALUES

        fixture = _benchmark_fixture(; sectors=sectors, pixel_size=pixel_size, plot_paving=plot_paving)
        cases = _backend_cases(;
            max_hits=max_hits,
            workgroupsize=workgroupsize,
            edge_accumulation=edge_accumulation,
            max_prechunk_instances=max_prechunk_instances,
        )
        _run_raycore_smokes(cases)
        _append_variant_results!(
            results,
            components,
            fixture,
            cases,
            plot_paving,
            max_hits,
            max_prechunk_instances,
            workgroupsize,
            edge_accumulation,
        )
    end

    println()
    println("ArchimedLight backend comparison")
    println(
        "workload=$(BENCH_WORKLOAD), warmups=$(BENCH_WARMUPS), samples=$(BENCH_SAMPLES), " *
        "sector_values=$(_configured_or_values(BENCH_SECTOR_VALUES, BENCH_HAS_SECTOR_OVERRIDE)), " *
        "pixel_size_values=$(_configured_or_values(BENCH_PIXEL_SIZE_VALUES, BENCH_HAS_PIXEL_SIZE_OVERRIDE)), " *
        "plot_paving_values=$(BENCH_WORKLOAD == "smoke" ? "n/a" : join(BENCH_PLOT_PAVING_VALUES, ',')), " *
        "max_hit_values=$(join(BENCH_MAX_HIT_VALUES, ',')), " *
        "max_prechunk_instance_values=$(_max_prechunk_instances_values_label(BENCH_MAX_PRECHUNK_INSTANCE_VALUES)), " *
        "workgroupsize_values=$(_workgroupsize_values_label(BENCH_WORKGROUPSIZE_VALUES)), " *
        "edge_accumulation_values=$(join(string.(BENCH_EDGE_ACCUMULATION_VALUES), ',')), " *
        "execution=$(BENCH_EXECUTION), " *
        "components=$(join(sort!(collect(components)), ','))",
    )
    println()
    _print_markdown_table(results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
