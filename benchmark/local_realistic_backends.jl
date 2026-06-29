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

using ArchimedLight
using GeometryBasics
using KernelAbstractions
using MultiScaleTreeGraph
using PlantGeom
using Printf
using Raycore
using Statistics

const BENCH_BUILD_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_BUILD_SAMPLES", "2"))
const BENCH_WARM_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_SAMPLES", "3"))
const BENCH_BACKENDS = split(lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_BACKENDS", "all")), ',')
const BENCH_WORKLOADS_RAW = get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS", "coffee,agripv")
const BENCH_WORKLOADS_EXPLICIT = haskey(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS")
const BENCH_WORKLOADS = [strip(lowercase(value)) for value in split(BENCH_WORKLOADS_RAW, ',') if !isempty(strip(value))]
const BENCH_BREAKDOWN = lowercase(get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN", "")) in ("1", "true", "yes", "on")
const BENCH_PREPARE_BREAKDOWN = lowercase(get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_PREPARE_BREAKDOWN", "")) in ("1", "true", "yes", "on")
const BENCH_MAX_PRECHUNK_INSTANCES = haskey(ENV, "ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES") ?
                                     parse(Int, ENV["ARCHIMEDLIGHT_BENCH_MAX_PRECHUNK_INSTANCES"]) :
                                     nothing

function _env_int(name::AbstractString, default::Int)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    return parse(Int, raw)
end

function _env_float(name::AbstractString, default::Float64)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    return parse(Float64, raw)
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

function _raycore_case(name, backend)
    kwargs = (
        backend=backend,
        max_hits_per_pixel=_env_int("ARCHIMEDLIGHT_BENCH_MAX_HITS", 32),
        workgroupsize=_env_int("ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", 256),
    )
    BENCH_MAX_PRECHUNK_INSTANCES !== nothing &&
        (kwargs = merge(kwargs, (; max_prechunk_instances=BENCH_MAX_PRECHUNK_INSTANCES)))
    interception = ArchimedLight.RaycoreInterceptionBackend(; kwargs...)
    return (
        name=name,
        interception_backend=interception,
        scattering_backend=ArchimedLight.RaycoreScatteringBackend(interception),
    )
end

function _cases()
    cases = [
        (name="normal_cpu", interception_backend=:raster_cpu, scattering_backend=nothing),
        _raycore_case("raycore_ka_cpu", KernelAbstractions.CPU()),
    ]
    metal = _optional_metal_backend()
    metal !== nothing && push!(cases, _raycore_case("raycore_metal_gpu", metal))
    cuda = _optional_cuda_backend()
    cuda !== nothing && push!(cases, _raycore_case("raycore_cuda_gpu", cuda))
    oneapi = _optional_oneapi_backend()
    oneapi !== nothing && push!(cases, _raycore_case("raycore_oneapi_gpu", oneapi))
    amdgpu = _optional_amdgpu_backend()
    amdgpu !== nothing && push!(cases, _raycore_case("raycore_amdgpu_gpu", amdgpu))

    selected = Set(backend for backend in strip.(BENCH_BACKENDS) if !isempty(backend))
    "all" in selected && return cases
    filtered = [case for case in cases if lowercase(case.name) in selected]
    isempty(filtered) && error(
        "ARCHIMEDLIGHT_BENCH_BACKENDS selected no backends. " *
        "Use all, normal_cpu, raycore_ka_cpu, raycore_metal_gpu, raycore_cuda_gpu, raycore_oneapi_gpu, or raycore_amdgpu_gpu. " *
        "GPU backends also require the matching ARCHIMEDLIGHT_BENCH_METAL, ARCHIMEDLIGHT_BENCH_CUDA, " *
        "ARCHIMEDLIGHT_BENCH_ONEAPI, or ARCHIMEDLIGHT_BENCH_AMDGPU flag.",
    )
    return filtered
end

function _coffee_workload()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    paving = _env_int("ARCHIMEDLIGHT_BENCH_COFFEE_PAVING", 2500)
    options, scene, meteo, models = ArchimedLight.read_config(config; plot_paving_override=paving)
    options = ArchimedLight.LightOptions(options; cache_radiation=true, meteo_range=nothing, scattering=true)
    return (
        name="coffee_multi_step",
        scene=scene,
        models=models,
        input=meteo,
        options=options,
        runner=(sim, input) -> ArchimedLight.run_light(sim, input),
        steps=length(ArchimedLight.prepare_meteo(meteo, options).rows),
    )
end

function _vpalm_modules()
    xpalm_root = get(ENV, "ARCHIMEDLIGHT_BENCH_XPALM_ROOT", "/Users/rvezy/Documents/dev/XPalm")
    isdir(xpalm_root) || error("XPalm checkout not found at $xpalm_root")
    package_dir = dirname(xpalm_root)

    old_load_path = copy(LOAD_PATH)
    try
        filter!(path -> path != xpalm_root && path != package_dir, LOAD_PATH)
        pushfirst!(LOAD_PATH, package_dir)
        pushfirst!(LOAD_PATH, xpalm_root)

        xpalm_mod = try
            Base.require(Main, :XPalm)
        catch err
            error("XPalm could not be loaded from $package_dir: $(sprint(showerror, err))")
        end

        vpalm_mod = try
            if !Base.invokelatest(isdefined, xpalm_mod, :VPalm)
                Base.invokelatest(Base.include, xpalm_mod, joinpath(xpalm_root, "src", "VPalm.jl"))
            end
            Base.invokelatest(getfield, xpalm_mod, :VPalm)
        catch err
            error("XPalm.VPalm could not be loaded: $(sprint(showerror, err))")
        end

        return (XPalm=xpalm_mod, VPalm=vpalm_mod, root=xpalm_root)
    finally
        empty!(LOAD_PATH)
        append!(LOAD_PATH, old_load_path)
    end
end

function _vpalm_merge_scale()
    raw = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_XPALM_MERGE_SCALE", "leaf"))
    raw == "none" && return :none
    raw == "leaflet" && return :leaflet
    raw == "leaf" && return :leaf
    raw == "plant" && return :plant
    error("Unsupported ARCHIMEDLIGHT_BENCH_XPALM_MERGE_SCALE=$raw. Use none, leaflet, leaf, or plant.")
end

function _xpalm_workload()
    modules = _vpalm_modules()
    VPalm = modules.VPalm
    parameter_file = joinpath(modules.root, "test", "references", "vpalm-parameter_file.yml")
    isfile(parameter_file) || error("VPalm parameter file not found at $parameter_file")

    read_parameters = Base.invokelatest(getproperty, VPalm, :read_parameters)
    parameters = Base.invokelatest(read_parameters, parameter_file)
    merge_scale = _vpalm_merge_scale()
    build_mockup = Base.invokelatest(getproperty, VPalm, :build_mockup)
    palm = Base.invokelatest(build_mockup, parameters; merge_scale=merge_scale, rng=nothing)

    spacing = _env_float("ARCHIMEDLIGHT_BENCH_XPALM_SPACING", 6.0)
    ground_n = _env_int("ARCHIMEDLIGHT_BENCH_XPALM_GROUND_N", 12)
    scene = PlantGeom.make_scene(domain=(-4.0, -4.0, spacing + 4.0, 4.0)) do s
        PlantGeom.add_plant!(s, palm; group="oil_palm", id=1, at=(0.0, 0.0, 0.0))
        PlantGeom.add_plant!(s, palm; group="oil_palm", id=2, at=(spacing, 0.0, 0.0), rotate=(z=180.0,), deg=true)
        PlantGeom.add_ground!(s; nx=ground_n, ny=ground_n, group="pavement", type="Cobblestone")
    end

    models = ArchimedLight.models_for(
        "oil_palm" => ("*" => ArchimedLight.translucent(par=0.15, nir=0.90),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
    options = ArchimedLight.LightOptions(
        turtle_sectors=_env_int("ARCHIMEDLIGHT_BENCH_XPALM_SECTORS", 8),
        pixel_size=_env_float("ARCHIMEDLIGHT_BENCH_XPALM_PIXEL_SIZE", 0.05),
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
        toricity=true,
        scattering_max_iter=_env_int("ARCHIMEDLIGHT_BENCH_XPALM_SCATTERING_MAX_ITER", 12),
    )
    skies = [
        ArchimedLight.SkyState(120.0, 35.0, 280.0, 180.0, 0.65, 0.35),
        ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
    ]
    dt = _env_float("ARCHIMEDLIGHT_BENCH_XPALM_STEP_SECONDS", 1800.0)

    return (
        name="xpalm_two_oil_palms",
        scene=scene,
        models=models,
        input=skies,
        options=options,
        runner=(sim, input) -> [ArchimedLight.run_light(sim, sky; step_duration_seconds=dt) for sky in input],
        steps=length(skies),
    )
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

function _agripv_workload()
    root = get(
        ENV,
        "ARCHIMEDLIGHT_BENCH_AGRIPV_ROOT",
        "/Users/rvezy/Documents/cirad/Articles_Rapports_Communications/Conférences/2026_fspm/coutellier_agripv",
    )
    opf_path = joinpath(root, "0_simulations", "archicrop", "wheat", "plant_1995-06-24.opf")
    isfile(opf_path) || error("Agrivoltaics wheat OPF not found at $opf_path")

    n_rows = _env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_ROWS", 2)
    plants_per_row = _env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_PLANTS_PER_ROW", 8)
    interrow = _env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_INTERROW", 0.25)
    intrarow = _env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_INTRAROW", 0.30)
    width = max(interrow * n_rows, interrow)
    scene_length = max(intrarow * plants_per_row, intrarow)

    wheat = PlantGeom.read_opf(opf_path, mtg_type=MultiScaleTreeGraph.NodeMTG)
    panel = _panel_mesh(
        width,
        _env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_PANEL_LENGTH", 3.0),
        _env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_PANEL_HEIGHT", 1.6),
        _env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_PANEL_INCLINATION", 25.0),
    )
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
        PlantGeom.add_ground!(
            s;
            nx=_env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_GROUND_NX", 18),
            ny=_env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_GROUND_NY", 24),
            group="pavement",
            type="Cobblestone",
        )
    end

    models = ArchimedLight.models_for(
        "wheat" => ("*" => ArchimedLight.translucent(par=0.15, nir=0.90),),
        "panel" => ("Panel" => ArchimedLight.translucent(par=0.0, nir=0.0),),
        "pavement" => ("Cobblestone" => ArchimedLight.translucent(par=0.12, nir=0.60),),
    )
    options = ArchimedLight.LightOptions(
        turtle_sectors=_env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_SECTORS", 8),
        pixel_size=_env_float("ARCHIMEDLIGHT_BENCH_AGRIPV_PIXEL_SIZE", 0.05),
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
        toricity=true,
        scattering_max_iter=_env_int("ARCHIMEDLIGHT_BENCH_AGRIPV_SCATTERING_MAX_ITER", 12),
    )
    skies = [
        ArchimedLight.SkyState(120.0, 35.0, 280.0, 180.0, 0.65, 0.35),
        ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
    ]

    return (
        name="agripv_wheat_panel",
        scene=scene,
        models=models,
        input=skies,
        options=options,
        runner=(sim, input) -> [ArchimedLight.run_light(sim, sky; step_duration_seconds=1800.0) for sky in input],
        steps=length(skies),
    )
end

function _selected_workloads()
    isempty(BENCH_WORKLOADS) && error(
        "ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS selected no workloads. Use coffee, agripv, xpalm, or all.",
    )

    selected = Set(BENCH_WORKLOADS)
    unknown = setdiff(selected, Set(["all", "coffee", "coffee_home", "home", "agripv", "agri", "wheat", "xpalm", "vpalm", "oil_palm", "oilpalm"]))
    isempty(unknown) || error(
        "Unsupported ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS value(s): $(join(sort!(collect(unknown)), ',')). " *
        "Use coffee, agripv, xpalm, or all.",
    )

    if "all" in selected
        selected = Set(["coffee", "agripv", "xpalm"])
    end

    workloads = NamedTuple[]
    if !isdisjoint(selected, Set(["coffee", "coffee_home", "home"]))
        push!(workloads, _coffee_workload())
    end

    if !isdisjoint(selected, Set(["agripv", "agri", "wheat"]))
        try
            push!(workloads, _agripv_workload())
        catch err
            if BENCH_WORKLOADS_EXPLICIT
                rethrow()
            end
            @warn "Skipping local agrivoltaics benchmark workload." reason=sprint(showerror, err)
        end
    end

    if !isdisjoint(selected, Set(["xpalm", "vpalm", "oil_palm", "oilpalm"]))
        try
            push!(workloads, _xpalm_workload())
        catch err
            if BENCH_WORKLOADS_EXPLICIT
                rethrow()
            end
            @warn "Skipping local XPalm benchmark workload." reason=sprint(showerror, err)
        end
    end

    isempty(workloads) && error(
        "No local realistic benchmark workloads are available. " *
        "Use ARCHIMEDLIGHT_LOCAL_BENCH_WORKLOADS=coffee to run only bundled fixtures.",
    )
    return workloads
end

function _totals(series)
    par = 0.0
    nir = 0.0
    for step in series
        totals = ArchimedLight._light_step_energy_totals(step)
        par += totals.incident_par_total
        nir += totals.incident_nir_total
    end
    return (incident_par=par, incident_nir=nir)
end

function _time_once(workload, case)
    sim = ArchimedLight.LightSimulation(
        workload.scene,
        workload.models;
        options=workload.options,
        interception_backend=case.interception_backend,
        scattering_backend=case.scattering_backend,
    )
    elapsed = @elapsed series = workload.runner(sim, workload.input)
    resolved_backend =
        sim.cache === nothing ? "unprepared" : string(nameof(typeof(sim.cache.resolved_interception_backend)))
    fallback_reason = sim.cache === nothing ? :unprepared : sim.cache.fallback_reason
    return (seconds=elapsed, series=series, sim=sim, resolved_backend=resolved_backend, fallback_reason=fallback_reason)
end

function _time_warm(workload, sim)
    times = Float64[]
    series = nothing
    for _ in 1:BENCH_WARM_SAMPLES
        elapsed = @elapsed series = workload.runner(sim, workload.input)
        push!(times, elapsed)
    end
    return (
        median=median(times),
        min=minimum(times),
        max=maximum(times),
        samples=copy(times),
        series=series,
    )
end

function _new_simulation(workload, case)
    return ArchimedLight.LightSimulation(
        workload.scene,
        workload.models;
        options=workload.options,
        interception_backend=case.interception_backend,
        scattering_backend=case.scattering_backend,
    )
end

function _time_build(workload, case)
    # One unreported run pays Julia method compilation and GPU kernel compilation.
    _time_once(workload, case)

    times = Float64[]
    series = nothing
    sim = nothing
    for _ in 1:BENCH_BUILD_SAMPLES
        run = _time_once(workload, case)
        push!(times, run.seconds)
        series = run.series
        sim = run.sim
    end
    return (
        median=median(times),
        min=minimum(times),
        max=maximum(times),
        samples=copy(times),
        series=series,
        sim=sim,
        resolved_backend=sim.cache === nothing ? "unprepared" : string(nameof(typeof(sim.cache.resolved_interception_backend))),
        fallback_reason=sim.cache === nothing ? :unprepared : sim.cache.fallback_reason,
        geometry_mode=sim.cache === nothing || sim.cache.raycore_data === nothing ? :none : sim.cache.raycore_data.geometry_mode,
    )
end

function _median_min_max(values::Vector{Float64})
    return (median=median(values), min=minimum(values), max=maximum(values))
end

_mib(bytes::Real) = Float64(bytes) / 1024.0^2

function _time_breakdown_once(workload, case)
    construct_timed = @timed _new_simulation(workload, case)
    sim = construct_timed.value
    prepare_timed = @timed ArchimedLight._ensure_light_cache!(sim)
    populate_timed = @timed workload.runner(sim, workload.input)
    reuse_timed = @timed workload.runner(sim, workload.input)
    summary = ArchimedLight.cache_summary(sim)
    resolved_backend =
        sim.cache === nothing ? "unprepared" : string(nameof(typeof(sim.cache.resolved_interception_backend)))
    fallback_reason = sim.cache === nothing ? :unprepared : sim.cache.fallback_reason
    geometry_mode = sim.cache === nothing || sim.cache.raycore_data === nothing ? :none : sim.cache.raycore_data.geometry_mode
    return (
        construct_seconds=construct_timed.time,
        prepare_seconds=prepare_timed.time,
        populate_seconds=populate_timed.time,
        reuse_seconds=reuse_timed.time,
        construct_bytes=construct_timed.bytes,
        prepare_bytes=prepare_timed.bytes,
        populate_bytes=populate_timed.bytes,
        reuse_bytes=reuse_timed.bytes,
        populate_series=populate_timed.value,
        reuse_series=reuse_timed.value,
        resolved_backend=resolved_backend,
        fallback_reason=fallback_reason,
        geometry_mode=geometry_mode,
        summary=summary,
    )
end

function _phase_seconds(timed)
    timed === nothing && return 0.0
    return timed.time
end

function _phase_mib(timed)
    timed === nothing && return 0.0
    return _mib(timed.bytes)
end

function _time_prepare_breakdown_once(workload, case)
    sim_timed = @timed _new_simulation(workload, case)
    sim = sim_timed.value

    check_timed = @timed ArchimedLight.check_simulation(sim)
    report = check_timed.value
    isempty(report.errors) || error("Invalid benchmark simulation: $(join(report.errors, "; "))")

    resolve_timed = @timed ArchimedLight._resolve_interception_backend(sim.interception_backend)
    ib = resolve_timed.value

    prepared_timed =
        if ib isa ArchimedLight.RasterCPUBackend || ib isa ArchimedLight.RaycoreInterceptionBackend
            @timed ArchimedLight._prepare_interception_data(
                sim.scene,
                sim.models,
                sim.options;
                include_budget_maps=true,
                include_raycore_instancing=ib isa ArchimedLight.RaycoreInterceptionBackend,
            )
        else
            nothing
        end
    prepared = prepared_timed === nothing ? nothing : prepared_timed.value

    render_timed =
        if prepared === nothing
            @timed ArchimedLight._light_render_geometry(
                ArchimedLight._scene_geometry_for_interception(sim.scene, sim.models, sim.options),
            )
        else
            @timed ArchimedLight._light_render_geometry(prepared.geometry)
        end

    prechunk_timed = nothing
    raycore_timed = nothing
    validation_timed = nothing
    retry_timed = nothing
    cache_mode_timed = nothing
    fallback_reason = :none
    resolved_backend = ib
    raycore_data = nothing
    raycore_was_chunked = false
    validation = nothing
    retry_used = false

    if ib isa ArchimedLight.RaycoreInterceptionBackend && prepared !== nothing
        prechunk_timed = @timed ArchimedLight._raycore_prechunk_instance_limit_status(
            prepared,
            ib.config;
            toricity=sim.options.toricity,
        )
        prechunk_status = prechunk_timed.value
        if prechunk_status.exceeded
            fallback_reason = :raycore_prechunk_instance_cap
            resolved_backend = ArchimedLight.RasterCPUBackend()
        else
            raycore_timed = @timed ArchimedLight._raycore_initial_scene_data(
                prepared,
                ib.config,
                sim.options;
                toricity=sim.options.toricity,
            )
            raycore_data, raycore_was_chunked = raycore_timed.value
        end
    end

    if raycore_data !== nothing
        validation_timed = @timed ArchimedLight._raycore_vertical_trace_validation(raycore_data, sim.options)
        validation = validation_timed.value
        if !validation.ok
            retry_timed = @timed begin
                raycore_was_chunked ?
                (nothing, validation) :
                ArchimedLight._raycore_retry_chunked_scene_data(
                    prepared,
                    ib.config,
                    sim.options,
                    validation;
                    toricity=sim.options.toricity,
                )
            end
            chunked_data, chunked_validation = retry_timed.value
            if chunked_data !== nothing
                raycore_data = chunked_data
                validation = chunked_validation
                retry_used = true
            else
                fallback_reason = :raycore_trace_validation
                resolved_backend = ArchimedLight.RasterCPUBackend()
                raycore_data = nothing
                validation = chunked_validation
            end
        end
    end

    limit = sim.memory_limit_bytes === nothing ?
            ArchimedLight._default_light_cache_memory_limit() :
            Int(sim.memory_limit_bytes)
    limit = max(limit, 1)
    cache_mode_timed = @timed ArchimedLight._cache_mode_for(prepared, sim.options, resolved_backend, limit)

    return (
        construct=sim_timed,
        check=check_timed,
        resolve=resolve_timed,
        prepared=prepared_timed,
        render=render_timed,
        prechunk=prechunk_timed,
        raycore=raycore_timed,
        validation=validation_timed,
        retry=retry_timed,
        cache_mode_timed=cache_mode_timed,
        resolved_backend=string(nameof(typeof(resolved_backend))),
        fallback_reason=fallback_reason,
        cache_mode=cache_mode_timed.value,
        raycore_instances=raycore_data === nothing ? 0 : length(raycore_data.tlas.instances),
        raycore_geometries=raycore_data === nothing ? 0 : Raycore.n_geometries(raycore_data.tlas),
        raycore_chunked=raycore_data !== nothing && raycore_data.chunked_tlas,
        geometry_mode=raycore_data === nothing ? :none : raycore_data.geometry_mode,
        validation_ratio=validation === nothing ? 1.0 : validation.ratio,
        validation_reference_pixels=validation === nothing ? 0 : validation.reference_pixels,
        validation_raycore_pixels=validation === nothing ? 0 : validation.raycore_pixels,
        retry_used=retry_used,
    )
end

function _time_prepare_breakdown(workload, case)
    _time_prepare_breakdown_once(workload, case)
    runs = [_time_prepare_breakdown_once(workload, case) for _ in 1:BENCH_BUILD_SAMPLES]
    median_phase(getphase) = (
        seconds=median([_phase_seconds(getphase(run)) for run in runs]),
        mib=median([_phase_mib(getphase(run)) for run in runs]),
    )
    last = runs[end]
    return (
        construct=median_phase(run -> run.construct),
        check=median_phase(run -> run.check),
        resolve=median_phase(run -> run.resolve),
        prepared=median_phase(run -> run.prepared),
        render=median_phase(run -> run.render),
        prechunk=median_phase(run -> run.prechunk),
        raycore=median_phase(run -> run.raycore),
        validation=median_phase(run -> run.validation),
        retry=median_phase(run -> run.retry),
        cache_mode_phase=median_phase(run -> run.cache_mode_timed),
        resolved_backend=last.resolved_backend,
        fallback_reason=last.fallback_reason,
        cache_mode=last.cache_mode,
        raycore_instances=last.raycore_instances,
        raycore_geometries=last.raycore_geometries,
        raycore_chunked=last.raycore_chunked,
        geometry_mode=last.geometry_mode,
        validation_ratio=last.validation_ratio,
        validation_reference_pixels=last.validation_reference_pixels,
        validation_raycore_pixels=last.validation_raycore_pixels,
        retry_used=last.retry_used,
    )
end

function _format_phase(phase)
    return @sprintf("%7.3f s/%8.1f MiB", phase.seconds, phase.mib)
end

function _print_prepare_breakdown_result(workload_name, case_name, timing)
    total_seconds =
        timing.construct.seconds +
        timing.check.seconds +
        timing.resolve.seconds +
        timing.prepared.seconds +
        timing.render.seconds +
        timing.prechunk.seconds +
        timing.raycore.seconds +
        timing.validation.seconds +
        timing.retry.seconds +
        timing.cache_mode_phase.seconds
    total_mib =
        timing.construct.mib +
        timing.check.mib +
        timing.resolve.mib +
        timing.prepared.mib +
        timing.render.mib +
        timing.prechunk.mib +
        timing.raycore.mib +
        timing.validation.mib +
        timing.retry.mib +
        timing.cache_mode_phase.mib
    @printf(
        "%-22s %-18s resolved=%-26s fallback=%-30s mode=%-20s prepared=%-22s raycore=%-22s validation=%-22s retry=%-22s total=%7.3f s/%8.1f MiB  instances=%5d  geometries=%5d  chunked=%5s  ratio=%6.3f  ref=%7d  gpu=%7d  retry=%s\n",
        workload_name,
        case_name,
        timing.resolved_backend,
        string(timing.fallback_reason),
        string(timing.geometry_mode),
        _format_phase(timing.prepared),
        _format_phase(timing.raycore),
        _format_phase(timing.validation),
        _format_phase(timing.retry),
        total_seconds,
        total_mib,
        timing.raycore_instances,
        timing.raycore_geometries,
        string(timing.raycore_chunked),
        timing.validation_ratio,
        timing.validation_reference_pixels,
        timing.validation_raycore_pixels,
        string(timing.retry_used),
    )
    return (
        workload=workload_name,
        backend=case_name,
        resolved_backend=timing.resolved_backend,
        fallback_reason=timing.fallback_reason,
        prepared_seconds=timing.prepared.seconds,
        raycore_seconds=timing.raycore.seconds,
        validation_seconds=timing.validation.seconds,
        retry_seconds=timing.retry.seconds,
        total_seconds=total_seconds,
        total_mib=total_mib,
        raycore_instances=timing.raycore_instances,
        raycore_geometries=timing.raycore_geometries,
        raycore_chunked=timing.raycore_chunked,
        geometry_mode=timing.geometry_mode,
        validation_ratio=timing.validation_ratio,
        retry_used=timing.retry_used,
    )
end

function main_prepare_breakdown()
    workloads = _selected_workloads()
    cases = _cases()
    results = NamedTuple[]

    println("Local realistic backend prepare breakdown")
    println("samples after warmup: ", BENCH_BUILD_SAMPLES)
    println("max hits: ", _env_int("ARCHIMEDLIGHT_BENCH_MAX_HITS", 32))
    println("workgroup size: ", _env_int("ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", 256))
    println("max prechunk instances: ", BENCH_MAX_PRECHUNK_INSTANCES === nothing ? "default" : string(BENCH_MAX_PRECHUNK_INSTANCES <= 0 ? "disabled" : BENCH_MAX_PRECHUNK_INSTANCES))
    println("reference instancing: ", get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING", "1"))
    println("phase cells report median seconds / median allocated MiB")
    println()
    @printf("%-22s %-18s %-35s %-39s %-25s %-22s %-22s %-22s %-22s %-21s %-13s %-14s %-9s %-10s %-9s %-9s %-8s\n", "workload", "backend", "resolved backend", "fallback", "geometry mode", "prepared", "raycore", "validation", "retry", "total", "instances", "geometries", "chunked", "ratio", "ref", "gpu", "retry")

    for workload in workloads
        @info "Workload" name=workload.name steps=workload.steps nodes=length(workload.scene.nodes)
        for case in cases
            push!(results, _print_prepare_breakdown_result(workload.name, case.name, _time_prepare_breakdown(workload, case)))
        end
    end
    return results
end

function _time_breakdown(workload, case)
    # One unreported run pays Julia method compilation and GPU kernel compilation.
    _time_breakdown_once(workload, case)

    construct_times = Float64[]
    prepare_times = Float64[]
    populate_times = Float64[]
    reuse_times = Float64[]
    construct_mib = Float64[]
    prepare_mib = Float64[]
    populate_mib = Float64[]
    reuse_mib = Float64[]
    last = nothing
    for _ in 1:BENCH_BUILD_SAMPLES
        run = _time_breakdown_once(workload, case)
        push!(construct_times, run.construct_seconds)
        push!(prepare_times, run.prepare_seconds)
        push!(populate_times, run.populate_seconds)
        push!(reuse_times, run.reuse_seconds)
        push!(construct_mib, _mib(run.construct_bytes))
        push!(prepare_mib, _mib(run.prepare_bytes))
        push!(populate_mib, _mib(run.populate_bytes))
        push!(reuse_mib, _mib(run.reuse_bytes))
        last = run
    end
    totals = _totals(last.reuse_series)
    return (
        construct=_median_min_max(construct_times),
        prepare=_median_min_max(prepare_times),
        populate=_median_min_max(populate_times),
        reuse=_median_min_max(reuse_times),
        construct_mib=_median_min_max(construct_mib),
        prepare_mib=_median_min_max(prepare_mib),
        populate_mib=_median_min_max(populate_mib),
        reuse_mib=_median_min_max(reuse_mib),
        resolved_backend=last.resolved_backend,
        fallback_reason=last.fallback_reason,
        geometry_mode=last.geometry_mode,
        summary=last.summary,
        incident_par=totals.incident_par,
        incident_nir=totals.incident_nir,
    )
end

function _print_result(workload_name, case_name, build, warm)
    build_totals = _totals(build.series)
    warm_totals = _totals(warm.series)
    @printf(
        "%-22s %-18s resolved=%-26s fallback=%-30s mode=%-20s build_median=%8.3f s  build_min=%8.3f s  cached_median=%8.3f s  PAR=%12.4e  NIR=%12.4e\n",
        workload_name,
        case_name,
        build.resolved_backend,
        string(build.fallback_reason),
        string(build.geometry_mode),
        build.median,
        build.min,
        warm.median,
        warm_totals.incident_par,
        warm_totals.incident_nir,
    )
    return (
        workload=workload_name,
        backend=case_name,
        resolved_backend=build.resolved_backend,
        fallback_reason=build.fallback_reason,
        geometry_mode=build.geometry_mode,
        build_median_seconds=build.median,
        build_min_seconds=build.min,
        build_max_seconds=build.max,
        cached_median_seconds=warm.median,
        cached_min_seconds=warm.min,
        cached_max_seconds=warm.max,
        incident_par=warm_totals.incident_par,
        incident_nir=warm_totals.incident_nir,
        build_incident_par=build_totals.incident_par,
        build_incident_nir=build_totals.incident_nir,
    )
end

function _print_breakdown_result(workload_name, case_name, timing)
    summary = timing.summary
    @printf(
        "%-22s %-18s resolved=%-26s fallback=%-30s mode=%-20s construct=%8.3f s/%8.1f MiB  prepare=%8.3f s/%8.1f MiB  populate=%8.3f s/%8.1f MiB  reuse=%8.3f s/%8.1f MiB  cached_turtles=%3d  full_sectors=%4d  PAR=%12.4e  NIR=%12.4e\n",
        workload_name,
        case_name,
        timing.resolved_backend,
        string(timing.fallback_reason),
        string(timing.geometry_mode),
        timing.construct.median,
        timing.construct_mib.median,
        timing.prepare.median,
        timing.prepare_mib.median,
        timing.populate.median,
        timing.populate_mib.median,
        timing.reuse.median,
        timing.reuse_mib.median,
        summary.cached_turtle_count,
        summary.cached_full_response_sector_count,
        timing.incident_par,
        timing.incident_nir,
    )
    return (
        workload=workload_name,
        backend=case_name,
        resolved_backend=timing.resolved_backend,
        fallback_reason=timing.fallback_reason,
        geometry_mode=timing.geometry_mode,
        construct_median_seconds=timing.construct.median,
        prepare_median_seconds=timing.prepare.median,
        populate_median_seconds=timing.populate.median,
        reuse_median_seconds=timing.reuse.median,
        construct_median_mib=timing.construct_mib.median,
        prepare_median_mib=timing.prepare_mib.median,
        populate_median_mib=timing.populate_mib.median,
        reuse_median_mib=timing.reuse_mib.median,
        cached_turtle_count=summary.cached_turtle_count,
        cached_full_response_sector_count=summary.cached_full_response_sector_count,
        incident_par=timing.incident_par,
        incident_nir=timing.incident_nir,
    )
end

function _relative_delta(candidate, reference)
    reference == 0 && return candidate == 0 ? 0.0 : Inf
    return abs(candidate - reference) / abs(reference)
end

function main_breakdown()
    workloads = _selected_workloads()
    cases = _cases()
    results = NamedTuple[]

    println("Local realistic backend cache breakdown")
    println("samples after warmup: ", BENCH_BUILD_SAMPLES)
    println("max hits: ", _env_int("ARCHIMEDLIGHT_BENCH_MAX_HITS", 32))
    println("workgroup size: ", _env_int("ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", 256))
    println("max prechunk instances: ", BENCH_MAX_PRECHUNK_INSTANCES === nothing ? "default" : string(BENCH_MAX_PRECHUNK_INSTANCES <= 0 ? "disabled" : BENCH_MAX_PRECHUNK_INSTANCES))
    println("reference instancing: ", get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING", "1"))
    println("phases: construct=LightSimulation constructor, prepare=_ensure_light_cache!, populate=first run after prepare, reuse=second run on same simulation")
    println("phase cells report median seconds / median allocated MiB")
    println()
    @printf("%-22s %-18s %-35s %-39s %-25s %-21s %-21s %-21s %-21s %-16s %-14s %-14s\n", "workload", "backend", "resolved backend", "fallback", "geometry mode", "construct", "prepare", "populate", "reuse", "cache", "PAR", "NIR")

    for workload in workloads
        @info "Workload" name=workload.name steps=workload.steps nodes=length(workload.scene.nodes)
        for case in cases
            push!(results, _print_breakdown_result(workload.name, case.name, _time_breakdown(workload, case)))
        end
    end

    println()
    println("Phase speedups are relative to normal_cpu within each workload.")
    for workload in workloads
        base_results = [r for r in results if r.workload == workload.name && r.backend == "normal_cpu"]
        if isempty(base_results)
            @printf("%-22s %-18s normal_cpu baseline was not selected; skipping phase speedups.\n", workload.name, "")
            continue
        end
        base = only(base_results)
        for result in (r for r in results if r.workload == workload.name)
            @printf(
                "%-22s %-18s prepare_speedup=%7.3fx  populate_speedup=%7.3fx  reuse_speedup=%7.3fx\n",
                result.workload,
                result.backend,
                base.prepare_median_seconds / result.prepare_median_seconds,
                base.populate_median_seconds / result.populate_median_seconds,
                base.reuse_median_seconds / result.reuse_median_seconds,
            )
        end
    end

    return results
end

function main()
    BENCH_PREPARE_BREAKDOWN && return main_prepare_breakdown()
    BENCH_BREAKDOWN && return main_breakdown()

    workloads = _selected_workloads()
    cases = _cases()
    results = NamedTuple[]

    println("Local realistic backend benchmark")
    println("build+run samples after warmup: ", BENCH_BUILD_SAMPLES)
    println("cached reuse samples: ", BENCH_WARM_SAMPLES)
    println("max hits: ", _env_int("ARCHIMEDLIGHT_BENCH_MAX_HITS", 32))
    println("workgroup size: ", _env_int("ARCHIMEDLIGHT_BENCH_WORKGROUPSIZE", 256))
    println("max prechunk instances: ", BENCH_MAX_PRECHUNK_INSTANCES === nothing ? "default" : string(BENCH_MAX_PRECHUNK_INSTANCES <= 0 ? "disabled" : BENCH_MAX_PRECHUNK_INSTANCES))
    println("reference instancing: ", get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCING", "1"))
    println()
    @printf("%-22s %-18s %-35s %-39s %-25s %-22s %-20s %-22s %-14s %-14s\n", "workload", "backend", "resolved backend", "fallback", "geometry mode", "build median", "build min", "cached median", "PAR", "NIR")

    for workload in workloads
        @info "Workload" name=workload.name steps=workload.steps nodes=length(workload.scene.nodes)
        for case in cases
            build = _time_build(workload, case)
            warm = _time_warm(workload, build.sim)
            push!(results, _print_result(workload.name, case.name, build, warm))
        end
    end

    println()
    println("Speedups are relative to normal_cpu warm median within each workload.")
    for workload in workloads
        base_results = [r for r in results if r.workload == workload.name && r.backend == "normal_cpu"]
        if isempty(base_results)
            @printf("%-22s %-18s normal_cpu baseline was not selected; skipping speedups.\n", workload.name, "")
            continue
        end
        base = only(base_results)
        for result in (r for r in results if r.workload == workload.name)
            build_speedup = base.build_median_seconds / result.build_median_seconds
            cached_speedup = base.cached_median_seconds / result.cached_median_seconds
            @printf(
                "%-22s %-18s resolved=%-16s fallback=%-30s mode=%-20s build_speedup=%7.3fx  cached_speedup=%7.3fx\n",
                result.workload,
                result.backend,
                result.resolved_backend,
                string(result.fallback_reason),
                string(result.geometry_mode),
                build_speedup,
                cached_speedup,
            )
        end
    end

    println()
    println("Energy deltas are relative to normal_cpu cached totals within each workload.")
    for workload in workloads
        base_results = [r for r in results if r.workload == workload.name && r.backend == "normal_cpu"]
        if isempty(base_results)
            @printf("%-22s %-18s normal_cpu baseline was not selected; skipping energy deltas.\n", workload.name, "")
            continue
        end
        base = only(base_results)
        for result in (r for r in results if r.workload == workload.name)
            @printf(
                "%-22s %-18s PAR_delta=%8.3f%%  NIR_delta=%8.3f%%\n",
                result.workload,
                result.backend,
                100 * _relative_delta(result.incident_par, base.incident_par),
                100 * _relative_delta(result.incident_nir, base.incident_nir),
            )
        end
    end

    return results
end

main()
