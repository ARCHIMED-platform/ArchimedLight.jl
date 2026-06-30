using GeometryBasics
using MultiScaleTreeGraph
using PlantGeom

const LOCAL_GPU_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_GPU_LOCAL", ""))
const LOCAL_GPU_REQUESTED = LOCAL_GPU_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const LOCAL_GPU_REQUIRED = LOCAL_GPU_FLAG in ("required", "force", "error")

const LOCAL_GPU_LARGE_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_GPU_LARGE", ""))
const LOCAL_GPU_LARGE_REQUESTED = LOCAL_GPU_LARGE_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const LOCAL_GPU_LARGE_REQUIRED = LOCAL_GPU_LARGE_FLAG in ("required", "force", "error")

function _local_skip_or_fail(message::AbstractString; required::Bool=LOCAL_GPU_REQUIRED)
    if required
        error(message)
    else
        @test_skip message
        return nothing
    end
end

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

function _step_totals(step)
    totals = ArchimedLight._light_step_energy_totals(step)
    return (
        incident_par=Float64(totals.incident_par_total),
        incident_nir=Float64(totals.incident_nir_total),
        absorbed_par=Float64(totals.absorbed_par_total),
        absorbed_nir=Float64(totals.absorbed_nir_total),
    )
end

function _assert_valid_step(step)
    totals = _step_totals(step)
    @test step.scattering !== nothing
    if step.scattering !== nothing && !step.scattering.converged
        @info "Local GPU test scattering reached the iteration limit" iterations=step.scattering.iterations
    end
    @test !isempty(step.budget.incident_energy.total.par)
    for value in values(totals)
        @test isfinite(value)
        @test value >= -sqrt(eps(Float64))
    end
    @test totals.incident_par > 0.0
end

function _assert_series_close(candidate, reference; label, rtol=5e-3, atol=1e-3)
    @test length(candidate) == length(reference)
    for i in eachindex(reference)
        cand = _step_totals(candidate[i])
        ref = _step_totals(reference[i])
        for key in keys(ref)
            @test isapprox(cand[key], ref[key]; rtol=rtol, atol=atol)
        end
    end
    @info "Local GPU parity" label steps=length(reference) rtol atol
end

function _run_meteo_series(name, scene, models, meteo, options; interception_backend=:raster_cpu, scattering_backend=nothing)
    sim = ArchimedLight.LightSimulation(
        scene,
        models;
        options=options,
        interception_backend=interception_backend,
        scattering_backend=scattering_backend,
    )
    elapsed = @elapsed series = ArchimedLight.run_light(sim, meteo)
    @info "Local meteo GPU test run" backend=name seconds=round(elapsed; digits=3) cache=ArchimedLight.cache_summary(sim)
    return (series=series, sim=sim)
end

function _run_sky_series(name, scene, models, skies, options; interception_backend=:raster_cpu, scattering_backend=nothing)
    sim = ArchimedLight.LightSimulation(
        scene,
        models;
        options=options,
        interception_backend=interception_backend,
        scattering_backend=scattering_backend,
    )
    dt = _env_float("ARCHIMEDLIGHT_TEST_GPU_LOCAL_STEP_SECONDS", 1800.0)
    elapsed = @elapsed series = [
        ArchimedLight.run_light(sim, sky; step_duration_seconds=dt)
        for sky in skies
    ]
    @info "Local sky GPU test run" backend=name seconds=round(elapsed; digits=3) cache=ArchimedLight.cache_summary(sim)
    return (series=series, sim=sim)
end

function _assert_resolved_interception_backend(run, backend_type)
    @test run.sim.cache !== nothing
    if run.sim.cache !== nothing
        @test run.sim.cache.resolved_interception_backend isa backend_type
    end
end

function _raycore_backends(; metal=false)
    backend = metal ? _metal_backend() : KernelAbstractions.CPU()
    backend === nothing && return nothing
    interception = ArchimedLight.RaycoreInterceptionBackend(
        backend=backend,
        max_hits_per_pixel=_env_int("ARCHIMEDLIGHT_TEST_GPU_MAX_HITS", 32),
        workgroupsize=_env_int("ARCHIMEDLIGHT_TEST_GPU_WORKGROUPSIZE", 256),
        validate=true,
    )
    return (
        interception=interception,
        scattering=ArchimedLight.RaycoreScatteringBackend(interception),
    )
end

function _coffee_fixture()
    config = joinpath(@__DIR__, "..", "..", "example_2", "config.yml")
    paving = _env_int("ARCHIMEDLIGHT_TEST_GPU_LOCAL_PAVING", 2500)
    options, scene, meteo, models = ArchimedLight.read_config(config; plot_paving_override=paving)
    options = ArchimedLight.LightOptions(
        options;
        cache_radiation=true,
        meteo_range=nothing,
        scattering=true,
    )
    return scene, models, meteo, options
end

function _vpalm_modules()
    xpalm_root = get(ENV, "ARCHIMEDLIGHT_TEST_XPALM_ROOT", "/Users/rvezy/Documents/dev/XPalm")
    isdir(xpalm_root) || return (nothing, "XPalm checkout not found at $xpalm_root")
    package_dir = dirname(xpalm_root)

    old_load_path = copy(LOAD_PATH)
    try
        filter!(path -> path != xpalm_root && path != package_dir, LOAD_PATH)
        pushfirst!(LOAD_PATH, package_dir)
        pushfirst!(LOAD_PATH, xpalm_root)

        xpalm_mod = try
            Base.require(Main, :XPalm)
        catch err
            return (nothing, "XPalm could not be loaded from $package_dir: $(sprint(showerror, err))")
        end

        vpalm_mod = try
            if !Base.invokelatest(isdefined, xpalm_mod, :VPalm)
                Base.invokelatest(Base.include, xpalm_mod, joinpath(xpalm_root, "src", "VPalm.jl"))
            end
            Base.invokelatest(getfield, xpalm_mod, :VPalm)
        catch err
            return (nothing, "XPalm.VPalm could not be loaded: $(sprint(showerror, err))")
        end

        return ((XPalm=xpalm_mod, VPalm=vpalm_mod, root=xpalm_root), nothing)
    finally
        empty!(LOAD_PATH)
        append!(LOAD_PATH, old_load_path)
    end
end

function _vpalm_merge_scale()
    raw = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_XPALM_MERGE_SCALE", "leaf"))
    raw == "none" && return :none
    raw == "leaflet" && return :leaflet
    raw == "leaf" && return :leaf
    raw == "plant" && return :plant
    error("Unsupported ARCHIMEDLIGHT_TEST_XPALM_MERGE_SCALE=$raw. Use none, leaflet, leaf, or plant.")
end

function _vpalm_two_palm_fixture()
    modules, err = _vpalm_modules()
    modules === nothing && return (nothing, err)

    VPalm = modules.VPalm
    parameter_file = joinpath(modules.root, "test", "references", "vpalm-parameter_file.yml")
    isfile(parameter_file) || return (nothing, "VPalm parameter file not found at $parameter_file")
    read_parameters = Base.invokelatest(getproperty, VPalm, :read_parameters)
    parameters = Base.invokelatest(read_parameters, parameter_file)
    merge_scale = _vpalm_merge_scale()
    build_mockup = Base.invokelatest(getproperty, VPalm, :build_mockup)
    palm = Base.invokelatest(build_mockup, parameters; merge_scale=merge_scale, rng=nothing)

    spacing = _env_float("ARCHIMEDLIGHT_TEST_XPALM_SPACING", 6.0)
    ground_n = _env_int("ARCHIMEDLIGHT_TEST_XPALM_GROUND_N", 12)
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
        turtle_sectors=_env_int("ARCHIMEDLIGHT_TEST_GPU_LARGE_SECTORS", 8),
        pixel_size=_env_float("ARCHIMEDLIGHT_TEST_GPU_LARGE_PIXEL_SIZE", 0.05),
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
        toricity=true,
        scattering_max_iter=_env_int("ARCHIMEDLIGHT_TEST_GPU_LARGE_SCATTERING_MAX_ITER", 12),
    )

    skies = [
        ArchimedLight.SkyState(120.0, 35.0, 280.0, 180.0, 0.65, 0.35),
        ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
    ]

    return ((scene=scene, models=models, options=options, skies=skies, source="XPalm VPalm", detail=merge_scale), nothing)
end

function _agripv_panel_mesh(width, length, height, inclination_deg)
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
    faces = GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)]
    return GeometryBasics.Mesh(points, faces)
end

function _agripv_wheat_fixture()
    root = get(
        ENV,
        "ARCHIMEDLIGHT_TEST_AGRIPV_ROOT",
        "/Users/rvezy/Documents/cirad/Articles_Rapports_Communications/Conférences/2026_fspm/coutellier_agripv",
    )
    opf_path = joinpath(root, "0_simulations", "archicrop", "wheat", "plant_1995-06-24.opf")
    isfile(opf_path) || return (nothing, "Agrivoltaics wheat OPF not found at $opf_path")

    n_rows = _env_int("ARCHIMEDLIGHT_TEST_AGRIPV_ROWS", 2)
    plants_per_row = _env_int("ARCHIMEDLIGHT_TEST_AGRIPV_PLANTS_PER_ROW", 8)
    interrow = _env_float("ARCHIMEDLIGHT_TEST_AGRIPV_INTERROW", 0.25)
    intrarow = _env_float("ARCHIMEDLIGHT_TEST_AGRIPV_INTRAROW", 0.30)
    ground_nx = _env_int("ARCHIMEDLIGHT_TEST_AGRIPV_GROUND_NX", 18)
    ground_ny = _env_int("ARCHIMEDLIGHT_TEST_AGRIPV_GROUND_NY", 24)

    width = max(interrow * n_rows, interrow)
    length = max(intrarow * plants_per_row, intrarow)
    wheat = PlantGeom.read_opf(opf_path, mtg_type=MultiScaleTreeGraph.NodeMTG)
    panel = _agripv_panel_mesh(
        width,
        _env_float("ARCHIMEDLIGHT_TEST_AGRIPV_PANEL_LENGTH", 3.0),
        _env_float("ARCHIMEDLIGHT_TEST_AGRIPV_PANEL_HEIGHT", 1.6),
        _env_float("ARCHIMEDLIGHT_TEST_AGRIPV_PANEL_INCLINATION", 25.0),
    )

    scene = PlantGeom.make_scene(domain=(0.0, 0.0, width, length)) do s
        PlantGeom.add_object!(s, panel; group="panel", type="Panel", id=1, at=(0.0, length / 2, 0.0))
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

    options = ArchimedLight.LightOptions(
        turtle_sectors=_env_int("ARCHIMEDLIGHT_TEST_GPU_LARGE_SECTORS", 8),
        pixel_size=_env_float("ARCHIMEDLIGHT_TEST_GPU_LARGE_PIXEL_SIZE", 0.05),
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
        toricity=true,
        scattering_max_iter=_env_int("ARCHIMEDLIGHT_TEST_GPU_LARGE_SCATTERING_MAX_ITER", 12),
    )

    skies = [
        ArchimedLight.SkyState(120.0, 35.0, 280.0, 180.0, 0.65, 0.35),
        ArchimedLight.SkyState(210.0, 52.0, 360.0, 260.0, 0.70, 0.30),
    ]

    detail = "$(n_rows)x$(plants_per_row) wheat plants"
    return ((scene=scene, models=models, options=options, skies=skies, source="agrivoltaics wheat", detail=detail), nothing)
end

@testset "Local cached coffee multi-step backends" begin
    if !LOCAL_GPU_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_GPU_LOCAL=1 to run local multi-step backend tests."
    else
        scene, models, meteo, options = _coffee_fixture()

        reference_run = _run_meteo_series("raster_cpu", scene, models, meteo, options)
        reference = reference_run.series
        _assert_resolved_interception_backend(reference_run, ArchimedLight.RasterCPUBackend)
        @test length(reference) >= 3
        foreach(_assert_valid_step, reference)

        ka = _raycore_backends()
        ka_run = _run_meteo_series(
            "raycore_ka_cpu",
            scene,
            models,
            meteo,
            options;
            interception_backend=ka.interception,
            scattering_backend=ka.scattering,
        )
        ka_series = ka_run.series
        _assert_resolved_interception_backend(ka_run, ArchimedLight.RaycoreInterceptionBackend)
        foreach(_assert_valid_step, ka_series)
        _assert_series_close(ka_series, reference; label="raycore_ka_cpu_vs_raster_cpu", rtol=2.5e-1, atol=1e-2)

        if METAL_REQUESTED
            metal = _raycore_backends(; metal=true)
            if metal !== nothing
                metal_run = _run_meteo_series(
                    "raycore_metal_gpu",
                    scene,
                    models,
                    meteo,
                    options;
                    interception_backend=metal.interception,
                    scattering_backend=metal.scattering,
                )
                metal_series = metal_run.series
                _assert_resolved_interception_backend(metal_run, ArchimedLight.RaycoreInterceptionBackend)
                foreach(_assert_valid_step, metal_series)
                _assert_series_close(metal_series, ka_series; label="raycore_metal_gpu_vs_ka_cpu", rtol=3e-2, atol=1e-2)
                _assert_series_close(metal_series, reference; label="raycore_metal_gpu_vs_raster_cpu", rtol=2.5e-1, atol=1e-2)
            end
        else
            @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to include Metal GPU in the local multi-step test."
        end
    end
end

@testset "Local realistic large scene" begin
    if !LOCAL_GPU_LARGE_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_GPU_LARGE=1 to run the local realistic large scene."
    else
        source = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_GPU_LARGE_SOURCE", "agripv"))
        load_path_before_fixture = copy(LOAD_PATH)
        fixture, err = if source in ("xpalm", "vpalm")
            _vpalm_two_palm_fixture()
        elseif source in ("auto", "any")
            vpalm_fixture, vpalm_err = _vpalm_two_palm_fixture()
            if vpalm_fixture === nothing
                @info "VPalm fixture unavailable; trying agrivoltaics fallback" reason=vpalm_err
                _agripv_wheat_fixture()
            else
                (vpalm_fixture, nothing)
            end
        elseif source in ("agripv", "agri", "wheat")
            _agripv_wheat_fixture()
        else
            (nothing, "Unsupported ARCHIMEDLIGHT_TEST_GPU_LARGE_SOURCE=$source. Use agripv, xpalm, or auto.")
        end
        if source in ("xpalm", "vpalm", "auto", "any")
            @test LOAD_PATH == load_path_before_fixture
        end
        if fixture === nothing
            _local_skip_or_fail(err; required=LOCAL_GPU_LARGE_REQUIRED)
        else
            @info "Local large fixture" source=fixture.source detail=fixture.detail nodes=length(fixture.scene.nodes) skies=length(fixture.skies)

            ka = _raycore_backends()
            ka_run = _run_sky_series(
                "raycore_ka_cpu_large",
                fixture.scene,
                fixture.models,
                fixture.skies,
                fixture.options;
                interception_backend=ka.interception,
                scattering_backend=ka.scattering,
            )
            ka_series = ka_run.series
            _assert_resolved_interception_backend(ka_run, ArchimedLight.RaycoreInterceptionBackend)
            @test length(ka_series) == length(fixture.skies)
            foreach(_assert_valid_step, ka_series)

            if lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_GPU_LARGE_REFERENCE", "")) in ("1", "true", "yes", "on")
                reference_run = _run_sky_series("raster_cpu_large", fixture.scene, fixture.models, fixture.skies, fixture.options)
                reference = reference_run.series
                _assert_resolved_interception_backend(reference_run, ArchimedLight.RasterCPUBackend)
                foreach(_assert_valid_step, reference)
                _assert_series_close(ka_series, reference; label="large_raycore_ka_cpu_vs_raster_cpu", rtol=2e-2, atol=1e-2)
            end

            if METAL_REQUESTED
                metal = _raycore_backends(; metal=true)
                if metal !== nothing
                    metal_run = try
                        _run_sky_series(
                            "raycore_metal_gpu_large",
                            fixture.scene,
                            fixture.models,
                            fixture.skies,
                            fixture.options;
                            interception_backend=metal.interception,
                            scattering_backend=metal.scattering,
                        )
                    catch err
                        @test err isa ArchimedLight.RaycoreValidationError
                        if !(err isa ArchimedLight.RaycoreValidationError)
                            rethrow()
                        end
                        @info "Metal large scene failed explicit Raycore validation." reason=err.reason stage=err.stage
                        nothing
                    end
                    if metal_run !== nothing
                        metal_series = metal_run.series
                        @test length(metal_series) == length(fixture.skies)
                        foreach(_assert_valid_step, metal_series)
                        resolved = metal_run.sim.cache === nothing ? nothing : metal_run.sim.cache.resolved_interception_backend
                        @test resolved isa ArchimedLight.RaycoreInterceptionBackend
                        _assert_series_close(metal_series, ka_series; label="large_raycore_metal_vs_ka_cpu", rtol=2e-2, atol=1e-2)
                    end
                end
            else
                @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to include Metal GPU in the large-scene test."
            end
        end
    end
end
