const _BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", ""))
const _BENCH_WANTS_METAL = _BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

include(joinpath(@__DIR__, "..", "scripts", "raycore_device_smoke.jl"))

using Dates
using Printf
using Statistics

const BENCH_WARMUPS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_WARMUPS", "2"))
const BENCH_SAMPLES = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_SAMPLES", "5"))
const BENCH_SECTORS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_SECTORS", "46"))
const BENCH_PIXEL_SIZE = parse(Float64, get(ENV, "ARCHIMEDLIGHT_BENCH_PIXEL_SIZE", "0.01"))
const BENCH_WORKLOAD = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_WORKLOAD", "smoke"))
const BENCH_PLOT_PAVING = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_PLOT_PAVING", "2500"))
const BENCH_BACKENDS = split(lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_BACKENDS", "all")), ',')
const BENCH_MAX_HITS = parse(Int, get(ENV, "ARCHIMEDLIGHT_BENCH_MAX_HITS", "32"))
const BENCH_EXECUTION = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_EXECUTION", "oneshot"))

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

function _backend_cases()
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
            interception_backend=ArchimedLight.RaycoreInterceptionBackend(
                backend=KernelAbstractions.CPU(),
                max_hits_per_pixel=BENCH_MAX_HITS,
            ),
            scattering_backend=nothing,
        ),
    ]
    cases[2] = merge(cases[2], (; scattering_backend=ArchimedLight.RaycoreScatteringBackend(cases[2].interception_backend)))

    metal_backend = _optional_metal_backend()
    if metal_backend !== nothing
        ib = ArchimedLight.RaycoreInterceptionBackend(
            backend=metal_backend,
            max_hits_per_pixel=BENCH_MAX_HITS,
        )
        push!(
            cases,
            (
                name="raycore_metal_gpu",
                backend=metal_backend,
                interception_backend=ib,
                scattering_backend=ArchimedLight.RaycoreScatteringBackend(ib),
            ),
        )
    end
    if BENCH_BACKENDS == ["all"]
        return cases
    end
    selected = Set(strip.(BENCH_BACKENDS))
    filtered = [case for case in cases if lowercase(case.name) in selected]
    isempty(filtered) && error(
        "ARCHIMEDLIGHT_BENCH_BACKENDS selected no backends. " *
        "Use all, normal_cpu, raycore_ka_cpu, or raycore_metal_gpu.",
    )
    return filtered
end

function _with_scattering(options::ArchimedLight.LightOptions, scattering::Bool)
    return ArchimedLight.LightOptions(options; scattering=scattering)
end

function _benchmark_fixture()
    if BENCH_WORKLOAD == "smoke"
        scene = _smoke_scene()
        models = _smoke_models()
        meteo_row = _benchmark_meteo_row(; ri_par_f=100.0, ri_nir_f=50.0)
        first_order_options = ArchimedLight.LightOptions(
            turtle_sectors=BENCH_SECTORS,
            all_in_turtle=false,
            scattering=false,
            pixel_size=BENCH_PIXEL_SIZE,
            toricity=false,
        )
        scattering_options = ArchimedLight.LightOptions(
            turtle_sectors=BENCH_SECTORS,
            all_in_turtle=false,
            scattering=true,
            pixel_size=BENCH_PIXEL_SIZE,
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
        options, scene, meteo, models = ArchimedLight.read_config(config_path; plot_paving_override=BENCH_PLOT_PAVING)
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
    return () -> ArchimedLight.run_light(sim, meteo_row)
end

function _first_order_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return () -> _first_order_call(scene, models, meteo_row, options, case)
    elseif BENCH_EXECUTION == "cached"
        return _cached_call(scene, models, meteo_row, options, case)
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _scattering_benchmark_call(scene, models, meteo_row, options, case)
    if BENCH_EXECUTION == "oneshot"
        return () -> _scattering_call(scene, models, meteo_row, options, case)
    elseif BENCH_EXECUTION == "cached"
        return _cached_call(scene, models, meteo_row, options, case; scattering_mode=:raycast)
    end
    error("Unsupported ARCHIMEDLIGHT_BENCH_EXECUTION=$BENCH_EXECUTION. Use oneshot or cached.")
end

function _print_markdown_table(results)
    println("| workload | backend | median ms | min ms | max ms | median alloc MiB |")
    println("|---|---:|---:|---:|---:|---:|")
    for row in results
        @printf(
            "| %s | %s | %.3f | %.3f | %.3f | %.3f |\n",
            row.workload,
            row.backend,
            row.median_ms,
            row.min_ms,
            row.max_ms,
            row.median_alloc_mb,
        )
    end
end

function main()
    fixture = _benchmark_fixture()
    scene = fixture.scene
    models = fixture.models
    meteo_row = fixture.meteo_row
    first_order_options = fixture.first_order_options
    scattering_options = fixture.scattering_options

    cases = _backend_cases()
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

    results = NamedTuple[]
    for case in cases
        push!(
            results,
            merge(
                _measure("first_order", _first_order_benchmark_call(scene, models, meteo_row, first_order_options, case)),
                (; workload="first_order", backend=case.name),
            ),
        )
        push!(
            results,
            merge(
                _measure("scattering", _scattering_benchmark_call(scene, models, meteo_row, scattering_options, case)),
                (; workload="scattering", backend=case.name),
            ),
        )
    end

    println()
    println("ArchimedLight backend comparison")
    println(
        "workload=$(fixture.name), warmups=$(BENCH_WARMUPS), samples=$(BENCH_SAMPLES), " *
        "sectors=$(first_order_options.turtle_sectors), pixel_size=$(first_order_options.pixel_size), " *
        "plot_paving=$(BENCH_WORKLOAD == "smoke" ? "n/a" : string(BENCH_PLOT_PAVING)), " *
        "max_hits=$(BENCH_MAX_HITS), execution=$(BENCH_EXECUTION)",
    )
    println()
    _print_markdown_table(results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
