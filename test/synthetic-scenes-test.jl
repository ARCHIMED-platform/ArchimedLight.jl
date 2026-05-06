@testmodule HelperModule begin
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    using GeometryBasics
    using Dates
    using ArchimedLight

    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))
end



@testitem "Synthetic case single_plate_direct" tags = [:synthetic, :fast, :single_plate_direct] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    run = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01))

    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_flux(run.budget), 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
end

@testitem "Synthetic case stacked_scattering" tags = [:synthetic, :fast, :stacked_scattering] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    run0 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01); models=models)
    @test isapprox(get(run0.first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_energy(run0.budget), 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
    budget = ArchimedLight.integrate_light(scene, models, first, scat, options; step_duration_seconds=1.0)

    @test get(scat.added_power.par, 2, 0.0) > 0.0
    @test get(HelperModule._incident_par_energy(budget), 2, 0.0) > 0.0
end

@testitem "Synthetic case toricity_wraparound" tags = [:synthetic, :fast, :toricity_wraparound] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1)])
    scene.scene_xy_bounds = (0.0, 0.0, 1.0, 1.0)
    sky = ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0)

    run0 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false))
    run1 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true))

    @test get(HelperModule._incident_par_initial_energy(run1.budget), 1, 0.0) > get(HelperModule._incident_par_initial_energy(run0.budget), 1, 0.0)
end

@testitem "Synthetic case virtual_sensor_transparency" tags = [:synthetic, :fast, :virtual_sensor_transparency] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="plant", type="plate", object_id=2),
    ])
    models = HelperModule._virtual_sensor_models()
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    run = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01); models=models)

    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case translucent stack follows Java visible area" tags = [:synthetic, :fast, :translucent_stack] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="leaf", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="leaf", type="plate", object_id=2),
    ])
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "leaf";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.25,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    run = HelperModule._run_direct(scene, sky, options; models=models)

    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(run.first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case run_light_step_matches_staged" tags = [:synthetic, :fast, :run_light_step_matches_staged] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    row = HelperModule._synthetic_meteo_row()

    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    budget = ArchimedLight.integrate_light(scene, models, first, nothing, options; meteo_row=row)
    step = ArchimedLight.run_light_step(scene, models, row, options)

    @test budget.incident_energy.total.par == step.budget.incident_energy.total.par
    @test first.projected_area_per_node == step.first_order.projected_area_per_node
end

@testitem "Synthetic case cache_radiation_parity" tags = [:synthetic, :fast, :cache_radiation_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_cached_series"))
    uncached = HelperModule._synthetic_options(cache_radiation=false)
    cached = HelperModule._synthetic_options(cache_radiation=true)

    series0 = ArchimedLight.run_light_series(scene, models, meteo, uncached)
    series1 = ArchimedLight.run_light_series(scene, models, meteo, cached)

    @test length(series0) == length(series1)
    for i in eachindex(series0)
        @test series0[i].budget.incident_flux.total.par == series1[i].budget.incident_flux.total.par
        @test series0[i].budget.incident_energy.total.par == series1[i].budget.incident_energy.total.par
    end
end

@testitem "Synthetic case PlantMeteo TimeStepTable input" tags = [:synthetic, :fast, :plantmeteo_timestep_table] setup = [HelperModule] begin
    using Dates

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=true)

    rows = [
        (
            date=DateTime(2020, 6, 21, 12, 0, 0),
            duration=Hour(1),
            Ri_PAR_f=120.0,
            Ri_NIR_f=80.0,
            clearness=0.6,
        ),
        (
            date=DateTime(2020, 6, 21, 13, 0, 0),
            duration=Hour(1),
            Ri_PAR_f=140.0,
            Ri_NIR_f=60.0,
            clearness=0.6,
        ),
    ]
    meteo = PM.TimeStepTable(rows, (latitude=15.0, source="synthetic_pm_namedtuple"))

    selected = ArchimedLight.prepare_meteo(meteo, options)
    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    step = ArchimedLight.run_light_step(scene, models, first(selected), options)
    step_raw = ArchimedLight.run_light_step(scene, models, first(meteo), options)

    @test selected isa PM.TimeStepTable
    @test length(selected) == 2
    @test length(series) == 2
    @test first(selected).Ri_PAR_f == 120.0
    @test step.budget.incident_energy.total.par == series[1].budget.incident_energy.total.par
    @test step_raw.budget.incident_energy.total.par == series[1].budget.incident_energy.total.par
end

@testitem "Synthetic case PlantMeteo Atmosphere input" tags = [:synthetic, :fast, :plantmeteo_atmosphere] setup = [HelperModule] begin
    using Dates

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    rows = PM.Atmosphere[
        PM.Atmosphere(date=DateTime(2020, 6, 21, 12, 0, 0), duration=Hour(1), T=25.0, Wind=1.0, Rh=0.6, Ri_PAR_f=120.0, Ri_NIR_f=80.0, clearness=0.6),
        PM.Atmosphere(date=DateTime(2020, 6, 21, 13, 0, 0), duration=Hour(1), T=26.0, Wind=1.0, Rh=0.6, Ri_PAR_f=100.0, Ri_NIR_f=50.0, clearness=0.5),
    ]
    meteo = PM.TimeStepTable(rows, (latitude=15.0, source="synthetic_pm_atmosphere"))

    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    sky = ArchimedLight.compute_sky(first(meteo), options)

    @test length(series) == 2
    @test isapprox(sky.ri_par_f, 120.0; atol=1e-9, rtol=1e-9)
    @test isapprox(sky.ri_nir_f, 80.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case generic table input" tags = [:synthetic, :fast, :generic_table_input] setup = [HelperModule] begin
    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    meteo = [
        (date="2020/06/21", hour_start="12:00:00", hour_end="13:00:00", latitude=15.0, T=25.0, Rh=0.60, Wind=1.0, Ri_SW_f=200.0, clearness=0.6, Cₐ=380.0),
        (date="2020/06/21", hour_start="13:00:00", hour_end="14:00:00", latitude=15.0, T=26.0, Rh=0.60, Wind=1.0, Ri_SW_f=240.0, clearness=0.6, Cₐ=380.0),
    ]

    selected = ArchimedLight.prepare_meteo(meteo, options)
    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    meteo_read = ArchimedLight.read_meteo(meteo)
    sky = ArchimedLight.compute_sky(first(selected), options)

    @test selected isa PM.TimeStepTable
    @test length(selected) == 2
    @test length(series) == 2
    @test meteo_read isa ArchimedLight.MeteoTable
    @test first(selected).Ri_SW_f == 200.0
    @test isapprox(sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case overlapping meteo steps option" tags = [:synthetic, :fast, :overlapping_meteo_steps] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()

    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 0), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 15), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_overlap"))

    strict_options = HelperModule._synthetic_options(cache_radiation=false)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.prepare_meteo(meteo, strict_options)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.run_light_series(scene, models, meteo, strict_options)

    permissive_options = ArchimedLight.LightOptions(strict_options; allow_overlapping_meteo_steps=true)
    selected = ArchimedLight.prepare_meteo(meteo, permissive_options)
    series = ArchimedLight.run_light_series(scene, models, meteo, permissive_options)

    @test length(selected.rows) == 2
    @test length(series) == 2

    config_path = tempname() * ".yml"
    try
        open(config_path, "w") do io
            write(io, "scene: dummy.ops\nmodels:\n  - dummy.yml\nmeteo: dummy.csv\nallowOverlappingMeteoSteps: true\ncomponent_variables:\n  sky_fraction: true\n")
        end
        parsed = ArchimedLight.read_options(config_path)
        @test parsed.allow_overlapping_meteo_steps
        @test parsed.include_sky_fraction
    finally
        rm(config_path; force=true)
    end
end

@testitem "Synthetic case cached_scattering_series_parity" tags = [:synthetic, :fast, :cached_scattering_series_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_cached_scattering_series"))
    uncached = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false)
    cached = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=true)

    series0 = ArchimedLight.run_light_series(scene, models, meteo, uncached)
    series1 = ArchimedLight.run_light_series(scene, models, meteo, cached)

    @test length(series0) == length(series1)
    for i in eachindex(series0)
        @test HelperModule._budgets_close(series0[i].budget, series1[i].budget; atol=1e-9, rtol=1e-9)
    end
end

@testitem "Synthetic case light_cache_manual_api" tags = [:synthetic, :fast, :light_cache_manual_api] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=50.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(14), duration_seconds=900.0, ri_par_f=90.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_manual_cache"))

    cache = ArchimedLight.prepare_light_cache(scene, models, options)
    summary0 = ArchimedLight.cache_summary(cache)
    @test summary0.mode == :full
    @test summary0.cached_turtle_count == 0

    series_cached = ArchimedLight.run_light_series(cache, meteo)
    series_uncached = ArchimedLight.run_light_series(scene, models, meteo, HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))

    @test length(series_cached) == length(series_uncached)
    for i in eachindex(series_cached)
        @test HelperModule._budgets_close(series_cached[i].budget, series_uncached[i].budget; atol=1e-9, rtol=1e-9)
    end

    step_cached = ArchimedLight.run_light_step(cache, rows[1])
    step_uncached = ArchimedLight.run_light_step(scene, models, rows[1], HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))
    @test HelperModule._budgets_close(step_cached.budget, step_uncached.budget; atol=1e-9, rtol=1e-9)

    summary1 = ArchimedLight.cache_summary(cache)
    @test summary1.cached_turtle_count == 1
    @test summary1.cached_full_response_sector_count > 0

    scene2 = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.5, group="middle", type="plate", object_id=3),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    cache2 = ArchimedLight.prepare_light_cache(scene2, models, options)
    step_rebuilt = ArchimedLight.run_light_step(cache2, rows[1])
    @test !HelperModule._budgets_close(step_cached.budget, step_rebuilt.budget; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case light_cache_extra_band_parity" tags = [:synthetic, :fast, :light_cache_extra_band_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    cached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    uncached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false)
    row0 = merge(HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0), (RI_UV_F=25.0,))
    row1 = merge(HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=600.0, ri_par_f=100.0, ri_nir_f=60.0), (RI_UV_F=10.0,))
    meteo = ArchimedLight.MeteoTable([row0, row1], (; source="synthetic_extra_band_cache"))

    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options)

    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test HelperModule._budgets_close(cached[i].budget, uncached[i].budget; atol=1e-9, rtol=1e-9)
    end
end

@testitem "Synthetic case light_cache_partial_lru" tags = [:synthetic, :fast, :light_cache_partial_lru] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=true)
    probe = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=10^9)

    row_a = HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(9), sun_azimut=120.0, sun_elevation=30.0)
    row_b = HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(15), sun_azimut=240.0, sun_elevation=35.0)

    summary0 = ArchimedLight.cache_summary(probe)
    @test summary0.mode == :partial

    step_probe = ArchimedLight.run_light_step(probe, row_a)
    probe_bytes = ArchimedLight.cache_summary(probe).resident_bytes
    @test probe_bytes > 0

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=probe_bytes + max(div(probe_bytes, 10), 1))

    step_a = ArchimedLight.run_light_step(cache, row_a)
    step_b = ArchimedLight.run_light_step(cache, row_b)
    uncached_a = ArchimedLight.run_light_step(scene, models, row_a, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false))
    uncached_b = ArchimedLight.run_light_step(scene, models, row_b, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false))

    @test HelperModule._budgets_close(step_a.budget, uncached_a.budget; atol=1e-9, rtol=1e-9)
    @test HelperModule._budgets_close(step_b.budget, uncached_b.budget; atol=1e-9, rtol=1e-9)

    summary = ArchimedLight.cache_summary(cache)
    @test summary.cached_turtle_count <= 1
    @test summary.resident_bytes <= cache.memory_limit_bytes
end

@testitem "Synthetic case light_cache_topology_fallback" tags = [:synthetic, :fast, :light_cache_topology_fallback] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    row = merge(HelperModule._synthetic_meteo_row(), (RI_UV_F=20.0,))

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=1)
    summary = ArchimedLight.cache_summary(cache)
    @test summary.mode == :topology_fallback

    cached = ArchimedLight.run_light_step(cache, row)
    uncached = ArchimedLight.run_light_step(scene, models, row, HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))
    @test HelperModule._budgets_close(cached.budget, uncached.budget; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case missing_models" tags = [:synthetic, :fast, :missing_models] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options()
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    @test_throws ErrorException ArchimedLight.compute_first_order(scene, ArchimedLight.LightModels(), turtle, fluxes, options)
end
