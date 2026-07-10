@testitem "Radiation input semantics preserve or daylight-average supplied energy" tags=[:core, :fast, :radiation_semantics] begin
    using Dates

    function meteo_row(hour_start, hour_end, duration_seconds)
        (
            date=Date(2020, 3, 20),
            hour_start=hour_start,
            hour_end=hour_end,
            step_duration=duration_seconds,
            latitude=0.0,
            RI_PAR_f=100.0,
            RI_NIR_f=50.0,
            RI_custom_f=40.0,
        )
    end

    base = ArchimedLight.LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        scattering=false,
        radiation_timestep_minutes=5.0,
    )
    interval_options = ArchimedLight.LightOptions(base; radiation_input_semantics=:interval_mean)
    sunlit_options = ArchimedLight.LightOptions(base; radiation_input_semantics=:sunlit_intensity)

    partial_row = meteo_row("05:00:00", "07:00:00", 7200.0)
    rise, set = ArchimedLight._java_sunrise_sunset_hours(0.0, dayofyear(partial_row.date))
    expected_daylight_fraction = clamp((min(7.0, set) - max(5.0, rise)) / 2.0, 0.0, 1.0)
    @test 0.0 < expected_daylight_fraction < 1.0

    interval_sky = ArchimedLight.compute_sky(partial_row, interval_options)
    sunlit_sky = ArchimedLight.compute_sky(partial_row, sunlit_options)
    @test interval_sky.ri_par_f == 100.0
    @test interval_sky.ri_nir_f == 50.0
    @test interval_sky.ri_sw_f == 150.0
    @test isapprox(sunlit_sky.ri_par_f, 100.0 * expected_daylight_fraction; atol=1e-10, rtol=1e-10)
    @test isapprox(sunlit_sky.ri_nir_f, 50.0 * expected_daylight_fraction; atol=1e-10, rtol=1e-10)
    @test isapprox(sunlit_sky.ri_sw_f, 150.0 * expected_daylight_fraction; atol=1e-10, rtol=1e-10)
    @test ArchimedLight._radiation_input_effective_scale(partial_row, interval_options) == 1.0
    @test isapprox(
        ArchimedLight._radiation_input_effective_scale(partial_row, sunlit_options),
        expected_daylight_fraction;
        atol=1e-10,
        rtol=1e-10,
    )

    interval_turtle = ArchimedLight.build_turtle(interval_options, interval_sky)
    sunlit_turtle = ArchimedLight.build_turtle(sunlit_options, sunlit_sky)
    interval_flux = ArchimedLight.compute_directional_fluxes(partial_row, interval_sky, interval_turtle, interval_options)
    sunlit_flux = ArchimedLight.compute_directional_fluxes(partial_row, sunlit_sky, sunlit_turtle, sunlit_options)
    @test isapprox(sum(interval_flux.par), interval_sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(interval_flux.nir), interval_sky.ri_nir_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(sunlit_flux.par), sunlit_sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(sunlit_flux.nir), sunlit_sky.ri_nir_f; atol=1e-9, rtol=1e-9)

    interval_custom = ArchimedLight._extra_band_irradiance(partial_row, interval_options)["CUSTOM"]
    sunlit_custom = ArchimedLight._extra_band_irradiance(partial_row, sunlit_options)["CUSTOM"]
    @test interval_custom == 40.0
    @test isapprox(sunlit_custom, 40.0 * expected_daylight_fraction; atol=1e-10, rtol=1e-10)
    interval_custom_flux = ArchimedLight._single_band_flux(
        interval_custom,
        partial_row,
        interval_sky,
        interval_turtle,
        interval_options,
    )
    sunlit_custom_flux = ArchimedLight._single_band_flux(
        sunlit_custom,
        partial_row,
        sunlit_sky,
        sunlit_turtle,
        sunlit_options,
    )
    @test isapprox(sum(interval_custom_flux), interval_custom; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(sunlit_custom_flux), sunlit_custom; atol=1e-9, rtol=1e-9)

    clearness_only_row = (
        date=partial_row.date,
        hour_start=partial_row.hour_start,
        hour_end=partial_row.hour_end,
        step_duration=partial_row.step_duration,
        latitude=partial_row.latitude,
        clearness=0.6,
        RI_custom_f=40.0,
    )
    for (semantics, expected_custom) in (
        (:interval_mean, 40.0),
        (:sunlit_intensity, 40.0 * expected_daylight_fraction),
    ), cache_radiation in (false, true)
        options = ArchimedLight.LightOptions(
            base;
            radiation_input_semantics=semantics,
            cache_radiation=cache_radiation,
        )
        sky = ArchimedLight.compute_sky(clearness_only_row, options)
        turtle = ArchimedLight.build_turtle(options, sky)
        custom = ArchimedLight._extra_band_irradiance(clearness_only_row, options)["CUSTOM"]
        custom_flux = ArchimedLight._single_band_flux(custom, clearness_only_row, sky, turtle, options)
        @test isapprox(custom, expected_custom; atol=1e-10, rtol=1e-10)
        @test isapprox(sum(custom_flux), custom; atol=1e-9, rtol=1e-9)
    end
    @test isapprox(sum(interval_flux.par) * 7200.0, 100.0 * 7200.0; atol=1e-6, rtol=1e-10)
    @test isapprox(
        sum(sunlit_flux.par) * 7200.0,
        100.0 * expected_daylight_fraction * 7200.0;
        atol=1e-6,
        rtol=1e-10,
    )

    night_row = meteo_row("00:00:00", "01:00:00", 3600.0)
    night_sunlit_sky = ArchimedLight.compute_sky(night_row, sunlit_options)
    @test night_sunlit_sky.ri_sw_f == 0.0
    @test night_sunlit_sky.ri_par_f == 0.0
    @test night_sunlit_sky.ri_nir_f == 0.0
    @test ArchimedLight._radiation_input_effective_scale(night_row, sunlit_options) == 0.0
    @test ArchimedLight._extra_band_irradiance(night_row, sunlit_options)["CUSTOM"] == 0.0
    night_sunlit_turtle = ArchimedLight.build_turtle(sunlit_options, night_sunlit_sky)
    night_sunlit_flux = ArchimedLight.compute_directional_fluxes(night_row, night_sunlit_sky, night_sunlit_turtle, sunlit_options)
    @test all(iszero, night_sunlit_flux.par)
    @test all(iszero, night_sunlit_flux.nir)
    @test all(
        iszero,
        ArchimedLight._single_band_flux(0.0, night_row, night_sunlit_sky, night_sunlit_turtle, sunlit_options),
    )

    night_interval_sky = ArchimedLight.compute_sky(night_row, interval_options)
    @test night_interval_sky.ri_par_f == 100.0
    @test night_interval_sky.ri_nir_f == 50.0
    @test ArchimedLight._extra_band_irradiance(night_row, interval_options)["CUSTOM"] == 40.0
    night_interval_turtle = ArchimedLight.build_turtle(interval_options, night_interval_sky)
    night_interval_flux = ArchimedLight.compute_directional_fluxes(night_row, night_interval_sky, night_interval_turtle, interval_options)
    @test isapprox(sum(night_interval_flux.par), 100.0; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(night_interval_flux.nir), 50.0; atol=1e-10, rtol=1e-10)
end

@testitem "Interval-mean direct energy survives an unrepresentable sun" tags=[:core, :fast, :radiation_semantics] begin
    using Dates

    options = ArchimedLight.LightOptions(
        turtle_sectors=46,
        all_in_turtle=true,
        scattering=false,
        radiation_timestep_minutes=5.0,
        radiation_input_semantics=:interval_mean,
    )

    night_row = (
        date=Date(2020, 3, 20),
        hour_start="00:00:00",
        hour_end="01:00:00",
        step_duration=3600.0,
        latitude=0.0,
        RI_PAR_f=100.0,
        RI_NIR_f=50.0,
        direct_fraction=1.0,
    )
    night_sky = ArchimedLight.compute_sky(night_row, options)
    @test night_sky.sun_elevation_deg <= 0.0
    @test night_sky.direct_fraction == 1.0
    night_turtle = ArchimedLight.build_turtle(options, night_sky)
    @test all(sector -> sector.source == :sky, night_turtle.sectors)
    night_flux = ArchimedLight.compute_directional_fluxes(night_row, night_sky, night_turtle, options)
    @test isapprox(sum(night_flux.par), night_sky.ri_par_f; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(night_flux.nir), night_sky.ri_nir_f; atol=1e-10, rtol=1e-10)

    dawn_row = (
        date=Date(2020, 3, 20),
        hour_start="05:00:00",
        hour_end="07:00:00",
        step_duration=7200.0,
        latitude=0.0,
        RI_PAR_f=100.0,
        RI_NIR_f=50.0,
        direct_fraction=1.0,
        sun_azimut=90.0,
        sun_elevation=-12.0,
    )
    dawn_sky = ArchimedLight.compute_sky(dawn_row, options)
    @test dawn_sky.sun_elevation_deg == -12.0
    dawn_turtle = ArchimedLight.build_turtle(options, dawn_sky)
    @test all(sector -> sector.source == :sky, dawn_turtle.sectors)
    dawn_flux = ArchimedLight.compute_directional_fluxes(dawn_row, dawn_sky, dawn_turtle, options)
    @test isapprox(sum(dawn_flux.par), dawn_sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(dawn_flux.nir), dawn_sky.ri_nir_f; atol=1e-9, rtol=1e-9)
end

@testitem "Explicit sun geometry and direct fraction propagate through substeps" tags=[:core, :fast, :radiation_semantics] begin
    using Dates

    row = (
        date=Date(2020, 6, 21),
        hour_start="12:00:00",
        hour_end="13:00:00",
        step_duration=3600.0,
        latitude=0.0,
        RI_PAR_f=100.0,
        RI_NIR_f=50.0,
        direct_fraction=1.0,
        sun_azimut=135.0,
        sun_elevation=35.0,
    )
    options = ArchimedLight.LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        scattering=false,
        radiation_input_semantics=:interval_mean,
    )

    sky = ArchimedLight.compute_sky(row, options)
    @test sky.sun_azimuth_deg == 135.0
    @test sky.sun_elevation_deg == 35.0
    @test sky.direct_fraction == 1.0
    @test sky.diffuse_fraction == 0.0

    turtle = ArchimedLight.build_turtle(options, sky)
    flux = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    sun_index = only(findall(sector -> sector.source == :sun, turtle.sectors))
    sky_indices = findall(sector -> sector.source == :sky, turtle.sectors)
    @test isapprox(flux.par[sun_index], sky.ri_par_f; atol=1e-10, rtol=1e-10)
    @test isapprox(flux.nir[sun_index], sky.ri_nir_f; atol=1e-10, rtol=1e-10)
    @test all(iszero, flux.par[sky_indices])
    @test all(iszero, flux.nir[sky_indices])
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-10, rtol=1e-10)
end

@testitem "Sunlit semantics honor explicit sun without interval metadata" tags=[:core, :fast, :radiation_semantics] begin
    options = ArchimedLight.LightOptions(
        radiation_input_semantics=:sunlit_intensity,
        all_in_turtle=false,
    )
    night_row = (
        RI_PAR_f=100.0,
        RI_NIR_f=50.0,
        sun_azimut=180.0,
        sun_elevation=-12.0,
        direct_fraction=1.0,
    )
    day_row = merge(night_row, (sun_elevation=35.0,))

    night = ArchimedLight.compute_sky(night_row, options)
    day = ArchimedLight.compute_sky(day_row, options)
    @test night.ri_sw_f == 0.0
    @test night.ri_par_f == 0.0
    @test night.ri_nir_f == 0.0
    @test day.ri_sw_f == 150.0
    @test day.ri_par_f == 100.0
    @test day.ri_nir_f == 50.0
end

@testitem "LightOptions validates radiation semantics and parses rotation" tags=[:core, :fast, :radiation_semantics] begin
    @test_throws ArgumentError ArchimedLight.LightOptions(radiation_input_semantics=:instantaneous)

    options = ArchimedLight.LightOptions(
        radiation_input_semantics=:sunlit_intensity,
        scene_rotation_deg=37.5,
    )
    copied = ArchimedLight.LightOptions(options; pixel_size=0.01)
    @test copied.radiation_input_semantics == :sunlit_intensity
    @test copied.scene_rotation_deg == 37.5

    legacy_names = filter(
        name -> name ∉ (:radiation_input_semantics, :scene_rotation_deg),
        fieldnames(ArchimedLight.LightOptions),
    )
    legacy = ArchimedLight.LightOptions(
        (getfield(ArchimedLight.LightOptions(), name) for name in legacy_names)...,
    )
    @test legacy.radiation_input_semantics == :interval_mean
    @test legacy.scene_rotation_deg == 0.0

    mktempdir() do dir
        path = joinpath(dir, "config.yml")
        write(path, "radiation_input_semantics: sunlit_intensity\nscene_rotation: 37.5\n")
        parsed = ArchimedLight.read_options(path)
        @test parsed.radiation_input_semantics == :sunlit_intensity
        @test parsed.scene_rotation_deg == 37.5
    end
end

@testitem "Scene rotation transforms the sun in a fixed local turtle" tags=[:core, :fast, :scene_rotation] begin
    azimuth(direction) = mod(rad2deg(atan(direction[1], direction[2])), 360.0)

    sky = ArchimedLight.SkyState(0.0, 30.0, 100.0, 50.0, 0.7, 0.3)
    options0 = ArchimedLight.LightOptions(turtle_sectors=6, all_in_turtle=false, scene_rotation_deg=0.0)
    options90 = ArchimedLight.LightOptions(options0; scene_rotation_deg=90.0)
    options360 = ArchimedLight.LightOptions(options0; scene_rotation_deg=360.0)

    turtle0 = ArchimedLight.build_turtle(options0, sky)
    original_directions = [copy(sector.direction) for sector in turtle0.sectors]
    turtle90 = ArchimedLight.build_turtle(options90, sky)
    turtle360 = ArchimedLight.build_turtle(options360, sky)

    @test [sector.direction for sector in turtle0.sectors] == original_directions
    @test length(turtle0.sectors) == length(turtle90.sectors) == length(turtle360.sectors)
    for i in eachindex(turtle0.sectors)
        if turtle0.sectors[i].source == :sky
            @test turtle90.sectors[i].direction == turtle0.sectors[i].direction
        else
            expected = ArchimedLight._scene_local_direction(turtle0.sectors[i].direction, options90)
            @test isapprox(turtle90.sectors[i].direction, expected; atol=1e-12, rtol=1e-12)
        end
        @test isapprox(turtle360.sectors[i].direction, turtle0.sectors[i].direction; atol=1e-12, rtol=1e-12)
    end

    sun0 = only(sector.direction for sector in turtle0.sectors if sector.source == :sun)
    sun90 = only(sector.direction for sector in turtle90.sectors if sector.source == :sun)
    @test isapprox(azimuth(sun0), 0.0; atol=1e-12)
    @test isapprox(azimuth(sun90), 90.0; atol=1e-12)

    sector_index = findfirst(sector -> sector.source == :sky && hypot(sector.direction[1], sector.direction[2]) > 0.1, turtle0.sectors)
    @test sector_index !== nothing
    sector_az0 = azimuth(turtle0.sectors[sector_index].direction)
    sector_az90 = azimuth(turtle90.sectors[sector_index].direction)
    @test sector_az90 == sector_az0

    all0 = ArchimedLight.LightOptions(options0; all_in_turtle=true)
    all90 = ArchimedLight.LightOptions(options90; all_in_turtle=true)
    all_turtle0 = ArchimedLight.build_turtle(all0, sky)
    all_turtle90 = ArchimedLight.build_turtle(all90, sky)
    flux0 = ArchimedLight.compute_directional_fluxes(sky, all_turtle0, all0)
    flux90 = ArchimedLight.compute_directional_fluxes(sky, all_turtle90, all90)
    @test !isapprox(flux90.par, flux0.par; atol=1e-10, rtol=1e-10)
    @test !isapprox(flux90.nir, flux0.nir; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(flux90.par), sum(flux0.par); atol=1e-10, rtol=1e-10)
    @test isapprox(sum(flux90.nir), sum(flux0.nir); atol=1e-10, rtol=1e-10)
    @test [sector.direction for sector in turtle0.sectors] == original_directions
end
