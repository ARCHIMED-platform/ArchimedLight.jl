@testitem "Corrected sunset hour-angle geometry" tags=[:sky, :solar_geometry, :fast] begin
    latitude = deg2rad(43.6)
    summer_declination = deg2rad(23.44)
    winter_declination = -summer_declination

    summer = ArchimedLight._java_sunset_hour_angle(latitude, summer_declination)
    winter = ArchimedLight._java_sunset_hour_angle(latitude, winter_declination)

    @test isapprox(summer, 2.0204829681571157; atol=1e-12, rtol=0.0)
    @test isapprox(winter, 1.1689930421877521; atol=1e-12, rtol=0.0)
    @test isapprox(
        ArchimedLight._java_sunset_hour_angle(-latitude, -summer_declination),
        summer;
        atol=1e-12,
        rtol=0.0,
    )

    legacy_x =
        ((sin(ArchimedLight._CENTER_OF_SOLAR_DISK_RAD) - sin(latitude) * sin(summer_declination)) /
         cos(latitude)) * cos(summer_declination)
    legacy = acos(clamp(legacy_x, -1.0, 1.0))
    @test !isapprox(summer, legacy; atol=1e-3, rtol=0.0)
end

@testitem "Polar and dawn forcing remain spectrally consistent" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    options = ArchimedLight.LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        radiation_timestep_minutes=1.2,
    )
    polar_day = (
        date=Date(2019, 6, 21),
        hour_start="00:00:00",
        step_duration=86400.0,
        latitude=68.0,
        clearness=0.6,
    )
    polar_night = merge(polar_day, (date=Date(2019, 12, 21),))

    day_sky = ArchimedLight.compute_sky(polar_day, options)
    day_turtle = ArchimedLight.build_turtle(options, day_sky)
    day_flux = ArchimedLight.compute_directional_fluxes(
        polar_day,
        day_sky,
        day_turtle,
        options,
    )
    @test isfinite(day_sky.ri_sw_f)
    @test day_sky.ri_sw_f > 0.0
    @test isapprox(day_sky.ri_par_f + day_sky.ri_nir_f, day_sky.ri_sw_f; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(day_flux.par), day_sky.ri_par_f; atol=1e-8, rtol=1e-8)

    night_sky = ArchimedLight.compute_sky(polar_night, options)
    night_turtle = ArchimedLight.build_turtle(options, night_sky)
    night_flux = ArchimedLight.compute_directional_fluxes(
        polar_night,
        night_sky,
        night_turtle,
        options,
    )
    @test night_sky.ri_sw_f == 0.0
    @test night_sky.ri_par_f == 0.0
    @test night_sky.ri_nir_f == 0.0
    @test all(iszero, night_flux.par)
    @test all(iszero, night_flux.nir)

    dawn = (
        date=Date(2016, 6, 12),
        hour_start="05:04:00",
        hour_end="05:34:00",
        step_duration=1800.0,
        latitude=15.0,
        clearness=0.75,
    )
    dawn_sky = ArchimedLight.compute_sky(dawn, options)
    @test dawn_sky.ri_sw_f > 0.0
    @test dawn_sky.ri_par_f > 0.0
    @test dawn_sky.ri_nir_f > 0.0
    @test isapprox(
        dawn_sky.ri_par_f + dawn_sky.ri_nir_f,
        dawn_sky.ri_sw_f;
        atol=1e-12,
        rtol=1e-12,
    )

    disabled_nir_options = ArchimedLight.LightOptions(options; nir_interception=false)
    numerical_dawn = merge(
        dawn,
        (
            hour_start=Time(5, 27, 57, 600),
            hour_end=Time(5, 33, 57, 600),
            step_duration=360.0,
        ),
    )
    disabled_nir_sky = ArchimedLight.compute_sky(numerical_dawn, disabled_nir_options)
    disabled_nir_turtle = ArchimedLight.build_turtle(disabled_nir_options, disabled_nir_sky)
    disabled_nir_flux = ArchimedLight.compute_directional_fluxes(
        numerical_dawn,
        disabled_nir_sky,
        disabled_nir_turtle,
        disabled_nir_options,
    )
    @test disabled_nir_sky.ri_sw_f > 0.0
    @test disabled_nir_sky.ri_par_f > 0.0
    @test disabled_nir_sky.ri_nir_f == 0.0
    @test all(iszero, disabled_nir_flux.nir)
end

@testitem "Polar sunrise and sunset regimes" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    arctic_latitude = deg2rad(68.0)
    solstice_declination = deg2rad(23.44)

    @test ArchimedLight._java_sunset_hour_angle(arctic_latitude, solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(arctic_latitude, -solstice_declination) == 0.0
    @test ArchimedLight._java_sunset_hour_angle(-arctic_latitude, -solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(-arctic_latitude, solstice_declination) == 0.0
    @test ArchimedLight._java_sunset_hour_angle(pi / 2, solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(pi / 2, -solstice_declination) == 0.0

    summer_date = Date(2019, 6, 21)
    winter_date = Date(2019, 12, 21)
    summer_doy = dayofyear(summer_date)
    winter_doy = dayofyear(winter_date)

    @test ArchimedLight._java_sunrise_sunset_hours(arctic_latitude, summer_doy) == (0.0, 24.0)
    @test ArchimedLight._java_sunrise_sunset_hours(arctic_latitude, winter_doy) == (0.0, 0.0)

    daylight_steps = ArchimedLight._java_substeps_v2(summer_date, 0.0, 24.0, 1.0, arctic_latitude)
    @test length(daylight_steps) == 24
    @test first(daylight_steps).start == 0.0
    @test isapprox(last(daylight_steps).stop, 24.0; atol=1e-12, rtol=0.0)
    @test isapprox(sum(step.duration for step in daylight_steps), 24.0; atol=1e-12, rtol=0.0)

    @test isempty(ArchimedLight._java_substeps_v2(winter_date, 0.0, 24.0, 1.0, arctic_latitude))
    @test ArchimedLight._java_extra_terrestrial_hourly_mj(arctic_latitude, winter_doy, 0.0, 24.0) == 0.0
    @test ArchimedLight._java_extra_terrestrial_hourly_mj(arctic_latitude, summer_doy, 0.0, 24.0) > 0.0
end

@testitem "Daylight clipping uses the corrected boundary" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    date = Date(2016, 6, 12)
    doy = dayofyear(date)
    latitude = deg2rad(15.0)
    rise, set = ArchimedLight._java_sunrise_sunset_hours(latitude, doy)

    @test isapprox(rise, 5.4937891731147355; atol=1e-12, rtol=0.0)
    @test isapprox(set, 18.499038173598567; atol=1e-12, rtol=0.0)

    steps = ArchimedLight._java_substeps_v2(date, 5.25, 5.75, 0.25, latitude)
    @test length(steps) == 2
    @test isapprox(first(steps).start, rise; atol=1e-12, rtol=0.0)
    @test isapprox(last(steps).stop, 5.75; atol=1e-12, rtol=0.0)
    @test isapprox(sum(step.duration for step in steps), 5.75 - rise; atol=1e-12, rtol=0.0)

    @test ArchimedLight._java_extra_terrestrial_hourly_mj(latitude, doy, 5.0, 5.5) == 0.0
    @test ArchimedLight._global_wm2_from_clearness(0.75, latitude, doy, 5.0, 5.5) == 0.0
end
