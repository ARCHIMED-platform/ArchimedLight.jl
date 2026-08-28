@testitem "Meteo resolver handles sentinels, aliases, and derivations" tags=[:meteo, :fast] begin
    using Dates

    options = ArchimedLight.LightOptions()
    base = (
        date=Date(2020, 6, 21),
        hour_start="12:00:00",
        hour_end="13:00:00",
        latitude=15.0,
    )

    mixed_aliases = merge(
        base,
        (
            Ri_SW_f=Inf,
            Ri_PAR_f=Inf,
            RI_PAR_f=120.0,
            Ri_NIR_f=Inf,
            RI_NIR_f=80.0,
            clearness=Inf,
        ),
    )
    mixed = ArchimedLight._resolve_meteo_step(mixed_aliases, options)
    @test isempty(mixed.errors)
    @test mixed.step !== nothing
    @test mixed.step.ri_sw_f == 200.0
    @test mixed.step.ri_par_f == 120.0
    @test mixed.step.ri_nir_f == 80.0
    @test mixed.step.sources[:ri_par_f] == :RI_PAR_f

    sky = ArchimedLight.compute_sky(mixed_aliases, options)
    @test sky.ri_sw_f == 200.0
    @test sky.ri_par_f == 120.0
    @test sky.ri_nir_f == 80.0

    fallback = ArchimedLight._resolve_meteo_step(
        merge(base, (Ri_SW_f=200.0, Ri_PAR_f=NaN, Ri_NIR_f=Inf)),
        options,
    )
    @test isempty(fallback.errors)
    @test fallback.step.ri_par_f == 96.0
    @test fallback.step.ri_nir_f == 104.0
    @test fallback.step.sources[:ri_par_f] == :derived_from_sw

    alias_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (Ri_PAR_f=100.0, RI_PAR_f=120.0, Ri_NIR_f=80.0)),
        options,
    )
    @test any(contains("conflicting aliases"), alias_conflict.errors)

    energy_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (Ri_SW_f=100.0, Ri_PAR_f=80.0, Ri_NIR_f=50.0)),
        options,
    )
    @test any(contains("conflicting total shortwave estimates"), energy_conflict.errors)

    missing = ArchimedLight._resolve_meteo_step(base, options)
    @test any(contains("cannot resolve incident radiation"), missing.errors)
end

@testitem "Resolved meteo drives partition and optional NIR paths" tags=[:meteo, :fast] begin
    explicit = (
        step_duration=3600.0,
        Ri_SW_f=200.0,
        Ri_PAR_f=96.0,
        Ri_NIR_f=104.0,
        sun_azimuth=180.0,
        sun_elevation=45.0,
        direct_fraction=1.0,
    )
    options = ArchimedLight.LightOptions(turtle_sectors=6, all_in_turtle=false)
    report = ArchimedLight.check_meteo(explicit; options=options)
    @test isempty(report.errors)

    resolved = ArchimedLight._resolved_meteo_step_or_error(explicit, options)
    sky = ArchimedLight.compute_sky(explicit, options; resolved_step=resolved)
    turtle = ArchimedLight.build_turtle(options, sky)
    flux = ArchimedLight.compute_directional_fluxes(
        explicit,
        sky,
        turtle,
        options;
        resolved_step=resolved,
    )
    sun_index = only(findall(sector -> sector.source == :sun, turtle.sectors))
    sky_indices = findall(sector -> sector.source == :sky, turtle.sectors)
    @test sky.direct_fraction == 1.0
    @test isapprox(flux.par[sun_index], sky.ri_par_f; atol=1e-10, rtol=1e-10)
    @test all(iszero, flux.par[sky_indices])

    small_angles = merge(explicit, (sun_azimuth=5.0, sun_elevation=5.0))
    small_angle_sky = ArchimedLight.compute_sky(small_angles, options)
    @test small_angle_sky.sun_azimuth_deg == 5.0
    @test small_angle_sky.sun_elevation_deg == 5.0

    components = merge(
        explicit,
        (
            Ri_SW_f=Inf,
            Ri_PAR_f=Inf,
            Ri_NIR_f=Inf,
            direct_fraction=Inf,
            Ri_SW_f_direct=40.0,
            Ri_SW_f_diffuse=60.0,
        ),
    )
    component_step = ArchimedLight._resolved_meteo_step_or_error(components, options)
    @test component_step.ri_sw_f == 100.0
    @test component_step.direct_fraction == 0.4
    component_summary = ArchimedLight.summarize_meteo(components; options=options)
    @test Set(component_summary.radiation_inputs) ==
          Set(["RI_SW_f_direct", "RI_SW_f_diffuse"])

    redundant_summary = ArchimedLight.summarize_meteo(explicit; options=options)
    @test Set(["RI_SW_f", "RI_PAR_f", "RI_NIR_f", "direct_fraction"]) <=
          Set(redundant_summary.radiation_inputs)

    no_nir_options = ArchimedLight.LightOptions(options; nir_interception=false)
    nir_only = merge(
        explicit,
        (Ri_SW_f=Inf, Ri_PAR_f=Inf, Ri_NIR_f=52.0, direct_fraction=0.5),
    )
    nir_disabled = ArchimedLight._resolved_meteo_step_or_error(nir_only, no_nir_options)
    @test nir_disabled.ri_sw_f == 100.0
    @test nir_disabled.ri_par_f == 48.0
    @test nir_disabled.ri_nir_f == 0.0

    unused_nir = merge(
        explicit,
        (Ri_SW_f=100.0, Ri_PAR_f=48.0, Ri_NIR_f=900.0, direct_fraction=0.5),
    )
    unused_nir_step =
        ArchimedLight._resolved_meteo_step_or_error(unused_nir, no_nir_options)
    @test unused_nir_step.ri_sw_f == 100.0
    @test unused_nir_step.ri_par_f == 48.0
    @test unused_nir_step.ri_nir_f == 0.0

    measured_bands = merge(
        explicit,
        (Ri_SW_f=Inf, Ri_PAR_f=80.0, Ri_NIR_f=20.0, direct_fraction=0.5),
    )
    measured_bands_step =
        ArchimedLight._resolved_meteo_step_or_error(measured_bands, no_nir_options)
    @test measured_bands_step.ri_sw_f == 100.0
    @test measured_bands_step.ri_par_f == 80.0
    @test measured_bands_step.ri_nir_f == 0.0

    thermal = merge(explicit, (Ri_TIR_f=450.0, RI_UV_f=12.0))
    thermal_step = ArchimedLight._resolved_meteo_step_or_error(thermal, options)
    @test !haskey(thermal_step.extra_irradiance, "TIR")
    @test thermal_step.extra_irradiance["UV"] == 12.0
    @test !haskey(
        ArchimedLight._extra_band_irradiance(thermal_step, thermal, options),
        "TIR",
    )

    custom_rows = [merge(explicit, (Ri_UV_f=Inf,))]
    custom_table = ArchimedLight.PlantMeteo.TimeStepTable(
        custom_rows,
        (RI_UV_f=12.0,),
    )
    metadata_custom = ArchimedLight._resolved_meteo_step_or_error(
        first(custom_table),
        options,
    )
    @test metadata_custom.extra_irradiance["UV"] == 12.0

    row_custom_table = ArchimedLight.PlantMeteo.TimeStepTable(
        [merge(explicit, (Ri_UV_f=15.0,))],
        (RI_UV_f=12.0,),
    )
    row_custom = ArchimedLight._resolved_meteo_step_or_error(
        first(row_custom_table),
        options,
    )
    @test row_custom.extra_irradiance["UV"] == 15.0
end

@testitem "Meteo resolver uses TimeStepTable metadata fallback" tags=[:meteo, :fast] begin
    using Dates

    rows = [(
        date=Date(2020, 6, 21),
        hour_start="12:00:00",
        hour_end="13:00:00",
        Ri_SW_f=Inf,
    )]
    table = ArchimedLight.PlantMeteo.TimeStepTable(
        rows,
        (latitude=15.0, Ri_SW_f=200.0),
    )
    result = ArchimedLight._resolve_meteo_step(first(table), ArchimedLight.LightOptions())
    @test isempty(result.errors)
    @test result.step.ri_sw_f == 200.0
    @test result.step.ri_par_f == 96.0
    @test result.step.sources[:ri_sw_f] == :Ri_SW_f
end

@testitem "Historical meteo use metadata is provenance only" tags=[:meteo, :fast] begin
    mktempdir() do dir
        path = joinpath(dir, "meteo.csv")
        write(
            path,
            """
            #' latitude: 15.0
            #' use: clearness, obsolete_source
            date;hour_start;hour_end;T;Wind;Rh;Ri_PAR_f;Ri_NIR_f
            2020/06/21;12:00:00;13:00:00;25.0;1.0;0.6;120.0;80.0
            """,
        )

        meteo = ArchimedLight.read_meteo(path)
        row = first(meteo)

        @test ArchimedLight.PlantMeteo.metadata(meteo, :use) ==
              "clearness, obsolete_source"
        @test isempty(ArchimedLight.check_meteo(meteo).errors)

        summary = ArchimedLight.summarize_meteo(meteo)
        @test Set(summary.radiation_inputs) == Set(["RI_PAR_f", "RI_NIR_f"])

        resolved = ArchimedLight._resolved_meteo_step_or_error(
            row,
            ArchimedLight.LightOptions(),
        )
        @test resolved.radiation_source == :spectral_sum
        @test resolved.ri_sw_f == 200.0
        @test resolved.ri_par_f == 120.0
        @test resolved.ri_nir_f == 80.0
    end
end

@testitem "Meteo file fallback is limited to incomplete Atmosphere schemas" tags=[:meteo, :fast] begin
    mktempdir() do dir
        radiation_only_path = joinpath(dir, "radiation-only.csv")
        write(
            radiation_only_path,
            """
            #' latitude: 15.0
            date;hour_start;hour_end;Ri_PAR_f;Ri_NIR_f
            2020/06/21;12:00:00;13:00:00;120.0;80.0
            """,
        )

        radiation_only = ArchimedLight.read_meteo(radiation_only_path)
        @test length(radiation_only) == 1
        @test first(radiation_only).Ri_PAR_f == 120.0
        @test first(radiation_only).Ri_NIR_f == 80.0

        invalid_humidity_path = joinpath(dir, "invalid-humidity.csv")
        write(
            invalid_humidity_path,
            """
            date;hour_start;hour_end;T;Wind;Rh;Ri_SW_f
            2020/06/21;12:00:00;13:00:00;25.0;1.0;60.0;200.0
            """,
        )
        humidity_error = try
            ArchimedLight.read_meteo(invalid_humidity_path)
            nothing
        catch err
            err
        end
        @test humidity_error isa ArgumentError
        @test occursin("Relative humidity", sprint(showerror, humidity_error))

        invalid_radiation_path = joinpath(dir, "invalid-radiation.csv")
        write(
            invalid_radiation_path,
            """
            date;hour_start;hour_end;T;Wind;Rh;Ri_SW_f
            2020/06/21;12:00:00;13:00:00;25.0;1.0;0.6;-1.0
            """,
        )
        radiation_error = try
            ArchimedLight.read_meteo(invalid_radiation_path)
            nothing
        catch err
            err
        end
        @test radiation_error isa ArgumentError
        @test occursin("Ri_SW_f", sprint(showerror, radiation_error))
        @test occursin("non-negative", sprint(showerror, radiation_error))

        invalid_date_path = joinpath(dir, "invalid-date.csv")
        write(
            invalid_date_path,
            """
            date;hour_start;hour_end;T;Wind;Rh;Ri_SW_f
            not-a-date;12:00:00;13:00:00;25.0;1.0;0.6;200.0
            """,
        )
        date_error = try
            ArchimedLight.read_meteo(invalid_date_path)
            nothing
        catch err
            err
        end
        @test date_error isa Exception
        @test occursin("date", lowercase(sprint(showerror, date_error)))
        @test occursin("cannot be parsed", sprint(showerror, date_error))
    end
end

@testitem "Meteo boundaries are configurable but conflicts remain mandatory" tags=[:meteo, :fast] begin
    using Dates

    row = (
        date=Date(2020, 6, 21),
        hour_start="12:00:00",
        hour_end="13:00:00",
        latitude=15.0,
        Ri_SW_f=100.0,
        Ri_PAR_f=-5.0,
        Ri_NIR_f=105.0,
    )
    options = ArchimedLight.LightOptions()

    strict = ArchimedLight.check_meteo(row; options=options)
    @test any(contains("must be non-negative"), strict.errors)

    permissive = ArchimedLight.check_meteo(
        row;
        options=options,
        check_boundaries=false,
    )
    @test isempty(permissive.errors)

    table = ArchimedLight.PlantMeteo.TimeStepTable([row])
    @test_throws ErrorException ArchimedLight.prepare_meteo(table, options)
    @test length(
        ArchimedLight.prepare_meteo(
            table,
            options;
            check_boundaries=false,
        ),
    ) == 1

    conflict = merge(row, (Ri_SW_f=90.0,))
    conflict_report = ArchimedLight.check_meteo(
        conflict;
        options=options,
        check_boundaries=false,
    )
    @test any(contains("conflicting total shortwave estimates"), conflict_report.errors)
end

@testitem "Meteo options parse and copy boundary policy" tags=[:meteo, :fast] begin
    options = ArchimedLight.LightOptions(check_meteo_boundaries=false)
    @test !options.check_meteo_boundaries
    @test !ArchimedLight.LightOptions(options; pixel_size=0.01).check_meteo_boundaries

    mktempdir() do dir
        path = joinpath(dir, "config.yml")
        write(path, "check_meteo_boundaries: false\n")
        @test !ArchimedLight.read_options(path).check_meteo_boundaries
    end


    fixture_config = joinpath(
        @__DIR__,
        "fast_fixtures",
        "sky_06_direct",
        "input",
        "config.yml",
    )
    sim, _ = ArchimedLight.read_simulation(fixture_config; check_boundaries=false)
    @test !sim.options.check_meteo_boundaries
end

@testitem "Meteo diagnostics preserve source rows after selection" tags=[:meteo, :fast] begin
    using Dates

    rows = [
        (
            date=Date(2020, 6, 21),
            hour_start="10:00:00",
            hour_end="11:00:00",
            latitude=15.0,
            Ri_SW_f=100.0,
            active=true,
        ),
        (
            date=Date(2020, 6, 21),
            hour_start="11:00:00",
            hour_end="12:00:00",
            latitude=15.0,
            Ri_SW_f=Inf,
            active=false,
        ),
        (
            date=Date(2020, 6, 21),
            hour_start="12:00:00",
            hour_end="13:00:00",
            latitude=15.0,
            Ri_SW_f=Inf,
            active=true,
        ),
    ]
    options = ArchimedLight.LightOptions(meteo_range="2, 3")
    report = ArchimedLight.check_meteo(rows; options=options)
    @test any(contains("Meteo row 3"), report.errors)
    @test !any(contains("Meteo row 2"), report.errors)

    bad_active_rows = [rows[1], merge(rows[2], (active="maybe",)), rows[3]]
    bad_active_report = ArchimedLight.check_meteo(bad_active_rows; options=options)
    @test any(contains("active value at meteo row 2"), bad_active_report.errors)

    datetime_range_rows = [
        merge(rows[1], (hour_start="bad", active=true)),
        merge(rows[3], (active=true,)),
    ]
    datetime_options = ArchimedLight.LightOptions(
        meteo_range="2020-06-21-12:00,2020-06-21-13:00",
    )
    datetime_report = ArchimedLight.check_meteo(
        datetime_range_rows;
        options=datetime_options,
    )
    @test any(contains("datetime interval at meteo row 1"), datetime_report.errors)
end
