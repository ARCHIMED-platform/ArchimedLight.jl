@testitem "Java light-source fixture uses physical Lambertian accounting" tags = [:emitter, :java_fixture, :slow] begin
    fixture_root = joinpath(
        @__DIR__,
        "..",
        "java_implementation",
        "archimed-lib-2018",
        "tests",
        "test-lightsource1",
    )
    config = joinpath(fixture_root, "config2.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config)
    row = first(ArchimedLight.prepare_meteo(meteo, options))
    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)

    emitted_par, emitted_nir = ArchimedLight._emitter_power_per_node(scene, models)
    @test length(emitted_par) == 1
    @test Set(keys(emitted_par)) == Set(keys(emitted_nir))
    source = only(keys(emitted_par))
    source_area = ArchimedLight._scene_area(scene, source, 0.0)

    # The historical fixture config deliberately uses gamma PAR=0.2 and
    # NIR=0.5. These are regression coefficients, not a pair to normalize.
    @test isapprox(emitted_par[source], pi * source_area * 100.0 * 0.2; atol=1e-12, rtol=1e-12)
    @test isapprox(emitted_nir[source], pi * source_area * 100.0 * 0.5; atol=1e-12, rtol=1e-12)
    @test isapprox(emitted_par[source] / emitted_nir[source], 0.4; atol=1e-14, rtol=1e-14)

    lamp_first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    received_par = sum(values(lamp_first.incident_power.par))
    received_nir = sum(values(lamp_first.incident_power.nir))
    escaped_par = sum(values(lamp_first.emitter_escaped_power.par))
    escaped_nir = sum(values(lamp_first.emitter_escaped_power.nir))

    # The fixture has no solar input and uses a toric receiver below the source.
    # It therefore provides a stable Java-scene topology regression while the
    # source magnitude follows the physical pi*A*L*gamma convention.
    @test isapprox(received_par + escaped_par, emitted_par[source]; atol=1e-8, rtol=1e-8)
    @test isapprox(received_nir + escaped_nir, emitted_nir[source]; atol=1e-8, rtol=1e-8)
    @test isapprox(escaped_par, 0.0; atol=1e-8)
    @test isapprox(escaped_nir, 0.0; atol=1e-8)

    # Java treated `radiance=100` as 100 W of component-total power, giving
    # 20 W PAR and 50 W NIR. The physical area/radiance result differs by pi*A;
    # dividing that factor out preserves the fixture's numerical regression.
    @test isapprox(received_par / (pi * source_area), 20.0; atol=1e-8, rtol=1e-8)
    @test isapprox(received_nir / (pi * source_area), 50.0; atol=1e-8, rtol=1e-8)

    # The fixture README's transfer oracle compares the one-plate solar scene
    # with the two-plate lamp scene: after removing the source magnitude, the
    # lower plate and ground must receive the same geometric shares.
    sun_options, sun_scene, sun_meteo, sun_models =
        ArchimedLight.read_config(joinpath(fixture_root, "config1.yml"))
    sun_row = first(ArchimedLight.prepare_meteo(sun_meteo, sun_options))
    sun_sky = ArchimedLight.compute_sky(sun_row, sun_options)
    sun_turtle = ArchimedLight.build_turtle(sun_options, sun_sky)
    sun_fluxes = ArchimedLight.compute_directional_fluxes(
        sun_row,
        sun_sky,
        sun_turtle,
        sun_options,
    )
    sun_first = ArchimedLight.compute_first_order(
        sun_scene,
        sun_models,
        sun_turtle,
        sun_fluxes,
        sun_options,
    )

    function group_power(first_order, fixture_scene, group, band)
        values = band == "PAR" ? first_order.incident_power.par : first_order.incident_power.nir
        return sum(
            (
                get(values, nid, 0.0) for nid in keys(fixture_scene.nodes) if
                ArchimedLight._scene_group(fixture_scene, nid, "") == group
            );
            init=0.0,
        )
    end

    xmin, ymin, xmax, ymax = sun_scene.scene_xy_bounds
    plot_area = (xmax - xmin) * (ymax - ymin)
    source_power = Dict("PAR" => emitted_par[source], "NIR" => emitted_nir[source])
    sun_power = Dict(
        "PAR" => sun_sky.ri_par_f * plot_area,
        "NIR" => sun_sky.ri_nir_f * plot_area,
    )
    # The two independently rasterized directional fields differ at the shared
    # plate boundary by roughly one pixel. Keep the oracle per component while
    # allowing that finite pixel-area error.
    for band in ("PAR", "NIR"), group in ("plate1", "pavement")
        lamp_share = group_power(lamp_first, scene, group, band) / source_power[band]
        sun_share = group_power(sun_first, sun_scene, group, band) / sun_power[band]
        @test isapprox(lamp_share, sun_share; atol=1e-3, rtol=1e-3)
    end
end
