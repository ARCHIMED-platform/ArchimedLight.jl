@testmodule HelperModule begin
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    using GeometryBasics
    using Dates
    using ArchimedLight
    import MultiScaleTreeGraph
    import PlantGeom

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

@testitem "Synthetic case translucent first-order and stack visibility" tags = [:synthetic, :fast, :translucent_stack] setup = [HelperModule] begin
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
    @test ArchimedLight._use_upper_hit_pixel_table(models, options)
    run = HelperModule._run_direct(scene, sky, options; models=models)

    @test run.first.projected_area_per_node[1] < 1.0
    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(run.first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 0.0; atol=1e-8, rtol=1e-8)

    stack_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01)
    stack_run = HelperModule._run_direct(scene, sky, stack_options; models=models)

    @test isapprox(get(stack_run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_run.first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case Lambertian emitter energy accounting" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=ArchimedLight.OpticalProperties(0.2, 0.5),
                    ),
                ),
            ),
        ),
    ])

    # Unequal solid-angle sectors retain their relative Omega*cos(theta)
    # weights, while the corrected hemispherical quadrature integrates to pi.
    quadrature_turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(1, HelperModule.SVector(0.0, 0.0, 1.0), 0.25, :sky),
        ArchimedLight.TurtleSector(2, HelperModule.SVector(sqrt(0.75), 0.0, 0.5), 0.75, :sky),
    ])
    projected_solid_angles = ArchimedLight._lambertian_projected_solid_angles(quadrature_turtle)
    @test isapprox(sum(projected_solid_angles), pi; atol=1e-14, rtol=1e-14)
    @test isapprox(projected_solid_angles[1] / projected_solid_angles[2], 2 / 3; atol=1e-14, rtol=1e-14)

    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )

    open_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=2.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
    ])
    open_emitted_par, open_emitted_nir = ArchimedLight._emitter_power_per_node(open_scene, models)
    open_first = ArchimedLight.compute_first_order(open_scene, models, turtle, zero_fluxes, options)
    @test isapprox(open_emitted_par[1], 4pi; atol=1e-12, rtol=1e-12)
    @test isapprox(open_emitted_nir[1], 10pi; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(values(open_first.incident_power.par)), 0.0; atol=1e-12)
    @test isapprox(open_first.emitter_escaped_power.par[1], open_emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(open_first.emitter_escaped_power.nir[1], open_emitted_nir[1]; atol=1e-10, rtol=1e-10)

    closed_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    emitted_par, emitted_nir = ArchimedLight._emitter_power_per_node(closed_scene, models)
    first = ArchimedLight.compute_first_order(closed_scene, models, turtle, zero_fluxes, options)
    prepared = ArchimedLight._prepare_interception_data(closed_scene, models, options; include_budget_maps=true)
    n_nodes = length(prepared.geometry.node_ids)
    n_sectors = options.turtle_sectors + (options.all_in_turtle ? 0 : 1)
    expected_cache_estimate =
        n_sectors * n_nodes * (sizeof(Float64) + sizeof(Int)) +
        2 * length(prepared.emitter_nodes) * n_nodes * 64 +
        4 * n_nodes * sizeof(Float64)
    @test ArchimedLight._estimate_light_cache_entry_bytes(prepared, options) ==
          expected_cache_estimate
    responses = ArchimedLight._build_sector_responses(prepared, closed_scene, models, turtle, options)
    cached_first = ArchimedLight._combine_sector_responses(responses, zero_fluxes)
    received_par = sum(values(first.incident_power.par))
    received_nir = sum(values(first.incident_power.nir))
    escaped_par = sum(values(first.emitter_escaped_power.par))
    escaped_nir = sum(values(first.emitter_escaped_power.nir))

    @test isapprox(emitted_par[1], 2pi; atol=1e-12, rtol=1e-12)
    @test isapprox(emitted_nir[1], 5pi; atol=1e-12, rtol=1e-12)
    @test isapprox(received_par + escaped_par, emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(received_nir + escaped_nir, emitted_nir[1]; atol=1e-10, rtol=1e-10)
    @test 0.0 < first.emitter_escaped_power.par[1] < emitted_par[1]
    @test isapprox(first.incident_power.par[2] + first.emitter_escaped_power.par[1], emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(first.incident_power.nir[2] + first.emitter_escaped_power.nir[1], emitted_nir[1]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.incident_power.par, first.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.incident_power.nir, first.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.emitter_escaped_power.par, first.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.emitter_escaped_power.nir, first.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)

    stream_options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=false,
        scattering=true,
        pixel_size=0.01,
        toricity=true,
    )
    stream_turtle = ArchimedLight.build_turtle(stream_options, sky)
    @test any(sector -> sector.source == :sun && sector.weight == 0.0, stream_turtle.sectors)
    stream_zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in stream_turtle.sectors],
        zeros(length(stream_turtle.sectors)),
        zeros(length(stream_turtle.sectors)),
    )
    stream_prepared = ArchimedLight._prepare_interception_data(
        closed_scene,
        models,
        stream_options;
        include_budget_maps=true,
    )
    stream_direct = ArchimedLight.compute_first_order(
        closed_scene,
        models,
        stream_turtle,
        stream_zero_fluxes,
        stream_options,
    )
    stream_first, _ = ArchimedLight._stream_first_order_with_scattering_topology(
        stream_prepared,
        closed_scene,
        models,
        stream_turtle,
        stream_zero_fluxes,
        stream_options,
    )
    @test HelperModule._dicts_close(stream_first.incident_power.par, stream_direct.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.incident_power.nir, stream_direct.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.emitter_escaped_power.par, stream_direct.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.emitter_escaped_power.nir, stream_direct.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.incident_power.par, first.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.incident_power.nir, first.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.emitter_escaped_power.par, first.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.emitter_escaped_power.nir, first.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)

    occluded_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    occluded_emitted_par, _ = ArchimedLight._emitter_power_per_node(occluded_scene, models)
    occluded_first = ArchimedLight.compute_first_order(occluded_scene, models, turtle, zero_fluxes, options)
    occluded_prepared = ArchimedLight._prepare_interception_data(
        occluded_scene,
        models,
        options;
        include_budget_maps=true,
    )
    occluded_geometry = occluded_prepared.geometry
    occluded_transfer = ArchimedLight._emitter_transfer_weights(
        occluded_geometry.vertices,
        occluded_geometry.faces,
        occluded_geometry.face2node,
        turtle,
        options,
        occluded_geometry.plotbox,
        occluded_prepared.emitter_nodes,
        occluded_prepared.emitter_node_mask,
        occluded_prepared.virtual_nodes,
        occluded_prepared.virtual_node_mask,
        occluded_geometry.node_ids,
        occluded_prepared.cache_ctx,
    )

    # A different emitter is an opaque first hit rather than being skipped;
    # duplicate triangles of the originating component never become receivers.
    @test get(occluded_transfer.received_fraction, (2, 1), 0.0) > 0.0
    @test get(occluded_transfer.received_fraction, (3, 2), 0.0) > 0.0
    @test all(to != source for ((to, source), _) in occluded_transfer.received_fraction)
    for source in keys(occluded_emitted_par)
        received_fraction = sum(
            (fraction for ((_, from), fraction) in occluded_transfer.received_fraction if from == source);
            init=0.0,
        )
        @test isapprox(
            received_fraction + occluded_transfer.escaped_fraction_per_node[source],
            1.0;
            atol=1e-12,
            rtol=1e-12,
        )
    end
    @test isapprox(sum(values(occluded_first.incident_power.par)) + sum(values(occluded_first.emitter_escaped_power.par)), sum(values(occluded_emitted_par)); atol=1e-10, rtol=1e-10)
end

@testitem "Emitter custom bands and virtual sensors preserve physical transfer" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    custom_gamma = ArchimedLight.OpticalProperties(
        0.2,
        0.5,
        HelperModule.OrderedDict{String,Any}("CUSTOM" => 0.3),
    )
    receiver_optics = ArchimedLight.OpticalProperties(
        0.1,
        0.2,
        HelperModule.OrderedDict{String,Any}("CUSTOM" => 0.4),
    )
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        optical_properties=receiver_optics,
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        optical_properties=receiver_optics,
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=custom_gamma,
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "sensor";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(model="VirtualSensor"),
                ),
            ),
        ),
    ])

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=true,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )
    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_budget_maps=true)
    geometry = prepared.geometry
    transfer = ArchimedLight._emitter_transfer_weights(
        geometry.vertices,
        geometry.faces,
        geometry.face2node,
        turtle,
        options,
        geometry.plotbox,
        prepared.emitter_nodes,
        prepared.emitter_node_mask,
        prepared.virtual_nodes,
        prepared.virtual_node_mask,
        geometry.node_ids,
        prepared.cache_ctx,
    )

    @test get(transfer.observed_fraction, (2, 1), 0.0) > 0.0
    @test !haskey(transfer.received_fraction, (2, 1))
    @test get(transfer.received_fraction, (3, 1), 0.0) > 0.0

    emitted_by_band = ArchimedLight._emitter_power_per_band_per_node(scene, models)
    custom_emitted = emitted_by_band["CUSTOM"][1]
    custom_first = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="CUSTOM",
        emitter_band_nir=nothing,
    )
    @test custom_first.incident_power.par[2] > 0.0
    @test custom_first.incident_power.par[3] > 0.0
    @test isapprox(
        custom_first.incident_power.par[3] + custom_first.emitter_escaped_power.par[1],
        custom_emitted;
        atol=1e-10,
        rtol=1e-10,
    )

    # An absent gamma must not inherit PAR emission in a custom-band pass.
    uv_first = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="UV",
        emitter_band_nir=nothing,
    )
    @test all(iszero, values(uv_first.incident_power.par))
    @test all(iszero, values(uv_first.emitter_escaped_power.par))

    # CUSTOM is emitted even without an RI_CUSTOM_f meteo column, while the
    # unrelated zero-irradiance UV band remains free of PAR emitter leakage.
    row = (RI_UV_f=0.0,)
    direct_0, direct_total, direct_irr, direct_escaped = ArchimedLight._compute_extra_band_light(
        scene,
        models,
        row,
        sky,
        turtle,
        options,
    )
    responses = ArchimedLight._build_sector_responses(prepared, scene, models, turtle, options)
    cached_0, cached_total, cached_irr, cached_escaped = ArchimedLight._compute_extra_band_light(
        scene,
        models,
        row,
        sky,
        turtle,
        options;
        responses_cache=responses,
    )
    @test direct_irr["CUSTOM"] == cached_irr["CUSTOM"] == 0.0
    @test direct_irr["UV"] == cached_irr["UV"] == 0.0
    @test direct_0["CUSTOM"][3] > 0.0
    @test all(iszero, values(direct_0["UV"]))
    @test HelperModule._dicts_close(cached_0["CUSTOM"], direct_0["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_total["CUSTOM"], direct_total["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_0["UV"], direct_0["UV"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_total["UV"], direct_total["UV"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_escaped["CUSTOM"], direct_escaped["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_escaped["UV"], direct_escaped["UV"]; atol=1e-10, rtol=1e-10)

    first_direct = ArchimedLight.compute_first_order(scene, models, turtle, zero_fluxes, options)
    first_cached = ArchimedLight._combine_sector_responses(responses, zero_fluxes)
    scat_direct = ArchimedLight.compute_scattering(scene, models, turtle, first_direct, options)
    scat_cached = ArchimedLight.compute_scattering(responses.scattering_topology, first_cached, options)
    @test HelperModule._dicts_close(scat_cached.added_power.par, scat_direct.added_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat_cached.added_power.nir, scat_direct.added_power.nir; atol=1e-10, rtol=1e-10)

    no_sensor_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    no_sensor_first =
        ArchimedLight.compute_first_order(no_sensor_scene, models, turtle, zero_fluxes, options)
    no_sensor_scat =
        ArchimedLight.compute_scattering(no_sensor_scene, models, turtle, no_sensor_first, options)
    for (with_sensor, without_sensor) in ((1, 1), (3, 2))
        @test isapprox(
            first_direct.incident_power.par[with_sensor],
            no_sensor_first.incident_power.par[without_sensor];
            atol=1e-10,
            rtol=1e-10,
        )
        @test isapprox(
            scat_direct.added_power.par[with_sensor],
            no_sensor_scat.added_power.par[without_sensor];
            atol=1e-10,
            rtol=1e-10,
        )
    end
    @test isapprox(
        first_direct.emitter_escaped_power.par[1],
        no_sensor_first.emitter_escaped_power.par[1];
        atol=1e-10,
        rtol=1e-10,
    )

    no_sensor_custom_first = ArchimedLight.compute_first_order(
        no_sensor_scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="CUSTOM",
        emitter_band_nir=nothing,
    )
    custom_scat = ArchimedLight.compute_scattering_band(
        scene,
        models,
        turtle,
        custom_first,
        options;
        band="CUSTOM",
    )
    no_sensor_custom_scat = ArchimedLight.compute_scattering_band(
        no_sensor_scene,
        models,
        turtle,
        no_sensor_custom_first,
        options;
        band="CUSTOM",
    )
    @test isapprox(
        custom_scat.added_power_per_node[3],
        no_sensor_custom_scat.added_power_per_node[2];
        atol=1e-10,
        rtol=1e-10,
    )

    meteo_row = HelperModule._synthetic_meteo_row(
        duration_seconds=600.0,
        ri_par_f=0.0,
        ri_nir_f=0.0,
        direct_fraction=0.0,
    )
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        [meteo_row],
        (; source="synthetic_emitter_cache_parity"),
    )
    uncached_options = ArchimedLight.LightOptions(options; cache_radiation=false)
    cached_options = ArchimedLight.LightOptions(options; cache_radiation=true)
    uncached_step = only(ArchimedLight.run_light_series(scene, models, meteo, uncached_options))
    cached_step = only(ArchimedLight.run_light_series(scene, models, meteo, cached_options))
    @test cached_step.extra_band_irradiance["CUSTOM"] == uncached_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test HelperModule._budgets_close(cached_step.budget, uncached_step.budget; atol=1e-9, rtol=1e-9)
    @test haskey(cached_step.budget.emitter_escaped_energy_per_band, "CUSTOM")
    @test isapprox(
        uncached_step.budget.extra_initial_energy_per_band["CUSTOM"][3] +
        uncached_step.budget.emitter_escaped_energy_per_band["CUSTOM"][1],
        custom_emitted * 600.0;
        atol=1e-7,
        rtol=1e-10,
    )
    @test HelperModule._dicts_close(
        cached_step.budget.extra_initial_energy_per_band["CUSTOM"],
        uncached_step.budget.extra_initial_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.budget.extra_energy_per_band["CUSTOM"],
        uncached_step.budget.extra_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.scattering.added_power.par,
        uncached_step.scattering.added_power.par;
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.scattering.added_power.nir,
        uncached_step.scattering.added_power.nir;
        atol=1e-9,
        rtol=1e-9,
    )

    uncached_sky_step = ArchimedLight.run_light(
        ArchimedLight.LightSimulation(scene, models; options=uncached_options),
        sky;
        step_duration_seconds=600.0,
    )
    cached_sky_step = ArchimedLight.run_light(
        ArchimedLight.LightSimulation(scene, models; options=cached_options),
        sky;
        step_duration_seconds=600.0,
    )
    @test uncached_sky_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test cached_sky_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test uncached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"][3] > 0.0
    @test HelperModule._dicts_close(
        cached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"],
        uncached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_sky_step.budget.extra_energy_per_band["CUSTOM"],
        uncached_sky_step.budget.extra_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
end

@testitem "Emitter wildcard resolution and same-node first hits" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    emitter_model(radiance) = ArchimedLight.TypeModel(
        interception=ArchimedLight.InterceptionModel(),
        light_emitter=ArchimedLight.EmitterModel(
            radiance=radiance,
            gamma=ArchimedLight.OpticalProperties(1.0, 0.0),
        ),
    )
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict("*" => emitter_model(1.0)),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "*" => emitter_model(2.0),
                "panel" => emitter_model(3.0),
            ),
        ),
    ])
    resolution_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=3.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="shade", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="other", type="panel", object_id=3),
    ])
    resolved_par, _ = ArchimedLight._emitter_power_per_node(resolution_scene, models)
    @test isapprox(resolved_par[1], 3pi; atol=1e-12, rtol=1e-12)
    @test isapprox(resolved_par[2], 2pi; atol=1e-12, rtol=1e-12)
    @test isapprox(resolved_par[3], pi; atol=1e-12, rtol=1e-12)

    # A shared-edge triangle duplicate is one emission origin. A second
    # surface of the same node at a genuinely lower depth remains both a
    # physical receiver and a separate emission origin.
    stack = ArchimedLight.SmallHitStack()
    push!(stack, (2.0, 1))
    push!(stack, (2.0 + eps(Float32), 1))
    push!(stack, (1.0, 1))
    push!(stack, (0.0, 2))
    projection = ArchimedLight.DirectionProjectionResult(
        Dict((1, 1) => stack),
        Dict{Int,Int}(),
        Dict{Int,Float64}(),
        Dict{Int,Float64}(),
    )
    edge_counts = Dict{UInt64,Int}()
    observed_edge_counts = Dict{UInt64,Int}()
    total_from = Dict{Int,Int}()
    ArchimedLight._accumulate_emitter_transfer_counts!(
        edge_counts,
        observed_edge_counts,
        total_from,
        projection,
        Set([1]);
        stacks_sorted=false,
    )
    @test total_from[1] == 2
    @test edge_counts[ArchimedLight._pack_emitter_edge(1, 1)] == 1
    @test edge_counts[ArchimedLight._pack_emitter_edge(2, 1)] == 1

    interleaved = ArchimedLight.SmallHitStack()
    push!(interleaved, (2.0 + eps(Float32), 1))
    push!(interleaved, (2.0, 2))
    push!(interleaved, (2.0 - eps(Float32), 1))
    push!(interleaved, (0.0, 3))
    interleaved_projection = ArchimedLight.DirectionProjectionResult(
        Dict((1, 1) => interleaved),
        Dict{Int,Int}(),
        Dict{Int,Float64}(),
        Dict{Int,Float64}(),
    )
    empty!(edge_counts)
    empty!(observed_edge_counts)
    empty!(total_from)
    ArchimedLight._accumulate_emitter_transfer_counts!(
        edge_counts,
        observed_edge_counts,
        total_from,
        interleaved_projection,
        Set([1, 2]);
        stacks_sorted=false,
    )
    @test total_from[1] == 1
    @test total_from[2] == 1

    # Re-map two separated quads to one emitting node. Equal-depth triangle
    # duplicates must be skipped, while the lower surface is a real first hit.
    base_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    nodes = copy(base_scene.nodes)
    nodes[1] = ArchimedLight.PlantGeom.SceneNodeData(
        nodes[1].area + nodes[2].area,
        nodes[1].barycenter,
        nodes[1].source_topology_id,
    )
    delete!(nodes, 2)
    folded_scene = ArchimedLight.PlantGeom.SceneGeometry(
        base_scene.mtg,
        base_scene.merged_mesh,
        [1, 1, 1, 1, 3, 3],
        nodes,
        base_scene.source_path,
        base_scene.scene_xy_bounds,
    )
    folded_models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict("panel" => emitter_model(1.0)),
        ),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )
    emitted_par, _ = ArchimedLight._emitter_power_per_node(folded_scene, folded_models)
    first = ArchimedLight.compute_first_order(folded_scene, folded_models, turtle, zero_fluxes, options)
    @test first.incident_power.par[1] > 0.0
    @test first.incident_power.par[3] > 0.0
    @test isapprox(
        first.incident_power.par[1] + first.incident_power.par[3] + first.emitter_escaped_power.par[1],
        emitted_par[1];
        atol=1e-10,
        rtol=1e-10,
    )
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
    @test isempty(first.emitter_escaped_power.par)
    @test isempty(first.emitter_escaped_power.nir)
    @test isempty(step.first_order.emitter_escaped_power.par)
    @test isempty(step.first_order.emitter_escaped_power.nir)
    @test isempty(step.budget.emitter_escaped_energy_per_band["PAR"])
    @test isempty(step.budget.emitter_escaped_energy_per_band["NIR"])
    @test step.component_metadata === nothing
    @test_throws ArgumentError ArchimedLight.component_values(step)

    fallback_row = merge(
        row,
        (
            hour_end="not-a-time",
            step_duration="not-a-duration",
            duration=1800.0,
        ),
    )
    fallback_budget = ArchimedLight.integrate_light(
        scene,
        models,
        first,
        nothing,
        options;
        meteo_row=fallback_row,
    )
    explicit_budget = ArchimedLight.integrate_light(
        scene,
        models,
        first,
        nothing,
        options;
        step_duration_seconds=1800.0,
    )
    @test fallback_budget.incident_energy.total.par ==
          explicit_budget.incident_energy.total.par
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
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_cached_series"))
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
    using Dates

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
    @test meteo_read isa PM.TimeStepTable
    @test first(selected).Ri_SW_f == 200.0
    @test isapprox(sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)

    datetime_only = [
        (date=DateTime(2020, 6, 21, 12), latitude=15.0, Ri_SW_f=200.0),
        (date=DateTime(2020, 6, 21, 13), latitude=15.0, Ri_SW_f=220.0),
    ]
    inferred = ArchimedLight.prepare_meteo(datetime_only, options)
    @test length(inferred) == 2
    @test ArchimedLight._resolved_meteo_step_or_error(first(inferred), options).duration_seconds == 3600.0
    datetime_range_options = ArchimedLight.LightOptions(
        options;
        meteo_range="2020-06-21-12:30,2020-06-21-12:45",
    )
    @test length(
        ArchimedLight.prepare_meteo(datetime_only, datetime_range_options),
    ) == 1

    datetime_with_inactive_bad_row = [
        (date="bad", latitude=15.0, Ri_SW_f=Inf, active=false),
        (date=DateTime(2020, 6, 21, 12), latitude=15.0, Ri_SW_f=200.0, active=true),
        (date=DateTime(2020, 6, 21, 13), latitude=15.0, Ri_SW_f=220.0, active=true),
    ]
    inferred_selected =
        ArchimedLight.prepare_meteo(datetime_with_inactive_bad_row, options)
    @test length(inferred_selected) == 2
    @test ArchimedLight._resolved_meteo_step_or_error(
        first(inferred_selected),
        options,
    ).duration_seconds == 3600.0
end

@testitem "Generic step_duration input and boundary execution overrides" tags = [:synthetic, :fast, :meteo] setup = [HelperModule] begin
    using Tables

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    row = (
        step_duration=3600.0,
        Ri_SW_f=100.0,
        Ri_PAR_f=-5.0,
        Ri_NIR_f=105.0,
        sun_azimuth=180.0,
        sun_elevation=45.0,
        direct_fraction=0.5,
    )
    generic = [row]
    column_table = Tables.columntable(generic)
    meteo = ArchimedLight.read_meteo(generic)
    @test meteo isa PM.TimeStepTable
    @test first(meteo).step_duration == 3600.0
    @test length(ArchimedLight.prepare_meteo(generic, options; check_boundaries=false)) == 1

    strict_sim = ArchimedLight.LightSimulation(scene, models; options=options)
    @test_throws ErrorException ArchimedLight.run_light(strict_sim, row)
    @test strict_sim.cache === nothing

    sim = ArchimedLight.LightSimulation(scene, models; options=options)
    cache = ArchimedLight.prepare_light_cache(
        scene,
        models,
        options;
        memory_limit_bytes=0,
    )
    steps = Any[
        ArchimedLight.run_light(sim, row; check_boundaries=false),
        only(ArchimedLight.run_light(sim, meteo; check_boundaries=false)),
        only(ArchimedLight.run_light(sim, column_table; check_boundaries=false)),
        ArchimedLight.run_light_step(scene, models, row, options; check_boundaries=false),
        ArchimedLight.run_light_step(cache, row; check_boundaries=false),
        only(ArchimedLight.run_light_series(scene, models, meteo, options; check_boundaries=false)),
        only(ArchimedLight.run_light_series(cache, meteo; check_boundaries=false)),
    ]
    @test all(step -> step.sky.ri_par_f == 0.0, steps)
    @test all(step -> step.sky.ri_nir_f == 105.0, steps)

    empty_meteo = ArchimedLight.read_meteo((
        step_duration=Float64[],
        Ri_SW_f=Float64[],
        sun_azimuth=Float64[],
        sun_elevation=Float64[],
        direct_fraction=Float64[],
    ))
    empty_sim = ArchimedLight.LightSimulation(scene, models; options=options)
    @test_throws ErrorException ArchimedLight.prepare_meteo(empty_meteo, options)
    @test_throws ErrorException ArchimedLight.run_light(empty_sim, empty_meteo)
    @test empty_sim.cache === nothing
end

@testitem "Synthetic case NIR disabled preserves PAR from shortwave" tags = [:synthetic, :fast, :nir_interception] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = ArchimedLight.LightOptions(
        turtle_sectors=46,
        all_in_turtle=true,
        pixel_size=0.01,
        toricity=true,
        cache_radiation=false,
        nir_interception=false,
    )
    meteo = [(
        date=Date("2020-06-21"),
        hour_start=Time(12),
        hour_end=Time(13),
        step_duration=3600.0,
        latitude=15.0,
        Ri_SW_f=200.0,
        direct_fraction=1.0,
        sun_azimut=180.0,
        sun_elevation=80.0,
    )]

    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    step = only(series)

    @test isapprox(step.sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)
    @test isapprox(step.sky.ri_par_f, 0.48 * 200.0; atol=1e-9, rtol=1e-9)
    @test isapprox(step.sky.ri_nir_f, 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(step.fluxes.par), step.sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(step.fluxes.nir), 0.0; atol=1e-12, rtol=1e-12)

    emitter_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    emitter_models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=ArchimedLight.OpticalProperties(0.2, 0.5),
                    ),
                ),
            ),
        ),
    ])
    emitter_meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        meteo,
        (; source="synthetic_nir_disabled_emitter"),
    )
    for cache_radiation in (false, true)
        emitter_options = ArchimedLight.LightOptions(
            options;
            scattering=false,
            cache_radiation=cache_radiation,
        )
        emitter_step = only(
            ArchimedLight.run_light_series(
                emitter_scene,
                emitter_models,
                emitter_meteo,
                emitter_options,
            ),
        )
        @test sum(values(emitter_step.first_order.incident_power.par)) > 0.0
        @test all(iszero, values(emitter_step.first_order.incident_power.nir))
        @test all(iszero, values(emitter_step.first_order.emitter_escaped_power.nir))
        @test all(iszero, values(emitter_step.budget.incident_energy.total.nir))
        @test all(
            iszero,
            values(emitter_step.budget.emitter_escaped_energy_per_band["NIR"]),
        )
    end
end

@testitem "Synthetic case overlapping meteo steps option" tags = [:synthetic, :fast, :overlapping_meteo_steps] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()

    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 0), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 15), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_overlap"))

    strict_options = HelperModule._synthetic_options(cache_radiation=false)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.prepare_meteo(meteo, strict_options)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.run_light_series(scene, models, meteo, strict_options)

    active_meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        [merge(rows[1], (active=true,)), merge(rows[2], (active=false,))],
        (; source="synthetic_overlap_inactive"),
    )
    @test length(ArchimedLight.prepare_meteo(active_meteo, strict_options)) == 1

    permissive_options = ArchimedLight.LightOptions(strict_options; allow_overlapping_meteo_steps=true)
    selected = ArchimedLight.prepare_meteo(meteo, permissive_options)
    series = ArchimedLight.run_light_series(scene, models, meteo, permissive_options)

    @test length(selected) == 2
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
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_cached_scattering_series"))
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
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_manual_cache"))

    cache = ArchimedLight.prepare_light_cache(scene, models, options)
    summary0 = ArchimedLight.cache_summary(cache)
    @test summary0.mode == :full
    @test summary0.cached_turtle_count == 0
    @test summary0.component_metadata_count == 0
    @test summary0.component_metadata_bytes == 0

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
    entry = only(values(cache.entries))
    @test entry.resident_bytes == ArchimedLight._turtle_cache_entry_retained_bytes(entry)
    @test summary1.resident_bytes == entry.resident_bytes

    scene2 = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.5, group="middle", type="plate", object_id=3),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    cache2 = ArchimedLight.prepare_light_cache(scene2, models, options)
    step_rebuilt = ArchimedLight.run_light_step(cache2, rows[1])
    @test !HelperModule._budgets_close(step_cached.budget, step_rebuilt.budget; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case plain light cache accounting" tags = [:synthetic, :fast, :plain_light_cache_accounting] setup = [HelperModule] begin
    using ArchimedLight

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        cache_radiation=true,
    )
    row = HelperModule._synthetic_meteo_row()

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=10^9)
    ArchimedLight.run_light_step(cache, row)
    entry = only(values(cache.entries))
    responses = entry.responses_cache
    @test responses.emitter_transfer === nothing
    @test responses.scattering_topology === nothing

    fast_retained = ArchimedLight._plain_turtle_cache_entry_retained_bytes(entry)
    walked_retained = ArchimedLight._summarysize_turtle_cache_entry_retained_bytes(entry)
    @test entry.resident_bytes == fast_retained
    @test fast_retained >= walked_retained
    @test ArchimedLight.cache_summary(cache).resident_bytes == fast_retained

    limited_bytes = max(cache.estimated_entry_bytes, fast_retained - 1)
    @test limited_bytes < fast_retained
    limited = ArchimedLight.prepare_light_cache(
        scene,
        models,
        options;
        memory_limit_bytes=limited_bytes,
    )
    @test ArchimedLight.cache_summary(limited).mode == :full
    ArchimedLight.run_light_step(limited, row)
    summary = ArchimedLight.cache_summary(limited)
    @test summary.mode == :topology_fallback
    @test summary.cached_turtle_count == 0
    @test summary.resident_bytes == 0
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
    meteo = ArchimedLight.PlantMeteo.TimeStepTable([row0, row1], (; source="synthetic_extra_band_cache"))

    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options)

    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test HelperModule._budgets_close(cached[i].budget, uncached[i].budget; atol=1e-9, rtol=1e-9)
    end

    probe = ArchimedLight.prepare_light_cache(
        scene,
        models,
        cached_options;
        memory_limit_bytes=10^9,
    )
    exact_limit = probe.estimated_entry_bytes
    limited = ArchimedLight.prepare_light_cache(
        scene,
        models,
        cached_options;
        memory_limit_bytes=exact_limit,
    )
    @test ArchimedLight.cache_summary(limited).mode == :full
    limited_step = ArchimedLight.run_light_step(limited, row0)
    uncached_step = ArchimedLight.run_light_step(scene, models, row0, uncached_options)
    @test HelperModule._budgets_close(
        limited_step.budget,
        uncached_step.budget;
        atol=1e-9,
        rtol=1e-9,
    )
    @test ArchimedLight.cache_summary(limited).resident_bytes <= exact_limit
    @test all(
        entry.resident_bytes == ArchimedLight._turtle_cache_entry_retained_bytes(entry)
        for entry in values(limited.entries)
    )
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
