@testitem "Voxel Beer-Lambert response and spectral energy balance" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (3, 2, 3), fill(0.4, 3, 2, 3))
    backend = VoxelCPUBackend(rays_per_voxel=4, traversal=:dda, boundary=:periodic)
    response = ArchimedLight.compute_voxel_direction_response(grid, (0, 0, -1), backend)

    expected_layers = [
        0.12150841951370239,
        0.14841067790985107,
        0.18126922845840454,
    ]
    for i in axes(grid.pad, 1), j in axes(grid.pad, 2), k in axes(grid.pad, 3)
        @test isapprox(response.interception_fraction[i, j, k], expected_layers[k]; atol=4e-8, rtol=0)
    end
    expected_column_interception = 1 - exp(-0.5 * 0.4 * 3.0)
    @test isapprox(sum(response.interception_fraction), 6expected_column_interception; atol=1e-12)
    @test isapprox(response.escaped_weight, 6exp(-0.5 * 0.4 * 3.0); atol=1e-12)
    @test isapprox(
        sum(response.interception_fraction) + response.escaped_weight,
        response.incoming_weight;
        atol=1e-12,
    )
    @test response.injected_weight == 0.0
    @test response.effective_ray_count == 4

    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1,
            ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0),
            1.0,
            :sky,
        ),
    ])
    directional_fluxes = ArchimedLight.DirectionalFluxes([1], [100.0], [50.0])
    cache = prepare_voxel_responses(grid, turtle, backend)
    result = compute_voxel_first_order(
        grid,
        turtle,
        directional_fluxes,
        backend;
        duration_seconds=3600.0,
        cache=cache,
    )
    @test isapprox(sum(result.par_energy) + result.par_escaped_energy, result.par_incoming_energy; atol=1e-8)
    @test isapprox(sum(result.nir_energy) + result.nir_escaped_energy, result.nir_incoming_energy; atol=1e-8)
    @test result.par_incoming_energy == 100.0 * 3 * 2 * 3600.0
    @test result.nir_incoming_energy == 50.0 * 3 * 2 * 3600.0
    @test all(isfinite, result.par_flux)
    @test all(result.par_flux .>= 0)

    @test_throws ArgumentError compute_voxel_first_order(
        grid,
        turtle,
        directional_fluxes,
        VoxelCPUBackend(rays_per_voxel=4, g=0.4, traversal=:dda, boundary=:periodic);
        duration_seconds=3600.0,
        cache=cache,
    )

    changed_turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1,
            ArchimedLight.StaticArrays.SVector(0.1, 0.0, sqrt(0.99)),
            1.0,
            :sky,
        ),
    ])
    @test_throws ArgumentError compute_voxel_first_order(
        grid,
        changed_turtle,
        directional_fluxes,
        backend;
        duration_seconds=3600.0,
        cache=cache,
    )

    empty_grid = VoxelGrid((0, 0, 0), (1, 1, 1), zeros(1, 1, 1))
    empty_result = compute_voxel_first_order(
        empty_grid,
        turtle,
        directional_fluxes,
        VoxelCPUBackend(rays_per_voxel=1);
        duration_seconds=3600.0,
    )
    @test empty_result.par_energy == zeros(1, 1, 1)
    @test empty_result.par_flux == zeros(1, 1, 1)
    @test empty_result.par_escaped_energy == empty_result.par_incoming_energy

    grid.pad[1, 1, 1] = 0.5
    @test_throws ArgumentError compute_voxel_first_order(
        grid,
        turtle,
        directional_fluxes,
        backend;
        duration_seconds=3600.0,
        cache=cache,
    )
end

@testitem "Voxel sampling and analytical Beer-Lambert invariants" tags=[:voxel, :fast] begin
    unequal = VoxelGrid((0, 0, 0), (4, 1, 1), fill(0.2, 2, 1, 1))
    offsets_1 = ArchimedLight._voxel_reference_ray_offsets(unequal, 1)
    offsets_3 = ArchimedLight._voxel_reference_ray_offsets(unequal, 3)
    @test offsets_1 == [(1.0, 0.5)]
    @test length(offsets_3) == 2
    @test all(0 < offset[1] < 2 && 0 < offset[2] < 1 for offset in offsets_3)

    pad = reshape([0.1, 0.2, 0.4], 1, 1, 3)
    layered = VoxelGrid((0, 0, 0), (1, 1, 3), pad)
    response = ArchimedLight.compute_voxel_direction_response(
        layered,
        (0, 0, -1),
        VoxelCPUBackend(rays_per_voxel=1, g=0.5),
    )
    remaining = let incoming = 1.0
        for k in 3:-1:1
            outgoing = incoming * exp(-0.5 * pad[1, 1, k])
            @test isapprox(
                response.interception_fraction[1, 1, k], incoming - outgoing; atol=1e-14
            )
            incoming = outgoing
        end
        incoming
    end
    @test isapprox(response.escaped_weight, remaining; atol=1e-14)

    coarse = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.4, 1, 1, 3))
    fine = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.4, 1, 1, 12))
    coarse_response = ArchimedLight.compute_voxel_direction_response(
        coarse, (0, 0, -1), VoxelCPUBackend(rays_per_voxel=1)
    )
    fine_response = ArchimedLight.compute_voxel_direction_response(
        fine, (0, 0, -1), VoxelCPUBackend(rays_per_voxel=1)
    )
    analytical = 1 - exp(-0.5 * 0.4 * 3)
    @test isapprox(sum(coarse_response.interception_fraction), analytical; atol=1e-14)
    @test isapprox(sum(fine_response.interception_fraction), analytical; atol=1e-14)

    transparent = ArchimedLight.compute_voxel_direction_response(
        coarse, (0, 0, -1), VoxelCPUBackend(rays_per_voxel=1, g=0.0)
    )
    @test iszero(sum(transparent.interception_fraction))
    @test transparent.escaped_weight == transparent.incoming_weight
end

@testitem "Voxel meteorological step and series reuse ArchimedLight sky fluxes" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=true,
        scattering=false,
        radiation_timestep_minutes=5.0,
    )
    row = (
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    backend = VoxelCPUBackend(rays_per_voxel=1)
    step = run_voxel_light_step(grid, row, options; backend=backend)
    @test step.duration_seconds == 3600.0
    @test sum(step.fluxes.par) == step.sky.ri_par_f
    @test sum(step.fluxes.nir) == step.sky.ri_nir_f
    @test isapprox(
        sum(step.first_order.par_energy) + step.first_order.par_escaped_energy,
        step.first_order.par_incoming_energy;
        atol=1e-8,
    )

    table = ArchimedLight._as_plantmeteo_table([row, merge(row, (Ri_PAR_f=80.0, Ri_NIR_f=40.0))])
    series = run_voxel_light_series(grid, table, options; backend=backend)
    @test length(series) == 2
    @test series[2].first_order.par_incoming_energy == 0.8 * series[1].first_order.par_incoming_energy

    @test_throws ArgumentError run_voxel_light_step(
        grid,
        row,
        LightOptions(options; scattering=true);
        backend=backend,
    )

    night = merge(row, (Ri_PAR_f=0.0, Ri_NIR_f=0.0, direct_fraction=0.0, sun_elevation=-10.0))
    night_step = run_voxel_light_step(grid, night, options; backend=backend)
    @test iszero(sum(night_step.first_order.par_energy))
    @test iszero(sum(night_step.first_order.nir_energy))
    @test night_step.first_order.par_incoming_energy == 0.0

    dynamic_options = LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        scattering=false,
        radiation_timestep_minutes=5.0,
    )
    dynamic_rows = [
        merge(row, (direct_fraction=1.0, sun_azimuth=120.0, sun_elevation=40.0)),
        merge(row, (direct_fraction=1.0, sun_azimuth=240.0, sun_elevation=40.0)),
    ]
    dynamic_series = run_voxel_light_series(
        grid, dynamic_rows, dynamic_options; backend=backend
    )
    @test length(dynamic_series) == 2
    @test last(dynamic_series[1].turtle.sectors).direction !=
          last(dynamic_series[2].turtle.sectors).direction
    @test all(
        isapprox(
            sum(item.first_order.par_energy) + item.first_order.par_escaped_energy,
            item.first_order.par_incoming_energy;
            atol=1e-8,
        ) for item in dynamic_series
    )
end

@testitem "Voxel ports the scientific intentions of Java TestVoxel1-TestVoxel4" tags=[:voxel, :java_fixture, :fast] begin
    function rotated_direction(azimuth_degrees; elevation_degrees=37.0)
        azimuth = deg2rad(azimuth_degrees)
        horizontal = cosd(elevation_degrees)
        return (
            horizontal * cos(azimuth),
            horizontal * sin(azimuth),
            -sind(elevation_degrees),
        )
    end

    # TestVoxel1: a horizontally uniform periodic canopy stays horizontally
    # uniform and its vertical interception follows one exponential profile.
    uniform = VoxelGrid((0, 0, 0), (7, 5, 6), fill(0.1, 7, 5, 6))
    combined = zeros(size(uniform))
    backend = VoxelCPUBackend(rays_per_voxel=4, boundary=:periodic)
    for azimuth in (0.0, 90.0, 180.0, 270.0)
        combined .+= ArchimedLight.compute_voxel_direction_response(
            uniform, rotated_direction(azimuth), backend
        ).interception_fraction
    end
    for k in axes(combined, 3)
        layer = @view combined[:, :, k]
        @test maximum(layer) - minimum(layer) <= 2e-14
    end
    logarithmic_layers = log.([sum(@view combined[:, :, k]) for k in axes(combined, 3)])
    layer_deltas = diff(logarithmic_layers)
    @test maximum(layer_deltas) - minimum(layer_deltas) <= 2e-13

    # TestVoxel2: rotating a direction at fixed elevation does not change the
    # total response of a periodic homogeneous ground layer.
    ground_pad = zeros(12, 12, 6)
    ground_pad[:, :, 1] .= 0.5
    ground = VoxelGrid((0, 0, 0), (6, 6, 3), ground_pad)
    ground_totals = [
        sum(
            ArchimedLight.compute_voxel_direction_response(
                ground, rotated_direction(azimuth), backend
            ).interception_fraction,
        ) for azimuth in 0.0:15.0:90.0
    ]
    @test maximum(ground_totals) - minimum(ground_totals) <= 2e-12

    # TestVoxel3: a square central pillar is stable under quarter-turns.
    pillar_pad = zeros(11, 11, 6)
    pillar_pad[5:7, 5:7, :] .= 0.5
    pillar = VoxelGrid((0, 0, 0), (5.5, 5.5, 3), pillar_pad)
    pillar_totals = [
        sum(
            ArchimedLight.compute_voxel_direction_response(
                pillar,
                rotated_direction(azimuth; elevation_degrees=87.0),
                VoxelCPUBackend(rays_per_voxel=9, boundary=:open),
            ).interception_fraction,
        ) for azimuth in (0.0, 90.0, 180.0, 270.0)
    ]
    @test maximum(pillar_totals) - minimum(pillar_totals) <= 2e-12

    # TestVoxel4: periodic and open boundaries agree for a vertical ray, while
    # an oblique open ray can escape through a side instead of wrapping.
    vertical_periodic = ArchimedLight.compute_voxel_direction_response(
        ground, (0, 0, -1), VoxelCPUBackend(rays_per_voxel=1, boundary=:periodic)
    )
    vertical_open = ArchimedLight.compute_voxel_direction_response(
        ground, (0, 0, -1), VoxelCPUBackend(rays_per_voxel=1, boundary=:open)
    )
    @test vertical_periodic.interception_fraction == vertical_open.interception_fraction
    oblique_periodic = ArchimedLight.compute_voxel_direction_response(
        ground, rotated_direction(0.0), VoxelCPUBackend(rays_per_voxel=1, boundary=:periodic)
    )
    oblique_open = ArchimedLight.compute_voxel_direction_response(
        ground, rotated_direction(0.0), VoxelCPUBackend(rays_per_voxel=1, boundary=:open)
    )
    @test sum(oblique_open.interception_fraction) < sum(oblique_periodic.interception_fraction)
    @test oblique_open.escaped_weight > oblique_periodic.escaped_weight
end
