@testitem "Voxel optical properties are explicit, validated, and band-specific" tags=[:voxel, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 1, 1), fill(0.5, 2, 1, 1))

    scalar = VoxelOpticalProperties(
        grid;
        par_reflectance=0.08,
        par_transmittance=0.07,
        nir_reflectance=0.45,
        nir_transmittance=0.40,
    )
    @test voxel_single_scattering_albedo(scalar, :par) ≈ 0.15
    @test voxel_single_scattering_albedo(scalar, :nir) ≈ 0.85
    @test voxel_absorptance(scalar, :par) == 0.85
    @test (@inferred voxel_single_scattering_albedo(scalar, :par)) ≈ 0.15
    @test scalar.par_reflectance isa Float64
    @test scalar.par_transmittance isa Float64
    @test VOXEL_PAR_WAVELENGTH_NM == (400.0, 700.0)
    @test VOXEL_NIR_WAVELENGTH_NM == (700.0, 2500.0)

    par_r = reshape([0.05, 0.10], size(grid))
    heterogeneous = VoxelOpticalProperties(
        grid;
        par_reflectance=par_r,
        par_transmittance=0.05,
        nir_reflectance=0.40,
        nir_transmittance=0.45,
    )
    @test voxel_single_scattering_albedo(heterogeneous, :par) == par_r .+ 0.05
    @test (@inferred voxel_single_scattering_albedo(heterogeneous, :par)) ≈ par_r .+ 0.05
    @test voxel_absorptance(heterogeneous, :par) ≈ 1 .- par_r .- 0.05
    @test voxel_single_scattering_albedo(heterogeneous, :par, CartesianIndex(2, 1, 1)) ≈ 0.15

    java = voxel_java_parity_optics(grid)
    generic = voxel_generic_green_leaf_optics(grid)
    @test voxel_single_scattering_albedo(java, :par) == 0.15
    @test voxel_single_scattering_albedo(java, :nir) == 0.90
    @test voxel_single_scattering_albedo(generic, :par) == 0.17
    @test voxel_single_scattering_albedo(generic, :nir) == 0.85

    @test_throws ArgumentError VoxelOpticalProperties(
        grid;
        par_reflectance=-0.01,
        par_transmittance=0.0,
        nir_reflectance=0.4,
        nir_transmittance=0.4,
    )
    @test_throws ArgumentError VoxelOpticalProperties(
        grid;
        par_reflectance=0.6,
        par_transmittance=0.5,
        nir_reflectance=0.4,
        nir_transmittance=0.4,
    )
    @test_throws ArgumentError VoxelOpticalProperties(
        grid;
        par_reflectance=NaN,
        par_transmittance=0.0,
        nir_reflectance=0.4,
        nir_transmittance=0.4,
    )
    @test_throws DimensionMismatch VoxelOpticalProperties(
        grid;
        par_reflectance=zeros(1, 1, 1),
        par_transmittance=0.0,
        nir_reflectance=0.4,
        nir_transmittance=0.4,
    )
    @test_throws ArgumentError VoxelOpticalProperties(
        grid;
        scattering_par=0.15,
        par_reflectance=0.1,
        par_transmittance=0.05,
        scattering_nir=0.9,
    )
    @test_throws ArgumentError voxel_single_scattering_albedo(scalar, :uv)

    ground = VoxelGroundOptics(grid; par_reflectance=0.1, nir_reflectance=fill(0.3, 2, 1))
    @test ground.par_reflectance == 0.1
    @test ground.nir_reflectance == fill(0.3, 2, 1)
    @test_throws DimensionMismatch VoxelGroundOptics(
        grid;
        par_reflectance=zeros(1, 1),
        nir_reflectance=0.3,
    )
end

@testitem "Voxel full-sphere quadrature is normalized and excludes the sun" tags=[:voxel, :scattering, :fast] begin
    sectors = [
        ArchimedLight.TurtleSector(
            1,
            ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0),
            0.25,
            :sky,
        ),
        ArchimedLight.TurtleSector(
            2,
            ArchimedLight.StaticArrays.SVector(0.6, 0.0, 0.8),
            0.75,
            :sky,
        ),
        ArchimedLight.TurtleSector(
            3,
            ArchimedLight.StaticArrays.SVector(0.0, 0.6, 0.8),
            1.0,
            :sun,
        ),
    ]
    quadrature = prepare_voxel_scattering_quadrature(ArchimedLight.TurtleGrid(sectors))
    @test length(quadrature.directions) == 4
    @test isapprox(sum(quadrature.weights), 1.0; atol=1e-15)
    @test all(quadrature.weights .>= 0)
    @test quadrature.directions[2] == .-quadrature.directions[1]
    @test quadrature.directions[4] == .-quadrature.directions[3]
    @test sum(quadrature.weights[1:2:end]) == 0.5
    @test sum(quadrature.weights[2:2:end]) == 0.5
    @test quadrature.weights == [0.125, 0.125, 0.375, 0.375]

    no_sky = ArchimedLight.TurtleGrid([sectors[end]])
    @test_throws ArgumentError prepare_voxel_scattering_quadrature(no_sky)
    invalid = VoxelScatteringQuadrature([(0.0, 0.0, 1.0)], [0.5])
    @test_throws ArgumentError ArchimedLight._validate_voxel_scattering_quadrature(invalid)
    horizontal = VoxelScatteringQuadrature([(1.0, 0.0, 0.0)], [1.0])
    @test_throws ArgumentError ArchimedLight._validate_voxel_scattering_quadrature(horizontal)
end

@testitem "Matrix-free internal voxel transport closes and matches limiting cases" tags=[:voxel, :scattering, :fast] begin
    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1,
            ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0),
            1.0,
            :sky,
        ),
    ])
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    grid = VoxelGrid((0, 0, 0), (1, 1, 1), fill(0.8, 1, 1, 1))
    backend = VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:periodic, traversal=:dda)
    emitted = fill(1.0, 1, 1, 1)
    result = apply_voxel_scattering_transport(grid, emitted, quadrature, backend)
    transmission = exp(-0.5 * 0.8 * 0.5)
    @test isapprox(sum(result.intercepted_energy), 1 - transmission; atol=1e-14)
    @test isapprox(result.escaped_top_energy, 0.5transmission; atol=1e-14)
    @test isapprox(result.escaped_bottom_energy, 0.5transmission; atol=1e-14)
    @test result.escaped_side_energy == 0.0
    @test abs(result.balance_residual) <= 1e-14
    apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                                     cache=prepare_voxel_scattering_transport(grid, quadrature, backend))
    allocation_cache = prepare_voxel_scattering_transport(grid, quadrature, backend)
    cached_bytes = @allocated apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        cache=allocation_cache,
    )
    @test cached_bytes <= 16_384

    reference = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:periodic, traversal=:reference),
    )
    @test reference.intercepted_energy ≈ result.intercepted_energy
    @test reference.escaped_top_energy ≈ result.escaped_top_energy
    @test reference.escaped_bottom_energy ≈ result.escaped_bottom_energy

    two_layers = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.8, 1, 1, 2))
    lower_source = zeros(1, 1, 2)
    lower_source[1, 1, 1] = 1.0
    upward_transfer = apply_voxel_scattering_transport(
        two_layers,
        lower_source,
        quadrature,
        backend,
    )
    @test upward_transfer.intercepted_energy[1, 1, 2] > 0
    @test upward_transfer.intercepted_energy[1, 1, 1] >
          upward_transfer.intercepted_energy[1, 1, 2]

    periodic_grid = VoxelGrid((0, 0, 0), (3, 1, 2), fill(0.6, 3, 1, 2))
    oblique = VoxelScatteringQuadrature(
        [(0.6, 0.0, 0.8), (-0.6, 0.0, -0.8)],
        [0.5, 0.5],
    )
    source_a = zeros(size(periodic_grid))
    source_b = zeros(size(periodic_grid))
    source_a[1, 1, 1] = 1.0
    source_b[2, 1, 1] = 1.0
    translated_a = apply_voxel_scattering_transport(
        periodic_grid,
        source_a,
        oblique,
        backend,
    )
    translated_b = apply_voxel_scattering_transport(
        periodic_grid,
        source_b,
        oblique,
        backend,
    )
    @test translated_b.intercepted_energy ≈ circshift(translated_a.intercepted_energy, (1, 0, 0))

    rotational_grid = VoxelGrid((0, 0, 0), (3, 3, 1), fill(0.6, 3, 3, 1))
    horizontal = 0.6
    vertical = 0.8
    rotational_quadrature = VoxelScatteringQuadrature(
        [
            (horizontal, 0.0, vertical),
            (-horizontal, 0.0, vertical),
            (0.0, horizontal, vertical),
            (0.0, -horizontal, vertical),
            (-horizontal, 0.0, -vertical),
            (horizontal, 0.0, -vertical),
            (0.0, -horizontal, -vertical),
            (0.0, horizontal, -vertical),
        ],
        fill(1 / 8, 8),
    )
    rotational_source = zeros(size(rotational_grid))
    rotational_source[2, 2, 1] = 1.0
    rotational_result = apply_voxel_scattering_transport(
        rotational_grid,
        rotational_source,
        rotational_quadrature,
        backend,
    )
    @test rotational_result.intercepted_energy ≈
          permutedims(rotational_result.intercepted_energy, (2, 1, 3))

    open_quadrature = VoxelScatteringQuadrature(
        [(sqrt(3) / 2, 0.0, 0.5), (-sqrt(3) / 2, 0.0, -0.5)],
        [0.5, 0.5],
    )
    open_result = apply_voxel_scattering_transport(
        grid,
        emitted,
        open_quadrature,
        VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:open, traversal=:dda),
    )
    @test open_result.escaped_side_energy > 0
    @test open_result.escaped_top_energy == 0
    @test open_result.escaped_bottom_energy == 0
    @test isapprox(
        sum(open_result.intercepted_energy) + open_result.escaped_side_energy,
        1.0;
        atol=1e-14,
    )

    open_vertical = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:open, traversal=:dda),
    )
    @test open_vertical.escaped_top_energy > 0
    @test open_vertical.escaped_bottom_energy > 0
    @test open_vertical.escaped_side_energy == 0

    dense_grid = VoxelGrid((0, 0, 0), (2, 1, 2), fill(0.4, 2, 1, 2))
    dense_emitted = reshape([0.3, 0.1, 0.4, 0.2], size(dense_grid))
    matrix = ArchimedLight._dense_voxel_scattering_exchange(dense_grid, quadrature, backend)
    matrix_free = apply_voxel_scattering_transport(
        dense_grid,
        dense_emitted,
        quadrature,
        backend,
    )
    @test vec(matrix_free.intercepted_energy) ≈ matrix * vec(dense_emitted)

    uniform = fill(1.0, size(dense_grid))
    uniform_result = apply_voxel_scattering_transport(
        dense_grid,
        uniform,
        quadrature,
        backend,
    )
    @test uniform_result.intercepted_energy[1, :, :] ≈
          uniform_result.intercepted_energy[2, :, :]

    black = VoxelGroundOptics(grid; par_reflectance=0.0, nir_reflectance=0.0)
    white = VoxelGroundOptics(grid; par_reflectance=1.0, nir_reflectance=1.0)
    no_foliage_source = zeros(1, 1, 1)
    ground_source = fill(1.0, 1, 1)
    black_result = apply_voxel_scattering_transport(
        grid,
        no_foliage_source,
        quadrature,
        backend;
        band=:par,
        ground=black,
        initial_ground_energy=ground_source,
    )
    @test black_result.ground_absorbed_energy == 1.0
    @test black_result.ground_reflected_energy == 0.0
    @test iszero(sum(black_result.intercepted_energy))

    foliage_source = fill(1.0, 1, 1, 1)
    without_ground = apply_voxel_scattering_transport(
        grid,
        foliage_source,
        quadrature,
        backend,
    )
    with_black_ground = apply_voxel_scattering_transport(
        grid,
        foliage_source,
        quadrature,
        backend;
        ground=black,
    )
    @test with_black_ground.intercepted_energy ≈ without_ground.intercepted_energy
    @test with_black_ground.escaped_top_energy ≈ without_ground.escaped_top_energy
    @test with_black_ground.ground_absorbed_energy ≈ without_ground.escaped_bottom_energy

    white_result = apply_voxel_scattering_transport(
        grid,
        no_foliage_source,
        quadrature,
        backend;
        band=:par,
        ground=white,
        initial_ground_energy=ground_source,
    )
    @test white_result.ground_absorbed_energy == 0.0
    @test white_result.ground_reflected_energy == 1.0
    @test isapprox(sum(white_result.intercepted_energy), 1 - exp(-0.4); atol=1e-14)
    @test isapprox(white_result.escaped_top_energy, exp(-0.4); atol=1e-14)
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid,
        no_foliage_source,
        quadrature,
        backend;
        initial_ground_energy=ground_source,
    )

    empty = VoxelGrid((0, 0, 0), (1, 1, 1), zeros(1, 1, 1))
    empty_result = apply_voxel_scattering_transport(
        empty,
        zeros(1, 1, 1),
        quadrature,
        backend,
    )
    @test iszero(sum(empty_result.intercepted_energy))
    @test_throws ArgumentError apply_voxel_scattering_transport(
        empty,
        ones(1, 1, 1),
        quadrature,
        backend,
    )

    cache = prepare_voxel_scattering_transport(grid, quadrature, backend)
    grid.pad[1, 1, 1] = 0.7
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        cache=cache,
    )
end

@testitem "Iterative voxel scattering is spectral, convergent, and energy-conserving" tags=[:voxel, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.7, 1, 1, 2))
    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1,
            ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0),
            1.0,
            :sky,
        ),
    ])
    fluxes = ArchimedLight.DirectionalFluxes([1], [100.0], [50.0])
    backend = VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:periodic, traversal=:dda)
    first = compute_voxel_first_order(
        grid,
        turtle,
        fluxes,
        backend;
        duration_seconds=3600.0,
    )

    zero_optics = VoxelOpticalProperties(
        grid;
        scattering_par=0.0,
        scattering_nir=0.0,
    )
    convergent_options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=true,
        scattering=true,
        scattering_max_iter=30,
        scattering_stop_ratio=1e-10,
    )
    zero = compute_voxel_scattering(
        grid,
        turtle,
        first,
        zero_optics,
        convergent_options,
        backend;
        duration_seconds=3600.0,
    )
    @test iszero(sum(zero.par.added_intercepted_energy))
    @test zero.par.absorbed_energy == first.par_energy
    @test zero.par.total_intercepted_energy == first.par_energy
    @test zero.par.converged
    @test zero.par.relative_balance_residual <= 1e-14

    high_optics = VoxelOpticalProperties(
        grid;
        scattering_par=0.95,
        scattering_nir=0.80,
    )
    one_iteration = LightOptions(
        convergent_options;
        scattering_max_iter=1,
        scattering_stop_ratio=0.0,
    )
    short = compute_voxel_scattering(
        grid,
        turtle,
        first,
        high_optics,
        one_iteration,
        backend;
        duration_seconds=3600.0,
    )
    @test short.par.iterations == 1
    @test !short.par.converged
    @test short.par.unresolved_energy > 0
    @test short.par.relative_balance_residual <= 1e-14

    five_iterations = LightOptions(one_iteration; scattering_max_iter=5)
    longer = compute_voxel_scattering(
        grid,
        turtle,
        first,
        high_optics,
        five_iterations,
        backend;
        duration_seconds=3600.0,
    )
    @test sum(longer.par.absorbed_energy) >= sum(short.par.absorbed_energy)
    @test longer.par.unresolved_energy <= short.par.unresolved_energy
    @test all(diff(longer.par.order_intercepted_energy) .<= 1e-8)

    spectral_optics = VoxelOpticalProperties(
        grid;
        scattering_par=0.0,
        scattering_nir=0.90,
    )
    spectral = compute_voxel_scattering(
        grid,
        turtle,
        first,
        spectral_optics,
        convergent_options,
        backend;
        duration_seconds=3600.0,
    )
    @test iszero(sum(spectral.par.added_intercepted_energy))
    @test sum(spectral.nir.added_intercepted_energy) > 0

    nir_disabled = compute_voxel_scattering(
        grid,
        turtle,
        first,
        spectral_optics,
        LightOptions(convergent_options; nir_scattering=false),
        backend;
        duration_seconds=3600.0,
    )
    @test !nir_disabled.nir.enabled
    @test nir_disabled.nir.iterations == 0
    @test nir_disabled.nir.unresolved_energy > 0

    ground = VoxelGroundOptics(grid; par_reflectance=0.2, nir_reflectance=0.5)
    grounded = compute_voxel_scattering(
        grid,
        turtle,
        first,
        high_optics,
        convergent_options,
        backend;
        duration_seconds=3600.0,
        ground=ground,
    )
    for (band, incoming, injected) in (
        (grounded.par, first.par_incoming_energy, first.par_injected_energy),
        (grounded.nir, first.nir_incoming_energy, first.nir_injected_energy),
    )
        accounted = sum(band.absorbed_energy) + band.ground_absorbed_energy +
                    band.escaped_top_energy + band.escaped_side_energy +
                    band.escaped_bottom_energy + band.unresolved_energy
        @test isapprox(accounted, incoming + injected; atol=1e-7, rtol=1e-12)
        @test band.relative_balance_residual <= 1e-12
    end

    java = compute_voxel_scattering(
        grid,
        turtle,
        first,
        voxel_java_parity_optics(grid),
        convergent_options,
        backend;
        duration_seconds=3600.0,
    )
    @test sum(java.nir.added_intercepted_energy) > sum(java.par.added_intercepted_energy)
end

@testitem "Voxel pipeline enables explicit scattering and preserves first order" tags=[:voxel, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    backend = VoxelCPUBackend(rays_per_voxel=1)
    row = (
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    disabled_options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=true,
        scattering=false,
        radiation_timestep_minutes=5.0,
    )
    optics = voxel_java_parity_optics(grid)
    ground = VoxelGroundOptics(grid; par_reflectance=0.1, nir_reflectance=0.3)
    plain = run_voxel_light_step(grid, row, disabled_options; backend=backend)
    ignored = run_voxel_light_step(
        grid,
        row,
        disabled_options;
        backend=backend,
        optics=optics,
        ground=ground,
    )
    @test plain.scattering === nothing
    @test ignored.scattering === nothing
    @test all(
        name -> isequal(getfield(plain.first_order, name), getfield(ignored.first_order, name)),
        fieldnames(typeof(plain.first_order)),
    )

    enabled_options = LightOptions(
        disabled_options;
        scattering=true,
        scattering_max_iter=20,
        scattering_stop_ratio=1e-8,
    )
    @test_throws ArgumentError run_voxel_light_step(
        grid,
        row,
        enabled_options;
        backend=backend,
    )
    enabled = run_voxel_light_step(
        grid,
        row,
        enabled_options;
        backend=backend,
        optics=optics,
        ground=ground,
    )
    @test all(
        name -> isequal(getfield(enabled.first_order, name), getfield(plain.first_order, name)),
        fieldnames(typeof(plain.first_order)),
    )
    @test enabled.scattering !== nothing
    @test enabled.scattering.par.enabled
    @test enabled.scattering.par.relative_balance_residual <= 1e-12

    rows = [row, merge(row, (Ri_PAR_f=80.0, Ri_NIR_f=40.0))]
    series = run_voxel_light_series(
        grid,
        rows,
        enabled_options;
        backend=backend,
        optics=optics,
        ground=ground,
    )
    @test length(series) == 2
    @test all(step -> step.scattering !== nothing, series)
    @test isapprox(
        sum(series[2].scattering.par.absorbed_energy),
        0.8sum(series[1].scattering.par.absorbed_energy);
        rtol=1e-12,
    )

    @test_throws ArgumentError run_voxel_light_step(
        grid,
        row,
        enabled_options;
        backend=VoxelCPUBackend(rays_per_voxel=1, traversal=:java_reference),
        optics=optics,
    )
    @test_throws ArgumentError run_voxel_light_step(
        grid,
        row,
        LightOptions(enabled_options; scattering_max_iter=0);
        backend=backend,
        optics=optics,
    )
end
