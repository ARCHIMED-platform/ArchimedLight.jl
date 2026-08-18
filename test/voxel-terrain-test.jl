@testitem "Voxel terrain geometry validates optics and exact intersections" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.45)
    @test soil.par_reflectance == 0.2
    @test soil.nir_reflectance == 0.45
    @test_throws ArgumentError SoilOpticalProperties(par_reflectance=-0.1, nir_reflectance=0.2)
    @test_throws ArgumentError SoilOpticalProperties(par_reflectance=0.1, nir_reflectance=Inf)
    @test_throws ArgumentError TerrainHit(1.0, (0, 0, 0), (0, 0, 2), 1, 1)
    @test_throws ArgumentError TerrainHit(1.0, (0, 0, 0), (0, 0, -1), 1, 1)

    grid = VoxelGrid((0, 0, 0), (2, 2, 3), zeros(2, 2, 3))
    plane = PlanarTerrain(grid, soil; elevation=0.75)
    hit = intersect_terrain(plane, (0.25, 0.25, 3.0), (0, 0, -1), 3.0)
    @test hit.distance == 2.25
    @test hit.point == (0.25, 0.25, 0.75)
    @test hit.normal == (0.0, 0.0, 1.0)
    @test hit.patch_id == 1
    @test hit.material_id == 1
    @test terrain_patch_count(plane) == 4
    @test terrain_patch_point(plane, 4) == (1.5, 1.5, 0.75)
    @test soil_optics(plane, hit).par_reflectance == 0.2
    bright = SoilOpticalProperties(par_reflectance=0.7, nir_reflectance=0.8)
    heterogeneous = PlanarTerrain(
        0.75,
        (0.0, 2.0),
        (0.0, 1.0),
        soil;
        cells=(2, 1),
        material_ids=reshape([1, 2], 2, 1),
        materials=Dict(1 => soil, 2 => bright),
    )
    bright_hit = intersect_terrain(heterogeneous, (1.5, 0.5, 2.0), (0, 0, -1))
    @test bright_hit.patch_id == 2
    @test soil_optics(heterogeneous, bright_hit, :par) == 0.7

    flat = HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        fill(0.75, 3, 3),
        soil,
    )
    for origin in ((0.2, 0.3, 3.0), (1.5, 0.4, 2.5), (0.5, 1.5, 2.0))
        planar_hit = intersect_terrain(plane, origin, (0.2, -0.1, -1.0), 5.0)
        height_hit = intersect_terrain(flat, origin, (0.2, -0.1, -1.0), 5.0)
        @test height_hit !== nothing
        @test height_hit.distance ≈ planar_hit.distance atol=1e-13
        @test all(isapprox.(height_hit.point, planar_hit.point; atol=1e-13))
        @test all(isapprox.(height_hit.normal, planar_hit.normal; atol=1e-13))
    end

    slope = HeightFieldTerrain(
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0 0.0; 1.0 1.0],
        soil,
    )
    slope_hit = intersect_terrain(slope, (0.25, 0.25, 2.0), (0, 0, -1), 3.0)
    @test all(isapprox.(slope_hit.point, (0.25, 0.25, 0.25); atol=1e-13))
    @test all(isapprox.(
        slope_hit.normal,
        (-inv(sqrt(2)), 0.0, inv(sqrt(2)));
        atol=1e-13,
    ))
    @test slope_hit.normal[3] > 0

    # The fixed SW-NE diagonal belongs deterministically to the lower patch id.
    diagonal_hits = [
        intersect_terrain(slope, (0.5, 0.5, 2.0), (0, 0, -1), 3.0).patch_id
        for _ in 1:20
    ]
    vertex_hits = [
        intersect_terrain(slope, (0.0, 0.0, 2.0), (0, 0, -1), 3.0).patch_id
        for _ in 1:20
    ]
    @test all(==(1), diagonal_hits)
    @test all(==(1), vertex_hits)

    interpolated = HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        [0.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.0],
        soil;
        normal_mode=:interpolated,
    )
    @test intersect_terrain(interpolated, (0.8, 0.8, 2.0), (0, 0, -1), 3.0).normal[3] > 0

    vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    mesh = TriangulatedTerrain(vertices, [(1, 2, 3)], soil)
    @test intersect_terrain(mesh, (0.2, 0.2, 1.0), (0, 0, -1), 2.0).patch_id == 1
    @test intersect_terrain(mesh, (0.2, 0.2, 1.0), (0, 0, -1)).distance == 1.0
    @test_throws ArgumentError TriangulatedTerrain(vertices, [(1, 2, 3)], soil; periodic=true)
end

@testitem "Terrain-aware DDA stops before subterranean voxels and classifies boundaries" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    pad = reshape([10.0, 0.0, 0.0], 1, 1, 3)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), pad)
    terrain = PlanarTerrain(grid, soil; elevation=1.2)

    dda = trace_voxel_ray(
        grid,
        (0.5, 0.5, 3.0),
        (0, 0, -1);
        boundary=:open,
        algorithm=:dda,
        terrain=terrain,
    )
    reference = trace_voxel_ray(
        grid,
        (0.5, 0.5, 3.0),
        (0, 0, -1);
        boundary=:open,
        algorithm=:reference,
        terrain=terrain,
    )
    @test dda.exit_reason == reference.exit_reason == :terrain
    @test dda.terrain_hit.point == reference.terrain_hit.point == (0.5, 0.5, 1.2)
    @test [segment.index for segment in dda] == [(1, 1, 3), (1, 1, 2)]
    @test [segment.length for segment in dda] ≈ [1.0, 0.8] atol=1e-14
    @test [segment.index for segment in reference] == [segment.index for segment in dda]
    @test [segment.length for segment in reference] ≈ [segment.length for segment in dda] atol=1e-14
    @test all(segment.index[3] != 1 for segment in dda)

    backend = VoxelCPUBackend(rays_per_voxel=4, g=0.5, boundary=:open, traversal=:dda)
    response = ArchimedLight.compute_voxel_direction_response(grid, (0, 0, -1), backend, terrain)
    @test response.interception_fraction[1, 1, 1] == 0.0
    @test sum(response.interception_fraction) == 0.0
    @test response.terrain_incident_weight ≈ 1.0
    @test response.escaped_bottom_weight == 0.0

    open_path = trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, -1); boundary=:open)
    @test open_path.exit_reason == :bottom
    @test open_path.terrain_hit === nothing
    @test [segment.index for segment in open_path] == [(1, 1, 3), (1, 1, 2), (1, 1, 1)]

    side_grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    side_terrain = PlanarTerrain(side_grid, soil; elevation=0.5)
    side = trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (1.0, 0.0, -0.1);
        boundary=:open,
        terrain=side_terrain,
    )
    downward = trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0.0, 0.0, -1.0);
        boundary=:open,
        terrain=side_terrain,
    )
    @test side.exit_reason == :side
    @test downward.exit_reason == :terrain

    incomplete = PlanarTerrain(0.5, (0.0, 0.8), (0.0, 1.0), soil)
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0, 0, -1);
        boundary=:open,
        terrain=incomplete,
    )
    mesh_gap = TriangulatedTerrain(
        [(0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.0, 1.0, 0.5)],
        [(1, 2, 3)],
        soil,
    )
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.9, 0.9, 2.0),
        (0, 0, -1);
        boundary=:open,
        terrain=mesh_gap,
    )

    nonperiodic = PlanarTerrain(side_grid, soil; elevation=0.5, periodic=false)
    periodic = PlanarTerrain(side_grid, soil; elevation=0.5, periodic=true)
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0, 0, -1);
        boundary=:periodic,
        terrain=nonperiodic,
    )
    @test trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (0.5, 0, -1);
        boundary=:periodic,
        terrain=periodic,
    ).exit_reason == :terrain

    tiled_height = HeightFieldTerrain(
        [0.0, 1.0],
        [0.0, 1.0],
        zeros(2, 2),
        soil;
        periodic=true,
    )
    @test trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (0.5, 0, -1);
        boundary=:periodic,
        terrain=tiled_height,
    ).exit_reason == :terrain
end

@testitem "Local-normal Lambertian soil uses only the fixed full-sphere quadrature" tags=[:voxel, :terrain, :scattering, :fast] begin
    inverse = inv(sqrt(3.0))
    directions = [
        (sx * inverse, sy * inverse, sz * inverse)
        for sx in (-1.0, 1.0), sy in (-1.0, 1.0), sz in (-1.0, 1.0)
    ] |> vec
    quadrature = VoxelScatteringQuadrature(directions, fill(1 / 8, 8))

    horizontal_indices, horizontal_weights = terrain_lambertian_weights(
        quadrature,
        (0.0, 0.0, 1.0),
    )
    @test length(horizontal_indices) == 4
    @test all(directions[index][3] > 0 for index in horizontal_indices)
    @test all(isapprox(weight, 0.25; atol=1e-15) for weight in horizontal_weights)

    normals = [
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (-0.3, 0.8, 0.5),
        (0.9, -0.1, 0.2),
        (-0.7, -0.4, 0.6),
    ]
    for normal in normals
        indices, weights = terrain_lambertian_weights(quadrature, normal)
        unit = normal ./ sqrt(sum(abs2, normal))
        @test all(weights .>= 0)
        @test isapprox(sum(weights), 1.0; atol=2e-15)
        @test all(sum(unit[axis] * directions[index][axis] for axis in 1:3) > 0 for index in indices)
        @test all(index in eachindex(directions) for index in indices)
    end

    tilted_indices, _ = terrain_lambertian_weights(quadrature, (0.95, 0.0, 0.31))
    @test any(directions[index][3] < 0 for index in tilted_indices)
    @test Set(directions[tilted_indices]) ⊆ Set(directions)
end

@testitem "Terrain soil absorption, reflection, cache invalidation, and energy close" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), reshape([0.0, 1.0, 0.0], 1, 1, 3))
    backend = VoxelCPUBackend(rays_per_voxel=1, g=0.5, boundary=:open, traversal=:dda)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)],
        [0.5, 0.5],
    )
    emitted = reshape([0.0, 1.0, 0.0], 1, 1, 3)
    black_optics = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    white_optics = SoilOpticalProperties(par_reflectance=1.0, nir_reflectance=1.0)
    terrain = PlanarTerrain(
        0.75,
        (0.0, 1.0),
        (0.0, 1.0),
        black_optics;
        materials=Dict(1 => black_optics, 2 => black_optics),
    )
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain)
    @test sum(cache.terrain_direction_weights[1][2]) == 1.0
    black = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        band=:par,
        terrain=terrain,
        cache=cache,
    )
    @test soil_absorbed_energy(black) > 0
    @test soil_reflected_energy(black) == 0.0
    @test black.escaped_bottom_energy == 0.0
    @test black.unresolved_energy == 0.0
    @test abs(black.balance_residual) <= 1e-14

    # Reflectance is optical state, not geometry: the same path cache remains valid.
    terrain.materials[1] = white_optics
    white = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        band=:par,
        terrain=terrain,
        cache=cache,
    )
    @test soil_absorbed_energy(white) == 0.0
    @test soil_reflected_energy(white) > 0
    @test white.escaped_top_energy > black.escaped_top_energy
    @test abs(white.balance_residual) <= 1e-14

    # Material assignment is geometry/scene state and rejects an old cache.
    terrain.material_ids[1, 1] = 2
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        terrain=terrain,
        cache=cache,
    )

    subterranean_grid = VoxelGrid((0, 0, 0), (1, 1, 2), reshape([1.0, 0.0], 1, 1, 2))
    subterranean_terrain = PlanarTerrain(
        subterranean_grid,
        black_optics;
        elevation=1.1,
    )
    subterranean_emission = reshape([1.0, 0.0], 1, 1, 2)
    @test_throws ArgumentError apply_voxel_scattering_transport(
        subterranean_grid,
        subterranean_emission,
        quadrature,
        backend;
        terrain=subterranean_terrain,
    )
end

@testitem "AAIGrid and AMAPVox ground_distance form a validated scene contract" tags=[:voxel, :terrain, :io, :fast] begin
    using LinearAlgebra: I

    soil = SoilOpticalProperties(par_reflectance=0.15, nir_reflectance=0.35)
    mktempdir() do directory
        voxel_path = joinpath(directory, "scene.vox")
        open(voxel_path, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 2 2 2")
            println(io, "#split: 2 2 2")
            println(io, "#type: TLS")
            println(io, "i j k ground_distance PadBVTotal")
            for i in 0:1, j in 0:1, k in 0:1
                println(io, "$i $j $k $(k + 0.5) 0.0")
            end
        end

        dtm_path = joinpath(directory, "terrain.asc")
        open(dtm_path, "w") do io
            println(io, "ncols 1")
            println(io, "nrows 1")
            println(io, "xllcenter 1")
            println(io, "yllcenter 1")
            println(io, "cellsize 2")
            println(io, "NODATA_value -9999")
            println(io, "0")
        end
        transform = TerrainCoordinateTransform(
            Matrix{Float64}(I, 4, 4);
            source_crs="EPSG:2154",
            target_crs="EPSG:2154",
            source_units="m",
            target_units="m",
        )
        @test transform.source_units == transform.target_units == "m"
        input = VoxelTerrainSceneInput(
            voxel_path,
            dtm_path;
            configuration_path=joinpath(directory, "amapvox.xml"),
            dtm_to_voxel=transform,
            dtm_used_by_amapvox=true,
        )
        scene = read_voxel_terrain_scene(input, soil)
        @test size(scene.grid) == (2, 2, 2)
        @test terrain_bounds(scene.terrain) == ((0.0, 0.0, 0.0), (2.0, 2.0, 0.0))
        @test scene.alignment.sample_count == 8
        @test scene.alignment.maximum_absolute_error == 0.0
        @test read_voxel_column(voxel_path, "GROUND_DISTANCE")[:, :, 2] == fill(1.5, 2, 2)
        @test read_voxel_column(voxel_path, "missing"; required=false) === nothing

        translated_matrix = Matrix{Float64}(I, 4, 4)
        translated_matrix[1, 4] = 3.0
        translated_matrix[2, 4] = -2.0
        translated_matrix[3, 4] = 0.5
        translated = height_field_terrain(
            read_aai_grid(dtm_path),
            soil;
            transform=TerrainCoordinateTransform(translated_matrix),
        )
        @test terrain_bounds(translated) == ((3.0, -2.0, 0.5), (5.0, 0.0, 0.5))

        rotated_matrix = Matrix{Float64}(I, 4, 4)
        rotated_matrix[1:2, 1:2] .= [0.0 -1.0; 1.0 0.0]
        rotated = height_field_terrain(
            read_aai_grid(dtm_path),
            soil;
            transform=TerrainCoordinateTransform(rotated_matrix),
        )
        @test rotated isa TriangulatedTerrain
        @test_throws ArgumentError height_field_terrain(
            read_aai_grid(dtm_path),
            soil;
            transform=TerrainCoordinateTransform(rotated_matrix),
            periodic=true,
        )

        ground_distance = read_voxel_column(voxel_path, "ground_distance")
        fallback = terrain_from_ground_distance(
            scene.grid,
            ground_distance,
            soil;
            dtm_used_by_amapvox=true,
        )
        @test terrain_bounds(fallback) == terrain_bounds(scene.terrain)
        @test_throws ArgumentError terrain_from_ground_distance(
            scene.grid,
            ground_distance,
            soil;
            dtm_used_by_amapvox=false,
        )

        shifted_path = joinpath(directory, "shifted.asc")
        open(shifted_path, "w") do io
            println(io, "ncols 1")
            println(io, "nrows 1")
            println(io, "xllcorner 0.25")
            println(io, "yllcorner 0")
            println(io, "cellsize 2")
            println(io, "NODATA_value -9999")
            println(io, "0")
        end
        @test_throws ArgumentError read_voxel_terrain_scene(
            VoxelTerrainSceneInput(voxel_path, shifted_path),
            soil,
        )

        nodata_path = joinpath(directory, "nodata.asc")
        open(nodata_path, "w") do io
            println(io, "ncols 2")
            println(io, "nrows 1")
            println(io, "xllcorner 0")
            println(io, "yllcorner 0")
            println(io, "cellsize 1")
            println(io, "NODATA_value -9999")
            println(io, "0 -9999")
        end
        @test_throws ArgumentError read_aai_grid(nodata_path)

        varying_path = joinpath(directory, "varying.asc")
        open(varying_path, "w") do io
            println(io, "ncols 2")
            println(io, "nrows 1")
            println(io, "xllcorner 0")
            println(io, "yllcorner 0")
            println(io, "cellsize 1")
            println(io, "NODATA_value -9999")
            println(io, "0 1")
        end
        varying = height_field_terrain(read_aai_grid(varying_path), soil)
        @test terrain_elevation_at(varying, 0.5, 0.5) == 0.0
        @test terrain_elevation_at(varying, 1.5, 0.5) == 1.0

        inconsistent = copy(ground_distance)
        inconsistent[1, 1, 2] += 0.1
        @test_throws ArgumentError validate_ground_distance(scene.grid, scene.terrain, inconsistent)
        @test_throws ArgumentError terrain_from_ground_distance(
            scene.grid,
            inconsistent,
            soil;
            dtm_used_by_amapvox=true,
        )
    end
end

@testitem "Terrain pipeline closes complete PAR and NIR balances" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    backend = VoxelCPUBackend(rays_per_voxel=4, g=0.5, boundary=:open, traversal=:dda)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.5)
    terrain = PlanarTerrain(grid, soil; elevation=0.0)
    optics = VoxelOpticalProperties(
        grid;
        par_reflectance=0.08,
        par_transmittance=0.07,
        nir_reflectance=0.45,
        nir_transmittance=0.40,
    )
    row = (
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=true,
        scattering=true,
        scattering_max_iter=80,
        scattering_stop_ratio=1e-10,
        nir_interception=true,
        nir_scattering=true,
        radiation_timestep_minutes=5.0,
    )
    result = run_voxel_light_step(
        grid,
        row,
        options;
        backend=backend,
        optics=optics,
        terrain=terrain,
    )
    @test result.first_order.par_escaped_bottom_energy == 0.0
    @test result.first_order.nir_escaped_bottom_energy == 0.0
    @test sum(result.first_order.par_terrain_energy) ==
          result.first_order.par_terrain_incident_energy
    @test sum(result.first_order.nir_terrain_energy) ==
          result.first_order.nir_terrain_incident_energy

    for (band, incoming, injected) in (
        (result.scattering.par, result.first_order.par_incoming_energy, result.first_order.par_injected_energy),
        (result.scattering.nir, result.first_order.nir_incoming_energy, result.first_order.nir_injected_energy),
    )
        accounted = sum(band.absorbed_energy) + soil_absorbed_energy(band) +
                    band.escaped_top_energy + band.escaped_side_energy +
                    band.escaped_bottom_energy + band.unresolved_energy
        @test isapprox(accounted, incoming + injected; atol=1e-7, rtol=1e-11)
        @test band.escaped_bottom_energy == 0.0
        @test band.absolute_balance_residual <= 1e-7
        @test soil_reflected_energy(band) >= 0.0
    end

    @test soil_reflected_energy(result.scattering.nir) >
          soil_reflected_energy(result.scattering.par)
end
