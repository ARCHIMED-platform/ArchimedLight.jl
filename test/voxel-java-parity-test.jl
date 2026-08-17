@testitem "Voxel Java baseline provenance and intentional traversal divergence" tags=[:voxel, :java_fixture, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "java_reference")
    @test isfile(joinpath(fixture, "source_sha256.txt"))
    @test isfile(joinpath(fixture, "baseline.tsv"))
    @test count(line -> !isempty(line), readlines(joinpath(fixture, "source_sha256.txt"))) == 9
    baseline = [split(line, '\t') for line in readlines(joinpath(fixture, "baseline.tsv"))[2:end]]

    function path_baseline(name)
        [
            (index=parse(Int, row[3]), length=parse(Float64, row[4])) for row in baseline if
            row[1] == "SEGMENT" && row[2] == name
        ]
    end

    function voxel_baseline(name)
        [
            (
                index=parse(Int, row[3]),
                coefficient=parse(Float64, row[4]),
                intercepted=parse(Float64, row[5]),
            ) for row in baseline if row[1] == "VOXEL" && row[2] == name
        ]
    end

    function total_baseline(name)
        row = only(filter(row -> row[1] == "TOTAL" && row[2] == name, baseline))
        return (intercepted=parse(Float64, row[3]), ray_count=parse(Int, row[4]))
    end

    function profile_baseline(name)
        [
            (
                layer=parse(Int, row[3]),
                coefficient=parse(Float64, row[4]),
                intercepted=parse(Float64, row[5]),
            ) for row in baseline if row[1] == "PROFILE" && row[2] == name
        ]
    end

    function layer_coefficients(response)
        [
            sum(
                response.interception_fraction[i, j, k] for
                i in axes(response.interception_fraction, 1),
                j in axes(response.interception_fraction, 2)
            ) for k in axes(response.interception_fraction, 3)
        ]
    end

    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))
    vertical = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:reference
    )
    @test [ArchimedLight.voxel_linear_index(grid, segment.index...) for segment in vertical] == [3, 2, 1]
    @test [s.length for s in vertical] == [1.0, 1.0, 1.0]

    java_vertical = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:java_reference
    )
    expected_vertical = path_baseline("vertical")
    @test [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in java_vertical] ==
          [item.index for item in expected_vertical]
    @test [segment.length for segment in java_vertical] ==
          [item.length for item in expected_vertical]

    # The Java diagonal fixture repeats voxel index 6. The corrected Julia
    # contract crosses an actual boundary at every segment transition.
    diagonal = trace_voxel_ray(
        grid, (0.25, 0.5, 3.0), (0.5, 0, -1); algorithm=:reference
    )
    julia_indices = [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in diagonal]
    @test julia_indices == [2, 1, 7, 6]
    @test julia_indices != [2, 7, 6, 6]

    java_diagonal = trace_voxel_ray(
        grid, (0.25, 0.5, 3.0), (0.5, 0, -1); algorithm=:java_reference
    )
    expected_diagonal = path_baseline("diagonal")
    @test [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in java_diagonal] ==
          [item.index for item in expected_diagonal]
    @test [segment.length for segment in java_diagonal] ==
          [item.length for item in expected_diagonal]

    java_corner = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (1, 1, -1); algorithm=:java_reference
    )
    expected_corner = path_baseline("corner")
    @test [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in java_corner] ==
          [item.index for item in expected_corner]
    @test [segment.length for segment in java_corner] ==
          [item.length for item in expected_corner]

    java_grid = VoxelGrid(
        (0, 0, 0),
        (3, 2, 3),
        fill(Float64(Float32(0.4)), 3, 2, 3),
    )
    java_vertical_response = ArchimedLight.compute_voxel_direction_response(
        java_grid,
        (0, 0, -1),
        VoxelCPUBackend(rays_per_voxel=4, traversal=:java_reference, boundary=:periodic),
    )
    java_toric = ArchimedLight.compute_voxel_direction_response(
        java_grid,
        (0.5, 0.2, -1),
        VoxelCPUBackend(rays_per_voxel=4, traversal=:java_reference, boundary=:periodic),
    )
    java_nontoric = ArchimedLight.compute_voxel_direction_response(
        java_grid,
        (0.5, 0.2, -1),
        VoxelCPUBackend(rays_per_voxel=4, traversal=:java_reference, boundary=:java_nontoric),
    )
    vertical_voxels = voxel_baseline("vertical_toric")
    toric_voxels = voxel_baseline("diagonal_toric")
    nontoric_voxels = voxel_baseline("diagonal_java_nontoric")
    vertical_total = total_baseline("vertical_toric")
    toric_total = total_baseline("diagonal_toric")
    nontoric_total = total_baseline("diagonal_java_nontoric")
    @test vec(permutedims(java_vertical_response.interception_fraction, (3, 2, 1))) ==
          [item.coefficient for item in vertical_voxels]
    @test vec(permutedims(java_toric.interception_fraction, (3, 2, 1))) ==
          [item.coefficient for item in toric_voxels]
    @test vec(permutedims(java_nontoric.interception_fraction, (3, 2, 1))) ==
          [item.coefficient for item in nontoric_voxels]
    @test sum(java_vertical_response.interception_fraction) * 2 == vertical_total.intercepted
    @test sum(java_toric.interception_fraction) * 2 == toric_total.intercepted
    @test sum(java_nontoric.interception_fraction) * 2 == nontoric_total.intercepted
    @test java_toric.effective_ray_count == toric_total.ray_count
    @test java_nontoric.effective_ray_count == nontoric_total.ray_count
    @test java_vertical_response.effective_ray_count == vertical_total.ray_count
    @test [java_toric.interception_fraction[1, 1, k] for k in 1:3] ==
          [0.1290108859539032, 0.16191202402114868, 0.20320379734039307]
    @test [java_nontoric.interception_fraction[1, 1, k] for k in 1:3] ==
          [0.20320382714271545, 0.16191202402114868, 0.20320379734039307]

    ground_pad = zeros(5, 5, 2)
    ground_pad[:, :, 1] .= 0.5
    ground = VoxelGrid((0, 0, 0), (5, 5, 2), ground_pad)
    java_ground = ArchimedLight.compute_voxel_direction_response(
        ground,
        (0.5, 0.2, -1),
        VoxelCPUBackend(rays_per_voxel=4, traversal=:java_reference, boundary=:periodic),
    )
    ground_profile = profile_baseline("ground_toric")
    ground_total = total_baseline("ground_toric")
    @test layer_coefficients(java_ground) == [item.coefficient for item in ground_profile]
    @test sum(java_ground.interception_fraction) * 2 == ground_total.intercepted
    @test java_ground.effective_ray_count == ground_total.ray_count

    pillar_pad = zeros(5, 5, 4)
    pillar_pad[3, 3, :] .= 0.5
    pillar = VoxelGrid((0, 0, 0), (5, 5, 4), pillar_pad)
    java_pillar = ArchimedLight.compute_voxel_direction_response(
        pillar,
        (0.05, 0.0, -1),
        VoxelCPUBackend(
            rays_per_voxel=9,
            traversal=:java_reference,
            boundary=:java_nontoric,
        ),
    )
    pillar_profile = profile_baseline("pillar_java_nontoric")
    pillar_total = total_baseline("pillar_java_nontoric")
    @test layer_coefficients(java_pillar) == [item.coefficient for item in pillar_profile]
    @test sum(java_pillar.interception_fraction) * 2 == pillar_total.intercepted
    @test java_pillar.effective_ray_count == pillar_total.ray_count
end
