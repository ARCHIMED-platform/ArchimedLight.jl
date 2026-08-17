@testitem "Voxel reference and DDA traversal agree" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))

    function signature(path)
        [(segment.index, segment.length) for segment in path]
    end

    vertical_reference = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:reference
    )
    vertical_dda = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:dda
    )
    @test [segment.index for segment in vertical_dda] == [(1, 1, 3), (1, 1, 2), (1, 1, 1)]
    @test [segment.length for segment in vertical_dda] == [1.0, 1.0, 1.0]
    @test signature(vertical_reference) == signature(vertical_dda)
    @test vertical_dda.exit_reason == :bottom

    cases = [
        ((0.25, 0.5, 3.0), (0.5, 0.0, -1.0)),
        ((0.5, 0.5, 3.0), (1.0, 1.0, -1.0)),
        ((1.5, 0.25, 3.0), (-0.7, 0.3, -0.4)),
        ((1.0, 1.0, 3.0), (-1.0, -1.0, -1.0)),
    ]
    for (origin, direction) in cases
        reference = trace_voxel_ray(grid, origin, direction; algorithm=:reference)
        dda = trace_voxel_ray(grid, origin, direction; algorithm=:dda)
        @test length(reference) == length(dda)
        @test [segment.index for segment in reference] == [segment.index for segment in dda]
        @test all(
            isapprox(a.length, b.length; atol=1e-12, rtol=1e-12)
            for (a, b) in zip(reference, dda)
        )
        @test all(
            reference.segments[index].index != reference.segments[index + 1].index
            for index in 1:(length(reference) - 1)
        )
    end

    diagonal = trace_voxel_ray(
        grid, (0.25, 0.5, 3.0), (0.5, 0.0, -1.0); algorithm=:dda
    )
    @test [segment.index for segment in diagonal] ==
          [(1, 1, 3), (1, 1, 2), (2, 1, 2), (2, 1, 1)]
    @test isapprox(sum(segment.length for segment in diagonal), 3sqrt(1.25); atol=1e-12)

    open_path = trace_voxel_ray(
        grid, (1.75, 0.5, 3.0), (1.0, 0.0, -1.0); boundary=:open
    )
    @test open_path.exit_reason == :side
    @test only(open_path.segments).index == (2, 1, 3)
    @test isapprox(only(open_path.segments).length, sqrt(0.125); atol=1e-12)

    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, 0))
    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, 1))
    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (1, 0, 0))
end

@testitem "Voxel DDA matches the reference on seeded rays" tags=[:voxel, :fast] begin
    using Random

    rng = MersenneTwister(0x564f58454c)
    grid = VoxelGrid((-1.5, 0.25, -2.0), (3.5, 4.25, 2.5), fill(0.2, 5, 4, 6))

    for boundary in (:periodic, :open), _ in 1:128
        origin = (
            grid.minimum[1] + rand(rng) * (grid.maximum[1] - grid.minimum[1]),
            grid.minimum[2] + rand(rng) * (grid.maximum[2] - grid.minimum[2]),
            grid.maximum[3],
        )
        direction = (4rand(rng) - 2, 4rand(rng) - 2, -(0.05 + 1.95rand(rng)))
        reference = trace_voxel_ray(
            grid, origin, direction; boundary=boundary, algorithm=:reference
        )
        dda = trace_voxel_ray(grid, origin, direction; boundary=boundary, algorithm=:dda)

        @test dda.exit_reason == reference.exit_reason
        @test [segment.index for segment in dda] == [segment.index for segment in reference]
        @test length(dda) == length(reference)
        @test all(
            isapprox(actual.length, expected.length; atol=2e-14, rtol=2e-13)
            for (actual, expected) in zip(dda, reference)
        )
        @test all(segment.length > 0 for segment in dda)
        @test all(
            dda.segments[index].index != dda.segments[index + 1].index
            for index in 1:(length(dda) - 1)
        )
    end
end
