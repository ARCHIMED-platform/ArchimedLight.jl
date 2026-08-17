@testitem "Voxel grid indexing, validation, and historical I/O" tags=[:voxel, :fast] begin
    pad = reshape(collect(0.0:0.1:1.1), 2, 2, 3)
    grid = VoxelGrid((10.0, 20.0, 2.0), (14.0, 22.0, 8.0), pad)

    @test size(grid) == (2, 2, 3)
    @test grid.voxel_size == (2.0, 1.0, 2.0)
    @test ArchimedLight.voxel_volume(grid) == 4.0
    @test ArchimedLight.voxel_horizontal_area(grid) == 2.0
    @test ArchimedLight.voxel_linear_index(grid, 1, 1, 1) == 1
    @test ArchimedLight.voxel_linear_index(grid, 1, 1, 3) == 3
    @test ArchimedLight.voxel_linear_index(grid, 1, 2, 1) == 4
    @test ArchimedLight.voxel_linear_index(grid, 2, 1, 1) == 7
    @test all(
        ArchimedLight.voxel_cartesian_index(grid, ArchimedLight.voxel_linear_index(grid, index.I...)) == index.I
        for index in CartesianIndices(grid.pad)
    )

    @test_throws ArgumentError VoxelGrid((0, 0, 0), (0, 1, 1), zeros(1, 1, 1))
    @test_throws ArgumentError VoxelGrid((0, 0, 0), (1, 1, 1), fill(-1.0, 1, 1, 1))
    @test_throws ArgumentError VoxelGrid((0, 0, 0), (1, 1, 1), fill(NaN, 1, 1, 1))

    mktempdir() do directory
        path = joinpath(directory, "roundtrip.vox")
        write_voxel_grid(path, grid)
        restored = read_voxel_grid(path)
        @test restored.minimum == grid.minimum
        @test restored.maximum == grid.maximum
        @test restored.voxel_size == grid.voxel_size
        @test restored.pad == grid.pad

        reordered = joinpath(directory, "reordered.vox")
        open(reordered, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 1 1 1")
            println(io, "#split: 1 1 1")
            println(io, "#type: TLS")
            println(io, "i j k ignored PadBVTotal")
            println(io, "0 0 0 99 0.75")
        end
        @test read_voxel_grid(reordered).pad[1, 1, 1] == 0.75

        duplicate = joinpath(directory, "duplicate.vox")
        open(duplicate, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 1 1 1")
            println(io, "#split: 1 1 1")
            println(io, "i j k PAD")
            println(io, "0 0 0 0.5")
            println(io, "0 0 0 0.5")
        end
        @test_throws ArgumentError read_voxel_grid(duplicate)

        missing = joinpath(directory, "missing.vox")
        open(missing, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 2 1 1")
            println(io, "#split: 2 1 1")
            println(io, "i j k PAD extra")
            println(io, "0 0 0 0.5 ignored")
        end
        @test_throws ArgumentError read_voxel_grid(missing)

        invalid_index = joinpath(directory, "invalid-index.vox")
        open(invalid_index, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 1 1 1")
            println(io, "#split: 1 1 1")
            println(io, "i j k PAD")
            println(io, "1 0 0 0.5")
        end
        @test_throws ArgumentError read_voxel_grid(invalid_index)
    end
end
