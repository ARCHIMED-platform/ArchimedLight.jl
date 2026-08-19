@testitem "Makie voxelplot renders terrain-aware voxel transport" tags=[:core, :voxel, :terrain, :fast] begin
    using CairoMakie

    pad = reshape(
        [
            0.0, 0.3, 0.6, 0.9,
            0.2, 0.5, 0.8, 1.1,
            0.0, 0.4, 0.7, 1.0,
            0.3, 0.6, 0.9, 1.2,
            0.1, 0.4, 0.8, 1.1,
            0.2, 0.5, 0.7, 1.0,
        ],
        2,
        3,
        4,
    )
    grid = VoxelGrid((0, 0, 0), (2, 3, 4), pad)
    original_pad = copy(grid.pad)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    terrain = HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.5, 3.0],
        [
            0.40 0.50 0.30
            0.60 0.90 0.50
            0.20 0.40 0.70
        ],
        soil,
    )
    upward = [
        (0.0, 0.0, 1.0),
        (0.8, 0.0, 0.6),
        (-0.8, 0.0, 0.6),
        (0.0, 0.8, 0.6),
        (0.0, -0.8, 0.6),
    ]
    directions = vcat(upward, [(-d[1], -d[2], -d[3]) for d in upward])
    quadrature = VoxelScatteringQuadrature(directions, fill(0.1, length(directions)))
    rays = [
        (origin=(0.5, 0.5, 4.0), direction=(0.05, 0.10, -1.0)),
        (origin=(1.5, 2.2, 4.0), direction=(-0.10, -0.05, -1.0)),
    ]

    figure, axis, plot = voxelplot(
        grid,
        terrain;
        rays=rays,
        quadrature=quadrature,
        cutaway=:front,
        show_normals=false,
    )

    @test axis isa LScene
    @test length(plot.plots) == 9
    @test plot.plots[1] isa Makie.Voxels
    @test plot.plots[2] isa Makie.Voxels
    @test plot.plots[3] isa Makie.LineSegments
    @test plot.plots[4] isa Makie.Mesh
    @test plot.plots[5] isa Makie.Wireframe
    @test plot.plots[6] isa Makie.LineSegments
    @test plot.plots[7] isa Makie.Scatter
    @test plot.plots[8] isa Makie.Arrows3D
    @test plot.plots[9] isa Makie.Arrows3D
    @test !Makie.to_value(plot.plots[1].transparency)
    @test !Makie.to_value(plot.plots[2].transparency)
    @test !Makie.to_value(plot.plots[4].transparency)
    @test !Makie.to_value(plot.plots[1].visible)
    @test Makie.to_value(plot.plots[2].visible)

    paths = Makie.to_value(plot[:voxel_ray_paths])
    hit_points = Makie.to_value(plot[:voxel_hit_points])
    reflected_weights = Makie.to_value(plot[:voxel_reflected_weights])
    @test length(paths) == length(rays)
    @test all(path -> path.exit_reason == :terrain, paths)
    @test length(hit_points) == length(rays)
    @test sum(reflected_weights) ≈ length(rays)
    @test grid.pad == original_pad

    Makie.update!(plot, cutaway=:none, show_normals=true)
    @test Makie.to_value(plot.plots[1].visible)
    @test !Makie.to_value(plot.plots[2].visible)
    @test Makie.to_value(plot.plots[9].visible)
    @test grid.pad == original_pad

    inplace_figure = Figure()
    inplace_axis = LScene(inplace_figure[1, 1])
    inplace_plot = voxelplot!(inplace_axis, grid, terrain; rays=rays)
    @test length(inplace_plot.plots) == 9
    @test length(Makie.to_value(inplace_plot[:voxel_ray_paths])) == length(rays)

    transparent_figure = Figure()
    transparent_axis = LScene(transparent_figure[1, 1])
    transparent_plot = voxelplot!(
        transparent_axis,
        grid,
        terrain;
        transparency=true,
    )
    @test Makie.to_value(transparent_plot.plots[1].transparency)
    @test Makie.to_value(transparent_plot.plots[2].transparency)
    @test Makie.to_value(transparent_plot.plots[4].transparency)

    no_terrain_figure, _, no_terrain_plot = voxelplot(
        grid;
        boundary=:periodic,
        rays=[(origin=(1.8, 1.5, 4.0), direction=(1.0, 0.0, -1.0))],
    )
    @test isempty(Makie.to_value(no_terrain_plot[:voxel_hit_points]))
    @test isempty(Makie.to_value(no_terrain_plot[:voxel_reflected_weights]))
    periodic_segments = Makie.to_value(no_terrain_plot[:voxel_ray_segments])
    @test all(
        point -> grid.minimum[1] - 1.0f-6 <= point[1] <= grid.maximum[1] + 1.0f-6,
        periodic_segments,
    )
    @test all(
        pair -> abs(periodic_segments[pair][1] - periodic_segments[pair + 1][1]) <=
                grid.voxel_size[1] + 1.0f-6,
        1:2:(length(periodic_segments) - 1),
    )

    @test_throws ArgumentError voxelplot(grid, terrain; cutaway=:side)
    @test_throws ArgumentError voxelplot(
        grid,
        terrain;
        rays=[(origin=(0.5, 0.5, 4.0), direction=(0.0, 0.0, 1.0))],
    )
    @test_throws ArgumentError voxelplot(
        grid,
        terrain;
        rays=[(origin=(3.0, 0.5, 4.0), direction=(0.0, 0.0, -1.0))],
    )

    mktemp() do path, io
        close(io)
        png_path = string(path, ".png")
        save(png_path, figure)
        @test isfile(png_path)
        @test filesize(png_path) > 0
    end

    CairoMakie.empty!(figure)
    CairoMakie.empty!(inplace_figure)
    CairoMakie.empty!(transparent_figure)
    CairoMakie.empty!(no_terrain_figure)
end
