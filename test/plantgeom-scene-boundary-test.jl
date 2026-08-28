@testitem "PlantGeom OPS boundary stays visualization-only" tags=[:core, :fast] begin
    using CairoMakie

    import GeometryBasics
    import MultiScaleTreeGraph
    import PlantGeom

    fixture_dir = joinpath(
        @__DIR__,
        "fast_fixtures",
        "simpleplant_16_notoric",
        "input",
    )
    options, scene, meteo, models = ArchimedLight.read_config(
        joinpath(fixture_dir, "config.yml"),
    )

    mesh_before = scene.merged_mesh
    vertices_before = collect(GeometryBasics.coordinates(scene.merged_mesh))
    faces_before = collect(GeometryBasics.faces(scene.merged_mesh))
    face2node_before = copy(scene.face2node)
    node_snapshot(scene) = Dict(
        id => (
            area=node.area,
            barycenter=node.barycenter,
            source_topology_id=node.source_topology_id,
            source_owner=node.source_owner,
        ) for (id, node) in scene.nodes
    )
    nodes_before = node_snapshot(scene)
    source_node_ids_before = MultiScaleTreeGraph.node_id.(
        MultiScaleTreeGraph.traverse(scene.mtg, identity),
    )
    areas_before = PlantGeom.node_areas(scene)
    summary_before = sprint(show, ArchimedLight.summarize_scene(scene))
    @test occursin("pavement", summary_before)
    render_before = ArchimedLight.light_render_geometry(scene, models, options)
    render_snapshot(render) = (
        vertices=copy(render.vertices),
        faces=copy(render.faces),
        face2node=copy(render.face2node),
    )
    rendered_before = render_snapshot(render_before)
    radiative_areas_before = ArchimedLight._interception_area_per_node_local(
        scene,
        models,
        options,
    )
    version_before = PlantGeom.scene_version(scene)

    boundary = @test_nowarn PlantGeom.scene_boundary_mesh(scene.mtg)
    @test length(GeometryBasics.faces(boundary)) == 2
    fig, ax, plot = PlantGeom.plantviz(
        scene.mtg;
        show_scene_boundary=true,
        cache=false,
    )
    @test plot.plots[end] isa Makie.Lines
    mktempdir() do tmp
        @test_nowarn CairoMakie.save(joinpath(tmp, "scene-boundary.png"), fig)
    end

    @test scene.merged_mesh === mesh_before
    @test collect(GeometryBasics.coordinates(scene.merged_mesh)) == vertices_before
    @test collect(GeometryBasics.faces(scene.merged_mesh)) == faces_before
    @test scene.face2node == face2node_before
    @test node_snapshot(scene) == nodes_before
    @test MultiScaleTreeGraph.node_id.(
        MultiScaleTreeGraph.traverse(scene.mtg, identity),
    ) == source_node_ids_before
    @test PlantGeom.node_areas(scene) == areas_before
    @test sprint(show, ArchimedLight.summarize_scene(scene)) == summary_before
    @test render_snapshot(
        ArchimedLight.light_render_geometry(scene, models, options),
    ) == rendered_before
    @test ArchimedLight._interception_area_per_node_local(scene, models, options) ==
          radiative_areas_before
    @test PlantGeom.scene_version(scene) == version_before
    @test !PlantGeom.has_geometry(scene.mtg)
end
