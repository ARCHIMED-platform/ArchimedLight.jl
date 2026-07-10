@testmodule QueryHelper begin
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    using GeometryBasics
    using Dates
    using ArchimedLight
    using JET
    import MultiScaleTreeGraph
    import PlantGeom

    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))
end

@testitem "Scene-aware light value queries" tags = [:query, :fast] setup = [QueryHelper] begin
    using ArchimedLight
    using JET
    import MultiScaleTreeGraph
    import Tables

    scene = QueryHelper._synthetic_horizontal_scene([
        (x0=0.0, x1=0.5, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Leaf", object_id=1, source_topology_id=101),
        (x0=0.5, x1=1.0, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Internode", object_id=1, source_topology_id=102),
        (x0=1.0, x1=1.5, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Leaf", object_id=2, source_topology_id=103),
        (x0=1.5, x1=2.0, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Internode", object_id=2, source_topology_id=104),
        (x0=2.0, x1=2.5, y0=0.0, y1=0.5, z=1.0, group="cacao", type="Leaf", object_id=3, source_topology_id=105),
        (x0=2.5, x1=3.0, y0=0.0, y1=0.5, z=0.0, group="pavement", type="Ground", object_id=-1, source_topology_id=106),
    ]; root_id=100)
    MultiScaleTreeGraph.attribute!(scene.mtg, :site, "field")
    MultiScaleTreeGraph.attribute!(MultiScaleTreeGraph.get_node(scene.mtg, 1), :rank, 1)

    options = LightOptions(
        QueryHelper._synthetic_options(
            sectors=1,
            all_in_turtle=false,
            scattering=false,
            toricity=false,
        );
        node_metadata_attributes=(:site, :rank),
    )
    models = QueryHelper._default_synthetic_models()
    sim = LightSimulation(scene, models; options=options)
    sky1 = SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    sky2 = SkyState(180.0, 90.0, 200.0, 0.0, 1.0, 0.0)
    step1 = run_light(sim, sky1; step_duration_seconds=60.0)
    step2 = run_light(sim, sky2; step_duration_seconds=60.0)
    series = [step1, step2]

    @test light_node_ids(sim; species="coffee") == [1, 2, 3, 4]
    @test light_node_ids(sim; group="coffee", object_id=1) == [1, 2]
    @test light_node_ids(sim; object_id=2, source_topology_id=103) == [3]
    @test light_node_ids(sim; species="coffee", object_id=[2]) == [3, 4]
    @test light_node_ids(sim; species="coffee", symbol=:Leaf) == [1, 3]
    @test light_node_ids(sim; species="coffee", type="Internode") == [2, 4]
    @test light_node_ids(sim; scale=1) == collect(1:6)
    @test light_node_ids(sim; node_ids=(4, 2)) == [2, 4]
    @test light_node_ids(sim; attributes=(rank=1,)) == [1]
    @test isempty(light_node_ids(sim; species="coffee", attributes=(site="field",)))
    @test light_node_ids(
        sim;
        species="coffee",
        attributes=(site="field",),
        inherit_attributes=true,
    ) == [1, 2, 3, 4]
    @test light_node_ids(
        sim;
        where=node -> isodd(MultiScaleTreeGraph.node_id(node)),
    ) == [1, 3, 5]

    @test_throws ArgumentError light_node_ids(sim; group="coffee", species="coffee")
    @test_throws ArgumentError light_node_ids(sim; node_ids=99)
    @test_throws ArgumentError light_node_ids(sim; attributes=(does_not_exist=1,))

    raw = light_metric_values(step1, :Ra_PAR_q)
    @test raw isa AbstractDict
    @test step1.node_metadata !== nothing
    @test step1.node_metadata === step2.node_metadata
    @test cache_summary(sim).node_metadata_count == 6
    @test cache_summary(sim).node_metadata_bytes > 0
    @test light_node_ids(step1; object_id=2, source_topology_id=103) == [3]

    leaves = light_metric_values(
        sim,
        step1,
        :Ra_PAR_q;
        species="coffee",
        symbol=:Leaf,
    )
    expected_columns = (
        :step_number,
        :node_id,
        :source_topology_id,
        :object_id,
        :item_id,
        :component_id,
        :group,
        :type,
        :symbol,
        :scale,
        :value,
    )
    @test propertynames(leaves) == expected_columns
    @test leaves.node_id == [1, 3]
    @test leaves.step_number == [1, 1]
    @test leaves.object_id == [1, 2]
    @test leaves.group == ["coffee", "coffee"]
    @test leaves.symbol == [:Leaf, :Leaf]
    @test leaves.value == [raw[1], raw[3]]
    @test eltype(leaves.value) == Float64
    @test Tables.istable(typeof(leaves))

    inherited = light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        species="coffee",
        attributes=(site="field",),
        inherit_attributes=true,
    )
    @test inherited.site == fill("field", 4)
    @test eltype(inherited.site) <: AbstractString

    leaf_series = light_metric_values(
        sim,
        series,
        :absorbed_par_energy;
        species="coffee",
        symbol=:Leaf,
    )
    @test leaf_series.node_id == [1, 3, 1, 3]
    @test leaf_series.step_number == [1, 1, 2, 2]
    @test leaf_series.value[3:4] ≈ 2 .* leaf_series.value[1:2]

    expected_plant_2 = sum(raw[nid] for nid in (3, 4))
    @test light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        species="coffee",
        object_id=2,
        reduce=sum,
    ) ≈ expected_plant_2

    totals = light_metric_values(
        sim,
        series,
        :absorbed_par_energy;
        species="coffee",
        reduce=sum,
    )
    @test totals[1] ≈ sum(raw[nid] for nid in 1:4)
    @test totals[2] ≈ 2 * totals[1]
    @test light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        species="missing",
        reduce=sum,
    ) == 0.0
    @test_throws ArgumentError light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        species="missing",
        reduce=maximum,
    )

    grouped = light_metric_values(
        sim,
        series,
        :absorbed_par_energy;
        species="coffee",
        reduce=sum,
        by=:object_id,
    )
    @test grouped.step_number == [1, 1, 2, 2]
    @test grouped.object_id == [1, 2, 1, 2]
    @test grouped.value[3:4] ≈ 2 .* grouped.value[1:2]
    @test Tables.istable(typeof(grouped))

    row_sink = light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        node_ids=1,
        sink=Tables.rowtable,
    )
    @test length(row_sink) == 1
    @test only(row_sink).node_id == 1

    @test_throws ArgumentError light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        by=:group,
    )
    @test_throws ArgumentError light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        reduce=sum,
        sink=Tables.rowtable,
    )

    dynamic_scene = QueryHelper._synthetic_horizontal_scene([
        (x0=0.0, x1=0.5, y0=0.0, y1=0.5, z=1.0, group="cacao", type="Leaf", object_id=3, source_topology_id=105),
        (x0=0.5, x1=1.0, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Leaf", object_id=2, source_topology_id=103),
        (x0=1.0, x1=1.5, y0=0.0, y1=0.5, z=1.0, group="coffee", type="Internode", object_id=2, source_topology_id=104),
        (x0=1.5, x1=2.0, y0=0.0, y1=0.5, z=0.0, group="pavement", type="Ground", object_id=-1, source_topology_id=106),
    ]; root_id=200)
    MultiScaleTreeGraph.attribute!(dynamic_scene.mtg, :site, "field")
    MultiScaleTreeGraph.attribute!(MultiScaleTreeGraph.get_node(dynamic_scene.mtg, 2), :rank, 1)
    update_scene!(sim, dynamic_scene)
    dynamic_step = run_light(sim, sky1; step_duration_seconds=60.0)
    @test dynamic_step.node_metadata !== step1.node_metadata
    @test light_node_ids(step1; species="coffee", object_id=2) == [3, 4]
    @test light_node_ids(dynamic_step; species="coffee", object_id=2) == [2, 3]

    dynamic_series = light_metric_values(
        sim,
        [step1, dynamic_step],
        :absorbed_par_energy;
        species="coffee",
        object_id=2,
    )
    @test dynamic_series.step_number == [1, 1, 2, 2]
    @test dynamic_series.node_id == [3, 4, 2, 3]
    @test dynamic_series.source_topology_id == [103, 104, 103, 104]
    @test dynamic_series.value[1:2] == [raw[3], raw[4]]
    dynamic_raw = light_metric_values(dynamic_step, :absorbed_par_energy)
    @test dynamic_series.value[3:4] == [dynamic_raw[2], dynamic_raw[3]]

    stable_leaf_series = light_metric_values(
        sim,
        [step1, dynamic_step],
        :absorbed_par_energy;
        object_id=2,
        source_topology_id=103,
    )
    @test stable_leaf_series.node_id == [3, 2]
    @test stable_leaf_series.source_topology_id == [103, 103]

    disabled_sim = LightSimulation(
        dynamic_scene,
        models;
        options=LightOptions(options; store_node_metadata=false),
    )
    disabled_step = run_light(disabled_sim, sky1; step_duration_seconds=60.0)
    @test disabled_step.node_metadata === nothing
    @test light_metric_values(
        disabled_sim,
        disabled_step,
        :absorbed_par_energy;
        species="coffee",
        reduce=sum,
    ) ≈ sum(dynamic_raw[nid] for nid in (2, 3))

    default_sim = LightSimulation(dynamic_scene, models; options=QueryHelper._synthetic_options())
    default_step = run_light(default_sim, sky1; step_duration_seconds=60.0)
    @test default_step.node_metadata !== nothing
    @test isempty(propertynames(default_step.node_metadata.attributes))

    JET.@test_call target_modules=(ArchimedLight,) light_node_ids(sim; species="coffee")
    JET.@test_call target_modules=(ArchimedLight,) light_metric_values(
        sim,
        step1,
        :absorbed_par_energy;
        node_ids=[1, 2],
    )
end
