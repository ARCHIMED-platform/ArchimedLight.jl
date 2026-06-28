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

@testitem "Synthetic case Raycore scene adapter metadata" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)

    data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )

    @test data.prepared.geometry.face2node_index == [1, 1]
    node_idx = ArchimedLight._raycore_closest_hit_node_index(data, (0.5, 0.5, 0.0), (0.0, 0.0, 1.0))
    @test node_idx == 1
    @test data.prepared.geometry.node_ids[node_idx] == 1

    traced = ArchimedLight._raycore_trace_top_hits(
        data,
        [(0.5, 0.5, 2.0), (2.0, 2.0, 2.0)],
        [(0.0, 0.0, -1.0), (0.0, 0.0, -1.0)];
        t_maxs=[3.0f0, 3.0f0],
    )
    @test traced.nodes == UInt32[1, 0]
    @test traced.heights[1] ≈ 1.0f0
    @test isinf(traced.distances[2])
end

@testitem "Synthetic case Raycore first-order top-hit" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first_cpu = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    first_raycore = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=:raycore_cpu)

    @test get(first_raycore.projected_area_per_node, 1, 0.0) ≈ get(first_cpu.projected_area_per_node, 1, 0.0)
    @test get(first_raycore.incident_power.par, 1, 0.0) ≈ get(first_cpu.incident_power.par, 1, 0.0)
    @test get(first_raycore.hits_per_node, 1, 0) == get(first_cpu.hits_per_node, 1, 0)

    stacked = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    stacked_first = ArchimedLight.compute_first_order(stacked, models, turtle, fluxes, options; backend=:raycore_cpu)
    @test isapprox(get(stacked_first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(stacked_first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(stacked_first.incident_power.par, 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stacked_first.incident_power.par, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
end

@testitem "Synthetic case Raycore full stack first-order" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    sensor_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="plant", type="plate", object_id=2),
    ])
    sensor_models = HelperModule._virtual_sensor_models()
    sensor_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    sensor_turtle = ArchimedLight.build_turtle(sensor_options, sky)
    sensor_fluxes = ArchimedLight.compute_directional_fluxes(sky, sensor_turtle, sensor_options)
    sensor_first = ArchimedLight.compute_first_order(sensor_scene, sensor_models, sensor_turtle, sensor_fluxes, sensor_options; backend=:raycore_cpu)
    @test isapprox(get(sensor_first.incident_power.par, 1, 0.0), 100.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(sensor_first.incident_power.par, 2, 0.0), 100.0; atol=1e-8, rtol=1e-8)

    leaf_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="leaf", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="leaf", type="plate", object_id=2),
    ])
    leaf_models = ArchimedLight.prepare_models([
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
    upper_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    upper_turtle = ArchimedLight.build_turtle(upper_options, sky)
    upper_fluxes = ArchimedLight.compute_directional_fluxes(sky, upper_turtle, upper_options)
    upper_first = ArchimedLight.compute_first_order(leaf_scene, leaf_models, upper_turtle, upper_fluxes, upper_options; backend=:raycore_cpu)
    @test isapprox(get(upper_first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(upper_first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)

    stack_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    stack_turtle = ArchimedLight.build_turtle(stack_options, sky)
    stack_fluxes = ArchimedLight.compute_directional_fluxes(sky, stack_turtle, stack_options)
    stack_data = ArchimedLight._prepare_raycore_interception_data(
        leaf_scene,
        leaf_models,
        stack_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    traced_stack = ArchimedLight._raycore_trace_stacks(
        stack_data,
        [(0.5, 0.5, 2.0)],
        [(0.0, 0.0, -1.0)];
        t_maxs=[3.0f0],
    )
    @test traced_stack.counts == Int32[2]
    @test traced_stack.nodes[1:2] == UInt32[1, 2]
    @test traced_stack.heights[1] ≈ 1.0f0
    @test traced_stack.heights[2] ≈ 0.1f0

    stack_first = ArchimedLight.compute_first_order(leaf_scene, leaf_models, stack_turtle, stack_fluxes, stack_options; backend=:raycore_cpu)
    @test isapprox(get(stack_first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_first.incident_power.par, 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(stack_first.incident_power.par, 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)

    err = @test_throws ErrorException ArchimedLight.compute_first_order(
        leaf_scene,
        leaf_models,
        stack_turtle,
        stack_fluxes,
        stack_options;
        backend=ArchimedLight.RaycoreInterceptionBackend(max_hits_per_pixel=1),
    )
    @test occursin("max_hits_per_pixel=1 exceeded", sprint(showerror, err.value))
end

@testitem "Synthetic case Raycore emitter transfer" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="bulb", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="target", type="plate", object_id=2),
    ])
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "bulb" => ArchimedLight.TypeModel(
                    light_emitter=ArchimedLight.EmitterModel(radiance=100.0, gamma=ArchimedLight.OpticalProperties(0.6, 0.4)),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "target";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=false, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 0.0, 1.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    cpu_first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    raycore_first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=:raycore_cpu)

    @test get(raycore_first.incident_power.par, 2, 0.0) > 0.0
    @test isapprox(get(raycore_first.incident_power.par, 2, 0.0), get(cpu_first.incident_power.par, 2, 0.0); atol=1e-8, rtol=1e-8)
    @test isapprox(get(raycore_first.incident_power.nir, 2, 0.0), get(cpu_first.incident_power.nir, 2, 0.0); atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case Raycore scattering topology" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    ib = ArchimedLight.RaycoreInterceptionBackend()
    sb = ArchimedLight.RaycoreScatteringBackend(ib)
    dense_sb = ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=:dense_atomic)
    sparse_sb = ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=:sparse_host_reduce)
    float32_sb = ArchimedLight.RaycoreScatteringBackend(ArchimedLight.RaycoreInterceptionBackend(scattering_eltype=Float32))

    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
    cpu_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=ArchimedLight.RaycastScatteringBackend())
    graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sb)
    dense_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=dense_sb)
    sparse_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sparse_sb)
    sparse_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:sparse_host_reduce),
    )
    sparse_sector_idx = findfirst(
        sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0,
        turtle.sectors,
    )
    @test sparse_sector_idx !== nothing
    sparse_traced = ArchimedLight._raycore_trace_direction_stack_nodes_device(
        sparse_data,
        turtle.sectors[sparse_sector_idx].direction,
        options,
    )
    @test !any(sparse_traced.overflow)
    sparse_edge_keys1 = ArchimedLight._raycore_scattering_edge_keys_from_device_traced_stacks(sparse_data, sparse_traced)
    sparse_edge_keys2 = ArchimedLight._raycore_scattering_edge_keys_from_device_traced_stacks(sparse_data, sparse_traced)
    @test sparse_edge_keys1.keys === sparse_data.edge_keys_host
    @test sparse_edge_keys1.counts === sparse_data.edge_key_counts_host
    @test sparse_edge_keys2.keys === sparse_edge_keys1.keys
    @test sparse_edge_keys2.counts === sparse_edge_keys1.counts
    @test sparse_edge_keys1.max_edges == 2 * (sparse_data.max_hits_per_pixel - 1)
    sparse_edge_counts = Dict{UInt64,Int}()
    compact_scratch = sparse_data.edge_compact_host
    ArchimedLight._merge_counted_packed_edge_keys!(
        sparse_edge_counts,
        sparse_edge_keys1.keys,
        sparse_edge_keys1.counts,
        sparse_edge_keys1.max_edges,
        sparse_data.edge_compact_host,
    )
    @test !isempty(sparse_edge_counts)
    @test sparse_data.edge_compact_host === compact_scratch
    dense_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:dense_atomic),
    )
    dense_traced = ArchimedLight._raycore_trace_direction_stack_nodes_device(
        dense_data,
        turtle.sectors[sparse_sector_idx].direction,
        options,
    )
    @test !any(dense_traced.overflow)
    dense_counts1 = ArchimedLight._raycore_scattering_dense_counts_from_device_traced_stacks(dense_data, dense_traced)
    dense_counts2 = ArchimedLight._raycore_scattering_dense_counts_from_device_traced_stacks(dense_data, dense_traced)
    @test dense_data.dense_edge_counts_dev !== nothing
    @test dense_counts1 isa Vector{Int32}
    @test dense_counts2 == dense_counts1
    @test any(!iszero, dense_counts1)
    @test graph.dense[] === nothing
    dense_initial_par = ArchimedLight._dense_initial_scattering_power(graph, first, nothing, "PAR")
    @test dense_initial_par === first.dense.incident_power.par
    dense_initial_nir32 = ArchimedLight._dense_initial_scattering_power(graph, first, nothing, "NIR", Float32)
    @test dense_initial_nir32 isa Vector{Float32}
    @test dense_initial_nir32 ≈ Float32.(first.dense.incident_power.nir)
    dense_initial_copy_check, _ = ArchimedLight._dense_scattering_band_arrays(
        graph,
        dense_initial_nir32,
        graph.coeff_nir_by_node,
        graph.default_coeff_nir,
        Float32,
    )
    @test dense_initial_copy_check === dense_initial_nir32
    override_initial = Dict(graph.node_ids[1] => 1.25)
    dense_override = ArchimedLight._dense_initial_scattering_power(graph, first, override_initial, "PAR", Float32)
    @test dense_override[1] == Float32(1.25)
    @test all(iszero, dense_override[2:end])
    cpu_scat = ArchimedLight.compute_scattering(cpu_graph, first, options; backend=ArchimedLight.RaycastScatteringBackend())
    scat = ArchimedLight.compute_scattering(graph, first, options; backend=sb)
    dense_cache = graph.dense[]
    @test dense_cache !== nothing
    @test cpu_scat.dense !== nothing
    @test cpu_scat.dense.node_ids == cpu_graph.node_ids
    @test cpu_scat.dense.added_power.par == [get(cpu_scat.added_power.par, nid, 0.0) for nid in cpu_scat.dense.node_ids]
    scat32 = ArchimedLight.compute_scattering(graph, first, options; backend=float32_sb)
    @test graph.dense[] === dense_cache

    @test length(graph.pair_counts) > 0
    @test graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test dense_graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test dense_graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test dense_graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test sparse_graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test sparse_graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test sparse_graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test graph.all_hits == cpu_graph.all_hits
    @test get(scat.added_power.par, 2, 0.0) > 0.0
    @test HelperModule._dicts_close(scat.added_power.par, cpu_scat.added_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat.added_power.nir, cpu_scat.added_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat32.added_power.par, cpu_scat.added_power.par; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(scat32.added_power.nir, cpu_scat.added_power.nir; atol=1e-5, rtol=1e-5)
    @test scat.iterations == cpu_scat.iterations
    @test scat.converged == cpu_scat.converged
    @test scat32.iterations == cpu_scat.iterations
    @test scat32.converged == cpu_scat.converged

    too_small_dense_sb = ArchimedLight.RaycoreScatteringBackend(
        ib;
        edge_accumulation=:dense_atomic,
        dense_edge_limit_bytes=1,
    )
    err = @test_throws ErrorException ArchimedLight.build_scattering_transfer_graph(
        scene,
        models,
        turtle,
        first,
        options;
        backend=too_small_dense_sb,
    )
    @test occursin("dense_edge_limit_bytes=1", sprint(showerror, err.value))
end

@testitem "Synthetic case Raycore toricity wraparound" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1)])
    scene.scene_xy_bounds = (0.0, 0.0, 1.0, 1.0)
    sky = ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0)
    models = HelperModule._default_synthetic_models()

    non_toric_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    non_toric_turtle = ArchimedLight.build_turtle(non_toric_options, sky)
    non_toric_fluxes = ArchimedLight.compute_directional_fluxes(sky, non_toric_turtle, non_toric_options)
    non_toric_first = ArchimedLight.compute_first_order(scene, models, non_toric_turtle, non_toric_fluxes, non_toric_options; backend=:raycore_cpu)

    toric_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true)
    toric_turtle = ArchimedLight.build_turtle(toric_options, sky)
    toric_fluxes = ArchimedLight.compute_directional_fluxes(sky, toric_turtle, toric_options)
    toric_cpu = ArchimedLight.compute_first_order(scene, models, toric_turtle, toric_fluxes, toric_options)
    toric_raycore = ArchimedLight.compute_first_order(scene, models, toric_turtle, toric_fluxes, toric_options; backend=:raycore_cpu)

    @test get(toric_raycore.hits_per_node, 1, 0) > get(non_toric_first.hits_per_node, 1, 0)
    @test get(toric_raycore.incident_power.par, 1, 0.0) > get(non_toric_first.incident_power.par, 1, 0.0)
    @test isapprox(get(toric_raycore.projected_area_per_node, 1, 0.0), get(toric_cpu.projected_area_per_node, 1, 0.0); atol=1e-8, rtol=1e-8)
    @test isapprox(get(toric_raycore.incident_power.par, 1, 0.0), get(toric_cpu.incident_power.par, 1, 0.0); atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case Raycore pipeline entry points" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    row = HelperModule._synthetic_meteo_row(; ri_par_f=100.0, ri_nir_f=30.0, direct_fraction=1.0)
    ib = ArchimedLight.RaycoreInterceptionBackend()
    sb = ArchimedLight.RaycoreScatteringBackend(ib)

    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
    graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sb)
    scat = ArchimedLight.compute_scattering(graph, first, options; backend=sb)
    budget = ArchimedLight.integrate_light(scene, models, first, scat, options; meteo_row=row)
    first_public_only = ArchimedLight.FirstOrderResult(first.projected_area_per_node, first.incident_power, first.hits_per_node)
    scat_public_only = ArchimedLight.ScatteringResult(scat.added_power, scat.iterations, scat.converged)
    budget_public_only = ArchimedLight.integrate_light(scene, models, first_public_only, scat_public_only, options; meteo_row=row)

    step = ArchimedLight.run_light_step(
        scene,
        models,
        row,
        options;
        interception_backend=ib,
        scattering_backend=sb,
    )
    par_only_step = ArchimedLight.run_light_step(
        scene,
        models,
        row,
        ArchimedLight.LightOptions(options; nir_scattering=false);
        interception_backend=ib,
        scattering_backend=sb,
    )

    @test step.first_order.projected_area_per_node == first.projected_area_per_node
    @test step.scattering.added_power.par == scat.added_power.par
    @test step.budget.incident_energy.total.par == budget.incident_energy.total.par
    @test first.dense !== nothing
    @test scat.dense !== nothing
    @test first.dense.node_ids == scat.dense.node_ids
    @test first.dense.incident_power.par == [get(first.incident_power.par, nid, 0.0) for nid in first.dense.node_ids]
    @test scat.dense.added_power.par == [get(scat.added_power.par, nid, 0.0) for nid in scat.dense.node_ids]
    @test HelperModule._budgets_close(budget, budget_public_only; atol=1e-10, rtol=1e-10)
    @test par_only_step.scattering.added_power.nir == Dict{Int,Float64}()
    @test par_only_step.scattering.dense !== nothing
    @test par_only_step.scattering.dense.node_ids == par_only_step.first_order.dense.node_ids
    @test par_only_step.scattering.dense.added_power.par ==
          [get(par_only_step.scattering.added_power.par, nid, 0.0) for nid in par_only_step.scattering.dense.node_ids]
    @test all(iszero, par_only_step.scattering.dense.added_power.nir)
end

@testitem "Synthetic case Raycore series response cache" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=90.0, ri_nir_f=70.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_raycore_series"))
    ib = ArchimedLight.RaycoreInterceptionBackend()

    uncached_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=true, scattering=false, pixel_size=0.01, cache_radiation=false, toricity=false)
    cached_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=true, scattering=false, pixel_size=0.01, cache_radiation=true, toricity=false)
    cached_cache = ArchimedLight.prepare_light_cache(scene, models, cached_options; interception_backend=ib)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options; interception_backend=ib)
    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options; interception_backend=ib)

    @test ArchimedLight.cache_summary(cached_cache).mode == :full
    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test cached[i].first_order.projected_area_per_node == uncached[i].first_order.projected_area_per_node
        @test cached[i].budget.incident_energy.total.par == uncached[i].budget.incident_energy.total.par
    end

    stacked = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    scatter_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true, toricity=false)
    scatter_uncached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false, toricity=false)
    sb = ArchimedLight.RaycoreScatteringBackend(ib)
    scatter_cache = ArchimedLight.prepare_light_cache(stacked, models, scatter_options; interception_backend=ib, scattering_backend=sb)
    scatter_sky = ArchimedLight.compute_sky(rows[1], scatter_options)
    scatter_turtle = ArchimedLight.build_turtle(scatter_options, scatter_sky)
    mesh_area_cache = scatter_cache.raycore_data.projected_mesh_area_cache
    area_ratio_cache = scatter_cache.raycore_data.area_ratio_cache
    mesh_area_keys = Set(ArchimedLight._raycore_direction_key(sector.direction) for sector in scatter_turtle.sectors)
    @test isempty(mesh_area_cache)
    @test isempty(area_ratio_cache)
    up_sector_idx = findfirst(sector -> Float32(sector.direction[3]) > 0.0f0, scatter_turtle.sectors)
    @test up_sector_idx !== nothing
    up_direction = scatter_turtle.sectors[up_sector_idx].direction
    top_projection = ArchimedLight._raycore_direction_projection_top_hit(scatter_cache.raycore_data, up_direction, scatter_options)
    top_ratio1 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :top_hit, up_direction, top_projection.projected_pixels_area)
    top_ratio2 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :top_hit, up_direction, top_projection.projected_pixels_area)
    full_projection = ArchimedLight._raycore_direction_projection_full_stack(scatter_cache.raycore_data, up_direction, scatter_options)
    full_ratio1 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :full_stack, up_direction, full_projection.projected_pixels_area)
    full_ratio2 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :full_stack, up_direction, full_projection.projected_pixels_area)
    @test top_ratio1 === top_ratio2
    @test full_ratio1 === full_ratio2
    @test top_ratio1 !== full_ratio1
    @test length(area_ratio_cache) == 2
    scatter_entry = ArchimedLight._get_turtle_cache_entry!(scatter_cache, scatter_turtle)
    @test Set(keys(mesh_area_cache)) == mesh_area_keys
    cached_mesh_area_count = length(mesh_area_cache)
    @test ArchimedLight._get_turtle_cache_entry!(scatter_cache, scatter_turtle) === scatter_entry
    @test length(mesh_area_cache) == cached_mesh_area_count
    topology = scatter_entry.responses_cache.scattering_topology
    @test topology !== nothing
    @test topology.dense_static[] === nothing
    unit_first = ArchimedLight._combine_sector_responses(
        scatter_entry.responses_cache,
        ArchimedLight._unit_directional_fluxes(scatter_entry.turtle, 1; par=1.0),
    )
    graph1 = ArchimedLight.build_scattering_transfer_graph(topology, unit_first, scatter_options, sb)
    @test topology.dense_static[] !== nothing
    graph2 = ArchimedLight.build_scattering_transfer_graph(topology, unit_first, scatter_options, sb)
    @test graph1.dense_static === topology.dense_static[]
    @test graph2.dense_static === topology.dense_static[]
    @test isempty(graph1.dense_static.device_cache)
    dev_static1 = ArchimedLight._scattering_static_edge_device_arrays(graph1.dense_static, sb.config.backend)
    dev_static2 = ArchimedLight._scattering_static_edge_device_arrays(graph2.dense_static, sb.config.backend)
    @test dev_static1 === dev_static2
    scatter_uncached = ArchimedLight.run_light_series(stacked, models, meteo, scatter_uncached_options; interception_backend=ib, scattering_backend=sb)
    scatter_cached = ArchimedLight.run_light_series(stacked, models, meteo, scatter_options; interception_backend=ib, scattering_backend=sb)

    @test ArchimedLight.cache_summary(scatter_cache).mode == :full
    @test length(scatter_cached) == length(scatter_uncached)
    for i in eachindex(scatter_cached)
        @test scatter_cached[i].first_order.projected_area_per_node == scatter_uncached[i].first_order.projected_area_per_node
        @test HelperModule._dicts_close(scatter_cached[i].scattering.added_power.par, scatter_uncached[i].scattering.added_power.par; atol=1e-10, rtol=1e-10)
        @test HelperModule._dicts_close(scatter_cached[i].budget.incident_energy.total.par, scatter_uncached[i].budget.incident_energy.total.par; atol=1e-8, rtol=1e-10)
    end

    extra_rows = [
        merge(rows[1], (RI_UV_F=25.0,)),
        merge(rows[2], (RI_UV_F=10.0,)),
    ]
    extra_meteo = ArchimedLight.MeteoTable(extra_rows, (; source="synthetic_raycore_extra_band_cache"))
    extra_uncached = ArchimedLight.run_light_series(stacked, models, extra_meteo, scatter_uncached_options; interception_backend=ib, scattering_backend=sb)
    extra_cached = ArchimedLight.run_light_series(stacked, models, extra_meteo, scatter_options; interception_backend=ib, scattering_backend=sb)

    @test length(extra_cached) == length(extra_uncached)
    for i in eachindex(extra_cached)
        @test HelperModule._budgets_close(extra_cached[i].budget, extra_uncached[i].budget; atol=1e-8, rtol=1e-10)
    end
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
    run = HelperModule._run_direct(scene, sky, options; models=models)

    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(run.first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 100.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 0.0; atol=1e-8, rtol=1e-8)

    stack_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01)
    stack_run = HelperModule._run_direct(scene, sky, stack_options; models=models)

    @test isapprox(get(stack_run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_run.first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)
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
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_cached_series"))
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
    @test meteo_read isa ArchimedLight.MeteoTable
    @test first(selected).Ri_SW_f == 200.0
    @test isapprox(sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case overlapping meteo steps option" tags = [:synthetic, :fast, :overlapping_meteo_steps] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()

    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 0), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 15), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_overlap"))

    strict_options = HelperModule._synthetic_options(cache_radiation=false)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.prepare_meteo(meteo, strict_options)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.run_light_series(scene, models, meteo, strict_options)

    permissive_options = ArchimedLight.LightOptions(strict_options; allow_overlapping_meteo_steps=true)
    selected = ArchimedLight.prepare_meteo(meteo, permissive_options)
    series = ArchimedLight.run_light_series(scene, models, meteo, permissive_options)

    @test length(selected.rows) == 2
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
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_cached_scattering_series"))
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
    meteo = ArchimedLight.MeteoTable(rows, (; source="synthetic_manual_cache"))

    cache = ArchimedLight.prepare_light_cache(scene, models, options)
    summary0 = ArchimedLight.cache_summary(cache)
    @test summary0.mode == :full
    @test summary0.cached_turtle_count == 0

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

    scene2 = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.5, group="middle", type="plate", object_id=3),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    cache2 = ArchimedLight.prepare_light_cache(scene2, models, options)
    step_rebuilt = ArchimedLight.run_light_step(cache2, rows[1])
    @test !HelperModule._budgets_close(step_cached.budget, step_rebuilt.budget; atol=1e-9, rtol=1e-9)
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
    meteo = ArchimedLight.MeteoTable([row0, row1], (; source="synthetic_extra_band_cache"))

    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options)

    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test HelperModule._budgets_close(cached[i].budget, uncached[i].budget; atol=1e-9, rtol=1e-9)
    end
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
