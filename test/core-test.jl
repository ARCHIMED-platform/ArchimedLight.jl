@testitem "Core smoke" tags=[:core, :fast] begin
    import PlantGeom

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    selected = ArchimedLight.prepare_meteo(fixture.meteo, fixture.options)
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)

    @test !isempty(selected.rows)
    @test !isempty(series)
    @test length(series) == length(selected.rows)
    @test !isempty(series[1].budget.incident_energy.total.par)
    @test length(series[1].turtle.sectors) == 17
    step_summary = sprint(show, series[1])
    step_pretty = sprint(show, MIME"text/plain"(), series[1])
    @test occursin("LightStepResult(", step_summary)
    @test occursin("sectors=17", step_summary)
    @test occursin("LightStepResult", step_pretty)
    @test occursin("sky        PAR", step_pretty)
    @test occursin("incident   PAR", step_pretty)
    @test occursin("absorbed   PAR", step_pretty)
    @test occursin("turtle     17 sectors", step_pretty)
    @test occursin("scattering off", step_pretty)

    options2, scene2, meteo2, models2 = ArchimedLight.read_config(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml"))
    @test options2.turtle_sectors == fixture.options.turtle_sectors
    @test length(meteo2.rows) == length(fixture.meteo.rows)
    @test collect(keys(models2.groups)) == collect(keys(fixture.models.groups))
    @test any(item -> item.group == "pavement", ArchimedLight.summarize_scene(scene2).group_types)
    @test haskey(scene2.mtg, :geometry)
    @test scene2.mtg[:geometry] === nothing

    _, scene3, _, _ = ArchimedLight.read_config(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml");
        plot_paving_override=25,
    )
    pavement = only(item for item in ArchimedLight.summarize_scene(scene3).group_types if item.group == "pavement")
    @test pavement.nodes == 25
    @test scene3.mtg[:geometry] === nothing

    raw_scene = ArchimedLight.read_scene(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "scene", "simple.ops"),
    )
    @test haskey(raw_scene.mtg, :geometry)
    @test raw_scene.mtg[:geometry] !== nothing
    PlantGeom.add_ground!(raw_scene; nx=2, ny=2)
    @test raw_scene.mtg[:geometry] === nothing
end

@testitem "RasterGPU backend matches CPU raster fixtures" tags=[:core, :fast, :raster_gpu] begin
    import KernelAbstractions

    include(joinpath(@__DIR__, "support.jl"))

    function first_order_pair(fixture_name; toricity)
        fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", fixture_name, "input"))
        options = ArchimedLight.LightOptions(
            fixture.options;
            scattering=false,
            cache_pixel_table=false,
            toricity=toricity,
        )
        row = first(fixture.meteo.rows)
        sky = ArchimedLight.compute_sky(row, options)
        turtle = ArchimedLight.build_turtle(options, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
        cpu = ArchimedLight.compute_first_order(
            fixture.scene,
            fixture.models,
            turtle,
            fluxes,
            options;
            backend=:raster_cpu,
        )
        gpu = ArchimedLight.compute_first_order(
            fixture.scene,
            fixture.models,
            turtle,
            fluxes,
            options;
            backend=ArchimedLight.RasterGPUBackend(
                backend=KernelAbstractions.CPU(),
                max_hits_per_pixel=512,
                tile_size=4,
                tile_face_capacity=4096,
            ),
        )
        return cpu.dense, gpu.dense
    end

    for (fixture_name, toricity) in (
        ("simpleplant_16_notoric", false),
        ("simpleplant_16_toric", true),
    )
        cpu, gpu = first_order_pair(fixture_name; toricity=toricity)
        @test cpu.node_ids == gpu.node_ids
        @test cpu.hits_per_node == gpu.hits_per_node
        @test cpu.projected_area_per_node ≈ gpu.projected_area_per_node atol = 1e-7 rtol = 1e-6
        @test cpu.incident_power.par ≈ gpu.incident_power.par atol = 1e-6 rtol = 1e-6
        @test cpu.incident_power.nir ≈ gpu.incident_power.nir atol = 1e-6 rtol = 1e-6
    end
end

@testitem "RasterGPU allocation preflight reports actionable diagnostics" tags=[:core, :fast, :raster_gpu] begin
    config = ArchimedLight.RasterGPUBackendConfig(
        max_hits_per_pixel=64,
        tile_size=1,
        tile_face_capacity=1024,
    )
    context = ArchimedLight._rastergpu_config_context(
        config;
        n_pixels=100,
        n_tiles=100,
        tile_size=1,
        tile_face_capacity=1024,
        max_hits=64,
        edge_key_capacity=12_600,
        stackless_top_hit=false,
    )

    err = try
        ArchimedLight._rastergpu_checked_buffer_count(
            "RasterGPU edge_keys_dev",
            UInt64,
            12_600,
            config;
            backend_limit=128,
            context=context,
        )
        nothing
    catch caught
        caught
    end
    @test err isa ErrorException
    msg = sprint(showerror, err)
    @test occursin("RasterGPU edge_keys_dev", msg)
    @test occursin("requested_count=12600", msg)
    @test occursin("eltype=UInt64", msg)
    @test occursin("requested_bytes=100800", msg)
    @test occursin("backend_buffer_limit_bytes=128", msg)
    @test occursin("max_hits_per_pixel=64", msg)
    @test occursin("tile_face_capacity=1024", msg)

    tile_err = try
        ArchimedLight._rastergpu_checked_buffer_count(
            "RasterGPU tile_faces_dev/tile_unwrapped_i_dev/tile_unwrapped_j_dev",
            Int32,
            1_000,
            config;
            backend_limit=512,
            context=context,
        )
        nothing
    catch caught
        caught
    end
    @test tile_err isa ErrorException
    tile_msg = sprint(showerror, tile_err)
    @test occursin("tile_faces_dev", tile_msg)
    @test occursin("tile_unwrapped_i_dev", tile_msg)
    @test occursin("tile_unwrapped_j_dev", tile_msg)

    overflow_err = try
        ArchimedLight._rastergpu_checked_buffer_count(
            "RasterGPU tile_faces_dev",
            Int32,
            Int128(typemax(Int)) + 1,
            config;
            backend_limit=nothing,
            context=context,
        )
        nothing
    catch caught
        caught
    end
    @test overflow_err isa ErrorException
    @test occursin("element count exceeds typemax(Int)", sprint(showerror, overflow_err))

    aggregate_err = try
        ArchimedLight._rastergpu_check_total_allocation_bytes(
            [
                "tile_faces_dev" => Int128(600),
                "tile_unwrapped_i_dev" => Int128(600),
                "nodes_dev" => Int128(400),
            ],
            config;
            memory_info=(working_set=Int128(2000), allocated=Int128(200), free=Int128(1800)),
            context=context,
        )
        nothing
    catch caught
        caught
    end
    @test aggregate_err isa ErrorException
    aggregate_msg = sprint(showerror, aggregate_err)
    @test occursin("aggregate device allocation preflight failed", aggregate_msg)
    @test occursin("requested_device_bytes=1600", aggregate_msg)
    @test occursin("metal_working_set_bytes=2000", aggregate_msg)
    @test occursin("tile_faces_dev=600", aggregate_msg)
end

@testitem "RasterGPU auto edge accumulation can exceed dense cap when it fits device budget" tags=[:core, :fast, :raster_gpu] begin
    config = ArchimedLight.RasterGPUBackendConfig(
        backend=:fake_gpu_backend,
        max_hits_per_pixel=128,
        tile_size=1,
        tile_face_capacity=256,
        edge_accumulation=:auto,
    )

    dense_pairs = Int128(22_373) * Int128(22_373)
    dense_bytes = dense_pairs * Int128(sizeof(Int32))
    @test dense_bytes > config.dense_edge_limit_bytes

    sparse_total_bytes = Int128(24_585_595_926)
    dense_total_bytes = Int128(14_400_000_000)
    memory_info = (
        working_set=Int128(30_150_672_384),
        allocated=Int128(1_785_856),
        free=Int128(30_148_886_528),
    )

    @test ArchimedLight._rastergpu_choose_dense_edge_accumulation(
        config;
        stackless_top_hit=false,
        dense_pairs=dense_pairs,
        dense_bytes=dense_bytes,
        sparse_total_bytes=sparse_total_bytes,
        dense_total_bytes=dense_total_bytes,
        memory_info=memory_info,
    )
    @test !ArchimedLight._rastergpu_choose_dense_edge_accumulation(
        config;
        stackless_top_hit=false,
        dense_pairs=dense_pairs,
        dense_bytes=dense_bytes,
        sparse_total_bytes=sparse_total_bytes,
        dense_total_bytes=dense_total_bytes,
        memory_info=nothing,
    )
    @test !ArchimedLight._rastergpu_choose_dense_edge_accumulation(
        config;
        stackless_top_hit=false,
        dense_pairs=dense_pairs,
        dense_bytes=dense_bytes,
        sparse_total_bytes=sparse_total_bytes,
        dense_total_bytes=Int128(25_000_000_000),
        memory_info=memory_info,
    )

    explicit_dense = ArchimedLight._rastergpu_config_with(config; edge_accumulation=:dense_atomic)
    @test !ArchimedLight._rastergpu_choose_dense_edge_accumulation(
        explicit_dense;
        stackless_top_hit=false,
        dense_pairs=dense_pairs,
        dense_bytes=dense_bytes,
        sparse_total_bytes=sparse_total_bytes,
        dense_total_bytes=dense_total_bytes,
        memory_info=memory_info,
    )
end

@testitem "RasterGPU top-hit scene data avoids full-stack scratch" tags=[:core, :fast, :raster_gpu] begin
    import KernelAbstractions

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    options = ArchimedLight.LightOptions(
        fixture.options;
        scattering=false,
        cache_pixel_table=false,
    )
    prepared = ArchimedLight._prepare_interception_data(
        fixture.scene,
        fixture.models,
        options;
        include_budget_maps=true,
    )
    @test prepared.upper_hit

    config = ArchimedLight.RasterGPUBackendConfig(
        backend=KernelAbstractions.CPU(),
        tile_face_capacity=64,
    )
    @test config.top_hit_tile_face_capacity == 64
    explicit_top_hit_config = ArchimedLight.RasterGPUBackendConfig(
        backend=KernelAbstractions.CPU(),
        tile_face_capacity=64,
        top_hit_tile_face_capacity=17,
    )
    @test explicit_top_hit_config.top_hit_tile_face_capacity == 17

    data = ArchimedLight._rastergpu_scene_data(prepared, config)
    @test data.tile_face_capacity == 64
    @test length(data.nodes_dev) == 1
    @test length(data.heights_dev) == 1
    @test isempty(data.nodes_host)
    @test isempty(data.heights_host)
    @test data.edge_keys_dev === nothing
    @test data.edge_key_counts_dev === nothing
    @test data.dense_edge_counts_dev === nothing
end

@testitem "RasterGPU auto edge accumulation avoids sparse scratch when dense is smaller" tags=[:core, :fast, :raster_gpu] begin
    import KernelAbstractions

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    options = ArchimedLight.LightOptions(
        fixture.options;
        scattering=true,
        cache_pixel_table=false,
        toricity=false,
    )
    prepared = ArchimedLight._prepare_interception_data(
        fixture.scene,
        fixture.models,
        options;
        include_budget_maps=true,
    )
    @test !prepared.upper_hit

    auto_config = ArchimedLight.RasterGPUBackendConfig(
        backend=KernelAbstractions.CPU(),
        max_hits_per_pixel=64,
        tile_size=1,
        tile_face_capacity=4096,
        edge_accumulation=:auto,
    )
    auto_data = ArchimedLight._rastergpu_scene_data(prepared, auto_config)
    @test !ArchimedLight._rastergpu_use_fused_dense_edges(auto_data)
    @test length(auto_data.nodes_dev) > 1
    @test length(auto_data.heights_dev) > 1
    @test auto_data.dense_edge_counts_dev !== nothing
    @test auto_data.dense_edge_counts_host !== nothing
    @test auto_data.edge_keys_dev === nothing
    @test auto_data.edge_key_counts_dev === nothing
    @test !ArchimedLight._rastergpu_fused_dense_supported(auto_config)
    @test ArchimedLight._rastergpu_fused_dense_supported(
        ArchimedLight.RasterGPUBackendConfig(
            backend=:fake_gpu_backend,
            max_hits_per_pixel=128,
            edge_accumulation=:auto,
        ),
    )
    @test !ArchimedLight._rastergpu_fused_dense_supported(
        ArchimedLight.RasterGPUBackendConfig(
            backend=:fake_gpu_backend,
            max_hits_per_pixel=129,
            edge_accumulation=:auto,
        ),
    )

    sparse_config = ArchimedLight.RasterGPUBackendConfig(
        backend=KernelAbstractions.CPU(),
        max_hits_per_pixel=64,
        tile_size=4,
        tile_face_capacity=4096,
        edge_accumulation=:sparse_host_reduce,
    )
    sparse_data = ArchimedLight._rastergpu_scene_data(prepared, sparse_config)
    @test sparse_data.dense_edge_counts_dev === nothing
    @test sparse_data.edge_keys_dev !== nothing
    @test sparse_data.edge_key_counts_dev !== nothing
end

@testitem "RasterGPU scattering backend matches CPU raycast fixture" tags=[:core, :fast, :raster_gpu] begin
    import KernelAbstractions

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    options = ArchimedLight.LightOptions(
        fixture.options;
        scattering=true,
        cache_pixel_table=false,
        toricity=false,
        scattering_max_iter=3,
    )
    row = first(fixture.meteo.rows)
    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)

    cpu_first = ArchimedLight.compute_first_order(
        fixture.scene,
        fixture.models,
        turtle,
        fluxes,
        options;
        backend=:raster_cpu,
    )
    cpu_scattering = ArchimedLight.compute_scattering(
        fixture.scene,
        fixture.models,
        turtle,
        cpu_first,
        options;
        backend=ArchimedLight.RaycastScatteringBackend(),
    )

    interception_backend = ArchimedLight.RasterGPUBackend(
        backend=KernelAbstractions.CPU(),
        max_hits_per_pixel=512,
        tile_size=4,
        tile_face_capacity=4096,
        edge_accumulation=:sparse_host_reduce,
    )
    scattering_backend = ArchimedLight.RasterGPUScatteringBackend(interception_backend)
    gpu_first = ArchimedLight.compute_first_order(
        fixture.scene,
        fixture.models,
        turtle,
        fluxes,
        options;
        backend=interception_backend,
    )
    gpu_scattering = ArchimedLight.compute_scattering(
        fixture.scene,
        fixture.models,
        turtle,
        gpu_first,
        options;
        backend=scattering_backend,
    )

    @test gpu_scattering.dense !== nothing
    @test gpu_scattering.iterations == cpu_scattering.iterations
    @test gpu_scattering.converged == cpu_scattering.converged
    for nid in cpu_first.dense.node_ids
        @test get(gpu_scattering.added_power.par, nid, 0.0) ≈
              get(cpu_scattering.added_power.par, nid, 0.0) atol = 1e-8 rtol = 1e-6
        @test get(gpu_scattering.added_power.nir, nid, 0.0) ≈
              get(cpu_scattering.added_power.nir, nid, 0.0) atol = 1e-8 rtol = 1e-6
    end
end

@testitem "Sky fraction can be stored and attached" tags=[:core, :fast] begin
    import MultiScaleTreeGraph

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    options = ArchimedLight.LightOptions(fixture.options; include_sky_fraction=true)
    @test options.include_sky_fraction
    meteo = ArchimedLight.prepare_meteo(fixture.meteo, options)
    series = ArchimedLight.run_light_series(
        fixture.scene,
        fixture.models,
        meteo,
        options,
    )

    @test length(series) == length(meteo.rows)
    @test series[1].sky_fraction !== nothing
    @test !isempty(series[1].sky_fraction)

    ArchimedLight.attach_light_step!(fixture.scene, series[1]; fields=[:area])

    scalar_area_found = Ref(false)
    MultiScaleTreeGraph.traverse!(fixture.scene.mtg) do node
        nid = Int(MultiScaleTreeGraph.node_id(node))
        if haskey(fixture.scene.nodes, nid) && haskey(node, :area)
            scalar_area_found[] = true
            @test node[:area] ≈ fixture.scene.nodes[nid].area
            return false
        end
        return true
    end
    @test scalar_area_found[]

    ArchimedLight.attach_light_series!(
        fixture.scene,
        series;
        fields=[:area, :absorbed_par_flux, :absorbed_nir_flux, :sky_fraction],
        names=Dict(:absorbed_nir_flux => :Ra_SW_f),
    )

    found = Ref(false)
    MultiScaleTreeGraph.traverse!(fixture.scene.mtg) do node
        if haskey(node, :sky_fraction)
            found[] = true
            @test node[:sky_fraction] isa Vector{Float64}
            @test length(node[:sky_fraction]) == length(series)
            @test node[:area] isa Vector{Float64}
            @test length(node[:area]) == length(series)
            expected_area = fixture.scene.nodes[Int(MultiScaleTreeGraph.node_id(node))].area
            @test all(v -> isapprox(v, expected_area), node[:area])
            @test haskey(node, :Ra_PAR_f)
            @test haskey(node, :Ra_SW_f)
            return false
        end
        return true
    end
    @test found[]
end

@testitem "SmallHitStack handles dense pixels" tags=[:core, :fast] begin
    stack = ArchimedLight.SmallHitStack()
    for i in 1:300
        push!(stack, (Float64(i), i))
    end
    @test length(stack) == 300
    @test stack[1] == (1.0, 1)
    @test stack[end] == (300.0, 300)
end
