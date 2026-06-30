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

@testitem "Raycore backend plumbing" tags=[:core, :fast, :raycore_backend] begin
    ib = ArchimedLight.RaycoreInterceptionBackend()

    @test ib isa ArchimedLight.InterceptionBackend
    @test ib.config.workgroupsize == 256
    @test ArchimedLight._raycore_default_workgroupsize(ArchimedLight.KernelAbstractions.CPU()) == 256
    @test ib.config.max_hits_per_pixel == 32
    @test ib.config.hit_epsilon == Float32(1.0f-4)
    raw_stack_config = ArchimedLight.RaycoreBackendConfig(hit_epsilon=-1.0f0)
    @test raw_stack_config.hit_epsilon == -1.0f0
    @test ib.config.edge_accumulation == :auto
    @test ib.config.dense_edge_limit_bytes == 512 * 1024^2
    @test ib.config.propagation_backend == :auto
    @test ib.config.max_prechunk_instances == ArchimedLight._raycore_default_max_prechunk_instances()
    @test ib.config.scattering_eltype == Float64
    @test !ib.config.validate
    no_cap_config = ArchimedLight.RaycoreBackendConfig(max_prechunk_instances=0)
    @test no_cap_config.max_prechunk_instances == typemax(Int)
    validating_config = ArchimedLight.RaycoreBackendConfig(validate=true)
    @test validating_config.validate

    resolved = ArchimedLight._resolve_interception_backend(:raycore_cpu)
    @test resolved isa ArchimedLight.RaycoreInterceptionBackend

    err = try
        ArchimedLight._resolve_interception_backend(:not_a_backend)
        nothing
    catch caught
        caught
    end
    @test err isa ErrorException
    @test occursin(":raycore_cpu", sprint(showerror, err))

    sb = ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=:sparse_host_reduce)
    @test sb isa ArchimedLight.ScatteringBackend
    @test sb.config.backend == ib.config.backend
    @test sb.config.edge_accumulation == :sparse_host_reduce
    @test sb.config.propagation_backend == :auto
    @test sb.config.max_prechunk_instances == ib.config.max_prechunk_instances
    @test sb.config.scattering_eltype == Float64
    @test sb.config.validate == ib.config.validate
    @test ArchimedLight._resolve_scattering_backend(:raycast, sb) === sb

    forced_cpu_sb = ArchimedLight.RaycoreScatteringBackend(ib; propagation_backend=:cpu)
    @test forced_cpu_sb.config.propagation_backend == :cpu
    validating_sb = ArchimedLight.RaycoreScatteringBackend(ib; validate=true)
    @test validating_sb.config.validate
    capped_sb = ArchimedLight.RaycoreScatteringBackend(ib; max_prechunk_instances=17)
    @test capped_sb.config.max_prechunk_instances == 17
    forced_device_sb = ArchimedLight.RaycoreScatteringBackend(propagation_backend=:device)
    @test forced_device_sb.config.propagation_backend == :device
    strict_error = ArchimedLight.RaycoreValidationError(
        :raycore_trace_validation,
        :unit_test,
        "diagnostic message",
        (; ok=false),
    )
    @test occursin("RaycoreValidationError(raycore_trace_validation, unit_test)", sprint(showerror, strict_error))
    @test occursin("diagnostic message", sprint(showerror, strict_error))

    raw_trace(counts, metadata, instances, distances; max_hits=3) = (
        counts=Int32.(counts),
        metadata=UInt32.(metadata),
        instance_indices=UInt32.(instances),
        distances=Float32.(distances),
        overflow=falses(length(counts)),
        max_hits=max_hits,
    )
    missing_ref = raw_trace([2], [1, 2, 0], [1, 1, 0], [1.0, 2.0, 0.0])
    missing_got = raw_trace([1], [1, 0, 0], [1, 0, 0], [1.0, 0.0, 0.0])
    missing = ArchimedLight._raycore_first_raw_stack_mismatch(missing_ref, missing_got)
    @test missing.kind == :count
    @test missing.diagnosis == :candidate_missing_hits
    @test missing.unmatched_reference_count == 1
    @test missing.unmatched_candidate_count == 0
    @test missing.first_unmatched_reference.metadata == UInt32(2)

    extra = ArchimedLight._raycore_first_raw_stack_mismatch(missing_got, missing_ref)
    @test extra.kind == :count
    @test extra.diagnosis == :candidate_extra_hits
    @test extra.unmatched_reference_count == 0
    @test extra.unmatched_candidate_count == 1
    @test extra.first_unmatched_candidate.metadata == UInt32(2)

    duplicate_suppression_ref = raw_trace([2], [1, 1, 0], [1, 2, 0], [1.0, 1.0, 0.0])
    duplicate_suppression_got = raw_trace([1], [1, 0, 0], [1, 0, 0], [1.0, 0.0, 0.0])
    duplicate_suppression = ArchimedLight._raycore_first_raw_stack_mismatch(
        duplicate_suppression_ref,
        duplicate_suppression_got,
    )
    @test duplicate_suppression.kind == :count
    @test duplicate_suppression.diagnosis == :candidate_missing_hits
    @test duplicate_suppression.unmatched_reference_count == 1
    @test duplicate_suppression.unmatched_candidate_count == 0
    @test duplicate_suppression.first_unmatched_reference.metadata == UInt32(1)
    @test duplicate_suppression.first_unmatched_reference.instance_index == UInt32(2)

    order_ref = raw_trace([2], [1, 2, 0], [1, 1, 0], [1.0, 2.0, 0.0])
    order_got = raw_trace([2], [2, 1, 0], [1, 1, 0], [2.0, 1.0, 0.0])
    order = ArchimedLight._raycore_first_raw_stack_mismatch(order_ref, order_got)
    @test order.kind == :metadata
    @test order.diagnosis == :raw_order_difference
    @test order.unmatched_reference_count == 0
    @test order.unmatched_candidate_count == 0

    instance_ref = raw_trace([2], [1, 1, 0], [1, 2, 0], [1.0, 2.0, 0.0])
    instance_got = raw_trace([2], [1, 1, 0], [1, 3, 0], [1.0, 2.0, 0.0])
    instance = ArchimedLight._raycore_first_raw_stack_mismatch(instance_ref, instance_got)
    @test instance.kind == :instance_index
    @test instance.diagnosis == :instance_index_difference
    @test instance.first_unmatched_reference.instance_index == UInt32(2)
    @test instance.first_unmatched_candidate.instance_index == UInt32(3)

    distance_ref = raw_trace([1], [1, 0, 0], [1, 0, 0], [1.0, 0.0, 0.0])
    distance_got = raw_trace([1], [1, 0, 0], [1, 0, 0], [1.01, 0.0, 0.0])
    distance = ArchimedLight._raycore_first_raw_stack_mismatch(distance_ref, distance_got; distance_atol=1.0f-4)
    @test distance.kind == :distance
    @test distance.diagnosis == :distance_difference

    ib32 = ArchimedLight.RaycoreInterceptionBackend(scattering_eltype=Float32)
    sb32 = ArchimedLight.RaycoreScatteringBackend(ib32)
    @test ib32.config.scattering_eltype == Float32
    @test sb32.config.scattering_eltype == Float32

    @test_throws ErrorException ArchimedLight.RaycoreBackendConfig(workgroupsize=0)
    @test_throws ErrorException ArchimedLight.RaycoreBackendConfig(edge_accumulation=:unsupported)
    @test_throws ErrorException ArchimedLight.RaycoreBackendConfig(propagation_backend=:unsupported)

    dense_counts = Int32[1]
    dense_area = Float32[1.0]
    @test !ArchimedLight._raycore_flat_stack_host_copy_required(false, false, dense_counts, dense_area, true)
    @test !ArchimedLight._raycore_flat_stack_host_copy_required(true, true, dense_counts, dense_area, true)
    @test ArchimedLight._raycore_flat_stack_host_copy_required(true, false, dense_counts, dense_area, true)
    @test ArchimedLight._raycore_flat_stack_host_copy_required(false, false, nothing, dense_area, true)
    @test ArchimedLight._raycore_flat_stack_host_copy_required(false, false, dense_counts, nothing, true)
    @test !ArchimedLight._raycore_flat_stack_host_copy_required(false, false, dense_counts, nothing, false)
    empty_shape = ArchimedLight._raycore_scene_shape_summary(nothing)
    @test empty_shape.geometry_mode == :none
    @test empty_shape.tlas_instances == 0
    @test empty_shape.tlas_geometries == 0
    @test empty_shape.raycore_traversal_stack_capacity == ArchimedLight._raycore_software_traversal_stack_capacity()
    @test empty_shape.max_blas_depth == 0
    @test !empty_shape.max_blas_depth_over_stack_capacity
    @test empty_shape.blas_node_count == 0
    @test empty_shape.expanded_face_count == 0
    @test empty_shape.dense_edge_pairs == 0
    @test empty_shape.edge_key_capacity == 0

    projected_area_by_sector = zeros(Float64, 4, 2)
    projected_area_active_indices = Int[]
    projected_area_active_offsets = [1, 0, 0]
    ArchimedLight._store_device_sector_area_and_active_indices!(
        projected_area_by_sector,
        projected_area_active_indices,
        projected_area_active_offsets,
        1,
        Float32[0.0, 2.0, 3.0, 0.0],
        nothing,
    )
    ArchimedLight._store_device_sector_area_and_active_indices!(
        projected_area_by_sector,
        projected_area_active_indices,
        projected_area_active_offsets,
        2,
        Float32[1.0, 0.0, 4.0, 0.0],
        [2.0, 2.0, 0.5, 2.0],
    )
    @test projected_area_by_sector == [
        0.0 2.0
        2.0 0.0
        3.0 2.0
        0.0 0.0
    ]
    @test projected_area_active_indices == [2, 3, 1, 3]
    @test projected_area_active_offsets == [1, 3, 5]

    backend = ArchimedLight.KernelAbstractions.CPU()
    node_counts = fill(Int32(-1), 4)
    counts = Int32[2, 0, 3]
    nodes = UInt32[
        1, 3, 0,
        0, 0, 0,
        2, 3, 3,
    ]
    clear_kernel = ArchimedLight._raycore_clear_node_counts_kernel!(backend, 4)
    clear_kernel(node_counts; ndrange=length(node_counts))
    count_kernel = ArchimedLight._raycore_stack_node_counts_kernel!(backend, 4)
    count_kernel(node_counts, counts, nodes, 3, length(node_counts); ndrange=length(counts))
    ArchimedLight.KernelAbstractions.synchronize(backend)
    @test node_counts == Int32[1, 1, 3, 0]

    sector_area = fill(-1.0f0, 4)
    area_counts = Int32[3, 2]
    area_nodes = UInt32[
        1, 2, 3,
        2, 4, 0,
    ]
    virtual_node_mask = Bool[false, true, false, false]
    node_transparency = Float32[0.5, 0.0, 0.0, 0.25]
    clear_area_kernel = ArchimedLight._raycore_clear_sector_area_kernel!(backend, 4)
    clear_area_kernel(sector_area; ndrange=length(sector_area))
    area_kernel = ArchimedLight._raycore_stack_sector_area_kernel!(backend, 4)
    area_kernel(
        sector_area,
        area_counts,
        area_nodes,
        virtual_node_mask,
        node_transparency,
        2.0f0,
        3,
        length(sector_area);
        ndrange=length(area_counts),
    )
    ArchimedLight.KernelAbstractions.synchronize(backend)
    @test sector_area ≈ Float32[1.0, 4.0, 2.0, 1.5]

    fused_counts = fill(Int32(-1), 4)
    fused_area = fill(-1.0f0, 4)
    clear_fused_kernel = ArchimedLight._raycore_clear_node_counts_and_sector_area_kernel!(backend, 4)
    clear_fused_kernel(fused_counts, fused_area; ndrange=length(fused_counts))
    fused_kernel = ArchimedLight._raycore_stack_node_counts_and_sector_area_kernel!(backend, 4)
    fused_kernel(
        fused_counts,
        fused_area,
        area_counts,
        area_nodes,
        virtual_node_mask,
        node_transparency,
        2.0f0,
        3,
        length(fused_counts);
        ndrange=length(area_counts),
    )
    ArchimedLight.KernelAbstractions.synchronize(backend)
    @test fused_counts == Int32[1, 2, 1, 1]
    @test fused_area ≈ sector_area
end

@testitem "Raycore raw stack diagnostic script loads" tags=[:core, :fast, :raycore_backend] begin
    script_path = joinpath(dirname(dirname(pathof(ArchimedLight))), "scripts", "raycore_raw_stack_compare.jl")
    @test isfile(script_path)
    include(script_path)
    @test isdefined(@__MODULE__, :main)
    @test isdefined(@__MODULE__, :_backend_from_env)
end

@testitem "Raycore raw stack diagnostic CPU smoke" tags=[:core, :fast, :raycore_backend] begin
    script_path = joinpath(dirname(dirname(pathof(ArchimedLight))), "scripts", "raycore_raw_stack_compare.jl")
    env = Dict(
        "ARCHIMEDLIGHT_RAW_STACK_BACKEND" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_BACKEND", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_WORKLOAD" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_WORKLOAD", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_SECTOR" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_SECTOR", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS", nothing),
        "ARCHIMEDLIGHT_RAW_STACK_OUTPUT" => get(ENV, "ARCHIMEDLIGHT_RAW_STACK_OUTPUT", nothing),
    )
    try
        ENV["ARCHIMEDLIGHT_RAW_STACK_BACKEND"] = "cpu"
        ENV["ARCHIMEDLIGHT_RAW_STACK_WORKLOAD"] = "simple"
        ENV["ARCHIMEDLIGHT_RAW_STACK_MAX_HITS_SWEEP"] = "8"
        ENV["ARCHIMEDLIGHT_RAW_STACK_SECTOR"] = "first_scattering"
        ENV["ARCHIMEDLIGHT_RAW_STACK_PIXEL_SIZE"] = "0.40"
        ENV["ARCHIMEDLIGHT_RAW_STACK_FACE_CHUNK_LIMITS"] = "auto,unchunked,2"
        ENV["ARCHIMEDLIGHT_RAW_STACK_OUTPUT"] = ""

        include(script_path)
        rows = main()
        @test length(rows) == 3
        @test [row.requested_face_chunk_limit for row in rows] == ["auto", "unchunked", "2"]
        @test all(row -> row.workload == "simple", rows)
        @test all(row -> row.backend == "cpu", rows)
        @test all(row -> row.max_hits == 8, rows)
        @test all(row -> row.ok, rows)
        @test all(row -> row.hit_ratio == 1.0, rows)
        @test all(row -> row.occupied_ratio == 1.0, rows)
        @test all(row -> row.reference_overflow == row.candidate_overflow == 0, rows)
        @test all(row -> row.mismatch === nothing, rows)
        @test all(row -> row.reference_matches_candidate_topology, rows)
        @test rows[1].reference_geometry_mode == rows[1].candidate_geometry_mode
        @test rows[2].reference_geometry_mode == rows[2].candidate_geometry_mode == :merged_mesh
        @test rows[2].reference_face_chunk_limit == "unchunked"
        @test !rows[2].candidate_chunked
        @test rows[3].reference_geometry_mode == rows[3].candidate_geometry_mode == :chunked_merged_mesh
        @test rows[3].reference_face_chunk_limit == "2"
        @test rows[3].candidate_chunked
        @test all(row -> row.raycore_traversal_stack_capacity == ArchimedLight._raycore_software_traversal_stack_capacity(), rows)
        @test all(row -> !row.reference_blas_depth_over_stack_capacity, rows)
        @test all(row -> !row.candidate_blas_depth_over_stack_capacity, rows)
        @test all(row -> row.reference_blas_depth == row.candidate_blas_depth, rows)
        @test all(row -> row.reference_blas_node_count == row.candidate_blas_node_count, rows)
        @test all(row -> row.reference_tlas_depth == row.candidate_tlas_depth, rows)
        @test all(row -> row.reference_instances == row.candidate_instances, rows)
    finally
        for (key, value) in env
            if value === nothing
                delete!(ENV, key)
            else
                ENV[key] = value
            end
        end
    end
end

@testitem "Local realistic benchmark driver selftest" tags=[:core, :fast, :raycore_backend] begin
    script_path = joinpath(dirname(dirname(pathof(ArchimedLight))), "benchmark", "local_realistic_backends.jl")
    @test isfile(script_path)
    env = Dict(
        "ARCHIMEDLIGHT_LOCAL_BENCH_SELFTEST" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_SELFTEST", nothing),
        "ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT", nothing),
        "ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN", nothing),
        "ARCHIMEDLIGHT_LOCAL_BENCH_PREPARE_BREAKDOWN" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_PREPARE_BREAKDOWN", nothing),
        "ARCHIMEDLIGHT_LOCAL_BENCH_STACK_PROFILE" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_STACK_PROFILE", nothing),
        "ARCHIMEDLIGHT_LOCAL_BENCH_STAGE_SPLIT" => get(ENV, "ARCHIMEDLIGHT_LOCAL_BENCH_STAGE_SPLIT", nothing),
        "ARCHIMEDLIGHT_BENCH_VALIDATE" => get(ENV, "ARCHIMEDLIGHT_BENCH_VALIDATE", nothing),
    )
    try
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_SELFTEST"] = "1"
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_OUTPUT"] = ""
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_BREAKDOWN"] = ""
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_PREPARE_BREAKDOWN"] = ""
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_STACK_PROFILE"] = ""
        ENV["ARCHIMEDLIGHT_LOCAL_BENCH_STAGE_SPLIT"] = ""
        ENV["ARCHIMEDLIGHT_BENCH_VALIDATE"] = "1"

        rows = include(script_path)
        @test length(rows) == 8
        @test all(row -> row.status == "ok", rows)
        strict_row = only(row for row in rows if row.selftest == "prepare_breakdown_raycore_strict")
        @test strict_row.resolved_backend == "RaycoreInterceptionBackend"
        @test strict_row.fallback_reason == :none
        @test strict_row.raycore_instances > 0
        @test strict_row.stack_validation_hit_ratio == 1.0
        @test strict_row.stack_validation_occupied_ratio == 1.0
        @test !strict_row.stack_retry_used
    finally
        for (key, value) in env
            if value === nothing
                delete!(ENV, key)
            else
                ENV[key] = value
            end
        end
    end
end

@testitem "Scene parameter benchmark script parses" tags=[:core, :fast, :raycore_backend] begin
    script_path = joinpath(dirname(dirname(pathof(ArchimedLight))), "benchmark", "scene_parameter_sweep.jl")
    @test isfile(script_path)
    @test Meta.parse("begin\n" * read(script_path, String) * "\nend") isa Expr
end

@testitem "Raycore optional device smoke" tags=[:core, :raycore_backend] begin
    using GeometryBasics
    import PlantGeom

    KA = ArchimedLight.KernelAbstractions

    function _optional_package(name::Symbol)
        Base.find_package(String(name)) === nothing &&
            error("ARCHIMEDLIGHT_TEST_$(uppercase(String(name)))=1 was set, but package $name is not available.")
        return Base.require(Main, name)
    end

    function _optional_backend(name::Symbol)
        if name == :CUDA
            mod = _optional_package(:CUDA)
            getproperty(mod, :functional)() ||
                error("ARCHIMEDLIGHT_TEST_CUDA=1 was set, but CUDA.functional() is false.")
            dev_array = Base.invokelatest(getproperty(mod, :CuArray), zeros(Float32, 1))
            return Base.invokelatest(KA.get_backend, dev_array)
        elseif name == :Metal
            mod = _optional_package(:Metal)
            dev_array = Base.invokelatest(getproperty(mod, :MtlArray), zeros(Float32, 1))
            return Base.invokelatest(KA.get_backend, dev_array)
        elseif name == :oneAPI
            mod = _optional_package(:oneAPI)
            dev_array = Base.invokelatest(getproperty(mod, :oneArray), zeros(Float32, 1))
            return Base.invokelatest(KA.get_backend, dev_array)
        elseif name == :AMDGPU
            mod = _optional_package(:AMDGPU)
            dev_array = Base.invokelatest(getproperty(mod, :ROCArray), zeros(Float32, 1))
            return Base.invokelatest(KA.get_backend, dev_array)
        end
        error("Unsupported optional backend $name")
    end

    function _smoke_scene()
        mesh = GeometryBasics.Mesh(
            GeometryBasics.Point3f[
                GeometryBasics.Point3f(0, 0, 1),
                GeometryBasics.Point3f(1, 0, 1),
                GeometryBasics.Point3f(1, 1, 1),
                GeometryBasics.Point3f(0, 1, 1),
            ],
            GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)],
        )
        return PlantGeom.make_scene(domain=(0.0, 0.0, 1.0, 1.0), source_path="raycore_optional_device_smoke") do builder
            PlantGeom.add_object!(
                builder,
                mesh;
                group="plate",
                type="plate",
                id=1,
                source_topology_id=1,
            )
        end
    end

    function _smoke_models()
        return ArchimedLight.prepare_models([
            ArchimedLight.GroupModel(
                "*";
                types=ArchimedLight.OrderedDict(
                    "*" => ArchimedLight.TypeModel(
                        interception=ArchimedLight.InterceptionModel(
                            model="Translucent",
                            transparency=0.0,
                            optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                        ),
                    ),
                ),
            ),
        ])
    end

    requested = Pair{String,Symbol}[
        "ARCHIMEDLIGHT_TEST_CUDA" => :CUDA,
        "ARCHIMEDLIGHT_TEST_METAL" => :Metal,
        "ARCHIMEDLIGHT_TEST_ONEAPI" => :oneAPI,
        "ARCHIMEDLIGHT_TEST_AMDGPU" => :AMDGPU,
    ]
    selected = [
        env => backend for (env, backend) in requested
        if lowercase(get(ENV, env, "")) in ("1", "true", "yes", "on", "required")
    ]

    if isempty(selected)
        @test true
    else
        scene = _smoke_scene()
        models = _smoke_models()
        options = ArchimedLight.LightOptions(
            turtle_sectors=1,
            all_in_turtle=false,
            scattering=false,
            pixel_size=0.01,
        )
        sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
        turtle = ArchimedLight.build_turtle(options, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
        raster = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
        raster_area = sum(values(raster.projected_area_per_node); init=0.0)
        raster_par = sum(values(raster.incident_power.par); init=0.0)

        for (env, name) in selected
            required = lowercase(get(ENV, env, "")) == "required"
            backend = try
                _optional_backend(name)
            catch err
                required && rethrow()
                @test_skip "Optional backend $name could not be initialized: $(sprint(showerror, err))"
                continue
            end
            raycore = try
                ArchimedLight.compute_first_order(
                    scene,
                    models,
                    turtle,
                    fluxes,
                    options;
                    backend=ArchimedLight.RaycoreInterceptionBackend(backend=backend),
                )
            catch err
                required && rethrow()
                @test_skip "Optional backend $name could not build Raycore data: $(sprint(showerror, err))"
                continue
            end
            raycore_area = sum(values(raycore.projected_area_per_node); init=0.0)
            raycore_par = sum(values(raycore.incident_power.par); init=0.0)
            @test raycore_area > 0.0
            @test isapprox(raycore_area, raster_area; atol=1e-10, rtol=1e-10)
            @test isapprox(raycore_par, raster_par; atol=1e-8, rtol=1e-8)

            data = ArchimedLight._prepare_raycore_interception_data(
                scene,
                models,
                options,
                ArchimedLight.RaycoreInterceptionBackend(
                    backend=backend,
                    max_hits_per_pixel=8,
                    edge_accumulation=:dense_atomic,
                ),
            )
            sector_idx = findfirst(sector -> sector.direction[3] > 0.0, turtle.sectors)
            @test sector_idx !== nothing
            traced_device = ArchimedLight._raycore_trace_direction_stack_nodes_device(
                data,
                turtle.sectors[sector_idx].direction,
                options,
            )
            @test !any(traced_device.overflow)
            device_counts = ArchimedLight._raycore_stack_node_counts_from_device_traced_stacks(data, traced_device)
            if device_counts === nothing
                @test_skip "Selected optional backend does not support device node-count reduction."
            else
                traced_host = ArchimedLight._raycore_copy_direction_stack_nodes_host!(data, traced_device)
                expected_counts = zeros(Int32, length(data.prepared.geometry.node_ids))
                @inbounds for pixel_idx in eachindex(traced_host.counts)
                    n_hits = Int(traced_host.counts[pixel_idx])
                    stack_base = (pixel_idx - 1) * traced_host.max_hits
                    for h in 1:n_hits
                        node_idx = Int(traced_host.nodes[stack_base+h])
                        expected_counts[node_idx] += Int32(1)
                    end
                end
                @test device_counts == expected_counts
                @test sum(device_counts) > 0
            end
        end
    end
end
