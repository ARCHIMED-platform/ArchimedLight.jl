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
    @test ib.config.edge_accumulation == :auto
    @test ib.config.dense_edge_limit_bytes == 512 * 1024^2
    @test ib.config.scattering_eltype == Float64

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
    @test sb.config.scattering_eltype == Float64
    @test ArchimedLight._resolve_scattering_backend(:raycast, sb) === sb

    ib32 = ArchimedLight.RaycoreInterceptionBackend(scattering_eltype=Float32)
    sb32 = ArchimedLight.RaycoreScatteringBackend(ib32)
    @test ib32.config.scattering_eltype == Float32
    @test sb32.config.scattering_eltype == Float32

    @test_throws ErrorException ArchimedLight.RaycoreBackendConfig(workgroupsize=0)
    @test_throws ErrorException ArchimedLight.RaycoreBackendConfig(edge_accumulation=:unsupported)
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
            return KA.get_backend(getproperty(mod, :CuArray)(zeros(Float32, 1)))
        elseif name == :Metal
            mod = _optional_package(:Metal)
            return KA.get_backend(getproperty(mod, :MtlArray)(zeros(Float32, 1)))
        elseif name == :oneAPI
            mod = _optional_package(:oneAPI)
            return KA.get_backend(getproperty(mod, :oneArray)(zeros(Float32, 1)))
        elseif name == :AMDGPU
            mod = _optional_package(:AMDGPU)
            return KA.get_backend(getproperty(mod, :ROCArray)(zeros(Float32, 1)))
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
    selected = [backend for (env, backend) in requested if lowercase(get(ENV, env, "")) in ("1", "true", "yes", "on")]

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

        for name in selected
            backend = _optional_backend(name)
            raycore = ArchimedLight.compute_first_order(
                scene,
                models,
                turtle,
                fluxes,
                options;
                backend=ArchimedLight.RaycoreInterceptionBackend(backend=backend),
            )
            raycore_area = sum(values(raycore.projected_area_per_node); init=0.0)
            raycore_par = sum(values(raycore.incident_power.par); init=0.0)
            @test raycore_area > 0.0
            @test isapprox(raycore_area, raster_area; atol=1e-10, rtol=1e-10)
            @test isapprox(raycore_par, raster_par; atol=1e-8, rtol=1e-8)
        end
    end
end
