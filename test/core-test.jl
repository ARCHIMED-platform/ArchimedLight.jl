@testitem "Core smoke" tags=[:core, :fast] begin
    import PlantGeom

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    selected = ArchimedLight.prepare_meteo(fixture.meteo, fixture.options)
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)

    @test !isempty(selected)
    @test !isempty(series)
    @test length(series) == length(selected)
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
    @test length(meteo2) == length(fixture.meteo)
    @test collect(keys(models2.groups)) == collect(keys(fixture.models.groups))
    @test any(item -> item.group == "pavement", ArchimedLight.summarize_scene(scene2).group_types)
    @test !PlantGeom.has_geometry(scene2.mtg)

    _, scene3, _, _ = ArchimedLight.read_config(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml");
        plot_paving_override=25,
    )
    pavement = only(item for item in ArchimedLight.summarize_scene(scene3).group_types if item.group == "pavement")
    @test pavement.nodes == 25
    @test !PlantGeom.has_geometry(scene3.mtg)

    raw_scene = ArchimedLight.read_scene(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "scene", "simple.ops"),
    )
    @test !PlantGeom.has_geometry(raw_scene.mtg)
    PlantGeom.add_ground!(raw_scene; nx=2, ny=2)
    @test !PlantGeom.has_geometry(raw_scene.mtg)
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

    @test length(series) == length(meteo)
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

    dense = ArchimedLight.DensePixelHits(ArchimedLight.SmallHitStack, 3)
    @test dense.stacks isa Vector{Union{Nothing,ArchimedLight.SmallHitStack}}
    @test all(isnothing, dense.stacks)
end
