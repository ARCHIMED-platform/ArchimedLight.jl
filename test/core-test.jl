@testitem "Core smoke" tags=[:core, :fast] begin
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
    @test !isempty([nid for (nid, node) in scene2.nodes if node.group == "pavement"])
    @test haskey(scene2.mtg, :geometry)
    @test scene2.mtg[:geometry] === nothing

    _, scene3, _, _ = ArchimedLight.read_config(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml");
        plot_paving_override=25,
    )
    @test length([nid for (nid, node) in scene3.nodes if node.group == "pavement"]) == 25
    @test scene3.mtg[:geometry] === nothing

    raw_scene = ArchimedLight.read_scene(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "scene", "simple.ops"),
    )
    @test haskey(raw_scene.mtg, :geometry)
    @test raw_scene.mtg[:geometry] !== nothing
    ArchimedLight.add_ground!(raw_scene; nx=2, ny=2)
    @test raw_scene.mtg[:geometry] === nothing
end

@testitem "Sky fraction can be stored and attached" tags=[:core, :fast] begin
    import MultiScaleTreeGraph

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    meteo = ArchimedLight.prepare_meteo(fixture.meteo, fixture.options)
    series = ArchimedLight.run_light_series(
        fixture.scene,
        fixture.models,
        meteo,
        fixture.options;
        include_sky_fraction=true,
    )

    @test length(series) == length(meteo.rows)
    @test series[1].sky_fraction !== nothing
    @test !isempty(series[1].sky_fraction)

    ArchimedLight.attach_light_series!(
        fixture.scene,
        series;
        fields=[:absorbed_par_flux, :absorbed_nir_flux, :sky_fraction],
        names=Dict(:absorbed_nir_flux => :Ra_SW_f),
    )

    found = Ref(false)
    MultiScaleTreeGraph.traverse!(fixture.scene.mtg) do node
        if haskey(node, :sky_fraction)
            found[] = true
            @test node[:sky_fraction] isa Vector{Float64}
            @test length(node[:sky_fraction]) == length(series)
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
