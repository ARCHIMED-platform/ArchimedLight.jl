@testitem "Makie extension lightplot colors and timestep updates" tags=[:core, :fast] begin
    using CairoMakie

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
    step = ArchimedLight.run_light_step(fixture.scene, fixture.models, row, fixture.options)

    node_ids = ArchimedLight.scene_node_ids(fixture.scene)
    metric_1 = Dict(nid => Float64(i) for (i, nid) in enumerate(node_ids))
    metric_2 = Dict(nid => 10.0 + Float64(i) for (i, nid) in enumerate(node_ids))
    metrics = [metric_1, metric_2]

    fig, ax, plt = ArchimedLight.lightplot(
        fixture.scene,
        fixture.models,
        fixture.options,
        metrics;
        color=:Ri_PAR_f,
        timestep=1,
        colorrange=(0.0, maximum(values(metric_2))),
    )

    child = only(plt.plots)
    observed_1 = copy(Makie.to_value(child.color))
    expected_1 = ArchimedLight.light_face_values(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_1;
        color=:Ri_PAR_f,
    )

    @test ax isa Axis3
    @test length(observed_1) == length(expected_1)
    @test all(isapprox.(observed_1, expected_1; atol=1f-6, rtol=1f-6))

    plt[:timestep][] = 2
    observed_2 = copy(Makie.to_value(child.color))
    expected_2 = ArchimedLight.light_face_values(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_2;
        color=:Ri_PAR_f,
    )

    @test all(isapprox.(observed_2, expected_2; atol=1f-6, rtol=1f-6))
    @test any(!isapprox(a, b) for (a, b) in zip(observed_1, observed_2))

    fig_s, _, plt_s = ArchimedLight.lightplot(
        fixture.scene,
        fixture.models,
        fixture.options,
        step;
        color=:Ri_PAR_f,
    )
    child_s = only(plt_s.plots)
    observed_s = copy(Makie.to_value(child_s.color))
    expected_s = ArchimedLight.light_face_values(
        fixture.scene,
        fixture.models,
        fixture.options,
        step;
        color=:Ri_PAR_f,
    )

    @test all(isapprox.(observed_s, expected_s; atol=1f-6, rtol=1f-6))

    fig_v, _, plt_v = ArchimedLight.lightplot(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_1;
        color=:Ri_PAR_f,
        interpolate=true,
    )
    child_v = only(plt_v.plots)
    observed_v = copy(Makie.to_value(child_v.color))
    expected_v = ArchimedLight.light_vertex_values(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_1;
        color=:Ri_PAR_f,
    )

    @test length(observed_v) == length(expected_v)
    @test all(isapprox.(observed_v, expected_v; atol=1f-6, rtol=1f-6))

    CairoMakie.empty!(fig)
    CairoMakie.empty!(fig_s)
    CairoMakie.empty!(fig_v)
end
