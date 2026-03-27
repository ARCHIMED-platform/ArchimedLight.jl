@testitem "Makie extension lightplot colors and timestep updates" tags=[:core, :fast] begin
    using CairoMakie

    include(joinpath(@__DIR__, "support.jl"))

    scale_map(values::Dict{Int,Float64}, factor::Float64) = Dict(nid => factor * value for (nid, value) in values)
    scale_spectral(values::ArchimedLight.SpectralNodeValues, factor::Float64) =
        ArchimedLight.SpectralNodeValues(scale_map(values.par, factor), scale_map(values.nir, factor))
    scale_initial_total(values::ArchimedLight.InitialTotalSpectralNodeValues, factor::Float64) =
        ArchimedLight.InitialTotalSpectralNodeValues(
            scale_spectral(values.initial, factor),
            scale_spectral(values.total, factor),
        )

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
    step = ArchimedLight.run_light_step(fixture.scene, fixture.models, row, fixture.options)
    render_geometry = ArchimedLight.light_render_geometry(step)

    budget2 = ArchimedLight.LightBudget(
        scale_initial_total(step.budget.incident_flux, 2.0),
        scale_initial_total(step.budget.incident_energy, 2.0),
        scale_initial_total(step.budget.absorbed_flux, 2.0),
        scale_initial_total(step.budget.absorbed_energy, 2.0),
        Dict{String,Dict{Int,Float64}}(band => scale_map(values, 2.0) for (band, values) in step.budget.extra_initial_energy_per_band),
        Dict{String,Dict{Int,Float64}}(band => scale_map(values, 2.0) for (band, values) in step.budget.extra_energy_per_band),
    )
    step2 = ArchimedLight.LightStepResult(
        step.sky,
        step.turtle,
        step.fluxes,
        step.first_order,
        step.scattering,
        budget2,
        step.extra_band_irradiance,
        render_geometry,
    )

    fig, ax, plt = ArchimedLight.lightplot(
        [step, step2];
        color=:Ri_PAR_f,
        timestep=1,
        colorrange=(0.0, 2.0 * maximum(values(step.budget.incident_flux.total.par))),
    )

    observed_1 = copy(Makie.to_value(plt[:light_color]))
    expected_1 = ArchimedLight.light_face_values(step; color=:Ri_PAR_f)

    @test ax isa LScene
    @test render_geometry === step.render_geometry
    @test length(observed_1) == length(expected_1)
    @test all(isapprox.(observed_1, expected_1; atol=1f-6, rtol=1f-6))
    expected_range = (
        min(minimum(expected_1), minimum(ArchimedLight.light_face_values(step2; color=:Ri_PAR_f))),
        max(maximum(expected_1), maximum(ArchimedLight.light_face_values(step2; color=:Ri_PAR_f))),
    )
    observed_range_1 = Makie.to_value(only(plt.plots).colorrange)
    @test all(isapprox.(collect(observed_range_1), collect(expected_range); atol=1f-6, rtol=1f-6))

    plt[:timestep][] = 2
    observed_2 = copy(Makie.to_value(plt[:light_color]))
    expected_2 = ArchimedLight.light_face_values([step, step2]; color=:Ri_PAR_f, timestep=2)

    @test all(isapprox.(observed_2, expected_2; atol=1f-6, rtol=1f-6))
    @test any(!isapprox(a, b) for (a, b) in zip(observed_1, observed_2))
    observed_range_2 = Makie.to_value(only(plt.plots).colorrange)
    @test all(isapprox.(collect(observed_range_2), collect(expected_range); atol=1f-6, rtol=1f-6))

    fig_s, _, plt_s = ArchimedLight.lightplot(step; color=:Ri_PAR_f)
    child_s = only(plt_s.plots)
    observed_s = copy(Makie.to_value(plt_s[:light_color]))
    expected_s = ArchimedLight.light_face_values(step; color=:Ri_PAR_f)

    @test all(isapprox.(observed_s, expected_s; atol=1f-6, rtol=1f-6))
    @test length(Makie.to_value(child_s[1]).color) == 3 * length(expected_s)

    tiled_geometry = ArchimedLight.tile_light_geometry(fixture.scene, step; nx=2, ny=2, centered=false)
    @test length(tiled_geometry.vertices) == 4 * length(render_geometry.vertices)
    @test length(tiled_geometry.faces) == 4 * length(render_geometry.faces)
    @test tiled_geometry.face2node == repeat(render_geometry.face2node, 4)

    fig_t, _, plt_t = ArchimedLight.lightplot(tiled_geometry, step; color=:Ri_PAR_f)
    tiled_expected = ArchimedLight.light_face_values(tiled_geometry, step; color=:Ri_PAR_f)
    @test all(isapprox.(Makie.to_value(plt_t[:light_color]), tiled_expected; atol=1f-6, rtol=1f-6))

    fig_v, _, plt_v = ArchimedLight.lightplot(step; color=:Ri_PAR_f, interpolate=true)
    observed_v = copy(Makie.to_value(plt_v[:light_color]))
    expected_v = ArchimedLight.light_vertex_values(step; color=:Ri_PAR_f)

    @test length(observed_v) == length(expected_v)
    @test all(isapprox.(observed_v, expected_v; atol=1f-6, rtol=1f-6))

    # Compatibility path for explicit metric maps still works when geometry is provided separately.
    node_ids = ArchimedLight.scene_node_ids(fixture.scene)
    metric_1 = Dict(nid => Float64(i) for (i, nid) in enumerate(node_ids))
    metric_fig, _, metric_plt = ArchimedLight.lightplot(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_1;
        color=:Ri_PAR_f,
    )
    metric_child = only(metric_plt.plots)
    metric_expected = ArchimedLight.light_face_values(
        fixture.scene,
        fixture.models,
        fixture.options,
        metric_1;
        color=:Ri_PAR_f,
    )
    @test all(isapprox.(Makie.to_value(metric_plt[:light_color]), metric_expected; atol=1f-6, rtol=1f-6))
    @test length(Makie.to_value(metric_child[1]).color) == 3 * length(metric_expected)

    mktemp() do path, io
        close(io)
        png_path = string(path, ".png")
        save(png_path, fig_s)
        @test isfile(png_path)
    end

    CairoMakie.empty!(fig)
    CairoMakie.empty!(fig_s)
    CairoMakie.empty!(fig_t)
    CairoMakie.empty!(fig_v)
    CairoMakie.empty!(metric_fig)
end
