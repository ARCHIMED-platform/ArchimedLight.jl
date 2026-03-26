@testitem "Sky fixture sky_06_direct" tags=[:fast_fixture, :fast, :sky_06_direct] begin
    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "sky_06_direct", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
    sky = ArchimedLight.compute_sky(row, fixture.options)
    turtle = ArchimedLight.build_turtle(fixture.options, sky)
    flux = ArchimedLight.compute_directional_fluxes(sky, turtle, fixture.options)

    @test fixture.options.turtle_sectors == 6
    @test fixture.options.all_in_turtle == false
    @test length(turtle.sectors) == 7
    @test count(s -> s.source == :sun, turtle.sectors) == 1
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
end

@testitem "Sky fixture sky_16_turtle" tags=[:fast_fixture, :fast, :sky_16_turtle] begin
    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "sky_16_turtle", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
    sky = ArchimedLight.compute_sky(row, fixture.options)
    turtle = ArchimedLight.build_turtle(fixture.options, sky)
    flux = ArchimedLight.compute_directional_fluxes(sky, turtle, fixture.options)

    @test fixture.options.turtle_sectors == 16
    @test fixture.options.all_in_turtle == true
    @test length(turtle.sectors) == 16
    @test count(s -> s.source == :sun, turtle.sectors) == 0
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
end

@testitem "Sky fixture sky_46_direct" tags=[:fast_fixture, :fast, :sky_46_direct] begin
    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "sky_46_direct", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
    sky = ArchimedLight.compute_sky(row, fixture.options)
    turtle = ArchimedLight.build_turtle(fixture.options, sky)
    flux = ArchimedLight.compute_directional_fluxes(sky, turtle, fixture.options)

    @test fixture.options.turtle_sectors == 46
    @test fixture.options.all_in_turtle == false
    @test length(turtle.sectors) == 47
    @test count(s -> s.source == :sun, turtle.sectors) == 1
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1e-6, rtol=1e-6)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
end

@testitem "Fast fixture simpleplant_16_notoric" tags=[:fast_fixture, :fast, :simpleplant_16_notoric] begin
    using CairoMakie
    using CSV
    using ReferenceTests
    using Tables

    include(joinpath(@__DIR__, "support.jl"))
    include(joinpath(@__DIR__, "fast_fixture_support.jl"))

    case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric")
    fixture = load_fixture_inputs(joinpath(case_root, "input"))
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)
    step = first(series)

    expected_csv = joinpath(case_root, "expected", "component_values.csv")
    expected_rows = collect(Tables.rowtable(CSV.File(expected_csv; delim=';', normalizenames=false)))
    observed_map = _source_topology_area_q(fixture.scene, fixture.models, fixture.options, step)

    @test length(observed_map) == length(expected_rows)
    for row in expected_rows
        key = (Int(row.source_topology_id), Int(row.object_id))
        @test haskey(observed_map, key)
        observed = observed_map[key]
        @test isapprox(observed.area, Float64(row.area); atol=1e-8, rtol=1e-6)
        @test isapprox(observed.Ri_PAR_0_q, Float64(row.Ri_PAR_0_q); atol=1e-4, rtol=1e-6)
    end

    fig = _render_ri_par_f_figure(fixture.scene, fixture.models, fixture.options, step; title="simpleplant_16_notoric | Ri_PAR_f")
    ref_png = joinpath(case_root, "expected", "ri_par_f_step0.png")
    @test_reference relpath(ref_png, @__DIR__) fig by = ReferenceTests.psnr_equality(35)
end

@testitem "Fast fixture simpleplant_16_toric" tags=[:fast_fixture, :fast, :simpleplant_16_toric] begin
    using CairoMakie
    using CSV
    using ReferenceTests
    using Tables

    include(joinpath(@__DIR__, "support.jl"))
    include(joinpath(@__DIR__, "fast_fixture_support.jl"))

    case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_toric")
    fixture = load_fixture_inputs(joinpath(case_root, "input"))
    series = ArchimedLight.run_light_series(fixture.scene, fixture.models, fixture.meteo, fixture.options)
    step = first(series)

    expected_csv = joinpath(case_root, "expected", "component_values.csv")
    expected_rows = collect(Tables.rowtable(CSV.File(expected_csv; delim=';', normalizenames=false)))
    observed_map = _source_topology_area_q(fixture.scene, fixture.models, fixture.options, step)

    @test length(observed_map) == length(expected_rows)
    for row in expected_rows
        key = (Int(row.source_topology_id), Int(row.object_id))
        @test haskey(observed_map, key)
        observed = observed_map[key]
        @test isapprox(observed.area, Float64(row.area); atol=1e-8, rtol=1e-6)
        @test isapprox(observed.Ri_PAR_0_q, Float64(row.Ri_PAR_0_q); atol=1e-4, rtol=1e-6)
    end

    fig = _render_ri_par_f_figure(fixture.scene, fixture.models, fixture.options, step; title="simpleplant_16_toric | Ri_PAR_f")
    ref_png = joinpath(case_root, "expected", "ri_par_f_step0.png")
    @test_reference relpath(ref_png, @__DIR__) fig by = ReferenceTests.psnr_equality(35)
end
