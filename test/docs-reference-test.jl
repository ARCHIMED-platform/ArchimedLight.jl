@testitem "Docs index coffee figure reference" tags=[:release, :slow, :docs, :docs_index_figure] begin
    using ReferenceTests

    include(joinpath(@__DIR__, "..", "scripts", "generate_home_figure.jl"))

    fig = build_home_figure()
    ref_png = joinpath(@__DIR__, "..", "docs", "src", "assets", "coffee_scene_light_interception.png")
    @test_reference ref_png fig by = ReferenceTests.psnr_equality(35)
end
