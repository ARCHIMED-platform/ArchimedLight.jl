@testitem "User-facing file example" tags=[:example, :fast] begin
    using ArchimedLight

    example_path = joinpath(@__DIR__, "..", "example_1", "full_featured_example.jl")
    example_module = Module(:ArchimedLightUserExample)
    Base.include(example_module, example_path)

    output_dir = mktempdir()
    step = example_module.run_example(; output_dir=output_dir, io=devnull)

    @test step isa LightStepResult
    @test isfile(joinpath(output_dir, "component_values.csv"))
    @test isfile(joinpath(output_dir, "scene_with_light.opf"))
end
