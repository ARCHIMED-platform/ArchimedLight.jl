@testitem "Release artifact builder output pruning" tags=[:release_artifact, :fast] begin
    include(joinpath(@__DIR__, "..", "scripts", "build_release_fixture_artifact.jl"))

    src = mktempdir()
    for name in ("output", "output2", "output_debug", "output-cache", "inputs")
        dir = joinpath(src, name)
        mkpath(dir)
        write(joinpath(dir, "marker.txt"), name)
    end
    write(joinpath(src, "output.txt"), "regular file")

    pruned = mktempdir()
    _copy_tree_filtered!(src, pruned; prune_output=true)

    @test isdir(joinpath(pruned, "inputs"))
    @test isfile(joinpath(pruned, "output.txt"))
    for name in ("output", "output2", "output_debug", "output-cache")
        @test !ispath(joinpath(pruned, name))
    end

    unpruned = mktempdir()
    _copy_tree_filtered!(src, unpruned; prune_output=false)

    for name in ("output", "output2", "output_debug", "output-cache", "inputs")
        @test isfile(joinpath(unpruned, name, "marker.txt"))
    end
end
