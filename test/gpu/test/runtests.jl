if !haskey(ENV, "ARCHIMEDLIGHT_TEST_METAL") &&
   Sys.isapple() &&
   Sys.ARCH == :aarch64 &&
   Base.find_package("Metal") !== nothing
    ENV["ARCHIMEDLIGHT_TEST_METAL"] = "1"
end

include(joinpath(@__DIR__, "..", "runtests.jl"))
