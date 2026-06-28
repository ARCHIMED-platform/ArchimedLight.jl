using KernelAbstractions
using Test

const METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_METAL", ""))
const METAL_REQUESTED = METAL_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const METAL_REQUIRED = METAL_FLAG in ("required", "force", "error")

function _skip_or_fail(message::AbstractString)
    if METAL_REQUIRED
        error(message)
    else
        @test_skip message
        return nothing
    end
end

function _metal_backend()
    if !Sys.isapple()
        return _skip_or_fail("Metal validation requires macOS.")
    end
    if Sys.ARCH != :aarch64
        return _skip_or_fail("Metal validation requires Apple Silicon (Sys.ARCH == :aarch64).")
    end
    if Base.find_package("Metal") === nothing
        return _skip_or_fail("Metal.jl is not available in this optional test environment. Run `julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'`.")
    end

    metal_mod = try
        Base.require(Main, :Metal)
    catch err
        return _skip_or_fail("Metal.jl could not be imported: $(sprint(showerror, err))")
    end

    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = try
        Base.invokelatest(array_type, zeros(Float32, 1))
    catch err
        return _skip_or_fail("Metal.jl is installed, but creating a Metal array failed: $(sprint(showerror, err))")
    end
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

@testset "Raycore Metal validation" begin
    if !METAL_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to run optional Metal validation."
    else
        backend = _metal_backend()
        if backend !== nothing
            include(joinpath(@__DIR__, "..", "..", "scripts", "raycore_device_smoke.jl"))
            result = run_raycore_device_smoke(;
                backend=backend,
                first_order_area_atol=1e-5,
                first_order_area_rtol=1e-5,
                first_order_power_atol=1e-4,
                first_order_power_rtol=1e-5,
                scattering_atol=1e-4,
                scattering_rtol=1e-4,
            )
            @test result.scattering_eltype == Float32
        end
    end
end
