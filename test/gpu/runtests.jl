pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..")))

using ArchimedLight
using KernelAbstractions
import Metal
using Test

const METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_METAL", ""))
const METAL_REQUESTED = METAL_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const METAL_REQUIRED = METAL_FLAG in ("required", "force", "error")

function _skip_or_fail(message::AbstractString)
    if METAL_REQUIRED
        error(message)
    end
    @test_skip message
    return nothing
end

function _metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        return nothing, nothing, "Metal tests require Apple Silicon."
    end
    metal = Metal
    Metal.functional() ||
        return nothing, metal, "Metal.functional() returned false."

    array_type = Metal.MtlArray
    backend = KernelAbstractions.get_backend(array_type(zeros(Float32, 1)))
    return backend, metal, nothing
end

@testset "RasterGPU Metal" begin
    if !METAL_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to run Metal validation."
    else
        backend, metal, failure = _metal_backend()
        if backend === nothing
            _skip_or_fail(failure)
        else
            @test Base.pkgversion(metal) >= v"1.10.3"
            atomic_support = KernelAbstractions.supports_atomics(backend)
            if atomic_support
                @test atomic_support
                config_path = joinpath(
                    @__DIR__,
                    "..",
                    "fast_fixtures",
                    "simpleplant_16_notoric",
                    "input",
                    "config.yml",
                )
                options, scene, meteo, models = ArchimedLight.read_config(config_path)
                options = ArchimedLight.LightOptions(
                    options;
                    scattering=true,
                    cache_radiation=false,
                    cache_pixel_table=false,
                    toricity=false,
                    turtle_sectors=6,
                    scattering_max_iter=3,
                )
                row = first(meteo.rows)
                sky = ArchimedLight.compute_sky(row, options)
                turtle = ArchimedLight.build_turtle(options, sky)
                fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)

                cpu_first = ArchimedLight.compute_first_order(
                    scene,
                    models,
                    turtle,
                    fluxes,
                    options;
                    backend=:raster_cpu,
                )
                interception_backend = ArchimedLight.RasterGPUBackend(
                    backend=backend,
                    max_hits_per_pixel=128,
                    tile_size=1,
                    tile_face_capacity=512,
                    edge_accumulation=:auto,
                    validate=true,
                )
                gpu_first = ArchimedLight.compute_first_order(
                    scene,
                    models,
                    turtle,
                    fluxes,
                    options;
                    backend=interception_backend,
                )

                @test gpu_first.dense !== nothing
                @test cpu_first.dense.node_ids == gpu_first.dense.node_ids
                @test cpu_first.dense.hits_per_node == gpu_first.dense.hits_per_node
                @test cpu_first.dense.projected_area_per_node ≈ gpu_first.dense.projected_area_per_node atol = 1e-5 rtol = 1e-5
                @test cpu_first.dense.incident_power.par ≈ gpu_first.dense.incident_power.par atol = 1e-4 rtol = 1e-5
                @test cpu_first.dense.incident_power.nir ≈ gpu_first.dense.incident_power.nir atol = 1e-4 rtol = 1e-5

                prepared = ArchimedLight._prepare_interception_data(
                    scene,
                    models,
                    options;
                    include_budget_maps=true,
                )
                gpu_data = ArchimedLight._rastergpu_scene_data(prepared, interception_backend.config)
                @test ArchimedLight._rastergpu_use_fused_dense_edges(gpu_data)

                cpu_scattering = ArchimedLight.compute_scattering(
                    scene,
                    models,
                    turtle,
                    cpu_first,
                    options;
                    backend=ArchimedLight.RaycastScatteringBackend(),
                )
                gpu_scattering = ArchimedLight.compute_scattering(
                    scene,
                    models,
                    turtle,
                    gpu_first,
                    options;
                    backend=ArchimedLight.RasterGPUScatteringBackend(interception_backend),
                )

                @test gpu_scattering.dense !== nothing
                @test cpu_scattering.iterations == gpu_scattering.iterations
                @test cpu_scattering.converged == gpu_scattering.converged
                @test cpu_scattering.dense.node_ids == gpu_scattering.dense.node_ids
                @test cpu_scattering.dense.added_power.par ≈ gpu_scattering.dense.added_power.par atol = 1e-4 rtol = 1e-4
                @test cpu_scattering.dense.added_power.nir ≈ gpu_scattering.dense.added_power.nir atol = 1e-4 rtol = 1e-4
            elseif METAL_REQUIRED
                @test atomic_support
            else
                @test_skip "The active Metal toolchain does not provide KernelAbstractions atomics (MSL 4.1 or newer is required)."
            end
        end
    end
end
