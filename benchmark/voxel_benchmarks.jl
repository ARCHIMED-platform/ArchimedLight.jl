# This file is included by `benchmarks.jl` after ArchimedLight,
# BenchmarkTools, and the top-level `SUITE` have been defined.

const VOXEL_BENCHMARK_METADATA = (
    julia_version=string(VERSION),
    boundary=:periodic,
    traversal=:dda,
    rays_per_voxel=4,
    turtle_sectors=16,
    repeated_steps=48,
    sizes=(small=(8, 8, 8), medium=(16, 16, 16), large=(32, 32, 32)),
)

function _voxel_benchmark_grid(dimensions::NTuple{3,Int})
    nx, ny, nz = dimensions
    pad = Array{Float64}(undef, dimensions)
    for k in 1:nz, j in 1:ny, i in 1:nx
        pad[i, j, k] = 0.2 + 0.3 * (i + 2j + 3k) / (nx + 2ny + 3nz)
    end
    return ArchimedLight.VoxelGrid((0.0, 0.0, 0.0), (Float64(nx), Float64(ny), Float64(nz)), pad)
end

function _voxel_apply_many_steps(
    grid,
    turtle,
    fluxes,
    backend,
    cache,
    count::Int,
)
    result = nothing
    for _ in 1:count
        result = ArchimedLight.compute_voxel_first_order(
            grid,
            turtle,
            fluxes,
            backend;
            duration_seconds=1800.0,
            cache=cache,
        )
    end
    return result
end

const VOXEL_SMALL = _voxel_benchmark_grid(VOXEL_BENCHMARK_METADATA.sizes.small)
const VOXEL_MEDIUM = _voxel_benchmark_grid(VOXEL_BENCHMARK_METADATA.sizes.medium)
const VOXEL_LARGE = _voxel_benchmark_grid(VOXEL_BENCHMARK_METADATA.sizes.large)
const VOXEL_BACKEND = ArchimedLight.VoxelCPUBackend(
    rays_per_voxel=VOXEL_BENCHMARK_METADATA.rays_per_voxel,
    boundary=VOXEL_BENCHMARK_METADATA.boundary,
    traversal=VOXEL_BENCHMARK_METADATA.traversal,
)
const VOXEL_ROW = (
    step_duration=1800.0,
    Ri_PAR_f=120.0,
    Ri_NIR_f=60.0,
    sun_azimuth=190.0,
    sun_elevation=40.0,
    direct_fraction=0.6,
)
const VOXEL_OPTIONS = ArchimedLight.LightOptions(
    turtle_sectors=VOXEL_BENCHMARK_METADATA.turtle_sectors,
    all_in_turtle=true,
    scattering=false,
    radiation_timestep_minutes=5.0,
)
const VOXEL_SKY = ArchimedLight.compute_sky(VOXEL_ROW, VOXEL_OPTIONS)
const VOXEL_TURTLE = ArchimedLight.build_turtle(VOXEL_OPTIONS, VOXEL_SKY)
const VOXEL_FLUXES = ArchimedLight.compute_directional_fluxes(
    VOXEL_ROW,
    VOXEL_SKY,
    VOXEL_TURTLE,
    VOXEL_OPTIONS,
)
const VOXEL_MEDIUM_CACHE = ArchimedLight.prepare_voxel_responses(
    VOXEL_MEDIUM,
    VOXEL_TURTLE,
    VOXEL_BACKEND,
)
const VOXEL_SERIES_ROWS = [
    merge(
        VOXEL_ROW,
        (
            Ri_PAR_f=120.0 - index,
            Ri_NIR_f=60.0 - index / 2,
            sun_azimuth=190.0 + index,
        ),
    ) for index in 1:24
]

SUITE["Voxel"] = BenchmarkGroup()
SUITE["Voxel"]["single ray"] = BenchmarkGroup()
SUITE["Voxel"]["direction construction"] = BenchmarkGroup()
SUITE["Voxel"]["cache and series"] = BenchmarkGroup()
SUITE["Voxel"]["scattering setup"] = BenchmarkGroup()
SUITE["Voxel"]["scattering apply"] = BenchmarkGroup()
SUITE["Voxel"]["scattering warmed"] = BenchmarkGroup()

const VOXEL_RAY_ORIGIN = (8.25, 7.5, 16.0)
const VOXEL_RAY_DIRECTION = (0.37, -0.21, -1.0)
SUITE["Voxel"]["single ray"]["16x16x16 reference"] =
    @benchmarkable ArchimedLight.trace_voxel_ray(
        $VOXEL_MEDIUM,
        $VOXEL_RAY_ORIGIN,
        $VOXEL_RAY_DIRECTION;
        boundary=:periodic,
        algorithm=:reference,
    )
SUITE["Voxel"]["single ray"]["16x16x16 DDA"] =
    @benchmarkable ArchimedLight.trace_voxel_ray(
        $VOXEL_MEDIUM,
        $VOXEL_RAY_ORIGIN,
        $VOXEL_RAY_DIRECTION;
        boundary=:periodic,
        algorithm=:dda,
    )

for (name, grid) in (
    "8x8x8 rays=4" => VOXEL_SMALL,
    "16x16x16 rays=4" => VOXEL_MEDIUM,
    "32x32x32 rays=4" => VOXEL_LARGE,
)
    SUITE["Voxel"]["direction construction"][name] =
        @benchmarkable ArchimedLight.compute_voxel_direction_response(
            $grid,
            (0.37, -0.21, -1.0),
            $VOXEL_BACKEND,
        ) evals = 1
end

for rays in (1, 4, 16)
    ray_backend = ArchimedLight.VoxelCPUBackend(
        rays_per_voxel=rays,
        boundary=:periodic,
        traversal=:dda,
    )
    SUITE["Voxel"]["direction construction"]["8x8x8 rays=$rays"] =
        @benchmarkable ArchimedLight.compute_voxel_direction_response(
            $VOXEL_SMALL,
            (0.37, -0.21, -1.0),
            $ray_backend,
        ) evals = 1
end

SUITE["Voxel"]["cache and series"]["build 16-direction cache 16x16x16"] =
    @benchmarkable ArchimedLight.prepare_voxel_responses(
        $VOXEL_MEDIUM,
        $VOXEL_TURTLE,
        $VOXEL_BACKEND,
    ) evals = 1

SUITE["Voxel"]["cache and series"]["apply cached response 48 steps 16x16x16"] =
    @benchmarkable _voxel_apply_many_steps(
        $VOXEL_MEDIUM,
        $VOXEL_TURTLE,
        $VOXEL_FLUXES,
        $VOXEL_BACKEND,
        $VOXEL_MEDIUM_CACHE,
        $(VOXEL_BENCHMARK_METADATA.repeated_steps),
    ) evals = 1

SUITE["Voxel"]["cache and series"]["full 24-step series 8x8x8"] =
    @benchmarkable ArchimedLight.run_voxel_light_series(
        $VOXEL_SMALL,
        $VOXEL_SERIES_ROWS,
        $VOXEL_OPTIONS;
        backend=$VOXEL_BACKEND,
    ) evals = 1

function _voxel_sparse_benchmark_grid(dimensions::NTuple{3,Int})
    pad = zeros(Float64, dimensions)
    for index in CartesianIndices(pad)
        i, j, k = index.I
        (i + 2j + 3k) % 5 == 0 && (pad[index] = 0.4)
    end
    return ArchimedLight.VoxelGrid(
        (0.0, 0.0, 0.0),
        Tuple(Float64.(dimensions)),
        pad,
    )
end

function _voxel_apply_many_scattering_steps(
    grid,
    turtle,
    first,
    optics,
    options,
    backend,
    cache,
    count::Int,
)
    result = nothing
    for _ in 1:count
        result = ArchimedLight.compute_voxel_scattering(
            grid,
            turtle,
            first,
            optics,
            options,
            backend;
            duration_seconds=1800.0,
            cache=cache,
        )
    end
    return result
end

const VOXEL_SCATTER_SMALL_DENSE = _voxel_benchmark_grid((4, 4, 4))
const VOXEL_SCATTER_MEDIUM_DENSE = _voxel_benchmark_grid((8, 8, 8))
const VOXEL_SCATTER_MEDIUM_SPARSE = _voxel_sparse_benchmark_grid((8, 8, 8))
const VOXEL_SCATTER_OPTIONS_6 = ArchimedLight.LightOptions(
    turtle_sectors=6,
    all_in_turtle=true,
    scattering=true,
    scattering_max_iter=20,
    scattering_stop_ratio=1e-6,
    radiation_timestep_minutes=5.0,
)
const VOXEL_SCATTER_OPTIONS_16 = ArchimedLight.LightOptions(
    VOXEL_SCATTER_OPTIONS_6;
    turtle_sectors=16,
)
const VOXEL_SCATTER_TURTLE_6 = ArchimedLight.build_turtle(
    VOXEL_SCATTER_OPTIONS_6,
    VOXEL_SKY,
)
const VOXEL_SCATTER_TURTLE_16 = VOXEL_TURTLE
const VOXEL_SCATTER_QUADRATURE_6 = ArchimedLight.prepare_voxel_scattering_quadrature(
    VOXEL_SCATTER_TURTLE_6,
)
const VOXEL_SCATTER_QUADRATURE_16 = ArchimedLight.prepare_voxel_scattering_quadrature(
    VOXEL_SCATTER_TURTLE_16,
)
const VOXEL_SCATTER_OPTICS_SMALL = ArchimedLight.voxel_generic_green_leaf_optics(
    VOXEL_SCATTER_SMALL_DENSE,
)
const VOXEL_SCATTER_OPTICS_MEDIUM = ArchimedLight.voxel_generic_green_leaf_optics(
    VOXEL_SCATTER_MEDIUM_DENSE,
)
const VOXEL_SCATTER_FLUXES_6 = ArchimedLight.compute_directional_fluxes(
    VOXEL_ROW,
    VOXEL_SKY,
    VOXEL_SCATTER_TURTLE_6,
    VOXEL_SCATTER_OPTIONS_6,
)
const VOXEL_SCATTER_FIRST_SMALL = ArchimedLight.compute_voxel_first_order(
    VOXEL_SCATTER_SMALL_DENSE,
    VOXEL_SCATTER_TURTLE_6,
    VOXEL_SCATTER_FLUXES_6,
    VOXEL_BACKEND;
    duration_seconds=1800.0,
)
const VOXEL_SCATTER_CACHE_SMALL_6 = ArchimedLight.prepare_voxel_scattering_transport(
    VOXEL_SCATTER_SMALL_DENSE,
    VOXEL_SCATTER_QUADRATURE_6,
    VOXEL_BACKEND,
)
const VOXEL_SCATTER_EMITTED_SMALL =
    0.85 .* VOXEL_SCATTER_FIRST_SMALL.nir_energy

for (name, grid, quadrature) in (
    (
        "dense 4x4x4 directions=12",
        VOXEL_SCATTER_SMALL_DENSE,
        VOXEL_SCATTER_QUADRATURE_6,
    ),
    (
        "dense 8x8x8 directions=12",
        VOXEL_SCATTER_MEDIUM_DENSE,
        VOXEL_SCATTER_QUADRATURE_6,
    ),
    (
        "sparse 8x8x8 directions=12",
        VOXEL_SCATTER_MEDIUM_SPARSE,
        VOXEL_SCATTER_QUADRATURE_6,
    ),
    (
        "dense 4x4x4 directions=32",
        VOXEL_SCATTER_SMALL_DENSE,
        VOXEL_SCATTER_QUADRATURE_16,
    ),
)
    SUITE["Voxel"]["scattering setup"][name] =
        @benchmarkable ArchimedLight.prepare_voxel_scattering_transport(
            $grid,
            $quadrature,
            $VOXEL_BACKEND,
        ) evals = 1
end

SUITE["Voxel"]["scattering apply"]["cached dense 4x4x4 directions=12"] =
    @benchmarkable ArchimedLight.apply_voxel_scattering_transport(
        $VOXEL_SCATTER_SMALL_DENSE,
        $VOXEL_SCATTER_EMITTED_SMALL,
        $VOXEL_SCATTER_QUADRATURE_6,
        $VOXEL_BACKEND;
        band=:nir,
        cache=$VOXEL_SCATTER_CACHE_SMALL_6,
    ) evals = 1

SUITE["Voxel"]["scattering warmed"]["8 repeated solves dense 4x4x4 directions=12"] =
    @benchmarkable _voxel_apply_many_scattering_steps(
        $VOXEL_SCATTER_SMALL_DENSE,
        $VOXEL_SCATTER_TURTLE_6,
        $VOXEL_SCATTER_FIRST_SMALL,
        $VOXEL_SCATTER_OPTICS_SMALL,
        $VOXEL_SCATTER_OPTIONS_6,
        $VOXEL_BACKEND,
        $VOXEL_SCATTER_CACHE_SMALL_6,
        8,
    ) evals = 1

SUITE["Voxel"]["scattering warmed"]["8 rebuilt solves dense 4x4x4 directions=12"] =
    @benchmarkable _voxel_apply_many_scattering_steps(
        $VOXEL_SCATTER_SMALL_DENSE,
        $VOXEL_SCATTER_TURTLE_6,
        $VOXEL_SCATTER_FIRST_SMALL,
        $VOXEL_SCATTER_OPTICS_SMALL,
        $VOXEL_SCATTER_OPTIONS_6,
        $VOXEL_BACKEND,
        nothing,
        8,
    ) evals = 1
