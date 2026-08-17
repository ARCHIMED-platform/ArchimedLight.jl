const VoxelScatteringCacheKey = Tuple{UInt,VoxelCPUBackend,UInt}

function _check_voxel_pipeline_options(
    grid::VoxelGrid,
    options::LightOptions,
    backend::VoxelCPUBackend,
    optics::Union{Nothing,VoxelOpticalProperties},
    ground::Union{Nothing,VoxelGroundOptics},
)
    options.scattering || return nothing
    optics === nothing && throw(ArgumentError(
        "voxel scattering requires explicit `optics=VoxelOpticalProperties(...)`; " *
        "the surface-model fallback coefficients are not applied to voxel grids",
    ))
    optics.grid_size == size(grid) || throw(DimensionMismatch("voxel optics do not match the grid"))
    ground === nothing || ground.grid_size == size(grid)[1:2] ||
        throw(DimensionMismatch("ground optics do not match the voxel grid"))
    backend.boundary in (:periodic, :open) || throw(ArgumentError(
        "voxel scattering supports only :periodic and :open boundaries",
    ))
    backend.traversal in (:dda, :reference) || throw(ArgumentError(
        "voxel scattering supports only :dda and :reference traversal",
    ))
    options.scattering_max_iter > 0 ||
        throw(ArgumentError("scattering_max_iter must be positive for voxel scattering"))
    isfinite(options.scattering_stop_ratio) && options.scattering_stop_ratio >= 0 ||
        throw(ArgumentError("scattering_stop_ratio must be finite and non-negative"))
    return nothing
end

function _voxel_scattering_cache_for_turtle!(
    cache_by_configuration::Dict{
        VoxelScatteringCacheKey,
        VoxelScatteringTransportCache{Float64},
    },
    grid::VoxelGrid,
    turtle::TurtleGrid,
    backend::VoxelCPUBackend,
)
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    key = (
        _voxel_grid_fingerprint(grid),
        backend,
        _voxel_scattering_quadrature_fingerprint(quadrature),
    )
    cache = get!(cache_by_configuration, key) do
        prepare_voxel_scattering_transport(grid, quadrature, backend)
    end
    return _validate_voxel_scattering_cache(grid, quadrature, backend, cache)
end

function _voxel_cache_for_turtle!(
    response_by_direction::Dict{NTuple{3,Float64},VoxelDirectionResponse{Float64}},
    grid::VoxelGrid,
    turtle::TurtleGrid,
    backend::VoxelCPUBackend,
)
    directions = NTuple{3,Float64}[]
    responses = VoxelDirectionResponse{Float64}[]
    for sector in turtle.sectors
        direction = _normalise_voxel_direction(ntuple(i -> -Float64(sector.direction[i]), 3))
        response = get!(response_by_direction, direction) do
            compute_voxel_direction_response(grid, direction, backend)
        end
        push!(directions, direction)
        push!(responses, response)
    end
    return VoxelResponseCache(_voxel_grid_fingerprint(grid), backend, directions, responses)
end

function _run_voxel_light_step_resolved(
    grid::VoxelGrid,
    meteo_row,
    options::LightOptions,
    backend::VoxelCPUBackend,
    resolved::ResolvedMeteoStep,
    response_by_direction::Dict{NTuple{3,Float64},VoxelDirectionResponse{Float64}};
    optics::Union{Nothing,VoxelOpticalProperties}=nothing,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    scattering_cache_by_configuration::Union{
        Nothing,
        Dict{VoxelScatteringCacheKey,VoxelScatteringTransportCache{Float64}},
    }=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    _check_voxel_pipeline_options(grid, options, backend, optics, ground)
    sky = compute_sky(
        meteo_row,
        options;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(
        meteo_row,
        sky,
        turtle,
        options;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
    cache = _voxel_cache_for_turtle!(response_by_direction, grid, turtle, backend)
    first_order = compute_voxel_first_order(
        grid,
        turtle,
        fluxes,
        backend;
        duration_seconds=resolved.duration_seconds,
        cache=cache,
    )
    scattering = if options.scattering
        transport_cache = _voxel_scattering_cache_for_turtle!(
            something(scattering_cache_by_configuration),
            grid,
            turtle,
            backend,
        )
        compute_voxel_scattering(
            grid,
            turtle,
            first_order,
            something(optics),
            options,
            backend;
            duration_seconds=resolved.duration_seconds,
            ground=ground,
            cache=transport_cache,
        )
    else
        nothing
    end
    return VoxelLightStepResult(
        sky,
        turtle,
        fluxes,
        first_order,
        scattering,
        Float64(resolved.duration_seconds),
    )
end

"""
    run_voxel_light_step(grid, meteo_row, options; backend=VoxelCPUBackend(),
                         optics=nothing, ground=nothing)

Run voxel interception for one meteorological row. When `options.scattering`
is true, explicit `VoxelOpticalProperties` must be supplied through `optics`.
An optional `VoxelGroundOptics` adds Lambertian ground reflection.
"""
function run_voxel_light_step(
    grid::VoxelGrid,
    meteo_row,
    options::LightOptions;
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    optics::Union{Nothing,VoxelOpticalProperties}=nothing,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    resolved = _resolved_meteo_step_or_error(
        meteo_row,
        options;
        check_boundaries=check_boundaries,
    )
    responses = Dict{NTuple{3,Float64},VoxelDirectionResponse{Float64}}()
    scattering_caches = options.scattering ?
                        Dict{
                            VoxelScatteringCacheKey,
                            VoxelScatteringTransportCache{Float64},
                        }() : nothing
    return _run_voxel_light_step_resolved(
        grid,
        meteo_row,
        options,
        backend,
        resolved,
        responses;
        optics=optics,
        ground=ground,
        scattering_cache_by_configuration=scattering_caches,
        check_boundaries=check_boundaries,
    )
end

function run_voxel_light_series(
    grid::VoxelGrid,
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions;
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    optics::Union{Nothing,VoxelOpticalProperties}=nothing,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    rows, resolved_steps = _prepare_meteo_for_series(
        meteo,
        options;
        check_boundaries=check_boundaries,
        return_resolved=true,
    )
    responses = Dict{NTuple{3,Float64},VoxelDirectionResponse{Float64}}()
    scattering_caches = options.scattering ?
                        Dict{
                            VoxelScatteringCacheKey,
                            VoxelScatteringTransportCache{Float64},
                        }() : nothing
    results = Vector{VoxelLightStepResult{Float64}}(undef, length(rows))
    for (index, row) in enumerate(rows)
        results[index] = _run_voxel_light_step_resolved(
            grid,
            row,
            options,
            backend,
            resolved_steps[index],
            responses;
            optics=optics,
            ground=ground,
            scattering_cache_by_configuration=scattering_caches,
            check_boundaries=check_boundaries,
        )
    end
    return results
end

"""
    run_voxel_light_series(grid, meteo, options; backend=VoxelCPUBackend(),
                           optics=nothing, ground=nothing)

Run a meteorological series while reusing first-order directional responses
and geometry-only internal scattering paths. Explicit voxel optics are required
when scattering is enabled.
"""
function run_voxel_light_series(
    grid::VoxelGrid,
    meteo,
    options::LightOptions;
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    optics::Union{Nothing,VoxelOpticalProperties}=nothing,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    Tables.istable(typeof(meteo)) ||
        throw(ArgumentError("meteo must be a PlantMeteo.TimeStepTable or Tables.jl-compatible table"))
    table = meteo isa PlantMeteo.TimeStepTable ? meteo : _as_plantmeteo_table(meteo)
    return run_voxel_light_series(
        grid,
        table,
        options;
        backend=backend,
        optics=optics,
        ground=ground,
        check_boundaries=check_boundaries,
    )
end
