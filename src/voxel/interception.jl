function _voxel_reference_ray_offsets(grid::VoxelGrid, requested::Integer)
    requested > 0 || throw(ArgumentError("requested ray count must be positive"))
    dx, dy, _ = grid.voxel_size
    spacing = sqrt(dx * dy / requested)
    count_x = max(1, floor(Int, dx / spacing))
    count_y = max(1, floor(Int, dy / spacing))
    step_x = dx / count_x
    step_y = dy / count_y
    offsets = NTuple{2,Float64}[]
    sizehint!(offsets, count_x * count_y)
    for j in 1:count_y, i in 1:count_x
        push!(offsets, ((i - 0.5) * step_x, (j - 0.5) * step_y))
    end
    return offsets
end

function _voxel_grid_fingerprint(grid::VoxelGrid)
    state = hash(size(grid), hash(grid.minimum, hash(grid.maximum)))
    for value in grid.pad
        state = hash(value, state)
    end
    return UInt(state)
end

function _voxel_path_for_backend(
    grid,
    origin,
    direction,
    backend::VoxelCPUBackend,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    return trace_voxel_ray(
        grid,
        origin,
        direction;
        boundary=backend.boundary,
        algorithm=backend.traversal,
        terrain=terrain,
    )
end

function _java_nontoric_wrap(previous, current)
    previous === nothing && return false
    return max(
        abs(previous[1] - current[1]),
        abs(previous[2] - current[2]),
        abs(previous[3] - current[3]),
    ) > 1
end

function _compute_java_direction_response(
    grid::VoxelGrid,
    propagation_direction,
    backend::VoxelCPUBackend,
)
    backend.boundary == :open &&
        throw(ArgumentError("the Java reference backend does not define true open boundaries"))
    direction = _java_float_direction(propagation_direction)
    offsets = _voxel_reference_ray_offsets(grid, backend.rays_per_voxel)
    effective_ray_count = length(offsets)
    ray_weight = Float32(1) / Float32(effective_ray_count)
    interception = zeros(Float32, size(grid))
    escaped_weight = 0.0f0
    escaped_top_weight = 0.0f0
    escaped_side_weight = 0.0f0
    escaped_bottom_weight = 0.0f0
    bottom_exit_fraction = zeros(Float32, size(grid)[1:2])
    injected_weight = 0.0f0
    nx, ny, _ = size(grid)
    top = grid.maximum[3]

    reference_paths = VoxelRayPath{Float64}[]
    for offset in offsets
        # Historical Java positions are local x/y offsets and do not add the
        # voxel-space minimum. Keeping that behavior is intentional here.
        origin = (offset[1], offset[2], top)
        push!(
            reference_paths,
            trace_voxel_ray(
                grid,
                origin,
                direction;
                boundary=:periodic,
                algorithm=:java_reference,
            ),
        )
    end

    for delta_i in (nx - 1):-1:0, delta_j in (ny - 1):-1:0
        previous_index = nothing
        for path in reference_paths
            incoming = ray_weight
            last_index = nothing
            for segment in path
                base_i, base_j, k = segment.index
                index = (mod1(base_i + delta_i, nx), mod1(base_j + delta_j, ny), k)
                hit_wall = _java_nontoric_wrap(previous_index, index)
                if backend.boundary == :java_nontoric && hit_wall
                    injected_weight += ray_weight - incoming
                    incoming = ray_weight
                end
                density = Float32(grid.pad[index...])
                if density != 0 && !isnan(density)
                    exponent = -Float64(Float32(backend.g) * density) * segment.length
                    transmission = Float32(exp(exponent))
                    outgoing = Float32(transmission * incoming)
                    interception[index...] = Float32(
                        interception[index...] + Float32(incoming - outgoing)
                    )
                    incoming = outgoing
                end
                previous_index = index
                last_index = index
            end
            escaped_weight += incoming
            if path.exit_reason == :bottom
                escaped_bottom_weight += incoming
                last_index === nothing || (bottom_exit_fraction[last_index[1], last_index[2]] += incoming)
            elseif path.exit_reason == :side
                escaped_side_weight += incoming
            elseif path.exit_reason == :top
                escaped_top_weight += incoming
            else
                error("unsupported voxel exit reason: $(path.exit_reason)")
            end
        end
    end

    values = Float64.(interception)
    incoming_weight = Float64(nx * ny)
    injected = Float64(injected_weight)
    escaped = Float64(escaped_weight)
    expected = incoming_weight + injected
    rounding_steps = nx * ny * effective_ray_count
    unit_roundoff = eps(Float32) / 2
    rounding_denominator = 1 - rounding_steps * unit_roundoff
    rounding_tolerance = if rounding_denominator > 0
        rounding_steps * unit_roundoff / rounding_denominator * max(expected, 1.0)
    else
        Inf
    end
    isapprox(sum(values) + escaped, expected; atol=rounding_tolerance, rtol=0) ||
        error("Java-compatible voxel directional energy balance does not close")
    return VoxelDirectionResponse(
        values,
        direction,
        incoming_weight,
        injected,
        escaped,
        Float64(escaped_top_weight),
        Float64(escaped_side_weight),
        Float64(escaped_bottom_weight),
        Float64.(bottom_exit_fraction),
        0.0,
        Float64[],
        effective_ray_count,
    )
end

"""
    compute_voxel_direction_response(grid, propagation_direction, backend)

Compute a unit-irradiance geometric response for one downward propagation
direction. Each horizontal voxel receives unit incident weight. The response
can subsequently be scaled independently for PAR and NIR.
"""
function compute_voxel_direction_response(
    grid::VoxelGrid,
    propagation_direction,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    _validate_terrain_for_grid(grid, terrain, backend.boundary)
    if backend.traversal == :java_reference
        terrain isa NoVoxelTerrain || throw(ArgumentError(
            "terrain is not supported by the Java reference traversal",
        ))
        return _compute_java_direction_response(grid, propagation_direction, backend)
    end
    direction = _normalise_voxel_direction(propagation_direction)
    offsets = _voxel_reference_ray_offsets(grid, backend.rays_per_voxel)
    effective_ray_count = length(offsets)
    ray_weight = 1.0 / effective_ray_count
    interception = zeros(Float64, size(grid))
    escaped_weight = 0.0
    escaped_top_weight = 0.0
    escaped_side_weight = 0.0
    escaped_bottom_weight = 0.0
    bottom_exit_fraction = zeros(Float64, size(grid)[1:2])
    terrain_incident_fraction = zeros(Float64, terrain_patch_count(terrain))
    terrain_incident_weight = 0.0
    injected_weight = 0.0
    nx, ny, _ = size(grid)
    dx, dy, _ = grid.voxel_size
    lo = grid.minimum
    top = grid.maximum[3]

    for start_j in 1:ny, start_i in 1:nx
        for offset in offsets
            origin = (
                lo[1] + (start_i - 1) * dx + offset[1],
                lo[2] + (start_j - 1) * dy + offset[2],
                top,
            )
            path = _voxel_path_for_backend(grid, origin, direction, backend, terrain)
            incoming = ray_weight
            previous_index = nothing
            last_index = nothing
            for segment in path
                if backend.boundary == :java_nontoric &&
                   _java_nontoric_wrap(previous_index, segment.index)
                    injected_weight += ray_weight - incoming
                    incoming = ray_weight
                end
                i, j, k = segment.index
                density = grid.pad[i, j, k]
                if density > 0 && backend.g > 0
                    transmission = exp(-backend.g * density * segment.length)
                    outgoing = transmission * incoming
                    interception[i, j, k] += incoming - outgoing
                    incoming = outgoing
                end
                previous_index = segment.index
                last_index = segment.index
            end
            if path.exit_reason == :terrain
                hit = something(path.terrain_hit)
                terrain_incident_fraction[hit.patch_id] += incoming
                terrain_incident_weight += incoming
            elseif path.exit_reason == :bottom
                escaped_weight += incoming
                escaped_bottom_weight += incoming
                last_index === nothing || (bottom_exit_fraction[last_index[1], last_index[2]] += incoming)
            elseif path.exit_reason == :side
                escaped_weight += incoming
                escaped_side_weight += incoming
            elseif path.exit_reason == :top
                escaped_weight += incoming
                escaped_top_weight += incoming
            else
                error("unsupported voxel exit reason: $(path.exit_reason)")
            end
        end
    end

    incoming_weight = Float64(nx * ny)
    intercepted_weight = sum(interception)
    expected = incoming_weight + injected_weight
    isapprox(
        intercepted_weight + escaped_weight + terrain_incident_weight,
        expected;
        atol=1e-11,
        rtol=1e-11,
    ) ||
        error("voxel directional energy balance does not close")
    return VoxelDirectionResponse(
        interception,
        direction,
        incoming_weight,
        injected_weight,
        escaped_weight,
        escaped_top_weight,
        escaped_side_weight,
        escaped_bottom_weight,
        bottom_exit_fraction,
        terrain_incident_weight,
        terrain_incident_fraction,
        effective_ray_count,
    )
end

"""
    prepare_voxel_responses(grid, turtle, backend=VoxelCPUBackend())

Precompute the band-independent voxel response for every upward source
direction in `turtle`. Directions are negated before ray propagation.
"""
function prepare_voxel_responses(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    _validate_terrain_for_grid(grid, terrain, backend.boundary)
    directions = NTuple{3,Float64}[]
    responses = VoxelDirectionResponse{Float64}[]
    sizehint!(directions, length(turtle.sectors))
    sizehint!(responses, length(turtle.sectors))
    for sector in turtle.sectors
        propagation = ntuple(i -> -Float64(sector.direction[i]), 3)
        direction = _normalise_voxel_direction(propagation)
        push!(directions, direction)
        push!(responses, compute_voxel_direction_response(grid, direction, backend, terrain))
    end
    return VoxelResponseCache(
        _voxel_grid_fingerprint(grid),
        _terrain_geometry_fingerprint(terrain),
        backend,
        directions,
        responses,
    )
end

function _validate_voxel_cache(
    grid,
    turtle,
    backend,
    terrain::AbstractVoxelTerrain,
    cache::VoxelResponseCache,
)
    cache.grid_fingerprint == _voxel_grid_fingerprint(grid) ||
        throw(ArgumentError("voxel response cache does not match the grid or PAD values"))
    cache.backend == backend || throw(ArgumentError("voxel response cache uses different backend settings"))
    cache.terrain_fingerprint == _terrain_geometry_fingerprint(terrain) || throw(ArgumentError(
        "voxel response cache does not match the terrain geometry or material assignment",
    ))
    length(cache.responses) == length(turtle.sectors) ||
        throw(ArgumentError("voxel response cache does not match turtle sector count"))
    for (index, sector) in enumerate(turtle.sectors)
        expected = _normalise_voxel_direction(ntuple(i -> -Float64(sector.direction[i]), 3))
        all(isapprox.(cache.directions[index], expected; atol=1e-13, rtol=1e-13)) ||
            throw(ArgumentError("voxel response cache direction mismatch at sector $index"))
    end
    return cache
end

function _voxel_flux_from_energy!(flux, energy, grid::VoxelGrid, duration_seconds)
    for index in CartesianIndices(grid.pad)
        plant_area = grid.pad[index] * voxel_volume(grid)
        flux[index] = plant_area == 0 || duration_seconds == 0 ? 0.0 :
                      energy[index] / (duration_seconds * plant_area)
    end
    return flux
end

"""
    compute_voxel_first_order(grid, turtle, fluxes, backend=VoxelCPUBackend();
                              duration_seconds, cache=nothing)

Apply directional PAR/NIR irradiances to precomputed unit responses. Returned
energy arrays are in joules per voxel and per step; flux arrays are in
`W m⁻² plant`.
"""
function compute_voxel_first_order(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    duration_seconds::Real,
    cache::Union{Nothing,VoxelResponseCache}=nothing,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    isfinite(duration_seconds) && duration_seconds >= 0 ||
        throw(ArgumentError("duration_seconds must be finite and non-negative"))
    sector_count = length(turtle.sectors)
    length(fluxes.sector_ids) == sector_count || throw(DimensionMismatch("flux sector count mismatch"))
    length(fluxes.par) == sector_count || throw(DimensionMismatch("PAR flux count mismatch"))
    length(fluxes.nir) == sector_count || throw(DimensionMismatch("NIR flux count mismatch"))
    [sector.id for sector in turtle.sectors] == fluxes.sector_ids ||
        throw(ArgumentError("flux sector ids do not match turtle sectors"))
    all(isfinite, fluxes.par) && all(>=(0), fluxes.par) ||
        throw(ArgumentError("PAR directional fluxes must be finite and non-negative"))
    all(isfinite, fluxes.nir) && all(>=(0), fluxes.nir) ||
        throw(ArgumentError("NIR directional fluxes must be finite and non-negative"))

    _validate_terrain_for_grid(grid, terrain, backend.boundary)
    responses = cache === nothing ? prepare_voxel_responses(grid, turtle, backend, terrain) :
                _validate_voxel_cache(grid, turtle, backend, terrain, cache)
    par_energy = zeros(Float64, size(grid))
    nir_energy = zeros(Float64, size(grid))
    area_time = voxel_horizontal_area(grid) * Float64(duration_seconds)
    par_incoming = 0.0
    nir_incoming = 0.0
    par_injected = 0.0
    nir_injected = 0.0
    par_escaped = 0.0
    nir_escaped = 0.0
    par_escaped_top = 0.0
    nir_escaped_top = 0.0
    par_escaped_side = 0.0
    nir_escaped_side = 0.0
    par_escaped_bottom = 0.0
    nir_escaped_bottom = 0.0
    par_bottom = zeros(Float64, size(grid)[1:2])
    nir_bottom = zeros(Float64, size(grid)[1:2])
    par_terrain = zeros(Float64, terrain_patch_count(terrain))
    nir_terrain = zeros(Float64, terrain_patch_count(terrain))
    par_terrain_incident = 0.0
    nir_terrain_incident = 0.0

    for sector_index in 1:sector_count
        response = responses.responses[sector_index]
        par_scale = fluxes.par[sector_index] * area_time
        nir_scale = fluxes.nir[sector_index] * area_time
        for index in eachindex(response.interception_fraction)
            fraction = response.interception_fraction[index]
            par_energy[index] += fraction * par_scale
            nir_energy[index] += fraction * nir_scale
        end
        par_incoming += response.incoming_weight * par_scale
        nir_incoming += response.incoming_weight * nir_scale
        par_injected += response.injected_weight * par_scale
        nir_injected += response.injected_weight * nir_scale
        par_escaped += response.escaped_weight * par_scale
        nir_escaped += response.escaped_weight * nir_scale
        par_escaped_top += response.escaped_top_weight * par_scale
        nir_escaped_top += response.escaped_top_weight * nir_scale
        par_escaped_side += response.escaped_side_weight * par_scale
        nir_escaped_side += response.escaped_side_weight * nir_scale
        par_escaped_bottom += response.escaped_bottom_weight * par_scale
        nir_escaped_bottom += response.escaped_bottom_weight * nir_scale
        par_terrain_incident += response.terrain_incident_weight * par_scale
        nir_terrain_incident += response.terrain_incident_weight * nir_scale
        for index in eachindex(response.bottom_exit_fraction)
            fraction = response.bottom_exit_fraction[index]
            par_bottom[index] += fraction * par_scale
            nir_bottom[index] += fraction * nir_scale
        end
        for index in eachindex(response.terrain_incident_fraction)
            fraction = response.terrain_incident_fraction[index]
            par_terrain[index] += fraction * par_scale
            nir_terrain[index] += fraction * nir_scale
        end
    end

    par_flux = similar(par_energy)
    nir_flux = similar(nir_energy)
    _voxel_flux_from_energy!(par_flux, par_energy, grid, Float64(duration_seconds))
    _voxel_flux_from_energy!(nir_flux, nir_energy, grid, Float64(duration_seconds))
    isapprox(
        sum(par_energy) + par_escaped + par_terrain_incident,
        par_incoming + par_injected;
        atol=1e-8,
        rtol=1e-10,
    ) ||
        error("voxel PAR energy balance does not close")
    isapprox(
        sum(nir_energy) + nir_escaped + nir_terrain_incident,
        nir_incoming + nir_injected;
        atol=1e-8,
        rtol=1e-10,
    ) ||
        error("voxel NIR energy balance does not close")

    return VoxelFirstOrderResult(
        par_energy,
        nir_energy,
        par_flux,
        nir_flux,
        par_incoming,
        nir_incoming,
        par_injected,
        nir_injected,
        par_escaped,
        nir_escaped,
        par_escaped_top,
        nir_escaped_top,
        par_escaped_side,
        nir_escaped_side,
        par_escaped_bottom,
        nir_escaped_bottom,
        par_bottom,
        nir_bottom,
        par_terrain_incident,
        nir_terrain_incident,
        par_terrain,
        nir_terrain,
    )
end
