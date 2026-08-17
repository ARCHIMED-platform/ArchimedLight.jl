function _voxel_scattering_quadrature_fingerprint(quadrature::VoxelScatteringQuadrature)
    state = hash(length(quadrature.directions))
    for index in eachindex(quadrature.directions)
        state = hash(quadrature.directions[index], state)
        state = hash(quadrature.weights[index], state)
    end
    return UInt(state)
end

function _validate_voxel_scattering_quadrature(quadrature::VoxelScatteringQuadrature)
    length(quadrature.directions) == length(quadrature.weights) ||
        throw(DimensionMismatch("voxel scattering directions and weights must have equal length"))
    isempty(quadrature.directions) &&
        throw(ArgumentError("voxel scattering quadrature cannot be empty"))
    total = 0.0
    for index in eachindex(quadrature.directions)
        direction = _normalise_voxel_transport_direction(quadrature.directions[index])
        all(isapprox.(direction, quadrature.directions[index]; atol=64eps(Float64), rtol=0)) ||
            throw(ArgumentError("voxel scattering quadrature directions must be normalized"))
        weight = Float64(quadrature.weights[index])
        isfinite(weight) && weight >= 0 || throw(ArgumentError(
            "voxel scattering quadrature weights must be finite and non-negative",
        ))
        total += weight
    end
    isapprox(total, 1.0; atol=64eps(Float64), rtol=0) || throw(ArgumentError(
        "voxel scattering quadrature weights must sum to one; got $total",
    ))
    return quadrature
end

"""
    prepare_voxel_scattering_quadrature(turtle)

Build a normalized full-sphere quadrature from the fixed diffuse-sky sectors
of a turtle. Explicit sun sectors are excluded. Every diffuse direction is
paired with its opposite direction, so isotropic foliage scattering sends
equal total energy into the upper and lower hemispheres.
"""
function prepare_voxel_scattering_quadrature(turtle::TurtleGrid)
    sky_sectors = [sector for sector in turtle.sectors if sector.source == :sky]
    isempty(sky_sectors) &&
        throw(ArgumentError("voxel scattering requires at least one diffuse sky sector"))

    hemisphere_weights = [max(Float64(sector.weight), 0.0) for sector in sky_sectors]
    weight_sum = sum(hemisphere_weights)
    if weight_sum <= 0
        fill!(hemisphere_weights, 1.0)
        weight_sum = length(hemisphere_weights)
    end
    hemisphere_weights ./= weight_sum

    directions = NTuple{3,Float64}[]
    weights = Float64[]
    sizehint!(directions, 2length(sky_sectors))
    sizehint!(weights, 2length(sky_sectors))
    for (sector, weight) in zip(sky_sectors, hemisphere_weights)
        upward = _normalise_voxel_transport_direction(sector.direction)
        upward[3] > 0 ||
            throw(ArgumentError("diffuse turtle directions must point into the upper hemisphere"))
        downward = ntuple(axis -> -upward[axis], 3)
        push!(directions, upward)
        push!(weights, 0.5weight)
        push!(directions, downward)
        push!(weights, 0.5weight)
    end
    isapprox(sum(weights), 1.0; atol=16eps(Float64), rtol=0) ||
        error("voxel scattering quadrature weights do not sum to one")
    return _validate_voxel_scattering_quadrature(
        VoxelScatteringQuadrature(directions, weights),
    )
end

function _voxel_center(grid::VoxelGrid, index::CartesianIndex{3})
    return ntuple(
        axis -> grid.minimum[axis] + (index[axis] - 0.5) * grid.voxel_size[axis],
        3,
    )
end

function _voxel_ground_center(grid::VoxelGrid, i::Int, j::Int)
    return (
        grid.minimum[1] + (i - 0.5) * grid.voxel_size[1],
        grid.minimum[2] + (j - 0.5) * grid.voxel_size[2],
        grid.minimum[3],
    )
end

"""
    prepare_voxel_scattering_transport(grid, quadrature, backend)

Cache internal source-to-boundary DDA paths for occupied voxel centres and
upward ground rays. Attenuation is deliberately not stored in the cache, so
the transport application remains an explicit Beer-Lambert calculation.
"""
function prepare_voxel_scattering_transport(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
)
    _validate_voxel_scattering_quadrature(quadrature)
    backend.boundary in (:periodic, :open) ||
        throw(ArgumentError("voxel scattering supports only :periodic and :open boundaries"))
    backend.traversal in (:dda, :reference) ||
        throw(ArgumentError("voxel scattering supports only :dda and :reference traversal"))

    source_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
    for index in CartesianIndices(grid.pad)
        grid.pad[index] > 0 || continue
        source_id = LinearIndices(grid.pad)[index]
        origin = _voxel_center(grid, index)
        for direction_index in eachindex(quadrature.directions)
            source_paths[(source_id, direction_index)] = _trace_voxel_transport_ray(
                grid,
                origin,
                quadrature.directions[direction_index];
                boundary=backend.boundary,
                algorithm=backend.traversal,
            )
        end
    end

    ground_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
    nx, ny, _ = size(grid)
    for j in 1:ny, i in 1:nx
        ground_id = LinearIndices((nx, ny))[i, j]
        origin = _voxel_ground_center(grid, i, j)
        for direction_index in eachindex(quadrature.directions)
            quadrature.directions[direction_index][3] > 0 || continue
            ground_paths[(ground_id, direction_index)] = _trace_voxel_transport_ray(
                grid,
                origin,
                quadrature.directions[direction_index];
                boundary=backend.boundary,
                algorithm=backend.traversal,
            )
        end
    end

    return VoxelScatteringTransportCache(
        _voxel_grid_fingerprint(grid),
        backend,
        _voxel_scattering_quadrature_fingerprint(quadrature),
        source_paths,
        ground_paths,
    )
end

function _validate_voxel_scattering_cache(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    cache::VoxelScatteringTransportCache,
)
    cache.grid_fingerprint == _voxel_grid_fingerprint(grid) ||
        throw(ArgumentError("voxel scattering cache does not match the grid or PAD values"))
    cache.backend == backend ||
        throw(ArgumentError("voxel scattering cache uses different backend settings"))
    cache.quadrature_fingerprint == _voxel_scattering_quadrature_fingerprint(quadrature) ||
        throw(ArgumentError("voxel scattering cache uses a different angular quadrature"))
    return cache
end

function _accumulate_voxel_path_transport!(
    intercepted,
    ground_incident,
    energy::Float64,
    path::VoxelRayPath,
    grid::VoxelGrid,
    backend::VoxelCPUBackend,
)
    incoming = energy
    last_index = nothing
    for segment in path
        index = CartesianIndex(segment.index)
        density = grid.pad[index]
        if density > 0 && backend.g > 0
            outgoing = incoming * exp(-backend.g * density * segment.length)
            intercepted[index] += incoming - outgoing
            incoming = outgoing
        end
        last_index = segment.index
    end
    if path.exit_reason == :bottom
        last_index === nothing && error("a bottom-exiting internal ray crossed no voxel")
        ground_incident[last_index[1], last_index[2]] += incoming
        return (top=0.0, side=0.0, bottom=0.0)
    elseif path.exit_reason == :top
        return (top=incoming, side=0.0, bottom=0.0)
    elseif path.exit_reason == :side
        return (top=0.0, side=incoming, bottom=0.0)
    end
    error("unsupported voxel exit reason: $(path.exit_reason)")
end

function _ground_band_field(ground::VoxelGroundOptics, band::Symbol)
    band == :par && return ground.par_reflectance
    band == :nir && return ground.nir_reflectance
    throw(ArgumentError("band must be :par or :nir"))
end

function _ground_lambertian_weights(quadrature::VoxelScatteringQuadrature)
    direction_indices = Int[]
    weights = Float64[]
    for index in eachindex(quadrature.directions)
        direction = quadrature.directions[index]
        direction[3] > 0 || continue
        push!(direction_indices, index)
        push!(weights, quadrature.weights[index] * direction[3])
    end
    total = sum(weights)
    total > 0 || error("voxel scattering quadrature has no upward Lambertian support")
    weights ./= total
    return direction_indices, weights
end

"""
    apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                                     band=:par, ground=nothing,
                                     initial_ground_energy=nothing, cache=nothing)

Apply one deterministic, matrix-free scattering transport step. Foliage
sources emit from voxel centres over the full sphere. A supplied ground is
Lambertian: radiation reaching it is split into absorbed and reflected energy,
and reflected energy is propagated upward immediately.
"""
function apply_voxel_scattering_transport(
    grid::VoxelGrid,
    emitted::AbstractArray{<:Real,3},
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    band::Symbol=:par,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    initial_ground_energy::Union{Nothing,AbstractArray{<:Real,2}}=nothing,
    cache::Union{Nothing,VoxelScatteringTransportCache}=nothing,
)
    size(emitted) == size(grid) || throw(DimensionMismatch("emitted energy must match the voxel grid"))
    ground === nothing || ground.grid_size == size(grid)[1:2] ||
        throw(DimensionMismatch("ground optics do not match the voxel grid"))
    initial_ground_energy === nothing || size(initial_ground_energy) == size(grid)[1:2] ||
        throw(DimensionMismatch("initial ground energy must match the horizontal voxel grid"))
    initial_ground_energy === nothing || ground !== nothing || throw(ArgumentError(
        "initial_ground_energy requires ground optical properties",
    ))
    band in (:par, :nir) || throw(ArgumentError("band must be :par or :nir"))

    transport_cache = cache === nothing ?
                      prepare_voxel_scattering_transport(grid, quadrature, backend) :
                      _validate_voxel_scattering_cache(grid, quadrature, backend, cache)
    intercepted = zeros(Float64, size(grid))
    ground_incident = zeros(Float64, size(grid)[1:2])
    if initial_ground_energy !== nothing
        for index in eachindex(initial_ground_energy)
            value = Float64(initial_ground_energy[index])
            isfinite(value) && value >= 0 ||
                throw(ArgumentError("initial ground energy must be finite and non-negative"))
            ground_incident[index] = value
        end
    end

    escaped_top = 0.0
    escaped_side = 0.0
    emitted_total = 0.0
    for index in CartesianIndices(grid.pad)
        energy = Float64(emitted[index])
        isfinite(energy) && energy >= 0 ||
            throw(ArgumentError("emitted voxel energy must be finite and non-negative"))
        energy == 0 && continue
        grid.pad[index] > 0 ||
            throw(ArgumentError("non-zero scattering energy was assigned to empty voxel $index"))
        emitted_total += energy
        source_id = LinearIndices(grid.pad)[index]
        for direction_index in eachindex(quadrature.directions)
            directional_energy = energy * quadrature.weights[direction_index]
            directional_energy == 0 && continue
            path = transport_cache.source_paths[(source_id, direction_index)]
            escaped = _accumulate_voxel_path_transport!(
                intercepted,
                ground_incident,
                directional_energy,
                path,
                grid,
                backend,
            )
            escaped_top += escaped.top
            escaped_side += escaped.side
        end
    end

    initial_ground_total = initial_ground_energy === nothing ? 0.0 : sum(Float64, initial_ground_energy)
    ground_absorbed = 0.0
    ground_reflected = 0.0
    escaped_bottom = 0.0
    if ground === nothing
        escaped_bottom = sum(ground_incident)
    else
        reflectance = _ground_band_field(ground, band)
        reflected_by_cell = zeros(Float64, size(ground_incident))
        for index in CartesianIndices(ground_incident)
            incident = ground_incident[index]
            coefficient = _voxel_optical_value(reflectance, index)
            reflected = coefficient * incident
            reflected_by_cell[index] = reflected
            ground_reflected += reflected
            ground_absorbed += incident - reflected
        end

        direction_indices, ground_weights = _ground_lambertian_weights(quadrature)
        secondary_ground = zeros(Float64, size(ground_incident))
        nx, ny = size(ground_incident)
        for index in CartesianIndices(reflected_by_cell)
            reflected = reflected_by_cell[index]
            reflected == 0 && continue
            ground_id = LinearIndices((nx, ny))[index]
            for (direction_index, weight) in zip(direction_indices, ground_weights)
                path = transport_cache.ground_paths[(ground_id, direction_index)]
                escaped = _accumulate_voxel_path_transport!(
                    intercepted,
                    secondary_ground,
                    reflected * weight,
                    path,
                    grid,
                    backend,
                )
                escaped_top += escaped.top
                escaped_side += escaped.side
            end
        end
        sum(secondary_ground) <= 64eps(Float64) * max(ground_reflected, 1.0) ||
            error("an upward ground-reflected ray unexpectedly returned to the ground")
    end

    expected = emitted_total + initial_ground_total
    accounted = sum(intercepted) + ground_absorbed + escaped_top + escaped_side + escaped_bottom
    residual = expected - accounted
    isapprox(accounted, expected; atol=1e-9 * max(expected, 1.0), rtol=1e-11) ||
        error("voxel scattering transport energy balance does not close: residual=$residual")
    return VoxelTransportResult(
        intercepted,
        ground_absorbed,
        ground_reflected,
        escaped_top,
        escaped_side,
        escaped_bottom,
        residual,
    )
end

function _voxel_band_albedo(optics::VoxelOpticalProperties, band::Symbol, index)
    return clamp(voxel_single_scattering_albedo(optics, band, index), 0.0, 1.0)
end

function _validate_voxel_first_order_for_scattering(
    grid::VoxelGrid,
    first::VoxelFirstOrderResult,
)
    for (band, energy, bottom, incoming, injected, escaped, top, side, lower) in (
        (
            :par,
            first.par_energy,
            first.par_bottom_energy,
            first.par_incoming_energy,
            first.par_injected_energy,
            first.par_escaped_energy,
            first.par_escaped_top_energy,
            first.par_escaped_side_energy,
            first.par_escaped_bottom_energy,
        ),
        (
            :nir,
            first.nir_energy,
            first.nir_bottom_energy,
            first.nir_incoming_energy,
            first.nir_injected_energy,
            first.nir_escaped_energy,
            first.nir_escaped_top_energy,
            first.nir_escaped_side_energy,
            first.nir_escaped_bottom_energy,
        ),
    )
        size(energy) == size(grid) ||
            throw(DimensionMismatch("first-order $band energy does not match the voxel grid"))
        size(bottom) == size(grid)[1:2] || throw(DimensionMismatch(
            "first-order $band bottom energy does not match the horizontal voxel grid",
        ))
        all(isfinite, energy) && all(>=(0), energy) || throw(ArgumentError(
            "first-order $band voxel energy must be finite and non-negative",
        ))
        all(isfinite, bottom) && all(>=(0), bottom) || throw(ArgumentError(
            "first-order $band bottom energy must be finite and non-negative",
        ))
        totals = (incoming, injected, escaped, top, side, lower)
        all(isfinite, totals) && all(>=(0), totals) || throw(ArgumentError(
            "first-order $band energy totals must be finite and non-negative",
        ))
        scale = max(incoming + injected, 1.0)
        isapprox(top + side + lower, escaped; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band boundary escape totals are inconsistent"))
        isapprox(sum(bottom), lower; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band bottom energy field is inconsistent"))
        isapprox(sum(energy) + escaped, incoming + injected; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band energy balance is inconsistent"))
    end
    return first
end

function _disabled_voxel_scattering_band(
    grid::VoxelGrid,
    first_energy,
    first_incoming,
    first_injected,
    first_top,
    first_side,
    first_bottom,
    first_bottom_by_cell,
    optics::VoxelOpticalProperties,
    band::Symbol,
    duration_seconds,
    ground,
)
    absorbed = zeros(Float64, size(grid))
    unresolved = 0.0
    for index in CartesianIndices(grid.pad)
        omega = _voxel_band_albedo(optics, band, index)
        absorbed[index] = (1 - omega) * first_energy[index]
        unresolved += omega * first_energy[index]
    end
    ground_absorbed = 0.0
    ground_reflected = 0.0
    escaped_bottom = first_bottom
    if ground !== nothing
        escaped_bottom = 0.0
        reflectance = _ground_band_field(ground, band)
        for index in CartesianIndices(first_bottom_by_cell)
            incident = first_bottom_by_cell[index]
            reflected = incident * _voxel_optical_value(reflectance, index)
            ground_reflected += reflected
            ground_absorbed += incident - reflected
        end
        unresolved += ground_reflected
    end
    absorbed_flux = similar(absorbed)
    _voxel_flux_from_energy!(absorbed_flux, absorbed, grid, duration_seconds)
    external = first_incoming + first_injected
    accounted = sum(absorbed) + ground_absorbed + first_top + first_side + escaped_bottom + unresolved
    residual = external - accounted
    relative = abs(residual) / max(external, eps(Float64))
    isapprox(accounted, external; atol=1e-8 * max(external, 1.0), rtol=1e-10) ||
        error("disabled voxel $band scattering energy balance does not close: residual=$residual")
    return VoxelBandScatteringResult(
        zeros(Float64, size(grid)),
        Float64.(first_energy),
        absorbed,
        absorbed_flux,
        unresolved - ground_reflected,
        ground_absorbed,
        ground_reflected,
        first_top,
        first_side,
        escaped_bottom,
        unresolved,
        0,
        unresolved <= 64eps(Float64) * max(external, 1.0),
        false,
        abs(residual),
        relative,
        Float64[sum(first_energy)],
    )
end

function _compute_voxel_scattering_band(
    grid::VoxelGrid,
    first_energy,
    first_incoming,
    first_injected,
    first_top,
    first_side,
    first_bottom,
    first_bottom_by_cell,
    optics::VoxelOpticalProperties,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    options::LightOptions,
    band::Symbol,
    duration_seconds;
    ground,
    cache,
    enabled::Bool=true,
)
    enabled || return _disabled_voxel_scattering_band(
        grid,
        first_energy,
        first_incoming,
        first_injected,
        first_top,
        first_side,
        first_bottom,
        first_bottom_by_cell,
        optics,
        band,
        duration_seconds,
        ground,
    )

    current = Float64.(first_energy)
    added = zeros(Float64, size(grid))
    total = copy(current)
    absorbed = zeros(Float64, size(grid))
    emitted = zeros(Float64, size(grid))
    initial_ground = ground === nothing ? nothing : Float64.(first_bottom_by_cell)
    first_order_reference = sum(current)
    source_reference = first_order_reference +
                       (initial_ground === nothing ? 0.0 : sum(initial_ground))
    threshold = options.scattering_stop_ratio * max(first_order_reference, eps(Float64))
    order_totals = Float64[sum(current)]
    scattered_total = 0.0
    ground_absorbed = 0.0
    ground_reflected = 0.0
    escaped_top = Float64(first_top)
    escaped_side = Float64(first_side)
    escaped_bottom = ground === nothing ? Float64(first_bottom) : 0.0
    iterations = 0
    converged = source_reference == 0.0

    for iteration in 1:options.scattering_max_iter
        source_reference == 0.0 && break
        iterations = iteration
        fill!(emitted, 0.0)
        for index in CartesianIndices(grid.pad)
            incoming = current[index]
            omega = _voxel_band_albedo(optics, band, index)
            absorbed[index] += (1 - omega) * incoming
            emitted[index] = omega * incoming
        end
        scattered_total += sum(emitted)
        transport = apply_voxel_scattering_transport(
            grid,
            emitted,
            quadrature,
            backend;
            band=band,
            ground=ground,
            initial_ground_energy=iteration == 1 ? initial_ground : nothing,
            cache=cache,
        )
        ground_absorbed += transport.ground_absorbed_energy
        ground_reflected += transport.ground_reflected_energy
        escaped_top += transport.escaped_top_energy
        escaped_side += transport.escaped_side_energy
        escaped_bottom += transport.escaped_bottom_energy
        current = transport.intercepted_energy
        added .+= current
        total .+= current
        current_total = sum(current)
        push!(order_totals, current_total)
        if current_total < threshold ||
           current_total <= 64eps(Float64) * max(source_reference, 1.0)
            converged = true
            break
        end
    end

    unresolved = sum(current)
    absorbed_flux = similar(absorbed)
    _voxel_flux_from_energy!(absorbed_flux, absorbed, grid, duration_seconds)
    external = Float64(first_incoming + first_injected)
    accounted = sum(absorbed) + ground_absorbed + escaped_top + escaped_side + escaped_bottom + unresolved
    residual = external - accounted
    relative = abs(residual) / max(external, eps(Float64))
    isapprox(accounted, external; atol=1e-8 * max(external, 1.0), rtol=1e-10) ||
        error("voxel $band multiple-scattering energy balance does not close: residual=$residual")
    return VoxelBandScatteringResult(
        added,
        total,
        absorbed,
        absorbed_flux,
        scattered_total,
        ground_absorbed,
        ground_reflected,
        escaped_top,
        escaped_side,
        escaped_bottom,
        unresolved,
        iterations,
        converged,
        true,
        abs(residual),
        relative,
        order_totals,
    )
end

"""
    compute_voxel_scattering(grid, turtle, first, optics, options, backend;
                             duration_seconds, ground=nothing, cache=nothing)

Compute deterministic isotropic multiple scattering independently for PAR and
NIR. The returned absorbed energy is distinct from intercepted energy, and the
unresolved residual is retained explicitly when the stopping criterion or
iteration limit terminates the series.
"""
function compute_voxel_scattering(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    first::VoxelFirstOrderResult,
    optics::VoxelOpticalProperties,
    options::LightOptions,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    duration_seconds::Real,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    cache::Union{Nothing,VoxelScatteringTransportCache}=nothing,
)
    options.scattering || throw(ArgumentError("options.scattering must be true"))
    options.scattering_max_iter > 0 ||
        throw(ArgumentError("scattering_max_iter must be positive"))
    isfinite(options.scattering_stop_ratio) && options.scattering_stop_ratio >= 0 ||
        throw(ArgumentError("scattering_stop_ratio must be finite and non-negative"))
    optics.grid_size == size(grid) || throw(DimensionMismatch("voxel optics do not match the grid"))
    ground === nothing || ground.grid_size == size(grid)[1:2] ||
        throw(DimensionMismatch("ground optics do not match the voxel grid"))
    isfinite(duration_seconds) && duration_seconds >= 0 ||
        throw(ArgumentError("duration_seconds must be finite and non-negative"))
    _validate_voxel_first_order_for_scattering(grid, first)
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    transport_cache = cache === nothing ?
                      prepare_voxel_scattering_transport(grid, quadrature, backend) :
                      _validate_voxel_scattering_cache(grid, quadrature, backend, cache)
    par = _compute_voxel_scattering_band(
        grid,
        first.par_energy,
        first.par_incoming_energy,
        first.par_injected_energy,
        first.par_escaped_top_energy,
        first.par_escaped_side_energy,
        first.par_escaped_bottom_energy,
        first.par_bottom_energy,
        optics,
        quadrature,
        backend,
        options,
        :par,
        Float64(duration_seconds);
        ground=ground,
        cache=transport_cache,
        enabled=true,
    )
    nir = _compute_voxel_scattering_band(
        grid,
        first.nir_energy,
        first.nir_incoming_energy,
        first.nir_injected_energy,
        first.nir_escaped_top_energy,
        first.nir_escaped_side_energy,
        first.nir_escaped_bottom_energy,
        first.nir_bottom_energy,
        optics,
        quadrature,
        backend,
        options,
        :nir,
        Float64(duration_seconds);
        ground=ground,
        cache=transport_cache,
        enabled=options.nir_interception && options.nir_scattering,
    )
    return VoxelScatteringResult(par, nir)
end

function _dense_voxel_scattering_exchange(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
)
    count = length(grid.pad)
    matrix = zeros(Float64, count, count)
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend)
    basis = zeros(Float64, size(grid))
    for source in eachindex(grid.pad)
        grid.pad[source] > 0 || continue
        basis[source] = 1.0
        result = apply_voxel_scattering_transport(
            grid,
            basis,
            quadrature,
            backend;
            cache=cache,
        )
        matrix[:, source] .= vec(result.intercepted_energy)
        basis[source] = 0.0
    end
    return matrix
end
