function _dda_first_boundary_time(origin, component, minimum, spacing, index)
    component == 0 && return Inf
    boundary = component > 0 ? minimum + index * spacing : minimum + (index - 1) * spacing
    time = (boundary - origin) / component
    delta = spacing / abs(component)
    tolerance = _voxel_tolerance(time)
    time <= tolerance && (time += delta)
    return time
end

function _trace_voxel_ray_dda(
    grid::VoxelGrid,
    origin,
    direction,
    boundary::Symbol;
    allow_upward::Bool=false,
)
    boundary in (:periodic, :open, :java_nontoric) ||
        throw(ArgumentError("unsupported voxel boundary mode: $boundary"))
    ray_direction = allow_upward ? _normalise_voxel_transport_direction(direction) :
                    _normalise_voxel_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, boundary)
    t_end, exit_reason = _voxel_exit_time(grid, ray_origin, ray_direction, boundary)
    tolerance = _voxel_tolerance(t_end)
    t_end <= tolerance && return VoxelRayPath(VoxelRaySegment{Float64}[], exit_reason)

    probe_time = min(t_end / 2, 256eps(Float64) * max(maximum(grid.maximum) - minimum(grid.minimum), 1.0))
    probe_time = max(probe_time, min(t_end / 2, tolerance))
    probe = _path_point(grid, ray_origin, ray_direction, probe_time, boundary)
    dimensions = size(grid)
    i = _voxel_axis_index(probe[1], grid.minimum[1], grid.voxel_size[1], dimensions[1])
    j = _voxel_axis_index(probe[2], grid.minimum[2], grid.voxel_size[2], dimensions[2])
    k = _voxel_axis_index(probe[3], grid.minimum[3], grid.voxel_size[3], dimensions[3])

    step_x = Int(sign(ray_direction[1]))
    step_y = Int(sign(ray_direction[2]))
    step_z = Int(sign(ray_direction[3]))
    delta_x = ray_direction[1] == 0 ? Inf : grid.voxel_size[1] / abs(ray_direction[1])
    delta_y = ray_direction[2] == 0 ? Inf : grid.voxel_size[2] / abs(ray_direction[2])
    delta_z = ray_direction[3] == 0 ? Inf : grid.voxel_size[3] / abs(ray_direction[3])
    next_x = _dda_first_boundary_time(
        ray_origin[1], ray_direction[1], grid.minimum[1], grid.voxel_size[1], i
    )
    next_y = _dda_first_boundary_time(
        ray_origin[2], ray_direction[2], grid.minimum[2], grid.voxel_size[2], j
    )
    next_z = _dda_first_boundary_time(
        ray_origin[3], ray_direction[3], grid.minimum[3], grid.voxel_size[3], k
    )

    crossing_rate = sum(
        abs(ray_direction[axis]) / grid.voxel_size[axis] for axis in 1:3
    )
    max_steps = max(8, ceil(Int, t_end * crossing_rate) + 8)
    segments = VoxelRaySegment{Float64}[]
    sizehint!(segments, max_steps)
    current_time = 0.0

    for _ in 1:max_steps
        next_time = min(next_x, next_y, next_z, t_end)
        next_time + tolerance >= current_time ||
            error("voxel DDA produced a decreasing traversal time")
        _append_voxel_segment!(segments, (i, j, k), current_time, next_time)
        if next_time >= t_end - tolerance
            return VoxelRayPath(segments, exit_reason)
        end

        hit_x = abs(next_x - next_time) <= tolerance
        hit_y = abs(next_y - next_time) <= tolerance
        hit_z = abs(next_z - next_time) <= tolerance
        (hit_x || hit_y || hit_z) || error("voxel DDA did not identify a crossed grid plane")
        if hit_x
            i += step_x
            next_x += delta_x
        end
        if hit_y
            j += step_y
            next_y += delta_y
        end
        if hit_z
            k += step_z
            next_z += delta_z
        end

        if k < 1
            return VoxelRayPath(segments, :bottom)
        end
        if k > dimensions[3]
            return VoxelRayPath(segments, :top)
        end

        if !(1 <= i <= dimensions[1])
            if _voxel_is_periodic(boundary)
                i = mod1(i, dimensions[1])
            else
                return VoxelRayPath(segments, :side)
            end
        end
        if !(1 <= j <= dimensions[2])
            if _voxel_is_periodic(boundary)
                j = mod1(j, dimensions[2])
            else
                return VoxelRayPath(segments, :side)
            end
        end
        current_time = next_time
    end
    error("voxel DDA exceeded its derived iteration bound")
end

function _trace_voxel_transport_ray(
    grid::VoxelGrid,
    origin,
    direction;
    boundary::Symbol=:periodic,
    algorithm::Symbol=:dda,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    path = if algorithm == :dda
        _trace_voxel_ray_dda(grid, origin, direction, boundary; allow_upward=true)
    elseif algorithm == :reference
        _trace_voxel_ray_reference(grid, origin, direction, boundary; allow_upward=true)
    else
        algorithm == :java_reference && throw(ArgumentError(
            "the Java reference traversal does not support internal upward rays",
        ))
        throw(ArgumentError("algorithm must be :dda or :reference for internal voxel transport"))
    end
    ray_direction = _normalise_voxel_transport_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, boundary)
    return _truncate_voxel_path_at_terrain(
        grid,
        path,
        ray_origin,
        ray_direction,
        terrain,
        boundary,
    )
end

"""
    trace_voxel_ray(grid, origin, direction;
                    boundary=:periodic, algorithm=:dda,
                    terrain=NoVoxelTerrain())

Trace a normalized propagation ray through `grid`. `direction` is normalized
internally and must point downward. `boundary=:periodic` wraps horizontally;
`:open` terminates at the first horizontal exit. `algorithm=:reference` uses a
correct sorted-plane implementation, while `:dda` uses the production streaming
DDA. `:java_reference` is a compatibility oracle that preserves known Java
wall-selection behavior.
"""
function trace_voxel_ray(
    grid::VoxelGrid,
    origin,
    direction;
    boundary::Symbol=:periodic,
    algorithm::Symbol=:dda,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    if algorithm == :java_reference
        terrain isa NoVoxelTerrain || throw(ArgumentError(
            "terrain is not supported by the Java reference traversal",
        ))
        return _trace_voxel_ray_java_reference(grid, origin, direction, boundary)
    end
    path = if algorithm == :dda
        _trace_voxel_ray_dda(grid, origin, direction, boundary)
    elseif algorithm == :reference
        _trace_voxel_ray_reference(grid, origin, direction, boundary)
    else
        throw(ArgumentError("algorithm must be :dda, :reference, or :java_reference"))
    end
    ray_direction = _normalise_voxel_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, boundary)
    return _truncate_voxel_path_at_terrain(
        grid,
        path,
        ray_origin,
        ray_direction,
        terrain,
        boundary,
    )
end
