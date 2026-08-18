const _SoilMaterialTable = Dict{Int,SoilOpticalProperties{Float64}}

function _normalise_terrain_vector(value, name::AbstractString)
    length(value) == 3 || throw(ArgumentError("$name must contain three coordinates"))
    vector = ntuple(axis -> Float64(value[axis]), 3)
    all(isfinite, vector) || throw(ArgumentError("$name coordinates must be finite"))
    magnitude = sqrt(sum(abs2, vector))
    magnitude > 0 || throw(ArgumentError("$name cannot be zero"))
    return ntuple(axis -> vector[axis] / magnitude, 3)
end

function _soil_material_table(materials)
    isempty(materials) && throw(ArgumentError("a terrain must define at least one soil material"))
    result = _SoilMaterialTable()
    for (identifier, optics) in pairs(materials)
        id = Int(identifier)
        id > 0 || throw(ArgumentError("soil material identifiers must be positive"))
        optics isa SoilOpticalProperties ||
            throw(ArgumentError("soil material $id is not a SoilOpticalProperties value"))
        result[id] = SoilOpticalProperties(
            par_reflectance=optics.par_reflectance,
            nir_reflectance=optics.nir_reflectance,
        )
    end
    return result
end

function _validate_material_ids(material_ids, materials, expected::Integer)
    length(material_ids) == expected || throw(DimensionMismatch(
        "terrain has $(length(material_ids)) material assignments; expected $expected",
    ))
    ids = Int.(vec(material_ids))
    for id in ids
        haskey(materials, id) || throw(ArgumentError("terrain references undefined soil material $id"))
    end
    return ids
end

"""Horizontal finite plane split into deterministic rectangular soil patches."""
struct PlanarTerrain{T<:AbstractFloat} <: AbstractVoxelTerrain
    elevation::T
    x_bounds::NTuple{2,T}
    y_bounds::NTuple{2,T}
    cells::NTuple{2,Int}
    material_ids::Matrix{Int}
    materials::_SoilMaterialTable
    periodic::Bool
end

function PlanarTerrain(
    elevation::Real,
    x_bounds,
    y_bounds,
    optics::SoilOpticalProperties;
    cells=(1, 1),
    material_id::Integer=1,
    material_ids=nothing,
    materials=Dict(Int(material_id) => optics),
    periodic::Bool=false,
)
    z = Float64(elevation)
    isfinite(z) || throw(ArgumentError("planar terrain elevation must be finite"))
    length(x_bounds) == 2 || throw(ArgumentError("x_bounds must contain two values"))
    length(y_bounds) == 2 || throw(ArgumentError("y_bounds must contain two values"))
    xb = (Float64(x_bounds[1]), Float64(x_bounds[2]))
    yb = (Float64(y_bounds[1]), Float64(y_bounds[2]))
    all(isfinite, xb) && xb[2] > xb[1] ||
        throw(ArgumentError("planar terrain x bounds must be finite and increasing"))
    all(isfinite, yb) && yb[2] > yb[1] ||
        throw(ArgumentError("planar terrain y bounds must be finite and increasing"))
    length(cells) == 2 || throw(ArgumentError("cells must contain two dimensions"))
    cell_count = (Int(cells[1]), Int(cells[2]))
    all(>(0), cell_count) || throw(ArgumentError("planar terrain cells must be positive"))
    table = _soil_material_table(materials)
    assignments = material_ids === nothing ? fill(Int(material_id), cell_count) : Int.(material_ids)
    size(assignments) == cell_count || throw(DimensionMismatch(
        "planar material_ids has size $(size(assignments)); expected $cell_count",
    ))
    _validate_material_ids(assignments, table, prod(cell_count))
    return PlanarTerrain(z, xb, yb, cell_count, Matrix(assignments), table, periodic)
end

function PlanarTerrain(
    grid::VoxelGrid,
    optics::SoilOpticalProperties;
    elevation::Real=grid.minimum[3],
    material_id::Integer=1,
    material_ids=nothing,
    materials=Dict(Int(material_id) => optics),
    periodic::Bool=false,
)
    return PlanarTerrain(
        elevation,
        (grid.minimum[1], grid.maximum[1]),
        (grid.minimum[2], grid.maximum[2]),
        optics;
        cells=size(grid)[1:2],
        material_id=material_id,
        material_ids=material_ids,
        materials=materials,
        periodic=periodic,
    )
end

# Separate method avoids clever tuple slicing in the public grid constructor.
function PlanarTerrain(
    grid::VoxelGrid;
    optics::SoilOpticalProperties,
    elevation::Real=grid.minimum[3],
    material_id::Integer=1,
    material_ids=nothing,
    materials=Dict(Int(material_id) => optics),
    periodic::Bool=false,
)
    return PlanarTerrain(
        elevation,
        (grid.minimum[1], grid.maximum[1]),
        (grid.minimum[2], grid.maximum[2]),
        optics;
        cells=size(grid)[1:2],
        material_id=material_id,
        material_ids=material_ids,
        materials=materials,
        periodic=periodic,
    )
end

"""Compact bounding-volume hierarchy used by triangulated terrain."""
struct TerrainBVH{T<:AbstractFloat}
    lower::Vector{NTuple{3,T}}
    upper::Vector{NTuple{3,T}}
    left::Vector{Int}
    right::Vector{Int}
    leaves::Vector{Vector{Int}}
end

"""General non-periodic triangle terrain with a cached BVH."""
struct TriangulatedTerrain{T<:AbstractFloat} <: AbstractVoxelTerrain
    vertices::Vector{NTuple{3,T}}
    triangles::Vector{NTuple{3,Int}}
    material_ids::Vector{Int}
    materials::_SoilMaterialTable
    facet_normals::Vector{NTuple{3,T}}
    vertex_normals::Vector{NTuple{3,T}}
    normal_mode::Symbol
    bvh::TerrainBVH{T}
end

function _triangle_geometry(vertices, triangle)
    a, b, c = (vertices[triangle[index]] for index in 1:3)
    edge1 = ntuple(axis -> b[axis] - a[axis], 3)
    edge2 = ntuple(axis -> c[axis] - a[axis], 3)
    raw = (
        edge1[2] * edge2[3] - edge1[3] * edge2[2],
        edge1[3] * edge2[1] - edge1[1] * edge2[3],
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
    )
    magnitude = sqrt(sum(abs2, raw))
    magnitude > 0 || throw(ArgumentError("terrain triangle has zero area"))
    normal = ntuple(axis -> raw[axis] / magnitude, 3)
    normal[3] >= 0 || (normal = ntuple(axis -> -normal[axis], 3))
    normal[3] > 128eps(Float64) || throw(ArgumentError(
        "terrain triangles must have an upward-facing component; vertical/overhanging facets are unsupported",
    ))
    return normal, magnitude
end

function _terrain_vertex_normals(vertices, triangles, facet_normals)
    accumulators = [(0.0, 0.0, 0.0) for _ in eachindex(vertices)]
    for (triangle, normal) in zip(triangles, facet_normals)
        _, area2 = _triangle_geometry(vertices, triangle)
        for vertex_id in triangle
            current = accumulators[vertex_id]
            accumulators[vertex_id] = ntuple(axis -> current[axis] + area2 * normal[axis], 3)
        end
    end
    result = NTuple{3,Float64}[]
    for normal in accumulators
        push!(result, _normalise_terrain_vector(normal, "terrain vertex normal"))
    end
    return result
end

function _triangle_bounds(vertices, triangle)
    points = (vertices[triangle[1]], vertices[triangle[2]], vertices[triangle[3]])
    lower = ntuple(axis -> min(points[1][axis], points[2][axis], points[3][axis]), 3)
    upper = ntuple(axis -> max(points[1][axis], points[2][axis], points[3][axis]), 3)
    centre = ntuple(axis -> (points[1][axis] + points[2][axis] + points[3][axis]) / 3, 3)
    return lower, upper, centre
end

function _build_terrain_bvh(vertices, triangles; leaf_size::Int=8)
    leaf_size > 0 || throw(ArgumentError("terrain BVH leaf_size must be positive"))
    triangle_bounds = [_triangle_bounds(vertices, triangle) for triangle in triangles]
    lowers = NTuple{3,Float64}[]
    uppers = NTuple{3,Float64}[]
    left = Int[]
    right = Int[]
    leaves = Vector{Int}[]

    function build!(ids::Vector{Int})
        node = length(lowers) + 1
        lower = ntuple(axis -> minimum(triangle_bounds[id][1][axis] for id in ids), 3)
        upper = ntuple(axis -> maximum(triangle_bounds[id][2][axis] for id in ids), 3)
        push!(lowers, lower)
        push!(uppers, upper)
        push!(left, 0)
        push!(right, 0)
        push!(leaves, Int[])
        if length(ids) <= leaf_size
            leaves[node] = sort(ids)
            return node
        end
        spans = ntuple(axis -> upper[axis] - lower[axis], 3)
        axis = argmax(spans)
        sort!(ids; by=id -> (triangle_bounds[id][3][axis], id))
        middle = length(ids) ÷ 2
        left[node] = build!(copy(ids[1:middle]))
        right[node] = build!(copy(ids[(middle + 1):end]))
        return node
    end

    build!(collect(eachindex(triangles)))
    return TerrainBVH(lowers, uppers, left, right, leaves)
end

function TriangulatedTerrain(
    vertices,
    triangles,
    optics::SoilOpticalProperties;
    material_id::Integer=1,
    material_ids=nothing,
    materials=Dict(Int(material_id) => optics),
    normal_mode::Symbol=:facet,
    periodic::Bool=false,
    leaf_size::Integer=8,
)
    periodic && throw(ArgumentError(
        "periodic generic triangle meshes are not supported; use HeightFieldTerrain",
    ))
    normal_mode in (:facet, :interpolated) || throw(ArgumentError(
        "terrain normal_mode must be :facet or :interpolated",
    ))
    points = NTuple{3,Float64}[]
    for (index, vertex) in enumerate(vertices)
        length(vertex) == 3 || throw(ArgumentError("terrain vertex $index must contain three coordinates"))
        point = ntuple(axis -> Float64(vertex[axis]), 3)
        all(isfinite, point) || throw(ArgumentError("terrain vertex $index is not finite"))
        push!(points, point)
    end
    length(points) >= 3 || throw(ArgumentError("a terrain mesh needs at least three vertices"))
    faces = NTuple{3,Int}[]
    normals = NTuple{3,Float64}[]
    for (index, triangle) in enumerate(triangles)
        length(triangle) == 3 || throw(ArgumentError("terrain triangle $index must have three vertices"))
        face = ntuple(axis -> Int(triangle[axis]), 3)
        all(vertex -> 1 <= vertex <= length(points), face) ||
            throw(ArgumentError("terrain triangle $index references an invalid vertex"))
        length(unique(face)) == 3 || throw(ArgumentError("terrain triangle $index repeats a vertex"))
        normal, _ = _triangle_geometry(points, face)
        push!(faces, face)
        push!(normals, normal)
    end
    isempty(faces) && throw(ArgumentError("a terrain mesh needs at least one triangle"))
    table = _soil_material_table(materials)
    assignments = material_ids === nothing ? fill(Int(material_id), length(faces)) : material_ids
    ids = _validate_material_ids(assignments, table, length(faces))
    vertex_normals = _terrain_vertex_normals(points, faces, normals)
    bvh = _build_terrain_bvh(points, faces; leaf_size=Int(leaf_size))
    return TriangulatedTerrain(
        points,
        faces,
        ids,
        table,
        normals,
        vertex_normals,
        normal_mode,
        bvh,
    )
end

"""
    HeightFieldTerrain(x, y, elevation, optics; diagonal=:southwest_northeast)

Triangulated height samples at increasing vertex coordinates `x` and `y`.
Every cell uses the fixed southwest-to-northeast diagonal. A BVH is built once
and reused. `periodic=true` tiles the height field horizontally.
"""
struct HeightFieldTerrain{T<:AbstractFloat} <: AbstractVoxelTerrain
    x::Vector{T}
    y::Vector{T}
    elevation::Matrix{T}
    diagonal::Symbol
    periodic::Bool
    mesh::TriangulatedTerrain{T}
end

function HeightFieldTerrain(
    x,
    y,
    elevation::AbstractMatrix,
    optics::SoilOpticalProperties;
    diagonal::Symbol=:southwest_northeast,
    periodic::Bool=false,
    material_id::Integer=1,
    material_ids=nothing,
    materials=Dict(Int(material_id) => optics),
    normal_mode::Symbol=:facet,
)
    diagonal == :southwest_northeast || throw(ArgumentError(
        "only diagonal=:southwest_northeast is supported",
    ))
    xv = Float64.(collect(x))
    yv = Float64.(collect(y))
    length(xv) >= 2 && all(diff(xv) .> 0) && all(isfinite, xv) ||
        throw(ArgumentError("height-field x coordinates must be finite and strictly increasing"))
    length(yv) >= 2 && all(diff(yv) .> 0) && all(isfinite, yv) ||
        throw(ArgumentError("height-field y coordinates must be finite and strictly increasing"))
    size(elevation) == (length(xv), length(yv)) || throw(DimensionMismatch(
        "height-field elevation has size $(size(elevation)); expected $((length(xv), length(yv)))",
    ))
    z = Float64.(elevation)
    all(isfinite, z) || throw(ArgumentError("height-field elevations must be finite (no NODATA)"))

    vertices = NTuple{3,Float64}[]
    for j in eachindex(yv), i in eachindex(xv)
        push!(vertices, (xv[i], yv[j], z[i, j]))
    end
    vertex_id(i, j) = (j - 1) * length(xv) + i
    triangles = NTuple{3,Int}[]
    for j in 1:(length(yv) - 1), i in 1:(length(xv) - 1)
        southwest = vertex_id(i, j)
        southeast = vertex_id(i + 1, j)
        northeast = vertex_id(i + 1, j + 1)
        northwest = vertex_id(i, j + 1)
        push!(triangles, (southwest, southeast, northeast))
        push!(triangles, (southwest, northeast, northwest))
    end
    assignments = if material_ids === nothing
        fill(Int(material_id), length(triangles))
    elseif material_ids isa AbstractMatrix &&
           size(material_ids) == (length(xv) - 1, length(yv) - 1)
        ids = Int[]
        for j in 1:(length(yv) - 1), i in 1:(length(xv) - 1)
            push!(ids, Int(material_ids[i, j]), Int(material_ids[i, j]))
        end
        ids
    else
        material_ids
    end
    mesh = TriangulatedTerrain(
        vertices,
        triangles,
        optics;
        material_id=material_id,
        material_ids=assignments,
        materials=materials,
        normal_mode=normal_mode,
        periodic=false,
    )
    return HeightFieldTerrain(xv, yv, z, diagonal, periodic, mesh)
end

terrain_bounds(::NoVoxelTerrain) = nothing
terrain_bounds(terrain::PlanarTerrain) = (
    (terrain.x_bounds[1], terrain.y_bounds[1], terrain.elevation),
    (terrain.x_bounds[2], terrain.y_bounds[2], terrain.elevation),
)

function terrain_bounds(terrain::TriangulatedTerrain)
    lower = ntuple(axis -> minimum(vertex[axis] for vertex in terrain.vertices), 3)
    upper = ntuple(axis -> maximum(vertex[axis] for vertex in terrain.vertices), 3)
    return lower, upper
end

terrain_bounds(terrain::HeightFieldTerrain) = terrain_bounds(terrain.mesh)

terrain_patch_count(::NoVoxelTerrain) = 0
terrain_patch_count(terrain::PlanarTerrain) = prod(terrain.cells)
terrain_patch_count(terrain::TriangulatedTerrain) = length(terrain.triangles)
terrain_patch_count(terrain::HeightFieldTerrain) = terrain_patch_count(terrain.mesh)

function _check_terrain_patch(terrain, patch_id::Integer)
    1 <= patch_id <= terrain_patch_count(terrain) || throw(BoundsError(terrain, patch_id))
    return Int(patch_id)
end

function terrain_patch_point(terrain::PlanarTerrain, patch_id::Integer)
    patch = _check_terrain_patch(terrain, patch_id)
    nx, ny = terrain.cells
    i = mod1(patch, nx)
    j = (patch - 1) ÷ nx + 1
    dx = (terrain.x_bounds[2] - terrain.x_bounds[1]) / nx
    dy = (terrain.y_bounds[2] - terrain.y_bounds[1]) / ny
    return (
        terrain.x_bounds[1] + (i - 0.5) * dx,
        terrain.y_bounds[1] + (j - 0.5) * dy,
        terrain.elevation,
    )
end

terrain_patch_normal(terrain::PlanarTerrain, patch_id::Integer) =
    (_check_terrain_patch(terrain, patch_id); (0.0, 0.0, 1.0))

function terrain_patch_material(terrain::PlanarTerrain, patch_id::Integer)
    patch = _check_terrain_patch(terrain, patch_id)
    i = mod1(patch, terrain.cells[1])
    j = (patch - 1) ÷ terrain.cells[1] + 1
    return terrain.material_ids[i, j]
end

function terrain_patch_point(terrain::TriangulatedTerrain, patch_id::Integer)
    patch = _check_terrain_patch(terrain, patch_id)
    triangle = terrain.triangles[patch]
    return ntuple(
        axis -> sum(terrain.vertices[vertex][axis] for vertex in triangle) / 3,
        3,
    )
end

terrain_patch_normal(terrain::TriangulatedTerrain, patch_id::Integer) =
    terrain.facet_normals[_check_terrain_patch(terrain, patch_id)]
terrain_patch_material(terrain::TriangulatedTerrain, patch_id::Integer) =
    terrain.material_ids[_check_terrain_patch(terrain, patch_id)]
terrain_patch_point(terrain::HeightFieldTerrain, patch_id::Integer) =
    terrain_patch_point(terrain.mesh, patch_id)
terrain_patch_normal(terrain::HeightFieldTerrain, patch_id::Integer) =
    terrain_patch_normal(terrain.mesh, patch_id)
terrain_patch_material(terrain::HeightFieldTerrain, patch_id::Integer) =
    terrain_patch_material(terrain.mesh, patch_id)

_terrain_materials(terrain::PlanarTerrain) = terrain.materials
_terrain_materials(terrain::TriangulatedTerrain) = terrain.materials
_terrain_materials(terrain::HeightFieldTerrain) = terrain.mesh.materials

_no_terrain_soil_optics_error() = throw(ArgumentError(
    "NoVoxelTerrain has no soil optical properties",
))

soil_optics(::NoVoxelTerrain, ::Integer) = _no_terrain_soil_optics_error()
soil_optics(::NoVoxelTerrain, ::TerrainHit) = _no_terrain_soil_optics_error()
soil_optics(::NoVoxelTerrain, ::TerrainHit, ::Symbol) = _no_terrain_soil_optics_error()

function soil_optics(terrain::AbstractVoxelTerrain, material_id::Integer)
    materials = _terrain_materials(terrain)
    haskey(materials, material_id) || throw(ArgumentError("undefined soil material $material_id"))
    return materials[Int(material_id)]
end

soil_optics(terrain::AbstractVoxelTerrain, hit::TerrainHit) =
    soil_optics(terrain, hit.material_id)

function _soil_band_reflectance(optics::SoilOpticalProperties, band::Symbol)
    band == :par && return optics.par_reflectance
    band == :nir && return optics.nir_reflectance
    throw(ArgumentError("band must be :par or :nir"))
end

soil_optics(terrain::AbstractVoxelTerrain, hit::TerrainHit, band::Symbol) =
    _soil_band_reflectance(soil_optics(terrain, hit), band)

function _soil_patch_reflectance(
    terrain::AbstractVoxelTerrain,
    patch_id::Integer,
    band::Symbol,
)
    material_id = terrain_patch_material(terrain, patch_id)
    return _soil_band_reflectance(soil_optics(terrain, material_id), band)
end

function _ray_aabb_intersects(origin, direction, lower, upper, maximum_distance)
    tmin = 0.0
    tmax = maximum_distance
    for axis in 1:3
        if abs(direction[axis]) <= eps(Float64)
            lower[axis] <= origin[axis] <= upper[axis] || return false
        else
            first = (lower[axis] - origin[axis]) / direction[axis]
            second = (upper[axis] - origin[axis]) / direction[axis]
            first > second && ((first, second) = (second, first))
            tmin = max(tmin, first)
            tmax = min(tmax, second)
            tmax + _voxel_tolerance(tmax) >= tmin || return false
        end
    end
    return tmax >= 0
end

function _ray_triangle_intersection(origin, direction, a, b, c, maximum_distance)
    edge1 = ntuple(axis -> b[axis] - a[axis], 3)
    edge2 = ntuple(axis -> c[axis] - a[axis], 3)
    pvec = (
        direction[2] * edge2[3] - direction[3] * edge2[2],
        direction[3] * edge2[1] - direction[1] * edge2[3],
        direction[1] * edge2[2] - direction[2] * edge2[1],
    )
    determinant = sum(edge1[axis] * pvec[axis] for axis in 1:3)
    scale = max(sqrt(sum(abs2, edge1)) * sqrt(sum(abs2, edge2)), 1.0)
    abs(determinant) > 128eps(Float64) * scale || return nothing
    inverse = inv(determinant)
    tvec = ntuple(axis -> origin[axis] - a[axis], 3)
    u = sum(tvec[axis] * pvec[axis] for axis in 1:3) * inverse
    tolerance = 256eps(Float64)
    -tolerance <= u <= 1 + tolerance || return nothing
    qvec = (
        tvec[2] * edge1[3] - tvec[3] * edge1[2],
        tvec[3] * edge1[1] - tvec[1] * edge1[3],
        tvec[1] * edge1[2] - tvec[2] * edge1[1],
    )
    v = sum(direction[axis] * qvec[axis] for axis in 1:3) * inverse
    -tolerance <= v && u + v <= 1 + tolerance || return nothing
    distance = sum(edge2[axis] * qvec[axis] for axis in 1:3) * inverse
    distance_scale = isfinite(maximum_distance) ? maximum_distance :
                     max(abs(distance), maximum(abs, origin), 1.0)
    ray_tolerance = 256eps(Float64) * max(distance_scale, 1.0)
    distance > ray_tolerance || return nothing
    isfinite(maximum_distance) && distance > maximum_distance + ray_tolerance && return nothing
    accepted_distance = isfinite(maximum_distance) ? min(distance, maximum_distance) : distance
    return (
        distance=accepted_distance,
        u=clamp(u, 0.0, 1.0),
        v=clamp(v, 0.0, 1.0),
    )
end

function _interpolated_triangle_normal(terrain, triangle, u, v)
    w = 1 - u - v
    values = ntuple(
        axis -> w * terrain.vertex_normals[triangle[1]][axis] +
                u * terrain.vertex_normals[triangle[2]][axis] +
                v * terrain.vertex_normals[triangle[3]][axis],
        3,
    )
    normal = _normalise_terrain_vector(values, "interpolated terrain normal")
    return normal[3] >= 0 ? normal : ntuple(axis -> -normal[axis], 3)
end

function _intersect_triangulated_base(terrain::TriangulatedTerrain, origin, direction, maximum_distance)
    best_distance = maximum_distance
    best_patch = 0
    best_uv = (0.0, 0.0)
    stack = Int[1]
    bounds = terrain_bounds(terrain)
    geometry_scale = maximum(ntuple(axis -> bounds[2][axis] - bounds[1][axis], 3))
    tolerance = 256eps(Float64) * max(geometry_scale, 1.0)
    while !isempty(stack)
        node = pop!(stack)
        _ray_aabb_intersects(
            origin,
            direction,
            terrain.bvh.lower[node],
            terrain.bvh.upper[node],
            best_distance,
        ) || continue
        if !isempty(terrain.bvh.leaves[node])
            for patch in terrain.bvh.leaves[node]
                triangle = terrain.triangles[patch]
                hit = _ray_triangle_intersection(
                    origin,
                    direction,
                    terrain.vertices[triangle[1]],
                    terrain.vertices[triangle[2]],
                    terrain.vertices[triangle[3]],
                    best_distance,
                )
                hit === nothing && continue
                if hit.distance < best_distance - tolerance ||
                   (abs(hit.distance - best_distance) <= tolerance && (best_patch == 0 || patch < best_patch))
                    best_distance = hit.distance
                    best_patch = patch
                    best_uv = (hit.u, hit.v)
                end
            end
        else
            # Push the higher index first so equal-distance traversal remains deterministic.
            terrain.bvh.right[node] > 0 && push!(stack, terrain.bvh.right[node])
            terrain.bvh.left[node] > 0 && push!(stack, terrain.bvh.left[node])
        end
    end
    best_patch == 0 && return nothing
    point = ntuple(axis -> origin[axis] + best_distance * direction[axis], 3)
    triangle = terrain.triangles[best_patch]
    normal = terrain.normal_mode == :facet ? terrain.facet_normals[best_patch] :
             _interpolated_triangle_normal(terrain, triangle, best_uv[1], best_uv[2])
    return TerrainHit(
        best_distance,
        point,
        normal,
        best_patch,
        terrain.material_ids[best_patch],
    )
end

function _periodic_tile_indices(lower, upper, origin, direction, maximum_distance)
    isfinite(maximum_distance) || throw(ArgumentError(
        "periodic terrain intersection requires a finite maximum_distance",
    ))
    endpoint = ntuple(axis -> origin[axis] + maximum_distance * direction[axis], 2)
    minimum_xy = (min(origin[1], endpoint[1]), min(origin[2], endpoint[2]))
    maximum_xy = (max(origin[1], endpoint[1]), max(origin[2], endpoint[2]))
    width = upper[1] - lower[1]
    height = upper[2] - lower[2]
    ix = floor(Int, (minimum_xy[1] - lower[1]) / width):floor(
        Int, (maximum_xy[1] - lower[1]) / width,
    )
    iy = floor(Int, (minimum_xy[2] - lower[2]) / height):floor(
        Int, (maximum_xy[2] - lower[2]) / height,
    )
    return ix, iy, width, height
end

intersect_terrain(::NoVoxelTerrain, origin, direction, maximum_distance::Real=Inf) = nothing

function intersect_terrain(
    terrain::PlanarTerrain,
    origin,
    direction,
    maximum_distance::Real=Inf,
)
    ray_origin = ntuple(axis -> Float64(origin[axis]), 3)
    ray_direction = _normalise_terrain_vector(direction, "terrain ray direction")
    limit = Float64(maximum_distance)
    isfinite(limit) && limit >= 0 || limit == Inf ||
        throw(ArgumentError("maximum terrain distance must be non-negative"))
    abs(ray_direction[3]) > eps(Float64) || return nothing
    distance = (terrain.elevation - ray_origin[3]) / ray_direction[3]
    tolerance = 256eps(Float64) * max(isfinite(limit) ? limit : abs(distance), 1.0)
    tolerance < distance <= limit + tolerance || return nothing
    x = ray_origin[1] + distance * ray_direction[1]
    y = ray_origin[2] + distance * ray_direction[2]
    if terrain.periodic
        x = terrain.x_bounds[1] + mod(x - terrain.x_bounds[1], terrain.x_bounds[2] - terrain.x_bounds[1])
        y = terrain.y_bounds[1] + mod(y - terrain.y_bounds[1], terrain.y_bounds[2] - terrain.y_bounds[1])
    else
        terrain.x_bounds[1] - tolerance <= x <= terrain.x_bounds[2] + tolerance || return nothing
        terrain.y_bounds[1] - tolerance <= y <= terrain.y_bounds[2] + tolerance || return nothing
        x = clamp(x, terrain.x_bounds[1], terrain.x_bounds[2])
        y = clamp(y, terrain.y_bounds[1], terrain.y_bounds[2])
    end
    nx, ny = terrain.cells
    i = clamp(floor(Int, (x - terrain.x_bounds[1]) /
                   (terrain.x_bounds[2] - terrain.x_bounds[1]) * nx) + 1, 1, nx)
    j = clamp(floor(Int, (y - terrain.y_bounds[1]) /
                   (terrain.y_bounds[2] - terrain.y_bounds[1]) * ny) + 1, 1, ny)
    patch = (j - 1) * nx + i
    return TerrainHit(
        max(distance, 0.0),
        (x, y, terrain.elevation),
        (0.0, 0.0, 1.0),
        patch,
        terrain.material_ids[i, j],
    )
end

function intersect_terrain(
    terrain::TriangulatedTerrain,
    origin,
    direction,
    maximum_distance::Real=Inf,
)
    ray_origin = ntuple(axis -> Float64(origin[axis]), 3)
    ray_direction = _normalise_terrain_vector(direction, "terrain ray direction")
    limit = Float64(maximum_distance)
    (isfinite(limit) && limit >= 0) || limit == Inf ||
        throw(ArgumentError("maximum terrain distance must be non-negative"))
    return _intersect_triangulated_base(terrain, ray_origin, ray_direction, limit)
end

function intersect_terrain(
    terrain::HeightFieldTerrain,
    origin,
    direction,
    maximum_distance::Real=Inf,
)
    ray_origin = ntuple(axis -> Float64(origin[axis]), 3)
    ray_direction = _normalise_terrain_vector(direction, "terrain ray direction")
    limit = Float64(maximum_distance)
    (isfinite(limit) && limit >= 0) || limit == Inf ||
        throw(ArgumentError("maximum terrain distance must be non-negative"))
    terrain.periodic || return _intersect_triangulated_base(
        terrain.mesh,
        ray_origin,
        ray_direction,
        limit,
    )

    lower, upper = terrain_bounds(terrain)
    ix, iy, width, height = _periodic_tile_indices(lower, upper, ray_origin, ray_direction, limit)
    best = nothing
    tolerance = _voxel_tolerance(limit)
    for tile_j in iy, tile_i in ix
        shifted_origin = (
            ray_origin[1] - tile_i * width,
            ray_origin[2] - tile_j * height,
            ray_origin[3],
        )
        hit = _intersect_triangulated_base(terrain.mesh, shifted_origin, ray_direction, limit)
        hit === nothing && continue
        if best === nothing || hit.distance < best.distance - tolerance ||
           (abs(hit.distance - best.distance) <= tolerance && hit.patch_id < best.patch_id)
            best = hit
        end
    end
    return best
end

function terrain_elevation_at(terrain::AbstractVoxelTerrain, x::Real, y::Real)
    bounds = terrain_bounds(terrain)
    bounds === nothing && throw(ArgumentError("open-bottom terrain has no elevation"))
    lower, upper = bounds
    margin = max(upper[3] - lower[3], 1.0)
    origin = (Float64(x), Float64(y), upper[3] + margin)
    hit = intersect_terrain(terrain, origin, (0.0, 0.0, -1.0), 2margin + (upper[3] - lower[3]))
    hit === nothing && throw(ArgumentError("terrain has no surface at ($x, $y)"))
    return hit.point[3]
end

function _terrain_periodic(terrain::AbstractVoxelTerrain)
    terrain isa PlanarTerrain && return terrain.periodic
    terrain isa HeightFieldTerrain && return terrain.periodic
    return false
end

function _validate_terrain_for_grid(
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain,
    boundary::Symbol,
)
    terrain isa NoVoxelTerrain && return terrain
    boundary in (:open, :periodic) || throw(ArgumentError(
        "terrain-aware voxel transport supports only :open and :periodic boundaries",
    ))
    lower, upper = terrain_bounds(terrain)
    scale = maximum(ntuple(axis -> grid.maximum[axis] - grid.minimum[axis], 3))
    tolerance = 256eps(Float64) * max(scale, 1.0)
    lower[3] >= grid.minimum[3] - tolerance || throw(ArgumentError(
        "terrain extends below the voxel-grid minimum and cannot intercept rays before bottom escape",
    ))
    upper[3] < grid.maximum[3] - tolerance || throw(ArgumentError(
        "terrain must remain below the voxel-grid top boundary",
    ))
    if boundary == :periodic
        _terrain_periodic(terrain) || throw(ArgumentError(
            "periodic voxel boundaries require terrain constructed with periodic=true",
        ))
        for axis in 1:2
            isapprox(lower[axis], grid.minimum[axis]; atol=tolerance, rtol=0) &&
            isapprox(upper[axis], grid.maximum[axis]; atol=tolerance, rtol=0) ||
                throw(ArgumentError("periodic terrain bounds must equal the voxel x/y bounds"))
        end
    else
        lower[1] <= grid.minimum[1] + tolerance && upper[1] >= grid.maximum[1] - tolerance ||
            throw(ArgumentError("terrain does not cover the full voxel x domain"))
        lower[2] <= grid.minimum[2] + tolerance && upper[2] >= grid.maximum[2] - tolerance ||
            throw(ArgumentError("terrain does not cover the full voxel y domain"))
    end
    return terrain
end

_terrain_geometry_fingerprint(::NoVoxelTerrain) = UInt(0)

function _terrain_geometry_fingerprint(terrain::PlanarTerrain)
    state = hash(:planar_terrain)
    for value in (terrain.elevation, terrain.x_bounds, terrain.y_bounds, terrain.cells, terrain.periodic)
        state = hash(value, state)
    end
    for id in terrain.material_ids
        state = hash(id, state)
    end
    return UInt(state)
end

function _terrain_geometry_fingerprint(terrain::TriangulatedTerrain)
    state = hash(:triangulated_terrain, hash(terrain.normal_mode))
    for vertex in terrain.vertices
        state = hash(vertex, state)
    end
    for triangle in terrain.triangles
        state = hash(triangle, state)
    end
    for id in terrain.material_ids
        state = hash(id, state)
    end
    return UInt(state)
end

function _terrain_geometry_fingerprint(terrain::HeightFieldTerrain)
    state = hash(:height_field_terrain, hash((terrain.diagonal, terrain.periodic)))
    state = hash(_terrain_geometry_fingerprint(terrain.mesh), state)
    return UInt(state)
end

function terrain_lambertian_weights(
    quadrature::VoxelScatteringQuadrature,
    normal,
)
    _validate_voxel_scattering_quadrature(quadrature)
    unit_normal = _normalise_terrain_vector(normal, "terrain normal")
    direction_indices = Int[]
    weights = Float64[]
    for index in eachindex(quadrature.directions)
        cosine = sum(unit_normal[axis] * quadrature.directions[index][axis] for axis in 1:3)
        cosine > 0 || continue
        push!(direction_indices, index)
        push!(weights, quadrature.weights[index] * cosine)
    end
    total = sum(weights)
    total > 0 || throw(ArgumentError(
        "voxel scattering quadrature has no support above the terrain normal",
    ))
    weights ./= total
    return direction_indices, weights
end

function _terrain_launch_point(grid::VoxelGrid, terrain, patch_id::Integer)
    point = terrain_patch_point(terrain, patch_id)
    normal = terrain_patch_normal(terrain, patch_id)
    scale = maximum(ntuple(axis -> grid.maximum[axis] - grid.minimum[axis], 3))
    offset = 1024eps(Float64) * max(scale, 1.0)
    return ntuple(axis -> point[axis] + offset * normal[axis], 3)
end

function _truncate_voxel_path_at_terrain(
    grid::VoxelGrid,
    path::VoxelRayPath,
    origin,
    direction,
    terrain::AbstractVoxelTerrain,
    boundary::Symbol,
)
    terrain isa NoVoxelTerrain && return path
    _validate_terrain_for_grid(grid, terrain, boundary)
    surface_at_origin = terrain_elevation_at(terrain, origin[1], origin[2])
    origin_tolerance = 512eps(Float64) * max(abs(surface_at_origin), abs(origin[3]), 1.0)
    origin[3] >= surface_at_origin - origin_tolerance || throw(ArgumentError(
        "voxel ray origin lies below the terrain surface",
    ))
    if isempty(path.segments)
        path.exit_reason == :bottom && throw(ArgumentError(
            "ray reached the voxel bottom without intersecting the supplied terrain",
        ))
        return path
    end
    maximum_distance = path.segments[end].t_exit
    hit = intersect_terrain(terrain, origin, direction, maximum_distance)
    if hit === nothing
        path.exit_reason == :bottom && throw(ArgumentError(
            "terrain has a horizontal coverage gap along a ray reaching the voxel bottom",
        ))
        return path
    end
    tolerance = _voxel_tolerance(maximum_distance)
    hit.distance <= maximum_distance + tolerance || return path

    segments = VoxelRaySegment{Float64}[]
    sizehint!(segments, length(path.segments))
    for segment in path.segments
        segment.t_enter < hit.distance - tolerance || break
        segment_exit = min(segment.t_exit, hit.distance)
        _append_voxel_segment!(segments, segment.index, segment.t_enter, segment_exit)
        segment.t_exit >= hit.distance - tolerance && break
    end
    return VoxelRayPath{Float64}(segments, :terrain, hit)
end
