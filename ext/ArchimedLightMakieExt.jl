module ArchimedLightMakieExt

using ArchimedLight
using GeometryBasics
using Makie
import PlantGeom
import Makie: Attributes, automatic, colormap_attributes!, generic_plot_attributes!, shading_attributes!

Makie.@recipe(LightPlot, payload) do scene
    attr = Attributes(
        color=:incident_par_flux,
        timestep=1,
        interpolate=false,
        fill_value=NaN,
        cycle=[],
        colormap=:viridis,
        colorscale=identity,
        colorrange=automatic,
        lowclip=automatic,
        highclip=automatic,
        nan_color=:transparent,
        alpha=1.0,
    )
    shading_attributes!(attr)
    generic_plot_attributes!(attr)
    return colormap_attributes!(attr, :viridis)
end

const _LIGHTPLOT_ATTRIBUTES = Makie.@DocumentedAttributes begin
    color = :incident_par_flux
    timestep = 1
    interpolate = false
    fill_value = NaN
    cycle = []
    colormap = :viridis
    colorscale = identity
    colorrange = automatic
    lowclip = automatic
    highclip = automatic
    nan_color = :transparent
    alpha = 1.0
    shading = true
    diffuse = 1.0
    specular = 0.2
    shininess = 32.0f0
    backlight = 0.0f0
    transformation = :automatic
    model = automatic
    visible = true
    transparency = false
    overdraw = false
    ssao = false
    inspectable = true
    depth_shift = 0.0f0
    space = :data
    fxaa = true
    inspector_label = automatic
    inspector_clear = automatic
    inspector_hover = automatic
    clip_planes = automatic
end

Makie.documented_attributes(::Type{<:LightPlot}) = _LIGHTPLOT_ATTRIBUTES

Makie.preferred_axis_type(::LightPlot) = Makie.LScene

_payload_geometry(
    payload::Union{LightStepResult,AbstractVector{<:LightStepResult}},
    timestep,
) = ArchimedLight._geometry_for_values(payload, Int(timestep))
_payload_geometry(payload::Tuple{LightRenderGeometry,Any}, _) = first(payload)
_payload_data(payload::Union{LightStepResult,AbstractVector{<:LightStepResult}}) = payload
_payload_data(payload::Tuple{LightRenderGeometry,Any}) = last(payload)

function Makie.plot!(plot::LightPlot)
    map!(plot.attributes, [:payload, :timestep], :light_geometry) do payload, timestep
        _payload_geometry(payload, timestep)
    end
    map!(plot.attributes, :light_geometry, :light_base_mesh) do geometry
        points = GeometryBasics.Point3d[
            GeometryBasics.Point3d(v[1], v[2], v[3]) for v in geometry.vertices
        ]
        GeometryBasics.Mesh(points, geometry.faces)
    end
    map!(
        plot.attributes,
        [:light_geometry, :payload, :color, :timestep, :interpolate, :fill_value],
        :light_color,
    ) do geometry, payload, color, timestep, interpolate, fill_value
        ArchimedLight._light_color_values(
            geometry,
            _payload_data(payload),
            color;
            timestep=Int(timestep),
            interpolate=Bool(interpolate),
            fill_value=Float64(fill_value),
        )
    end
    map!(
        plot.attributes,
        [:payload, :color, :fill_value, :colorrange],
        :light_colorrange,
    ) do payload, color, fill_value, colorrange
        colorrange === automatic || return colorrange
        default_colorrange = ArchimedLight._automatic_light_colorrange(
            _payload_data(payload),
            color;
            fill_value=Float64(fill_value),
        )
        isnothing(default_colorrange) ? automatic : default_colorrange
    end
    map!(
        plot.attributes,
        [:light_base_mesh, :light_color, :interpolate],
        :light_mesh,
    ) do mesh, color, interpolate
        float_color = Float32.(color)
        mesh_color = Bool(interpolate) ? float_color : GeometryBasics.per_face(float_color, GeometryBasics.faces(mesh))
        GeometryBasics.mesh(mesh, color=mesh_color)
    end

    mesh!(
        plot,
        plot[:light_mesh];
        interpolate=plot[:interpolate],
        shading=plot[:shading],
        colormap=plot[:colormap],
        colorrange=plot[:light_colorrange],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        nan_color=plot[:nan_color],
        alpha=plot[:alpha],
    )
    return plot
end

Makie.@recipe(VoxelPlot, payload) do scene
    attr = Attributes(
        rays=(),
        quadrature=nothing,
        boundary=:open,
        pad_threshold=0.0,
        cutaway=:none,
        voxel_alpha=0.35,
        voxel_gap=0.08,
        voxel_colormap=:YlGn,
        voxel_colorrange=automatic,
        terrain_alpha=0.92,
        terrain_colormap=:copper,
        terrain_colorrange=(0.0, 1.0),
        lattice_color=(:slategray, 0.38),
        lattice_linewidth=0.7,
        terrain_wireframe_color=(:black, 0.30),
        terrain_wireframe_linewidth=0.55,
        incoming_color=:gold,
        incoming_linewidth=3.0,
        hit_color=:orangered,
        reflected_colormap=:plasma,
        normal_color=:dodgerblue,
        normal_length=0.12,
        reflected_length=0.30,
        arrow_markerscale=0.025,
        show_voxels=true,
        show_lattice=true,
        show_terrain=true,
        show_rays=true,
        show_normals=false,
        cycle=[],
    )
    shading_attributes!(attr)
    generic_plot_attributes!(attr)
    return attr
end

const _VOXELPLOT_ATTRIBUTES = Makie.@DocumentedAttributes begin
    rays = ()
    quadrature = nothing
    boundary = :open
    pad_threshold = 0.0
    cutaway = :none
    voxel_alpha = 0.35
    voxel_gap = 0.08
    voxel_colormap = :YlGn
    voxel_colorrange = automatic
    terrain_alpha = 0.92
    terrain_colormap = :copper
    terrain_colorrange = (0.0, 1.0)
    lattice_color = (:slategray, 0.38)
    lattice_linewidth = 0.7
    terrain_wireframe_color = (:black, 0.30)
    terrain_wireframe_linewidth = 0.55
    incoming_color = :gold
    incoming_linewidth = 3.0
    hit_color = :orangered
    reflected_colormap = :plasma
    normal_color = :dodgerblue
    normal_length = 0.12
    reflected_length = 0.30
    arrow_markerscale = 0.025
    show_voxels = true
    show_lattice = true
    show_terrain = true
    show_rays = true
    show_normals = false
    cycle = []
    shading = true
    diffuse = 1.0
    specular = 0.2
    shininess = 32.0f0
    backlight = 0.0f0
    transformation = :automatic
    model = automatic
    visible = true
    transparency = false
    overdraw = false
    ssao = false
    inspectable = true
    depth_shift = 0.0f0
    space = :data
    inspector_label = automatic
    inspector_clear = automatic
    inspector_hover = automatic
    clip_planes = automatic
end

Makie.documented_attributes(::Type{<:VoxelPlot}) = _VOXELPLOT_ATTRIBUTES

Makie.preferred_axis_type(::VoxelPlot) = Makie.LScene

function _voxelplot_grid_terrain(payload)
    payload isa Tuple && length(payload) == 2 || throw(ArgumentError(
        "voxelplot payload must contain a VoxelGrid and an AbstractVoxelTerrain",
    ))
    grid, terrain = payload
    grid isa VoxelGrid || throw(ArgumentError("voxelplot requires a VoxelGrid"))
    terrain isa AbstractVoxelTerrain || throw(ArgumentError(
        "voxelplot terrain must be an AbstractVoxelTerrain",
    ))
    return grid, terrain
end

function _voxelplot_pad(grid::VoxelGrid, cutaway)
    cutaway in (:none, :front) || throw(ArgumentError(
        "voxelplot cutaway must be :none or :front",
    ))
    pad = Float32.(grid.pad)
    if cutaway == :front
        pad[:, 1:cld(size(pad, 2), 2), :] .= 0.0f0
    end
    return pad
end

function _voxelplot_lattice(grid::VoxelGrid)
    nx, ny, nz = size(grid)
    xs = range(grid.minimum[1], grid.maximum[1]; length=nx + 1)
    ys = range(grid.minimum[2], grid.maximum[2]; length=ny + 1)
    zs = range(grid.minimum[3], grid.maximum[3]; length=nz + 1)
    points = GeometryBasics.Point3f[]
    for y in ys, z in zs
        push!(points, GeometryBasics.Point3f(first(xs), y, z))
        push!(points, GeometryBasics.Point3f(last(xs), y, z))
    end
    for x in xs, z in zs
        push!(points, GeometryBasics.Point3f(x, first(ys), z))
        push!(points, GeometryBasics.Point3f(x, last(ys), z))
    end
    for x in xs, y in ys
        push!(points, GeometryBasics.Point3f(x, y, first(zs)))
        push!(points, GeometryBasics.Point3f(x, y, last(zs)))
    end
    return points
end

_empty_voxel_terrain_mesh() = GeometryBasics.Mesh(
    GeometryBasics.Point3f[],
    GeometryBasics.TriangleFace{Int}[],
)

function _planar_terrain_mesh(terrain::PlanarTerrain)
    nx, ny = terrain.cells
    dx = (terrain.x_bounds[2] - terrain.x_bounds[1]) / nx
    dy = (terrain.y_bounds[2] - terrain.y_bounds[1]) / ny
    points = GeometryBasics.Point3f[]
    faces = GeometryBasics.TriangleFace{Int}[]
    reflectance = Float32[]
    for j in 1:ny, i in 1:nx
        x0 = terrain.x_bounds[1] + (i - 1) * dx
        x1 = x0 + dx
        y0 = terrain.y_bounds[1] + (j - 1) * dy
        y1 = y0 + dy
        base = length(points)
        push!(points, GeometryBasics.Point3f(x0, y0, terrain.elevation))
        push!(points, GeometryBasics.Point3f(x1, y0, terrain.elevation))
        push!(points, GeometryBasics.Point3f(x1, y1, terrain.elevation))
        push!(points, GeometryBasics.Point3f(x0, y1, terrain.elevation))
        push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 2, base + 3))
        push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 3, base + 4))
        patch_id = (j - 1) * nx + i
        value = soil_optics(terrain, terrain_patch_material(terrain, patch_id)).par_reflectance
        push!(reflectance, Float32(value), Float32(value))
    end
    return GeometryBasics.Mesh(points, faces), reflectance
end

function _triangulated_terrain_mesh(terrain::TriangulatedTerrain)
    points = GeometryBasics.Point3f[
        GeometryBasics.Point3f(point[1], point[2], point[3]) for point in terrain.vertices
    ]
    faces = GeometryBasics.TriangleFace{Int}[
        GeometryBasics.TriangleFace{Int}(face[1], face[2], face[3]) for face in terrain.triangles
    ]
    reflectance = Float32[
        soil_optics(terrain, material_id).par_reflectance for material_id in terrain.material_ids
    ]
    return GeometryBasics.Mesh(points, faces), reflectance
end

_voxelplot_terrain_mesh(::NoVoxelTerrain) = (_empty_voxel_terrain_mesh(), Float32[])
_voxelplot_terrain_mesh(terrain::PlanarTerrain) = _planar_terrain_mesh(terrain)
_voxelplot_terrain_mesh(terrain::TriangulatedTerrain) = _triangulated_terrain_mesh(terrain)
_voxelplot_terrain_mesh(terrain::HeightFieldTerrain) = _triangulated_terrain_mesh(terrain.mesh)

function _voxelplot_ray(ray, index::Integer, grid::VoxelGrid)
    hasproperty(ray, :origin) || throw(ArgumentError(
        "voxelplot ray $index must define an origin field",
    ))
    hasproperty(ray, :direction) || throw(ArgumentError(
        "voxelplot ray $index must define a direction field",
    ))
    raw_origin = getproperty(ray, :origin)
    raw_direction = getproperty(ray, :direction)
    length(raw_origin) == 3 || throw(ArgumentError(
        "voxelplot ray $index origin must contain three coordinates",
    ))
    length(raw_direction) == 3 || throw(ArgumentError(
        "voxelplot ray $index direction must contain three coordinates",
    ))
    origin = ntuple(axis -> Float64(raw_origin[axis]), 3)
    direction = ntuple(axis -> Float64(raw_direction[axis]), 3)
    all(isfinite, origin) || throw(ArgumentError("voxelplot ray $index origin must be finite"))
    all(isfinite, direction) || throw(ArgumentError("voxelplot ray $index direction must be finite"))
    all(axis -> grid.minimum[axis] <= origin[axis] <= grid.maximum[axis], 1:3) ||
        throw(ArgumentError("voxelplot ray $index origin must lie inside the voxel grid"))
    magnitude = sqrt(sum(abs2, direction))
    magnitude > 0 || throw(ArgumentError("voxelplot ray $index direction cannot be zero"))
    unit_direction = ntuple(axis -> direction[axis] / magnitude, 3)
    unit_direction[3] < 0 || throw(ArgumentError(
        "voxelplot ray $index direction must point downward",
    ))
    return origin, unit_direction
end

function _voxelplot_segment_points(grid, origin, direction, segment, boundary)
    start = [origin[axis] + segment.t_enter * direction[axis] for axis in 1:3]
    stop = [origin[axis] + segment.t_exit * direction[axis] for axis in 1:3]
    if boundary == :periodic
        midpoint_t = (segment.t_enter + segment.t_exit) / 2
        for axis in 1:2
            width = grid.maximum[axis] - grid.minimum[axis]
            midpoint = origin[axis] + midpoint_t * direction[axis]
            wrapped_midpoint = mod(midpoint - grid.minimum[axis], width) + grid.minimum[axis]
            shift = wrapped_midpoint - midpoint
            start[axis] += shift
            stop[axis] += shift
        end
    end
    return GeometryBasics.Point3f(start...), GeometryBasics.Point3f(stop...)
end

function _voxelplot_ray_geometry(
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain,
    rays,
    boundary,
    quadrature,
    normal_length,
    reflected_length,
)
    boundary in (:open, :periodic) || throw(ArgumentError(
        "voxelplot boundary must be :open or :periodic",
    ))
    quadrature === nothing || quadrature isa VoxelScatteringQuadrature ||
        throw(ArgumentError("voxelplot quadrature must be a VoxelScatteringQuadrature or nothing"))
    diagonal = sqrt(sum((grid.maximum[axis] - grid.minimum[axis])^2 for axis in 1:3))
    ray_paths = VoxelRayPath[]
    ray_segments = GeometryBasics.Point3f[]
    hit_points = GeometryBasics.Point3f[]
    normal_origins = GeometryBasics.Point3f[]
    normal_directions = GeometryBasics.Vec3f[]
    reflected_origins = GeometryBasics.Point3f[]
    reflected_directions = GeometryBasics.Vec3f[]
    reflected_weights = Float32[]

    for (index, ray) in enumerate(rays)
        origin, direction = _voxelplot_ray(ray, index, grid)
        path = trace_voxel_ray(
            grid,
            origin,
            direction;
            boundary=boundary,
            terrain=terrain,
        )
        push!(ray_paths, path)
        for segment in path.segments
            start, stop = _voxelplot_segment_points(
                grid,
                origin,
                direction,
                segment,
                boundary,
            )
            push!(ray_segments, start, stop)
        end
        hit = path.terrain_hit
        hit === nothing && continue
        hit_point = GeometryBasics.Point3f(hit.point...)
        push!(hit_points, hit_point)
        push!(normal_origins, hit_point)
        push!(normal_directions, GeometryBasics.Vec3f(
            (diagonal * Float64(normal_length) * hit.normal[axis] for axis in 1:3)...,
        ))
        quadrature === nothing && continue
        direction_ids, weights = terrain_lambertian_weights(quadrature, hit.normal)
        maximum_weight = maximum(weights)
        for (direction_id, weight) in zip(direction_ids, weights)
            reflected = quadrature.directions[direction_id]
            scale = diagonal * Float64(reflected_length) * weight / maximum_weight
            push!(reflected_origins, hit_point)
            push!(reflected_directions, GeometryBasics.Vec3f(
                (scale * reflected[axis] for axis in 1:3)...,
            ))
            push!(reflected_weights, Float32(weight))
        end
    end
    return (
        ray_paths,
        ray_segments,
        hit_points,
        normal_origins,
        normal_directions,
        reflected_origins,
        reflected_directions,
        reflected_weights,
    )
end

function Makie.plot!(plot::VoxelPlot)
    grid, _ = _voxelplot_grid_terrain(Makie.to_value(plot[:payload]))
    x = (Float32(grid.minimum[1]), Float32(grid.maximum[1]))
    y = (Float32(grid.minimum[2]), Float32(grid.maximum[2]))
    z = (Float32(grid.minimum[3]), Float32(grid.maximum[3]))

    map!(plot.attributes, :payload, [:voxel_pad_full, :voxel_pad_front]) do payload
        current_grid, _ = _voxelplot_grid_terrain(payload)
        (_voxelplot_pad(current_grid, :none), _voxelplot_pad(current_grid, :front))
    end
    map!(
        plot.attributes,
        [:show_voxels, :cutaway],
        [:voxel_full_visible, :voxel_front_visible],
    ) do show_voxels, cutaway
        cutaway in (:none, :front) || throw(ArgumentError(
            "voxelplot cutaway must be :none or :front",
        ))
        (Bool(show_voxels) && cutaway == :none, Bool(show_voxels) && cutaway == :front)
    end
    map!(plot.attributes, :pad_threshold, :voxel_is_air) do threshold
        limit = Float32(threshold)
        isfinite(limit) && limit >= 0 || throw(ArgumentError(
            "voxelplot pad_threshold must be finite and non-negative",
        ))
        value -> !isfinite(value) || value <= limit
    end
    map!(plot.attributes, :payload, :voxel_lattice) do payload
        current_grid, _ = _voxelplot_grid_terrain(payload)
        _voxelplot_lattice(current_grid)
    end
    map!(plot.attributes, :payload, [:voxel_terrain_mesh, :voxel_terrain_values]) do payload
        _, current_terrain = _voxelplot_grid_terrain(payload)
        _voxelplot_terrain_mesh(current_terrain)
    end
    map!(
        plot.attributes,
        [:voxel_terrain_mesh, :voxel_terrain_values],
        :voxel_colored_terrain_mesh,
    ) do terrain_mesh, values
        face_values = GeometryBasics.per_face(values, GeometryBasics.faces(terrain_mesh))
        GeometryBasics.mesh(terrain_mesh; color=face_values)
    end
    map!(
        plot.attributes,
        [:payload, :rays, :boundary, :quadrature, :normal_length, :reflected_length],
        [
            :voxel_ray_paths,
            :voxel_ray_segments,
            :voxel_hit_points,
            :voxel_normal_origins,
            :voxel_normal_directions,
            :voxel_reflected_origins,
            :voxel_reflected_directions,
            :voxel_reflected_weights,
        ],
    ) do payload, rays, boundary, quadrature, normal_length, reflected_length
        current_grid, current_terrain = _voxelplot_grid_terrain(payload)
        _voxelplot_ray_geometry(
            current_grid,
            current_terrain,
            rays,
            boundary,
            quadrature,
            normal_length,
            reflected_length,
        )
    end
    map!(plot.attributes, [:payload, :arrow_markerscale], :voxel_arrow_markerscale) do payload, fraction
        current_grid, _ = _voxelplot_grid_terrain(payload)
        diagonal = sqrt(sum(
            (current_grid.maximum[axis] - current_grid.minimum[axis])^2 for axis in 1:3
        ))
        diagonal * Float64(fraction)
    end

    voxels!(
        plot,
        x,
        y,
        z,
        plot[:voxel_pad_full];
        is_air=plot[:voxel_is_air],
        colormap=plot[:voxel_colormap],
        colorrange=plot[:voxel_colorrange],
        alpha=plot[:voxel_alpha],
        gap=plot[:voxel_gap],
        depthsorting=true,
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:voxel_full_visible],
    )
    voxels!(
        plot,
        x,
        y,
        z,
        plot[:voxel_pad_front];
        is_air=plot[:voxel_is_air],
        colormap=plot[:voxel_colormap],
        colorrange=plot[:voxel_colorrange],
        alpha=plot[:voxel_alpha],
        gap=plot[:voxel_gap],
        depthsorting=true,
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:voxel_front_visible],
    )
    linesegments!(
        plot,
        plot[:voxel_lattice];
        color=plot[:lattice_color],
        linewidth=plot[:lattice_linewidth],
        transparency=true,
        visible=plot[:show_lattice],
    )
    mesh!(
        plot,
        plot[:voxel_colored_terrain_mesh];
        colormap=plot[:terrain_colormap],
        colorrange=plot[:terrain_colorrange],
        alpha=plot[:terrain_alpha],
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:show_terrain],
    )
    wireframe!(
        plot,
        plot[:voxel_terrain_mesh];
        color=plot[:terrain_wireframe_color],
        linewidth=plot[:terrain_wireframe_linewidth],
        transparency=true,
        visible=plot[:show_terrain],
    )
    linesegments!(
        plot,
        plot[:voxel_ray_segments];
        color=plot[:incoming_color],
        linewidth=plot[:incoming_linewidth],
        transparency=true,
        visible=plot[:show_rays],
    )
    scatter!(
        plot,
        plot[:voxel_hit_points];
        color=plot[:hit_color],
        markersize=13,
        strokecolor=:white,
        strokewidth=1,
        visible=plot[:show_rays],
    )
    arrows3d!(
        plot,
        plot[:voxel_reflected_origins],
        plot[:voxel_reflected_directions];
        color=plot[:voxel_reflected_weights],
        colormap=plot[:reflected_colormap],
        colorrange=(0.0, 1.0),
        markerscale=plot[:voxel_arrow_markerscale],
        transparency=true,
        visible=plot[:show_rays],
    )
    arrows3d!(
        plot,
        plot[:voxel_normal_origins],
        plot[:voxel_normal_directions];
        color=plot[:normal_color],
        markerscale=plot[:voxel_arrow_markerscale],
        visible=plot[:show_normals],
    )
    return plot
end

_lightplot_kwargs(kwargs) = _lightplot_kwargs(NamedTuple(kwargs))
_lightplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

_voxelplot_kwargs(kwargs) = _voxelplot_kwargs(NamedTuple(kwargs))
_voxelplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

function _validate_voxelplot_kwargs(grid::VoxelGrid, kwargs::NamedTuple)
    _voxelplot_pad(grid, get(kwargs, :cutaway, :none))
    boundary = get(kwargs, :boundary, :open)
    boundary in (:open, :periodic) || throw(ArgumentError(
        "voxelplot boundary must be :open or :periodic",
    ))
    quadrature = get(kwargs, :quadrature, nothing)
    quadrature === nothing || quadrature isa VoxelScatteringQuadrature || throw(ArgumentError(
        "voxelplot quadrature must be a VoxelScatteringQuadrature or nothing",
    ))
    for (index, ray) in enumerate(get(kwargs, :rays, ()))
        _voxelplot_ray(ray, index, grid)
    end
    for (name, value, allow_zero) in (
        (:pad_threshold, get(kwargs, :pad_threshold, 0.0), true),
        (:normal_length, get(kwargs, :normal_length, 0.12), false),
        (:reflected_length, get(kwargs, :reflected_length, 0.30), false),
        (:arrow_markerscale, get(kwargs, :arrow_markerscale, 0.025), false),
    )
        number = Float64(value)
        valid = isfinite(number) && (allow_zero ? number >= 0 : number > 0)
        valid || throw(ArgumentError(
            "voxelplot $name must be finite and $(allow_zero ? "non-negative" : "positive")",
        ))
    end
    return kwargs
end

function ArchimedLight.lightplot(data::Union{LightStepResult,AbstractVector{<:LightStepResult}}; kwargs...)
    lightplot(data; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot(geometry::LightRenderGeometry, data; kwargs...)
    lightplot((geometry, data); _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, data::Union{LightStepResult,AbstractVector{<:LightStepResult}}; kwargs...)
    lightplot!(axis, data; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, geometry::LightRenderGeometry, data; kwargs...)
    lightplot!(axis, (geometry, data); _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions, data; kwargs...)
    payload = (ArchimedLight.light_render_geometry(scene, models, options), data)
    lightplot(payload; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions, data; kwargs...)
    payload = (ArchimedLight.light_render_geometry(scene, models, options), data)
    lightplot!(axis, payload; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.voxelplot(
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    kwargs...,
)
    options = _validate_voxelplot_kwargs(grid, _voxelplot_kwargs(kwargs))
    voxelplot((grid, terrain); options...)
end

function ArchimedLight.voxelplot!(
    parent,
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    kwargs...,
)
    options = _validate_voxelplot_kwargs(grid, _voxelplot_kwargs(kwargs))
    voxelplot!(parent, (grid, terrain); options...)
end

end
