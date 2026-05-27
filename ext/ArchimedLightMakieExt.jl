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

_payload_geometry(payload::Union{LightStepResult,AbstractVector{<:LightStepResult}}) = ArchimedLight.light_render_geometry(payload)
_payload_geometry(payload::Tuple{LightRenderGeometry,Any}) = first(payload)
_payload_data(payload::Union{LightStepResult,AbstractVector{<:LightStepResult}}) = payload
_payload_data(payload::Tuple{LightRenderGeometry,Any}) = last(payload)

function Makie.plot!(plot::LightPlot)
    map!(plot.attributes, :payload, :light_geometry) do payload
        _payload_geometry(payload)
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

_lightplot_kwargs(kwargs) = _lightplot_kwargs(NamedTuple(kwargs))
_lightplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

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

end
