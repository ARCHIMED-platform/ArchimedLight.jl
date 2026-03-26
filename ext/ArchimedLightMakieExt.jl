module ArchimedLightMakieExt

using ArchimedLight
using GeometryBasics
using Makie
import Makie: Attributes, automatic, colormap_attributes!, generic_plot_attributes!, shading_attributes!

Makie.@recipe(LightPlot, scene, models, options, data) do scene
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

Makie.preferred_axis_type(::LightPlot) = Makie.Axis3

function Makie.plot!(plot::LightPlot)
    map!(plot.attributes, [:scene, :models, :options], :light_geometry) do scene, models, options
        ArchimedLight._scene_geometry_for_interception(scene, models, options)
    end
    map!(plot.attributes, :light_geometry, :light_mesh) do geometry
        points = GeometryBasics.Point3d[
            GeometryBasics.Point3d(v[1], v[2], v[3]) for v in geometry.vertices
        ]
        GeometryBasics.Mesh(points, geometry.faces)
    end
    map!(
        plot.attributes,
        [:light_geometry, :data, :color, :timestep, :interpolate, :fill_value],
        :light_color,
    ) do geometry, data, color, timestep, interpolate, fill_value
        ArchimedLight._light_color_values(
            geometry,
            data,
            color;
            timestep=Int(timestep),
            interpolate=Bool(interpolate),
            fill_value=Float64(fill_value),
        )
    end

    mesh!(
        plot,
        plot[:light_mesh];
        color=plot[:light_color],
        interpolate=plot[:interpolate],
        shading=plot[:shading],
        colormap=plot[:colormap],
        colorrange=plot[:colorrange],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        nan_color=plot[:nan_color],
        alpha=plot[:alpha],
    )
    return plot
end

_lightplot_kwargs(kwargs) = _lightplot_kwargs(NamedTuple(kwargs))
_lightplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

function ArchimedLight.lightplot(scene, models, options, data; kwargs...)
    lightplot(scene, models, options, data; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, scene, models, options, data; kwargs...)
    lightplot!(axis, scene, models, options, data; _lightplot_kwargs(kwargs)...)
end

end
