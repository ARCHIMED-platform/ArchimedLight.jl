const _ATTR_TO_BUDGET_FIELD = Dict{Symbol,Symbol}(attr => field for (field, attr) in _DEFAULT_BUDGET_ATTRS)

function _normalize_light_selector(selector::Symbol)
    haskey(_DEFAULT_BUDGET_ATTRS, selector) && return selector
    haskey(_ATTR_TO_BUDGET_FIELD, selector) && return _ATTR_TO_BUDGET_FIELD[selector]
    error(
        "Unknown light selector `$selector`. Use a LightBudget field like " *
        "`:incident_par_flux` or an ARCHIMED attribute like `:Ri_PAR_f`.",
    )
end

function _series_index_error(timestep::Integer, n::Integer)
    throw(BoundsError(1:n, timestep))
end

"""
    light_metric_values(step, selector)
    light_metric_values(steps, selector; timestep=1)

Return the per-node metric map selected by `selector`.

`selector` may be either a runtime budget field such as `:incident_par_flux` or
the corresponding ARCHIMED attribute name such as `:Ri_PAR_f`.
"""
function light_metric_values(step::LightStepResult, selector::Symbol)
    _budget_node_field(step, _normalize_light_selector(selector))
end

function light_metric_values(
    steps::AbstractVector{<:LightStepResult},
    selector::Symbol;
    timestep::Integer=1,
)
    1 <= timestep <= length(steps) || _series_index_error(timestep, length(steps))
    light_metric_values(steps[timestep], selector)
end

function _metric_map(
    data::LightStepResult,
    selector::Symbol;
    timestep::Integer=1,
)
    timestep == 1 || _series_index_error(timestep, 1)
    light_metric_values(data, selector)
end

function _metric_map(
    data::AbstractVector{<:LightStepResult},
    selector::Symbol;
    timestep::Integer=1,
)
    light_metric_values(data, selector; timestep=timestep)
end

function _metric_map(
    data::AbstractDict{<:Integer,<:Real},
    ::Symbol;
    timestep::Integer=1,
)
    timestep == 1 || _series_index_error(timestep, 1)
    data
end

function _metric_map(
    data::AbstractVector{<:AbstractDict{<:Integer,<:Real}},
    ::Symbol;
    timestep::Integer=1,
)
    1 <= timestep <= length(data) || _series_index_error(timestep, length(data))
    data[timestep]
end

function _face_values_from_metric_map(
    face2node::Vector{Int},
    metric_map::AbstractDict{<:Integer,<:Real};
    fill_value::Float64=NaN,
)
    Float64[Float64(get(metric_map, nid, fill_value)) for nid in face2node]
end

function _vertex_values_from_metric_map(
    vertices,
    faces,
    face2node::Vector{Int},
    metric_map::AbstractDict{<:Integer,<:Real};
    fill_value::Float64=NaN,
)
    v_sum = zeros(Float64, length(vertices))
    v_count = zeros(Int, length(vertices))
    for i in eachindex(faces)
        f = faces[i]
        v = Float64(get(metric_map, face2node[i], fill_value))
        isfinite(v) || continue
        for vid in (Int(f[1]), Int(f[2]), Int(f[3]))
            v_sum[vid] += v
            v_count[vid] += 1
        end
    end
    Float64[v_count[i] > 0 ? (v_sum[i] / v_count[i]) : fill_value for i in eachindex(vertices)]
end

function _light_color_values(
    geometry,
    data,
    color;
    timestep::Integer=1,
    interpolate::Bool=false,
    fill_value::Float64=NaN,
)
    if color isa Symbol
        metric_map = _metric_map(data, color; timestep=timestep)
        return interpolate ?
            _vertex_values_from_metric_map(geometry.vertices, geometry.faces, geometry.face2node, metric_map; fill_value=fill_value) :
            _face_values_from_metric_map(geometry.face2node, metric_map; fill_value=fill_value)
    elseif color isa AbstractDict{<:Integer,<:Real}
        return interpolate ?
            _vertex_values_from_metric_map(geometry.vertices, geometry.faces, geometry.face2node, color; fill_value=fill_value) :
            _face_values_from_metric_map(geometry.face2node, color; fill_value=fill_value)
    elseif color isa AbstractVector{<:Real}
        values = Float64.(color)
        expected = interpolate ? length(geometry.vertices) : length(geometry.faces)
        length(values) == expected || error(
            "Explicit color vector has length $(length(values)), expected $(expected) " *
            "for interpolate=$(interpolate).",
        )
        return values
    else
        error(
            "Unsupported `color` input of type $(typeof(color)). Use a selector Symbol, " *
            "a Dict keyed by node id, or an explicit numeric color vector.",
        )
    end
end

"""
    light_face_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)

Return one scalar color value per rendered face for direct Makie mesh coloring.
"""
function light_face_values(
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    geometry = _scene_geometry_for_interception(scene, models, options)
    _light_color_values(geometry, data, color; timestep=timestep, interpolate=false, fill_value=fill_value)
end

"""
    light_vertex_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)

Return one scalar color value per rendered vertex, obtained by averaging the
values of adjacent faces.
"""
function light_vertex_values(
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    geometry = _scene_geometry_for_interception(scene, models, options)
    _light_color_values(geometry, data, color; timestep=timestep, interpolate=true, fill_value=fill_value)
end

"""
    lightplot

Makie plotting entry point provided by the `Makie` package extension.
Load `Makie` or `CairoMakie` before calling it.
"""
function lightplot end

"""
    lightplot!

In-place Makie plotting entry point provided by the `Makie` package extension.
Load `Makie` or `CairoMakie` before calling it.
"""
function lightplot! end
