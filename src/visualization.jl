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

function _missing_render_geometry_error()
    error(
        "This light result does not carry render geometry. Use a LightStepResult returned by " *
        "`run_light_step`/`run_light_series`, or construct one with `render_geometry=`.",
    )
end

"""
    light_render_geometry(scene, models, options)
    light_render_geometry(step)
    light_render_geometry(steps)

Return the render-ready geometry associated with a light simulation.

For a scene, this applies the same model-dependent filtering as the solver. For
results returned by `run_light_step` or `run_light_series`, it returns the
stored render geometry directly.
"""
function light_render_geometry(
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions,
)
    _light_render_geometry(_scene_geometry_for_interception(scene, models, options))
end

function light_render_geometry(step::LightStepResult)
    isnothing(step.render_geometry) && _missing_render_geometry_error()
    step.render_geometry
end

function light_render_geometry(steps::AbstractVector{<:LightStepResult})
    isempty(steps) && error("Cannot infer render geometry from an empty step series.")
    light_render_geometry(first(steps))
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

_is_light_series_data(::AbstractVector{<:LightStepResult}) = true
_is_light_series_data(::AbstractVector{<:AbstractDict{<:Integer,<:Real}}) = true
_is_light_series_data(data) = false

function _finite_extrema(values)
    lo = Inf
    hi = -Inf
    found = false
    for value in values
        v = Float64(value)
        isfinite(v) || continue
        lo = min(lo, v)
        hi = max(hi, v)
        found = true
    end
    return found ? (lo, hi) : nothing
end

_merge_extrema(a::Nothing, b::Nothing) = nothing
_merge_extrema(a::Tuple{Float64,Float64}, b::Nothing) = a
_merge_extrema(a::Nothing, b::Tuple{Float64,Float64}) = b
_merge_extrema(a::Tuple{Float64,Float64}, b::Tuple{Float64,Float64}) = (min(a[1], b[1]), max(a[2], b[2]))

function _light_color_extrema(
    data::LightStepResult,
    color::Symbol;
    timestep::Integer=1,
)
    timestep == 1 || _series_index_error(timestep, 1)
    _finite_extrema(values(light_metric_values(data, color)))
end

function _light_color_extrema(
    data::AbstractVector{<:LightStepResult},
    color::Symbol;
    timestep::Integer=1,
)
    acc = nothing
    for step in data
        acc = _merge_extrema(acc, _finite_extrema(values(light_metric_values(step, color))))
    end
    return acc
end

function _light_color_extrema(
    data::AbstractDict{<:Integer,<:Real},
    ::Symbol;
    timestep::Integer=1,
)
    timestep == 1 || _series_index_error(timestep, 1)
    _finite_extrema(values(data))
end

function _light_color_extrema(
    data::AbstractVector{<:AbstractDict{<:Integer,<:Real}},
    ::Symbol;
    timestep::Integer=1,
)
    acc = nothing
    for metric_map in data
        acc = _merge_extrema(acc, _finite_extrema(values(metric_map)))
    end
    return acc
end

function _light_color_extrema(
    data,
    color::AbstractVector{<:Real};
    timestep::Integer=1,
)
    timestep == 1 || _series_index_error(timestep, 1)
    _finite_extrema(color)
end

function _automatic_light_colorrange(
    data,
    color;
    fill_value::Float64=NaN,
)
    _is_light_series_data(data) || return nothing
    return _light_color_extrema(data, color; timestep=1)
end

function _geometry_for_values(
    data::LightStepResult,
    ::Integer=1,
)
    light_render_geometry(data)
end

function _geometry_for_values(
    data::AbstractVector{<:LightStepResult},
    timestep::Integer=1,
)
    1 <= timestep <= length(data) || _series_index_error(timestep, length(data))
    light_render_geometry(data[timestep])
end

"""
    light_face_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
    light_face_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)

Return one scalar color value per rendered face for direct Makie mesh coloring.

When `data` is a `LightStepResult` or a series of steps returned by
`run_light_step` / `run_light_series`, the stored render geometry is used and no
scene or model inputs are needed.
"""
function light_face_values(
    geometry::LightRenderGeometry,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    _light_color_values(geometry, data, color; timestep=timestep, interpolate=false, fill_value=fill_value)
end

function light_face_values(
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    light_face_values(light_render_geometry(scene, models, options), data; color=color, timestep=timestep, fill_value=fill_value)
end

function light_face_values(
    data::Union{LightStepResult,AbstractVector{<:LightStepResult}};
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    geometry = _geometry_for_values(data, timestep)
    light_face_values(geometry, data; color=color, timestep=timestep, fill_value=fill_value)
end

"""
    light_vertex_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
    light_vertex_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)

Return one scalar color value per rendered vertex, obtained by averaging the
values of adjacent faces.

When `data` is a `LightStepResult` or a series of steps returned by
`run_light_step` / `run_light_series`, the stored render geometry is used and no
scene or model inputs are needed.
"""
function light_vertex_values(
    geometry::LightRenderGeometry,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    _light_color_values(geometry, data, color; timestep=timestep, interpolate=true, fill_value=fill_value)
end

function light_vertex_values(
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions,
    data;
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    light_vertex_values(light_render_geometry(scene, models, options), data; color=color, timestep=timestep, fill_value=fill_value)
end

function light_vertex_values(
    data::Union{LightStepResult,AbstractVector{<:LightStepResult}};
    color::Union{Symbol,AbstractDict,AbstractVector}=:incident_par_flux,
    timestep::Integer=1,
    fill_value::Float64=NaN,
)
    geometry = _geometry_for_values(data, timestep)
    light_vertex_values(geometry, data; color=color, timestep=timestep, fill_value=fill_value)
end

"""
    lightplot

Makie plotting entry point provided by the `Makie` package extension.
Load `Makie` or `CairoMakie` before calling it. The preferred API is
`lightplot(step; ...)` or `lightplot(steps; ...)`.
"""
function lightplot end

"""
    lightplot!

In-place Makie plotting entry point provided by the `Makie` package extension.
Load `Makie` or `CairoMakie` before calling it. The preferred API is
`lightplot!(axis, step; ...)` or `lightplot!(axis, steps; ...)`.
"""
function lightplot! end
