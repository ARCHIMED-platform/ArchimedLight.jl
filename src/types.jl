import OrderedCollections: OrderedDict

abstract type InterceptionBackend end

struct RasterCPUBackend <: InterceptionBackend end

abstract type ScatteringBackend end

struct RaycastScatteringBackend <: ScatteringBackend end

struct LinksScatteringBackend <: ScatteringBackend end

mutable struct OpticalProperties
    par::Float64
    nir::Float64
    extras::OrderedDict{String,Any}
end

OpticalProperties(par::Real=0.0, nir::Real=0.0) = OpticalProperties(Float64(par), Float64(nir), OrderedDict{String,Any}())

mutable struct InterceptionModel
    use::Union{Nothing,String}
    model::String
    sensor::Bool
    transparency::Float64
    optical_properties::Union{Nothing,OpticalProperties}
    variants::OrderedDict{String,OrderedDict{String,Any}}
    extras::OrderedDict{String,Any}
end

function InterceptionModel(;
    model::AbstractString="Translucent",
    sensor::Bool=false,
    transparency::Real=0.0,
    optical_properties::Union{Nothing,OpticalProperties}=nothing,
    use::Union{Nothing,AbstractString}=nothing,
    variants::OrderedDict{String,OrderedDict{String,Any}}=OrderedDict{String,OrderedDict{String,Any}}(),
    extras::OrderedDict{String,Any}=OrderedDict{String,Any}(),
)
    InterceptionModel(
        use === nothing ? nothing : String(use),
        String(model),
        sensor,
        Float64(transparency),
        optical_properties,
        variants,
        extras,
    )
end

mutable struct EmitterModel
    model::String
    radiance::Float64
    gamma::OpticalProperties
    extras::OrderedDict{String,Any}
end

function EmitterModel(;
    model::AbstractString="LambertianEmitter",
    radiance::Real=0.0,
    gamma::OpticalProperties=OpticalProperties(0.48, 0.52),
    extras::OrderedDict{String,Any}=OrderedDict{String,Any}(),
)
    EmitterModel(String(model), Float64(radiance), gamma, extras)
end

mutable struct TypeModel
    interception::Union{Nothing,InterceptionModel}
    light_emitter::Union{Nothing,EmitterModel}
    extras::OrderedDict{String,Any}
end

function TypeModel(;
    interception::Union{Nothing,InterceptionModel}=nothing,
    light_emitter::Union{Nothing,EmitterModel}=nothing,
    extras::OrderedDict{String,Any}=OrderedDict{String,Any}(),
)
    TypeModel(interception, light_emitter, extras)
end

mutable struct GroupModel
    group::String
    types::OrderedDict{String,TypeModel}
    extras::OrderedDict{String,Any}
end

function GroupModel(
    group::AbstractString;
    types::OrderedDict{String,TypeModel}=OrderedDict{String,TypeModel}(),
    extras::OrderedDict{String,Any}=OrderedDict{String,Any}(),
)
    GroupModel(String(group), types, extras)
end

struct LightModels
    groups::OrderedDict{String,GroupModel}
end

LightModels() = LightModels(OrderedDict{String,GroupModel}())

function LightModels(groups::AbstractVector{<:GroupModel})
    ordered = OrderedDict{String,GroupModel}()
    for group in groups
        ordered[group.group] = group
    end
    LightModels(ordered)
end

Base.values(models::LightModels) = values(models.groups)
Base.keys(models::LightModels) = keys(models.groups)
Base.haskey(models::LightModels, key) = haskey(models.groups, key)
Base.getindex(models::LightModels, key) = models.groups[key]

"""
    LightOptions

Runtime controls for interception, scattering, and caching.

`pixel_hit_stack_mode` selects how per-pixel hit stacks are stored during raster
projection:
- `"auto"`: current default optimized path
- `"small"`: force inline small stacks with spillover allocation
- `"vector"`: force the legacy `Vector` stack representation
"""
Base.@kwdef struct LightOptions
    all_in_turtle::Bool = false
    turtle_sectors::Int = 46
    pixel_size::Float64 = 0.25 / 100.0
    area_ratio::Bool = true
    scattering::Bool = true
    scattering_max_iter::Int = 20
    scattering_stop_ratio::Float64 = 0.01
    scattering_coeff_par::Float64 = 0.15
    scattering_coeff_nir::Float64 = 0.30
    cache_radiation::Bool = false
    cache_pixel_table::Bool = false
    pixel_hit_stack_mode::String = "auto"
    toricity::Bool = true
    radiation_timestep_minutes::Float64 = 15.0
    nir_interception::Bool = true
    nir_scattering::Bool = true
    java_logged_turtle_dirs::Bool = false
    meteo_range::Union{Nothing,String} = nothing
    debug::Bool = false
    log_debug::Bool = false
    debug_drop_leading_hit::Union{Nothing,NamedTuple{(:node_id, :x, :y),Tuple{Int,Int,Int}}} = nothing
end

function LightOptions(old::LightOptions; kwargs...)
    params = Dict{Symbol,Any}(
        :all_in_turtle => old.all_in_turtle,
        :turtle_sectors => old.turtle_sectors,
        :pixel_size => old.pixel_size,
        :area_ratio => old.area_ratio,
        :scattering => old.scattering,
        :scattering_max_iter => old.scattering_max_iter,
        :scattering_stop_ratio => old.scattering_stop_ratio,
        :scattering_coeff_par => old.scattering_coeff_par,
        :scattering_coeff_nir => old.scattering_coeff_nir,
        :cache_radiation => old.cache_radiation,
        :cache_pixel_table => old.cache_pixel_table,
        :pixel_hit_stack_mode => old.pixel_hit_stack_mode,
        :toricity => old.toricity,
        :radiation_timestep_minutes => old.radiation_timestep_minutes,
        :nir_interception => old.nir_interception,
        :nir_scattering => old.nir_scattering,
        :java_logged_turtle_dirs => old.java_logged_turtle_dirs,
        :meteo_range => old.meteo_range,
        :debug => old.debug,
        :log_debug => old.log_debug,
        :debug_drop_leading_hit => old.debug_drop_leading_hit,
    )
    for (k, v) in kwargs
        params[k] = v
    end
    return LightOptions(; params...)
end

struct MeteoTable
    rows::Vector{NamedTuple}
    metadata::NamedTuple
end

struct SceneNodeData{T}
    area::T
    barycenter::NTuple{3,T}
    group::String
    type::String
    source_topology_id::Int
    object_id::Int
end

mutable struct SceneGeometry{MTG,Mesh,T}
    mtg::MTG
    merged_mesh::Mesh
    face2node::Vector{Int}
    nodes::Dict{Int,SceneNodeData{T}}
    source_path::String
    scene_xy_bounds::Union{Nothing,NTuple{4,T}}
end

"""
    scene_node(scene, node_id)

Return the [`SceneNodeData`](@ref) entry for `node_id`, or `nothing` when the
node is absent from `scene`.
"""
scene_node(scene::SceneGeometry, node_id::Integer) = get(scene.nodes, Int(node_id), nothing)

"""
    scene_node_ids(scene)

Return the sorted geometry node ids present in `scene`.
"""
scene_node_ids(scene::SceneGeometry) = sort!(collect(keys(scene.nodes)))

"""
    node_areas(scene)

Return a `Dict{Int,Float64}` mapping each geometry node id to its component area
in the prepared scene.
"""
node_areas(scene::SceneGeometry) = Dict(nid => node.area for (nid, node) in scene.nodes)

"""
    node_barycenters(scene)

Return a `Dict` mapping each geometry node id to its `(x, y, z)` barycenter in
scene coordinates.
"""
node_barycenters(scene::SceneGeometry) = Dict(nid => node.barycenter for (nid, node) in scene.nodes)

function _scene_node_field(scene::SceneGeometry, node_id::Integer, field::Symbol, default)
    node = scene_node(scene, node_id)
    node === nothing ? default : getfield(node, field)
end

_scene_area(scene::SceneGeometry, node_id::Integer, default=0.0) = _scene_node_field(scene, node_id, :area, default)
_scene_barycenter(scene::SceneGeometry, node_id::Integer, default=(NaN, NaN, NaN)) =
    _scene_node_field(scene, node_id, :barycenter, default)
_scene_group(scene::SceneGeometry, node_id::Integer, default="") = _scene_node_field(scene, node_id, :group, default)
_scene_type(scene::SceneGeometry, node_id::Integer, default="") = _scene_node_field(scene, node_id, :type, default)
_scene_source_topology_id(scene::SceneGeometry, node_id::Integer, default=Int(node_id)) =
    _scene_node_field(scene, node_id, :source_topology_id, default)
_scene_object_id(scene::SceneGeometry, node_id::Integer, default=-1) =
    _scene_node_field(scene, node_id, :object_id, default)

struct SkyState
    sun_azimuth_deg::Float64
    sun_elevation_deg::Float64
    ri_sw_f::Float64
    ri_par_f::Float64
    ri_nir_f::Float64
    direct_fraction::Float64
    diffuse_fraction::Float64
end

SkyState(
    sun_azimuth_deg::Float64,
    sun_elevation_deg::Float64,
    ri_par_f::Float64,
    ri_nir_f::Float64,
    direct_fraction::Float64,
    diffuse_fraction::Float64,
) = SkyState(
    sun_azimuth_deg,
    sun_elevation_deg,
    ri_par_f + ri_nir_f,
    ri_par_f,
    ri_nir_f,
    direct_fraction,
    diffuse_fraction,
)

struct TurtleSector
    id::Int
    direction
    weight::Float64
    source::Symbol
end

struct TurtleGrid
    sectors::Vector{TurtleSector}
end

struct DirectionalFluxes
    sector_ids::Vector{Int}
    par::Vector{Float64}
    nir::Vector{Float64}
end

struct SpectralNodeValues
    par::Dict{Int,Float64}
    nir::Dict{Int,Float64}
end

struct InitialTotalSpectralNodeValues
    initial::SpectralNodeValues
    total::SpectralNodeValues
end

struct FirstOrderResult
    projected_area_per_node::Dict{Int,Float64}
    incident_power::SpectralNodeValues
    hits_per_node::Dict{Int,Int}
end

struct ScatteringResult
    added_power::SpectralNodeValues
    iterations::Int
    converged::Bool
end

"""
    ScatteringPairCounts

Compact transfer-edge storage for scattering graphs.

The hot topology builder accumulates counts with packed integer keys, then materializes this
container once so downstream code can still iterate edges as `((to, from), count)` pairs
without keeping tuple-key dictionaries in the graph.
"""
struct ScatteringPairCounts
    to_nodes::Vector{Int}
    from_nodes::Vector{Int}
    counts::Vector{Int}
end

Base.length(pair_counts::ScatteringPairCounts) = length(pair_counts.counts)
Base.isempty(pair_counts::ScatteringPairCounts) = isempty(pair_counts.counts)
Base.eltype(::Type{ScatteringPairCounts}) = Pair{Tuple{Int,Int},Int}

function Base.iterate(pair_counts::ScatteringPairCounts, state::Int=1)
    state > length(pair_counts.counts) && return nothing
    item = ((pair_counts.to_nodes[state], pair_counts.from_nodes[state]), pair_counts.counts[state])
    return item, state + 1
end

function ScatteringPairCounts(pair_counts::Dict{Tuple{Int,Int},Int})
    to_nodes = Int[]
    from_nodes = Int[]
    counts = Int[]
    sizehint!(to_nodes, length(pair_counts))
    sizehint!(from_nodes, length(pair_counts))
    sizehint!(counts, length(pair_counts))
    for ((to, from), count) in pair_counts
        push!(to_nodes, to)
        push!(from_nodes, from)
        push!(counts, count)
    end
    return ScatteringPairCounts(to_nodes, from_nodes, counts)
end

struct ScatteringTransferGraph
    pair_counts::ScatteringPairCounts
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
    coeff_par_by_node::Dict{Int,Float64}
    coeff_nir_by_node::Dict{Int,Float64}
    default_coeff_par::Float64
    default_coeff_nir::Float64
end

function ScatteringTransferGraph(
    pair_counts::Dict{Tuple{Int,Int},Int},
    all_hits::Dict{Int,Int},
    node_ids::Vector{Int},
    node_group::Dict{Int,String},
    node_type::Dict{Int,String},
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}},
    coeff_par_by_node::Dict{Int,Float64},
    coeff_nir_by_node::Dict{Int,Float64},
    default_coeff_par::Float64,
    default_coeff_nir::Float64,
)
    return ScatteringTransferGraph(
        ScatteringPairCounts(pair_counts),
        all_hits,
        node_ids,
        node_group,
        node_type,
        group_type_coeffs,
        coeff_par_by_node,
        coeff_nir_by_node,
        default_coeff_par,
        default_coeff_nir,
    )
end

struct LightBudget
    incident_flux::InitialTotalSpectralNodeValues
    incident_energy::InitialTotalSpectralNodeValues
    absorbed_flux::InitialTotalSpectralNodeValues
    absorbed_energy::InitialTotalSpectralNodeValues
    extra_initial_energy_per_band::Dict{String,Dict{Int,Float64}}
    extra_energy_per_band::Dict{String,Dict{Int,Float64}}
end

struct LightStepResult
    sky::SkyState
    turtle::TurtleGrid
    fluxes::DirectionalFluxes
    first_order::FirstOrderResult
    scattering::Union{Nothing,ScatteringResult}
    budget::LightBudget
    extra_band_irradiance::Dict{String,Float64}
end

function _format_decimal(value::Real; digits::Int=3)
    x = round(Float64(value); digits=digits)
    s = string(x)
    occursin('e', s) && return s
    occursin('E', s) && return s
    s = replace(s, r"0+$" => "")
    s = replace(s, r"\.$" => "")
    return s
end

function _format_scaled_quantity(value::Real, unit::AbstractString)
    x = Float64(value)
    ax = abs(x)
    if unit == "J"
        if ax >= 1.0e6
            return _format_decimal(x / 1.0e6) * " MJ"
        elseif ax >= 1.0e3
            return _format_decimal(x / 1.0e3) * " kJ"
        end
    elseif unit == "W m^-2"
        if ax >= 1.0e3
            return _format_decimal(x / 1.0e3) * " kW m^-2"
        end
    end
    return _format_decimal(x) * " " * unit
end

@inline _format_degrees(value::Real) = _format_decimal(value; digits=1) * "°"
@inline _format_percent(value::Real) = _format_decimal(100.0 * Float64(value); digits=1) * "%"

function _print_light_step_rule(io::IO)
    printstyled(io, "  ------------------------------------------------------------"; color=:light_black)
    println(io)
end

function _print_light_step_label(io::IO, label::AbstractString)
    printstyled(io, label; color=:light_blue, bold=true)
end

function _print_light_step_band(io::IO, band::AbstractString)
    color = uppercase(band) == "PAR" ? :green : uppercase(band) == "NIR" ? :yellow : :magenta
    printstyled(io, band; color=color, bold=true)
end

function _print_light_step_row(io::IO, key::AbstractString, value::AbstractString)
    print(io, "  ")
    _print_light_step_label(io, rpad(key, 11))
    println(io, value)
end

function _print_light_step_energy_row(io::IO, key::AbstractString, par_value::AbstractString, nir_value::AbstractString)
    print(io, "  ")
    _print_light_step_label(io, rpad(key, 11))
    _print_light_step_band(io, "PAR")
    print(io, " ", par_value, "  |  ")
    _print_light_step_band(io, "NIR")
    println(io, " ", nir_value)
end

function _light_step_energy_totals(step::LightStepResult)
    budget = step.budget
    incident_par_total = sum(values(budget.incident_energy.total.par))
    incident_nir_total = sum(values(budget.incident_energy.total.nir))
    incident_par_initial = sum(values(budget.incident_energy.initial.par))
    incident_nir_initial = sum(values(budget.incident_energy.initial.nir))
    absorbed_par_total = sum(values(budget.absorbed_energy.total.par))
    absorbed_nir_total = sum(values(budget.absorbed_energy.total.nir))
    (
        incident_par_total=incident_par_total,
        incident_nir_total=incident_nir_total,
        incident_par_added=incident_par_total - incident_par_initial,
        incident_nir_added=incident_nir_total - incident_nir_initial,
        absorbed_par_total=absorbed_par_total,
        absorbed_nir_total=absorbed_nir_total,
    )
end

function _light_step_sector_counts(step::LightStepResult)
    sky_count = count(sector -> sector.source == :sky, step.turtle.sectors)
    sun_count = count(sector -> sector.source == :sun, step.turtle.sectors)
    other_count = length(step.turtle.sectors) - sky_count - sun_count
    return (sky=sky_count, sun=sun_count, other=other_count)
end

function _light_step_extra_bands_summary(step::LightStepResult)
    isempty(step.extra_band_irradiance) && return ""
    parts = String[]
    for band in sort!(collect(keys(step.extra_band_irradiance)))
        push!(parts, uppercase(band) * " " * _format_scaled_quantity(step.extra_band_irradiance[band], "W m^-2"))
    end
    return join(parts, " | ")
end

function Base.show(io::IO, step::LightStepResult)
    totals = _light_step_energy_totals(step)
    sector_counts = _light_step_sector_counts(step)
    scattering_summary =
        if step.scattering === nothing
            "off"
        else
            iter = step.scattering.iterations
            status = step.scattering.converged ? "converged" : "maxiter"
            "$(iter) iters, $(status)"
        end
    print(
        io,
        "LightStepResult(",
        "PAR=", _format_scaled_quantity(totals.incident_par_total, "J"),
        ", NIR=", _format_scaled_quantity(totals.incident_nir_total, "J"),
        ", sky=", _format_scaled_quantity(step.sky.ri_par_f, "W m^-2"), " PAR",
        " / ", _format_scaled_quantity(step.sky.ri_nir_f, "W m^-2"), " NIR",
        ", sectors=", length(step.turtle.sectors),
        " [sky=", sector_counts.sky, ", sun=", sector_counts.sun,
        sector_counts.other > 0 ? ", other=$(sector_counts.other)" : "",
        "]",
        ", scattering=", scattering_summary,
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", step::LightStepResult)
    if get(io, :compact, false)
        show(io, step)
        return
    end
    totals = _light_step_energy_totals(step)
    sector_counts = _light_step_sector_counts(step)
    printstyled(io, "LightStepResult"; color=:cyan, bold=true)
    println(io)
    _print_light_step_rule(io)
    _print_light_step_energy_row(
        io,
        "sky",
        _format_scaled_quantity(step.sky.ri_par_f, "W m^-2"),
        _format_scaled_quantity(step.sky.ri_nir_f, "W m^-2"),
    )
    _print_light_step_row(
        io,
        "sun",
        "azimuth " * _format_degrees(step.sky.sun_azimuth_deg) *
        "  |  elevation " * _format_degrees(step.sky.sun_elevation_deg),
    )
    _print_light_step_row(
        io,
        "mix",
        "direct " * _format_percent(step.sky.direct_fraction) *
        "  |  diffuse " * _format_percent(step.sky.diffuse_fraction) *
        "  |  SW " * _format_scaled_quantity(step.sky.ri_sw_f, "W m^-2"),
    )
    _print_light_step_row(
        io,
        "turtle",
        string(length(step.turtle.sectors), " sectors (sky=", sector_counts.sky, ", sun=", sector_counts.sun,
            sector_counts.other > 0 ? ", other=$(sector_counts.other)" : "", ")"),
    )
    _print_light_step_rule(io)
    _print_light_step_energy_row(
        io,
        "incident",
        _format_scaled_quantity(totals.incident_par_total, "J"),
        _format_scaled_quantity(totals.incident_nir_total, "J"),
    )
    _print_light_step_energy_row(
        io,
        "absorbed",
        _format_scaled_quantity(totals.absorbed_par_total, "J"),
        _format_scaled_quantity(totals.absorbed_nir_total, "J"),
    )

    if step.scattering === nothing
        _print_light_step_row(io, "scattering", "off  |  iterations 0  |  converged n/a")
    else
        print(io, "  ")
        _print_light_step_label(io, rpad("scattering", 11))
        printstyled(io, "on"; color=:magenta, bold=true)
        print(io, "  |  iterations ", step.scattering.iterations)
        print(io, "  |  converged ")
        printstyled(io, string(step.scattering.converged); color=(step.scattering.converged ? :green : :yellow), bold=true)
        println(io)
        _print_light_step_energy_row(
            io,
            "added",
            _format_scaled_quantity(totals.incident_par_added, "J"),
            _format_scaled_quantity(totals.incident_nir_added, "J"),
        )
    end
    extra_summary = _light_step_extra_bands_summary(step)
    if !isempty(extra_summary)
        _print_light_step_row(io, "extra bands", extra_summary)
    end
    _print_light_step_rule(io)
end
