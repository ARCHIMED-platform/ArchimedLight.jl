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

scene_node(scene::SceneGeometry, node_id::Integer) = get(scene.nodes, Int(node_id), nothing)
scene_node_ids(scene::SceneGeometry) = sort!(collect(keys(scene.nodes)))
node_areas(scene::SceneGeometry) = Dict(nid => node.area for (nid, node) in scene.nodes)
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

struct ScatteringTransferGraph
    pair_counts::Dict{Tuple{Int,Int},Int}
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},Dict{String,Float64}}
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
