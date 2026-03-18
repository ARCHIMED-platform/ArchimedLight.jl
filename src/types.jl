import OrderedCollections: OrderedDict

abstract type InterceptionBackend end

struct RasterCPUBackend <: InterceptionBackend end

abstract type ScatteringBackend end

struct RaycastScatteringBackend <: ScatteringBackend end

struct LinksScatteringBackend <: ScatteringBackend end

mutable struct LightConfigPaths
    config::String
    scene::String
    meteo::String
    models::OrderedDict{String,String}
    base_dir::String
end

mutable struct OpticalProperties
    par::Float64
    nir::Float64
    extras::OrderedDict{String,Any}
end

mutable struct InterceptionModelConfig
    use::Union{Nothing,String}
    model::String
    transparency::Float64
    optical_properties::Union{Nothing,OpticalProperties}
    variants::OrderedDict{String,OrderedDict{String,Any}}
    extras::OrderedDict{String,Any}
end

mutable struct LightEmitterConfig
    model::String
    radiance::Float64
    gamma::OpticalProperties
    extras::OrderedDict{String,Any}
end

mutable struct TypeModelConfig
    interception::Union{Nothing,InterceptionModelConfig}
    light_emitter::Union{Nothing,LightEmitterConfig}
    plot_paving::Int
    extras::OrderedDict{String,Any}
end

mutable struct GroupModelConfig
    group::String
    types::OrderedDict{String,TypeModelConfig}
    extras::OrderedDict{String,Any}
end

mutable struct LightGeneralConfig
    all_in_turtle::Bool
    turtle_sectors::Int
    pixel_size::Float64
    area_ratio::Bool
    scattering::Bool
    scattering_max_iter::Int
    scattering_stop_ratio::Float64
    scattering_coeff_par::Float64
    scattering_coeff_nir::Float64
    cache_radiation::Bool
    cache_pixel_table::Bool
    toricity::Bool
    radiation_timestep_minutes::Float64
    nir_interception::Bool
    nir_scattering::Bool
    java_logged_turtle_dirs::Bool
    meteo_range::Union{Nothing,String}
    debug::Bool
    log_debug::Bool
    debug_drop_leading_hit::Union{Nothing,NamedTuple{(:node_id, :x, :y),Tuple{Int,Int,Int}}}
end

mutable struct LightOutputsConfig
    output_directory::String
    simulation_directory::String
    write_summary::Bool
    export_ops::Any
    component_variables::OrderedDict{String,Bool}
    scene_variables::OrderedDict{String,Bool}
    opf_variables::OrderedDict{String,Bool}
    opf_overwrite_variables::Bool
end

mutable struct LightConfig
    source_path::String
    paths::LightConfigPaths
    general::LightGeneralConfig
    outputs::LightOutputsConfig
    models::OrderedDict{String,GroupModelConfig}
    extras::OrderedDict{String,Any}
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

struct SceneGeometry{MTG,Mesh,T}
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
