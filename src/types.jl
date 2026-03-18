import OrderedCollections: OrderedDict

abstract type InterceptionBackend end

struct RasterCPUBackend <: InterceptionBackend end

abstract type ScatteringBackend end

struct RaycastScatteringBackend <: ScatteringBackend end

struct LinksScatteringBackend <: ScatteringBackend end

mutable struct LightConfigSourceFiles
    config::String
    scene::String
    meteo::String
    models::OrderedDict{String,String}
    constants::Union{Nothing,String}
    base_dir::String
end

mutable struct LightOutputVariables
    component::OrderedDict{String,Bool}
    scene::OrderedDict{String,Bool}
    opf::OrderedDict{String,Bool}
end

mutable struct LightOutputsConfig
    directory::String
    simulation_directory::String
    write_summary::Bool
    export_ops
    opf_overwrite_variables::Bool
    variables::LightOutputVariables
end

mutable struct LightConfig
    source_path::String
    source_files::LightConfigSourceFiles
    general::OrderedDict{String,Any}
    models::OrderedDict{String,OrderedDict{String,Any}}
    outputs::LightOutputsConfig
    constants::OrderedDict{String,Any}
end

struct MeteoTable
    rows::Vector{NamedTuple}
    metadata::NamedTuple
end

struct SceneGeometry{MTG,MeshT,T<:Real}
    mtg::MTG
    merged_mesh::MeshT
    face2node::Vector{Int}
    total_area_per_node::Dict{Int,T}
    barycenter_per_node::Dict{Int,NTuple{3,T}}
    source_topology_id_per_node::Dict{Int,Int}
    object_id_per_node::Dict{Int,Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    source_path::String
    scene_xy_bounds::Union{Nothing,NTuple{4,T}}
end

struct SkyState{T<:Real}
    sun_azimuth_deg::T
    sun_elevation_deg::T
    ri_sw_f::T
    ri_par_f::T
    ri_nir_f::T
    direct_fraction::T
    diffuse_fraction::T
end

SkyState(
    sun_azimuth_deg::Real,
    sun_elevation_deg::Real,
    ri_sw_f::Real,
    ri_par_f::Real,
    ri_nir_f::Real,
    direct_fraction::Real,
    diffuse_fraction::Real,
) = begin
    values = promote(sun_azimuth_deg, sun_elevation_deg, ri_sw_f, ri_par_f, ri_nir_f, direct_fraction, diffuse_fraction)
    SkyState(values...)
end

SkyState(
    sun_azimuth_deg::Real,
    sun_elevation_deg::Real,
    ri_par_f::Real,
    ri_nir_f::Real,
    direct_fraction::Real,
    diffuse_fraction::Real,
) = SkyState(
    sun_azimuth_deg,
    sun_elevation_deg,
    ri_par_f + ri_nir_f,
    ri_par_f,
    ri_nir_f,
    direct_fraction,
    diffuse_fraction,
)

struct TurtleSector{DirectionT,WeightT<:Real}
    id::Int
    direction::DirectionT
    weight::WeightT
    source::Symbol
end

struct TurtleGrid{SectorT<:TurtleSector}
    sectors::Vector{SectorT}
end

struct DirectionalFluxes{T<:Real}
    sector_ids::Vector{Int}
    par::Vector{T}
    nir::Vector{T}
end

struct FirstOrderResult{T<:Real}
    projected_area_per_node::Dict{Int,T}
    incident_par_power_per_node::Dict{Int,T}
    incident_nir_power_per_node::Dict{Int,T}
    hits_per_node::Dict{Int,Int}
end

struct ScatteringResult{T<:Real}
    added_par_power_per_node::Dict{Int,T}
    added_nir_power_per_node::Dict{Int,T}
    iterations::Int
    converged::Bool
end

struct OpticalCoefficients{T<:Real}
    par::T
    nir::T
end

struct ScatteringTransferGraph{T<:Real}
    pair_counts::Dict{Tuple{Int,Int},Int}
    all_hits::Dict{Int,Int}
    node_ids::Vector{Int}
    node_group::Dict{Int,String}
    node_type::Dict{Int,String}
    group_type_coeffs::Dict{Tuple{String,String},OpticalCoefficients{T}}
end

struct BandLightBudget{T<:Real}
    initial_flux_per_node::Dict{Int,T}
    flux_per_node::Dict{Int,T}
    initial_energy_per_node::Dict{Int,T}
    energy_per_node::Dict{Int,T}
end

struct SpectralLightBudget{T<:Real}
    par::BandLightBudget{T}
    nir::BandLightBudget{T}
end

struct LightBudget{T<:Real}
    incident::SpectralLightBudget{T}
    absorbed::SpectralLightBudget{T}
    extra_initial_energy_per_band::Dict{Symbol,Dict{Int,T}}
    extra_energy_per_band::Dict{Symbol,Dict{Int,T}}
end

struct LightStepResult{SkyT<:SkyState,TurtleT<:TurtleGrid,FluxT<:DirectionalFluxes,FirstT<:FirstOrderResult,ScatteringT,BudgetT<:LightBudget,ExtraT<:AbstractDict}
    sky::SkyT
    turtle::TurtleT
    fluxes::FluxT
    first_order::FirstT
    scattering::ScatteringT
    budget::BudgetT
    extra_band_irradiance::ExtraT
end
