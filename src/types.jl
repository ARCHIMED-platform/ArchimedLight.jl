struct LightConfig
    scene::String
    meteo::String
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
    raw::Dict{String,Any}
end

struct MeteoTable
    rows::Vector{NamedTuple}
    metadata::NamedTuple
end

struct SceneGeometry
    mtg
    merged_mesh
    face2node::Vector{Int}
    total_area_per_node::Dict{Int,Float64}
    node_group::Dict{Int,String}
    source_path::String
    scene_xy_bounds::Union{Nothing,NTuple{4,Float64}}
end

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

struct FirstOrderResult
    projected_area_per_node::Dict{Int,Float64}
    incident_par_power_per_node::Dict{Int,Float64}
    incident_nir_power_per_node::Dict{Int,Float64}
    hits_per_node::Dict{Int,Int}
end

struct ScatteringResult
    added_par_power_per_node::Dict{Int,Float64}
    added_nir_power_per_node::Dict{Int,Float64}
    iterations::Int
    converged::Bool
end

struct LightBudget
    ri_par_f_per_node::Dict{Int,Float64}
    ri_nir_f_per_node::Dict{Int,Float64}
    ri_par_0_q_per_node::Dict{Int,Float64}
    ri_nir_0_q_per_node::Dict{Int,Float64}
    ri_par_q_per_node::Dict{Int,Float64}
    ri_nir_q_per_node::Dict{Int,Float64}
    extra_0_q_per_band::Dict{String,Dict{Int,Float64}}
    extra_q_per_band::Dict{String,Dict{Int,Float64}}
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
