"""
    VoxelGrid(minimum, maximum, pad)

A regular Cartesian voxel grid. `minimum` and `maximum` are three-coordinate
corners in metres. `pad[i, j, k]` is plant area density in
`m² plant m⁻³`. Julia indices are 1-based; historical `.vox` files are
converted explicitly at the I/O boundary.
"""
struct VoxelGrid{T<:AbstractFloat}
    minimum::NTuple{3,T}
    maximum::NTuple{3,T}
    pad::Array{T,3}
    voxel_size::NTuple{3,T}
end

"""Abstract interface for a surface that bounds voxel transport from below."""
abstract type AbstractVoxelTerrain end

"""Open lower boundary: rays may leave the voxel grid through `:bottom`."""
struct NoVoxelTerrain <: AbstractVoxelTerrain end

"""
    SoilOpticalProperties(; par_reflectance, nir_reflectance)

Lambertian soil reflectance in the operational PAR and NIR bands. Soil is
opaque in the initial terrain model, so absorptance is `1 - reflectance` and
transmittance is zero.
"""
struct SoilOpticalProperties{T<:AbstractFloat}
    par_reflectance::T
    nir_reflectance::T
end

function SoilOpticalProperties(; par_reflectance::Real, nir_reflectance::Real)
    par = Float64(par_reflectance)
    nir = Float64(nir_reflectance)
    isfinite(par) && 0 <= par <= 1 ||
        throw(ArgumentError("soil PAR reflectance must be finite and lie in [0, 1]"))
    isfinite(nir) && 0 <= nir <= 1 ||
        throw(ArgumentError("soil NIR reflectance must be finite and lie in [0, 1]"))
    return SoilOpticalProperties(par, nir)
end

"""The nearest positive terrain intersection along a normalized ray."""
struct TerrainHit{T<:AbstractFloat}
    distance::T
    point::NTuple{3,T}
    normal::NTuple{3,T}
    patch_id::Int
    material_id::Int

    function TerrainHit(
        distance::T,
        point::NTuple{3,T},
        normal::NTuple{3,T},
        patch_id::Int,
        material_id::Int,
    ) where {T<:AbstractFloat}
        isfinite(distance) && distance > 0 ||
            throw(ArgumentError("terrain-hit distance must be finite and positive"))
        all(isfinite, point) || throw(ArgumentError("terrain-hit point must be finite"))
        all(isfinite, normal) || throw(ArgumentError("terrain-hit normal must be finite"))
        magnitude = sqrt(sum(abs2, normal))
        isapprox(magnitude, one(T); atol=128eps(T), rtol=128eps(T)) ||
            throw(ArgumentError("terrain-hit normal must be normalized"))
        normal[3] > 0 || throw(ArgumentError("terrain-hit normal must point from soil toward air"))
        patch_id > 0 || throw(ArgumentError("terrain-hit patch_id must be positive"))
        material_id > 0 || throw(ArgumentError("terrain-hit material_id must be positive"))
        return new{T}(distance, point, normal, patch_id, material_id)
    end
end

function TerrainHit(distance::Real, point, normal, patch_id::Integer, material_id::Integer)
    values = promote(
        Float64(distance),
        (Float64(point[1]), Float64(point[2]), Float64(point[3]))...,
        (Float64(normal[1]), Float64(normal[2]), Float64(normal[3]))...,
    )
    return TerrainHit(
        values[1],
        (values[2], values[3], values[4]),
        (values[5], values[6], values[7]),
        Int(patch_id),
        Int(material_id),
    )
end

function VoxelGrid(minimum, maximum, pad::AbstractArray{<:Real,3})
    length(minimum) == 3 || throw(ArgumentError("minimum must contain three coordinates"))
    length(maximum) == 3 || throw(ArgumentError("maximum must contain three coordinates"))
    all(>(0), size(pad)) || throw(ArgumentError("a voxel grid cannot have an empty dimension"))

    lo = ntuple(i -> Float64(minimum[i]), 3)
    hi = ntuple(i -> Float64(maximum[i]), 3)
    all(isfinite, lo) || throw(ArgumentError("minimum coordinates must be finite"))
    all(isfinite, hi) || throw(ArgumentError("maximum coordinates must be finite"))
    all(i -> hi[i] > lo[i], 1:3) ||
        throw(ArgumentError("each maximum coordinate must be greater than minimum"))

    values = Array{Float64,3}(undef, size(pad))
    for index in eachindex(pad)
        value = Float64(pad[index])
        isfinite(value) || throw(ArgumentError("PAD values must be finite"))
        value >= 0 || throw(ArgumentError("PAD values must be non-negative"))
        values[index] = value
    end

    dimensions = size(values)
    spacing = ntuple(i -> (hi[i] - lo[i]) / dimensions[i], 3)
    return VoxelGrid{Float64}(lo, hi, values, spacing)
end

Base.size(grid::VoxelGrid) = size(grid.pad)
Base.axes(grid::VoxelGrid) = axes(grid.pad)
Base.getindex(grid::VoxelGrid, i::Integer, j::Integer, k::Integer) = grid.pad[i, j, k]

"""Return the voxel volume in `m³`."""
voxel_volume(grid::VoxelGrid) = prod(grid.voxel_size)

"""Return the horizontal area `dx * dy` of one voxel in `m²`."""
voxel_horizontal_area(grid::VoxelGrid) = grid.voxel_size[1] * grid.voxel_size[2]

"""Return the plant area represented by voxel `(i, j, k)` in `m² plant`."""
voxel_plant_area(grid::VoxelGrid, i::Integer, j::Integer, k::Integer) =
    grid.pad[i, j, k] * voxel_volume(grid)

"""
    voxel_linear_index(grid, i, j, k)

Return a 1-based linear identifier matching the historical Java order after
adding one: `((i - 1) * ny + (j - 1)) * nz + k`.
"""
function voxel_linear_index(grid::VoxelGrid, i::Integer, j::Integer, k::Integer)
    nx, ny, nz = size(grid)
    checkbounds(Bool, grid.pad, i, j, k) || throw(BoundsError(grid, (i, j, k)))
    return ((Int(i) - 1) * ny + (Int(j) - 1)) * nz + Int(k)
end

"""Convert a 1-based historical-order identifier to `(i, j, k)`."""
function voxel_cartesian_index(grid::VoxelGrid, index::Integer)
    nx, ny, nz = size(grid)
    1 <= index <= nx * ny * nz || throw(BoundsError(grid, index))
    zero_based = Int(index) - 1
    i = zero_based ÷ (ny * nz) + 1
    remainder = zero_based % (ny * nz)
    j = remainder ÷ nz + 1
    k = remainder % nz + 1
    return (i, j, k)
end

"""One ordered ray segment inside a voxel."""
struct VoxelRaySegment{T<:AbstractFloat}
    index::NTuple{3,Int}
    length::T
    t_enter::T
    t_exit::T
end

"""
An ordered ray path through a voxel grid. `exit_reason` is `:terrain`, `:top`,
`:bottom`, or `:side`. A terrain path carries its nearest `terrain_hit`;
no-terrain paths retain `nothing`. Public first-order rays point downward;
internal scattering rays may leave through either vertical boundary.
"""
struct VoxelRayPath{T<:AbstractFloat}
    segments::Vector{VoxelRaySegment{T}}
    exit_reason::Symbol
    terrain_hit::Union{Nothing,TerrainHit{T}}
end

VoxelRayPath(segments::Vector{VoxelRaySegment{T}}, exit_reason::Symbol) where {T<:AbstractFloat} =
    VoxelRayPath{T}(segments, exit_reason, nothing)

Base.length(path::VoxelRayPath) = length(path.segments)
Base.iterate(path::VoxelRayPath, state...) = iterate(path.segments, state...)

const VoxelOpticalField3D = Union{Float64,Array{Float64,3}}
const VoxelOpticalField2D = Union{Float64,Array{Float64,2}}
"""Operational voxel PAR wavelength interval in nanometres."""
const VOXEL_PAR_WAVELENGTH_NM = (400.0, 700.0)
"""Operational voxel NIR wavelength interval in nanometres."""
const VOXEL_NIR_WAVELENGTH_NM = (700.0, 2500.0)

function _validated_voxel_optical_field(value, dimensions, name::AbstractString)
    if value isa Real
        scalar = Float64(value)
        isfinite(scalar) || throw(ArgumentError("$name must be finite"))
        0.0 <= scalar <= 1.0 || throw(ArgumentError("$name must lie in [0, 1]"))
        return scalar
    end
    value isa AbstractArray ||
        throw(ArgumentError("$name must be a scalar or an array aligned with the voxel grid"))
    size(value) == dimensions ||
        throw(DimensionMismatch("$name has size $(size(value)); expected $dimensions"))
    result = Array{Float64}(undef, dimensions)
    for index in eachindex(value)
        coefficient = Float64(value[index])
        isfinite(coefficient) || throw(ArgumentError("$name values must be finite"))
        0.0 <= coefficient <= 1.0 ||
            throw(ArgumentError("$name values must lie in [0, 1]"))
        result[index] = coefficient
    end
    return result
end

@inline _voxel_optical_value(field::Float64, index) = field
@inline _voxel_optical_value(field::AbstractArray, index) = field[index]

"""
    VoxelOpticalProperties(grid; par_reflectance, par_transmittance,
                           nir_reflectance, nir_transmittance)

Leaf optical properties for voxel multiple scattering. Each coefficient can be
a scalar or a three-dimensional array aligned with `grid.pad`. The
single-scattering albedo in one band is `reflectance + transmittance`; geometric
gaps remain represented by PAD and Beer-Lambert extinction.

For compatibility with historical ARCHIMED material files, pass
`scattering_par` and `scattering_nir` instead. Their total scattering factors
are split equally between reflection and transmission, matching the historical
bilambertian approximation.
"""
struct VoxelOpticalProperties{
    PR<:VoxelOpticalField3D,
    PT<:VoxelOpticalField3D,
    NR<:VoxelOpticalField3D,
    NT<:VoxelOpticalField3D,
}
    grid_size::NTuple{3,Int}
    par_reflectance::PR
    par_transmittance::PT
    nir_reflectance::NR
    nir_transmittance::NT
end

function VoxelOpticalProperties(
    grid::VoxelGrid;
    par_reflectance=nothing,
    par_transmittance=nothing,
    nir_reflectance=nothing,
    nir_transmittance=nothing,
    scattering_par=nothing,
    scattering_nir=nothing,
)
    if scattering_par !== nothing
        (par_reflectance === nothing && par_transmittance === nothing) ||
            throw(ArgumentError("provide either PAR reflectance/transmittance or scattering_par"))
        par_reflectance = Float64(scattering_par) / 2
        par_transmittance = Float64(scattering_par) / 2
    end
    if scattering_nir !== nothing
        (nir_reflectance === nothing && nir_transmittance === nothing) ||
            throw(ArgumentError("provide either NIR reflectance/transmittance or scattering_nir"))
        nir_reflectance = Float64(scattering_nir) / 2
        nir_transmittance = Float64(scattering_nir) / 2
    end
    any(isnothing, (par_reflectance, par_transmittance, nir_reflectance, nir_transmittance)) &&
        throw(ArgumentError("PAR and NIR reflectance and transmittance are required"))

    dimensions = size(grid)
    par_r = _validated_voxel_optical_field(par_reflectance, dimensions, "PAR reflectance")
    par_t = _validated_voxel_optical_field(par_transmittance, dimensions, "PAR transmittance")
    nir_r = _validated_voxel_optical_field(nir_reflectance, dimensions, "NIR reflectance")
    nir_t = _validated_voxel_optical_field(nir_transmittance, dimensions, "NIR transmittance")
    for index in CartesianIndices(grid.pad)
        _voxel_optical_value(par_r, index) + _voxel_optical_value(par_t, index) <= 1 + 16eps(Float64) ||
            throw(ArgumentError("PAR reflectance + transmittance exceeds one at $index"))
        _voxel_optical_value(nir_r, index) + _voxel_optical_value(nir_t, index) <= 1 + 16eps(Float64) ||
            throw(ArgumentError("NIR reflectance + transmittance exceeds one at $index"))
    end
    return VoxelOpticalProperties(dimensions, par_r, par_t, nir_r, nir_t)
end

"""Historical ARCHIMED example optics: total scattering 0.15 PAR and 0.90 NIR."""
voxel_java_parity_optics(grid::VoxelGrid) =
    VoxelOpticalProperties(grid; scattering_par=0.15, scattering_nir=0.90)

"""Documented generic-green-leaf example optics: total scattering 0.17 PAR and 0.85 NIR."""
voxel_generic_green_leaf_optics(grid::VoxelGrid) =
    VoxelOpticalProperties(grid; scattering_par=0.17, scattering_nir=0.85)

function _voxel_band_optical_fields(optics::VoxelOpticalProperties, band::Symbol)
    band == :par && return (optics.par_reflectance, optics.par_transmittance)
    band == :nir && return (optics.nir_reflectance, optics.nir_transmittance)
    throw(ArgumentError("band must be :par or :nir"))
end

"""Return the dimensionless single-scattering albedo `reflectance + transmittance`."""
function voxel_single_scattering_albedo(
    optics::VoxelOpticalProperties,
    band::Symbol,
    index::CartesianIndex{3},
)
    reflectance, transmittance = _voxel_band_optical_fields(optics, band)
    return _voxel_optical_value(reflectance, index) +
           _voxel_optical_value(transmittance, index)
end

function voxel_single_scattering_albedo(
    optics::VoxelOpticalProperties{Float64,Float64,Float64,Float64},
    band::Symbol,
)
    reflectance, transmittance = _voxel_band_optical_fields(optics, band)
    return reflectance + transmittance
end

function voxel_single_scattering_albedo(optics::VoxelOpticalProperties, band::Symbol)
    reflectance, transmittance = _voxel_band_optical_fields(optics, band)
    return [
        _voxel_optical_value(reflectance, index) +
        _voxel_optical_value(transmittance, index)
        for index in CartesianIndices(optics.grid_size)
    ]
end

"""Return foliage absorptance `1 - reflectance - transmittance`."""
voxel_absorptance(optics::VoxelOpticalProperties, band::Symbol, index::CartesianIndex{3}) =
    1 - voxel_single_scattering_albedo(optics, band, index)

voxel_absorptance(
    optics::VoxelOpticalProperties{Float64,Float64,Float64,Float64},
    band::Symbol,
) = 1 - voxel_single_scattering_albedo(optics, band)

voxel_absorptance(optics::VoxelOpticalProperties, band::Symbol) =
    1 .- voxel_single_scattering_albedo(optics, band)

"""
    VoxelGroundOptics(grid; par_reflectance, nir_reflectance)

Optional Lambertian ground reflectance. Coefficients can be scalars or arrays
aligned with the horizontal `(x, y)` grid. Energy not reflected is absorbed by
the ground.
"""
struct VoxelGroundOptics{P<:VoxelOpticalField2D,N<:VoxelOpticalField2D}
    grid_size::NTuple{2,Int}
    par_reflectance::P
    nir_reflectance::N
end

function VoxelGroundOptics(grid::VoxelGrid; par_reflectance, nir_reflectance)
    dimensions = size(grid)[1:2]
    par = _validated_voxel_optical_field(par_reflectance, dimensions, "ground PAR reflectance")
    nir = _validated_voxel_optical_field(nir_reflectance, dimensions, "ground NIR reflectance")
    return VoxelGroundOptics(dimensions, par, nir)
end

"""
    VoxelCPUBackend(; rays_per_voxel=1000, g=0.5,
                      boundary=:periodic, traversal=:dda, lidf=:spherical)

CPU settings for first-order voxel interception. `boundary` accepts
`:periodic`, `:open`, or the historical diagnostic mode `:java_nontoric`.
`traversal` accepts `:dda`, `:reference`, or the compatibility-only
`:java_reference`. Only spherical LIDF is currently supported.
"""
struct VoxelCPUBackend <: InterceptionBackend
    rays_per_voxel::Int
    g::Float64
    boundary::Symbol
    traversal::Symbol
    lidf::Symbol
end

function VoxelCPUBackend(;
    rays_per_voxel::Integer=1000,
    g::Real=0.5,
    boundary::Symbol=:periodic,
    traversal::Symbol=:dda,
    lidf::Symbol=:spherical,
)
    rays_per_voxel > 0 || throw(ArgumentError("rays_per_voxel must be positive"))
    isfinite(g) && g >= 0 || throw(ArgumentError("g must be finite and non-negative"))
    boundary in (:periodic, :open, :java_nontoric) ||
        throw(ArgumentError("boundary must be :periodic, :open, or :java_nontoric"))
    traversal in (:dda, :reference, :java_reference) ||
        throw(ArgumentError("traversal must be :dda, :reference, or :java_reference"))
    lidf == :spherical || throw(ArgumentError("only lidf=:spherical is supported"))
    return VoxelCPUBackend(Int(rays_per_voxel), Float64(g), boundary, traversal, lidf)
end

"""Unit-irradiance response of a voxel grid for one propagation direction."""
struct VoxelDirectionResponse{T<:AbstractFloat}
    interception_fraction::Array{T,3}
    direction::NTuple{3,T}
    incoming_weight::T
    injected_weight::T
    escaped_weight::T
    escaped_top_weight::T
    escaped_side_weight::T
    escaped_bottom_weight::T
    bottom_exit_fraction::Array{T,2}
    terrain_incident_weight::T
    terrain_incident_fraction::Vector{T}
    effective_ray_count::Int
end

"""Reusable directional responses aligned with a `TurtleGrid`."""
struct VoxelResponseCache{T<:AbstractFloat}
    grid_fingerprint::UInt
    terrain_fingerprint::UInt
    backend::VoxelCPUBackend
    directions::Vector{NTuple{3,T}}
    responses::Vector{VoxelDirectionResponse{T}}
end

"""First-order spectral interception values for one voxel light step."""
struct VoxelFirstOrderResult{T<:AbstractFloat}
    par_energy::Array{T,3}
    nir_energy::Array{T,3}
    par_flux::Array{T,3}
    nir_flux::Array{T,3}
    par_incoming_energy::T
    nir_incoming_energy::T
    par_injected_energy::T
    nir_injected_energy::T
    par_escaped_energy::T
    nir_escaped_energy::T
    par_escaped_top_energy::T
    nir_escaped_top_energy::T
    par_escaped_side_energy::T
    nir_escaped_side_energy::T
    par_escaped_bottom_energy::T
    nir_escaped_bottom_energy::T
    par_bottom_energy::Array{T,2}
    nir_bottom_energy::Array{T,2}
    par_terrain_incident_energy::T
    nir_terrain_incident_energy::T
    par_terrain_energy::Vector{T}
    nir_terrain_energy::Vector{T}
end

"""Normalized full-sphere angular quadrature for isotropic voxel scattering."""
struct VoxelScatteringQuadrature{T<:AbstractFloat}
    directions::Vector{NTuple{3,T}}
    weights::Vector{T}
end

"""Geometry-dependent internal-ray paths reused by the matrix-free transport operator."""
struct VoxelScatteringTransportCache{T<:AbstractFloat}
    grid_fingerprint::UInt
    terrain_fingerprint::UInt
    backend::VoxelCPUBackend
    quadrature_fingerprint::UInt
    source_paths::Dict{Tuple{Int,Int},VoxelRayPath{T}}
    ground_paths::Dict{Tuple{Int,Int},VoxelRayPath{T}}
    terrain_paths::Dict{Tuple{Int,Int},VoxelRayPath{T}}
    terrain_direction_weights::Dict{Int,Tuple{Vector{Int},Vector{T}}}
    terrain_patch_count::Int
end

"""One conservative application of foliage and optional soil transport."""
struct VoxelTransportResult{T<:AbstractFloat}
    intercepted_energy::Array{T,3}
    ground_absorbed_energy::T
    ground_reflected_energy::T
    escaped_top_energy::T
    escaped_side_energy::T
    escaped_bottom_energy::T
    unresolved_energy::T
    balance_residual::T
end

"""Complete iterative scattering result for one spectral band."""
struct VoxelBandScatteringResult{T<:AbstractFloat}
    added_intercepted_energy::Array{T,3}
    total_intercepted_energy::Array{T,3}
    absorbed_energy::Array{T,3}
    absorbed_flux::Array{T,3}
    scattered_energy::T
    ground_absorbed_energy::T
    ground_reflected_energy::T
    escaped_top_energy::T
    escaped_side_energy::T
    escaped_bottom_energy::T
    unresolved_energy::T
    iterations::Int
    converged::Bool
    enabled::Bool
    absolute_balance_residual::T
    relative_balance_residual::T
    order_intercepted_energy::Vector{T}
end

"""Soil absorption diagnostic (the stored `ground_*` field name is legacy-compatible)."""
soil_absorbed_energy(result::Union{VoxelTransportResult,VoxelBandScatteringResult}) =
    result.ground_absorbed_energy

"""Soil reflection diagnostic (an internal transfer, not an external source)."""
soil_reflected_energy(result::Union{VoxelTransportResult,VoxelBandScatteringResult}) =
    result.ground_reflected_energy

"""PAR and NIR voxel multiple-scattering results."""
struct VoxelScatteringResult{T<:AbstractFloat}
    par::VoxelBandScatteringResult{T}
    nir::VoxelBandScatteringResult{T}
end

"""Complete voxel result for one meteorological row."""
struct VoxelLightStepResult{T<:AbstractFloat}
    sky::SkyState
    turtle::TurtleGrid
    fluxes::DirectionalFluxes
    first_order::VoxelFirstOrderResult{T}
    scattering::Union{Nothing,VoxelScatteringResult{T}}
    duration_seconds::Float64
end

VoxelLightStepResult(sky, turtle, fluxes, first_order::VoxelFirstOrderResult{T}, duration_seconds) where {T} =
    VoxelLightStepResult(sky, turtle, fluxes, first_order, nothing, Float64(duration_seconds))
