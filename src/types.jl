import OrderedCollections: OrderedDict

"""
    InterceptionBackend

Abstract supertype for first-order interception backends.
"""
abstract type InterceptionBackend end

"""
    RasterCPUBackend()

Reference interception backend based on CPU raster projection.
"""
struct RasterCPUBackend <: InterceptionBackend end

"""
    ScatteringBackend

Abstract supertype for multiple-scattering backends.
"""
abstract type ScatteringBackend end

"""
    RaycastScatteringBackend()

Scattering backend that reconstructs transfer topology from directional
ray-visibility stacks.
"""
struct RaycastScatteringBackend <: ScatteringBackend end

"""
    LinksScatteringBackend()

Scattering backend that uses precomputed link-style transfer relationships.
"""
struct LinksScatteringBackend <: ScatteringBackend end

"""
    OpticalProperties(par=0.0, nir=0.0)

Per-waveband scattering coefficients for a component or emitter.

`par` and `nir` are the built-in ARCHIMED bands. Additional coefficients can be
stored in `extras`.

For plant components, these coefficients are typically used as the scattered
fraction in each waveband. Expected values are usually in `[0, 1]`.

Typical values:

- leaves in PAR: often around `0.05-0.25`
- leaves in NIR: often around `0.4-0.9`
- emitters: `gamma` often uses `PAR=0.48`, `NIR=0.52` to split total radiance
  into the built-in bands

Examples:

- `OpticalProperties(0.15, 0.90)` for a strongly NIR-scattering leaf
- `OpticalProperties(0.48, 0.52)` for a neutral PAR/NIR emitter split
"""
mutable struct OpticalProperties
    par::Float64
    nir::Float64
    extras::OrderedDict{String,Any}
end

OpticalProperties(par::Real=0.0, nir::Real=0.0) = OpticalProperties(Float64(par), Float64(nir), OrderedDict{String,Any}())

"""
    InterceptionModel(; model="Translucent", sensor=false, transparency=0.0, optical_properties=nothing, use=nothing, variants=..., extras=...)

Radiative behavior attached to one scene component type.

This stores the interception model name, transparency, optical coefficients, and
optional named variants from historical ARCHIMED model files.

Main fields:

- `model`: historical ARCHIMED interception model name. The common runtime case
  is `"Translucent"`.
- `sensor`: marks the type as a virtual sensor. Virtual sensors receive light
  diagnostics but remain transparent in the transfer logic.
- `transparency`: first-order transmitted fraction. Typical values are in
  `[0, 1]`, with `0.0` for fully intercepting surfaces and larger values for
  partially transmitting components.
- `optical_properties`: waveband-dependent scattering coefficients used by the
  scattering stage.
- `use` and `variants`: preserved support for historical YAML files that define
  several named parameter sets under one `Interception` block.

Typical pattern for ordinary canopy elements:

`InterceptionModel(model="Translucent", transparency=0.0, optical_properties=OpticalProperties(0.15, 0.90))`
"""
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

"""
    EmitterModel(; model="LambertianEmitter", radiance=0.0, gamma=OpticalProperties(0.48, 0.52), extras=...)

Emission model for artificial or diagnostic light sources.

`radiance` is the total emitted radiance-like magnitude attached to the source
type, and `gamma` splits that magnitude across wavebands. In the common PAR/NIR
case, `gamma=OpticalProperties(0.48, 0.52)` means "48% of the emitted energy is
treated as PAR and 52% as NIR".

Most canopy simulations do not need explicit emitters and rely only on sky and
sun forcing, but emitters are useful for artificial lighting setups, synthetic
tests, or debugging scenes.
"""
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

"""
    TypeModel(; interception=nothing, light_emitter=nothing, extras=...)

Model definition for one component type inside a functional group.

This is the level that corresponds to one entry under the `Type:` block in a
historical YAML model file.

Typical contents:

- `interception`: how this component intercepts, transmits, and scatters light
- `light_emitter`: optional emission model for artificial sources
- `extras`: preserved non-light metadata from input files

Examples:

- a leaf type with only interception behavior
- a lamp type with only `light_emitter`
- a diagnostic sensor type with `interception=InterceptionModel(sensor=true, ...)`
"""
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

"""
    GroupModel(group; types=OrderedDict(), extras=...)

Model definition for one functional group, keyed by component type.

`group` should match the functional-group name carried by the prepared scene
geometry, for example `"coffee"` or `"pavement"`.

`types` maps scene type names such as `"Leaf"`, `"Metamer"`, or
`"Cobblestone"` to [`TypeModel`](@ref) values.

Wildcard support:

- an exact type key such as `"Leaf"` applies only to that type
- the wildcard key `"*"` acts as a fallback for any type in the group

This allows compact interactive setups such as:

```julia
GroupModel(
    "coffee";
    types=OrderedDict(
        "*" => TypeModel(
            interception=InterceptionModel(
                model="Translucent",
                optical_properties=OpticalProperties(0.15, 0.90),
            ),
        ),
        "Stem" => TypeModel(
            interception=InterceptionModel(
                model="Translucent",
                transparency=0.2,
                optical_properties=OpticalProperties(0.10, 0.50),
            ),
        ),
    ),
)
```

In that example, all coffee components use the wildcard model except `"Stem"`,
which is overridden explicitly.
"""
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

"""
    LightModels(groups)

Top-level collection of all functional-group models used by one simulation.

This is the object passed to `run_light_step`, `run_light_series`, and the lower
level interception/scattering functions.

`groups` maps functional-group names to [`GroupModel`](@ref) values. The
special group key `"*"` is also supported as a global fallback.

When the solver resolves a scene node `(group, type)`, it matches models with
the following precedence:

1. exact group and exact type
2. exact group and wildcard type `"*"`
3. wildcard group `"*"` and exact type
4. wildcard group `"*"` and wildcard type `"*"`

This makes `LightModels` convenient for:

- precise canopy parameterizations with several plant groups
- synthetic scenes where one default model should cover many node types
- gradual interactive workflows where you start with fallback models and then
  add explicit overrides

Examples:

- `LightModels([GroupModel("coffee"; ...), GroupModel("pavement"; ...)])`
- `LightModels([GroupModel("*"; types=OrderedDict("*" => default_type_model))])`
"""
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

Fields:

- `all_in_turtle`: if `false`, keep the direct beam as a separate sun sector;
  if `true`, redistribute it into the turtle sectors.
- `turtle_sectors`: number of diffuse sky sectors. Common values are `16`, `46`,
  or denser grids for smoother angular resolution.
- `pixel_size`: raster pixel size in meters. Smaller values improve geometric
  fidelity but increase runtime and memory use.
- `area_ratio`: enable the ARCHIMED projected-area correction used for parity
  with the historical implementation.
- `scattering`: enable multiple scattering after first-order interception.
- `scattering_max_iter`: maximum number of scattering iterations.
- `scattering_stop_ratio`: stop scattering when the current scattered energy is
  below this fraction of the initial intercepted energy.
- `scattering_coeff_par`: fallback PAR scattering coefficient when the model
  has no PAR [`OpticalProperties`](@ref).
- `scattering_coeff_nir`: fallback NIR scattering coefficient when the model
  has no NIR [`OpticalProperties`](@ref).
- `cache_radiation`: reuse directional responses across series steps when
  possible.
- `include_sky_fraction`: store the per-node `sky_fraction` map in each
  [`LightStepResult`](@ref). Leave `false` unless downstream code needs it.
  When options are read from a config file, this is enabled by requesting
  `sky_fraction` in `component_variables` or `opf_variables`.
- `cache_pixel_table`: cache raster pixel tables for repeated projections.
- `pixel_hit_stack_mode`: storage mode for per-pixel hit stacks. Supported
  values are `"auto"`, `"small"`, and `"vector"`.
- `toricity`: enable horizontal periodic wrapping of the simulated plot.
- `radiation_timestep_minutes`: internal radiative substep used when a meteo
  row covers a coarser interval.
- `allow_overlapping_meteo_steps`: keep overlapping meteo intervals instead of
  rejecting the series during meteo preparation.
- `nir_interception`: include NIR in directional fluxes and first-order
  interception.
- `nir_scattering`: include NIR in the multiple-scattering stage. This has no
  effect if `nir_interception=false`.
- `java_logged_turtle_dirs`: use the Java-compatibility turtle direction path
  used in parity/debug workflows.
- `meteo_range`: optional historical range selector applied during meteo
  preparation, for example `"2, 5"` or a datetime range.
- `debug`: enable debug-only compatibility hooks.
- `log_debug`: emit additional debug logging where implemented.
- `debug_drop_leading_hit`: optional `(node_id, x, y)` hook used to remove a
  leading raster hit at one pixel for parity debugging.

Typical starting point for simple runs:

`LightOptions(turtle_sectors=46, pixel_size=0.0025, scattering=true)`
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
    include_sky_fraction::Bool = false
    cache_pixel_table::Bool = false
    pixel_hit_stack_mode::String = "auto"
    toricity::Bool = true
    radiation_timestep_minutes::Float64 = 15.0
    allow_overlapping_meteo_steps::Bool = false
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
        :include_sky_fraction => old.include_sky_fraction,
        :cache_pixel_table => old.cache_pixel_table,
        :pixel_hit_stack_mode => old.pixel_hit_stack_mode,
        :toricity => old.toricity,
        :radiation_timestep_minutes => old.radiation_timestep_minutes,
        :allow_overlapping_meteo_steps => old.allow_overlapping_meteo_steps,
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

"""
    MeteoTable(rows, metadata)

Normalized in-memory meteorological forcing table.

`rows` is a vector of named tuples, and `metadata` stores site-level context
such as latitude, altitude, and source file information.

This is the interactive equivalent of a parsed meteorological CSV file.

Typical row fields include:

- time and date information such as `date`, `hour_start`, `hour_end`
- radiative forcing such as `RI_SW_f`, `RI_PAR_f`, `RI_NIR_f`
- optional sun diagnostics such as `sun_azimuth`, `sun_elevation`
- optional partitioning inputs such as `direct_fraction` or `clearness`

The table is usually passed through `prepare_meteo(meteo, options)` before
running a series simulation.
"""
struct MeteoTable
    rows::Vector{NamedTuple}
    metadata::NamedTuple
end

"""
    SceneNodeData

Per-node geometric and semantic information extracted into a prepared
[`SceneGeometry`](@ref).

Fields:

- `area`: component area in scene units after geometry preparation
- `barycenter`: component barycenter in scene coordinates
- `group`: functional group used for model matching
- `type`: type name used for model matching
- `source_topology_id`: original topology identifier from the source MTG or file
- `object_id`: object identifier used to group components belonging to the same
  plant or object when available
"""
struct SceneNodeData{T}
    area::T
    barycenter::NTuple{3,T}
    group::String
    type::String
    source_topology_id::Int
    object_id::Int
end

"""
    SceneGeometry

Prepared scene representation used by the light solver.

It stores the merged mesh, face-to-node map, per-node metadata, source path,
and optional xy bounds used during raster projection.

This is the geometry object returned by `prepare_scene` or `read_scene`. It is
the scene-side input consumed by the solver.

Main fields:

- `mtg`: prepared scene graph
- `merged_mesh`: merged triangle mesh used by the projection engine
- `face2node`: mapping from merged mesh faces back to scene node ids
- `nodes`: per-node metadata and geometry summaries
- `source_path`: origin of the scene when it came from disk
- `scene_xy_bounds`: optional explicit xy limits for raster projection

In most workflows you do not construct this type manually; instead you build or
load a scene and then call `prepare_scene(...)`.
"""
mutable struct SceneGeometry{MTG,Mesh,T}
    mtg::MTG
    merged_mesh::Mesh
    face2node::Vector{Int}
    nodes::Dict{Int,SceneNodeData{T}}
    source_path::String
    scene_xy_bounds::Union{Nothing,NTuple{4,T}}
end

"""
    LightRenderGeometry

Render-ready geometry used to visualize one simulated light result.

This stores the exact face subset and face-to-node mapping used by the light
solver after model-dependent filtering has been applied.
"""
struct LightRenderGeometry
    vertices
    faces
    face2node::Vector{Int}
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

"""
    SkyState(sun_azimuth_deg, sun_elevation_deg, ri_par_f, ri_nir_f, direct_fraction, diffuse_fraction)

Instantaneous radiative forcing state used to build the turtle and directional
fluxes for one light step.

Arguments:

- `sun_azimuth_deg`: sun azimuth in degrees on the horizontal plane. The package
  uses `0°` along `+y`, `90°` along `+x`, `180°` along `-y`, and `270°` along
  `-x`. If your scene uses the common `x=east`, `y=north` convention, this means
  `0°=north`, `90°=east`, `180°=south`, `270°=west`.
- `sun_elevation_deg`: sun elevation in degrees above the horizon. Typical
  daylight values are between `5` and `80`. Use `90` for a sun at zenith; in
  simple examples, `89` is often a good practical choice when you want a nearly
  vertical beam.
- `ri_par_f`: incident PAR irradiance in `W m^-2` on a horizontal plane.
  Typical daylight values are often in the `100-600` range. Clear midday
  conditions are often around `300-500`.
- `ri_nir_f`: incident NIR irradiance in `W m^-2` on a horizontal plane.
  Typical daylight values are often in the `100-700` range, commonly of the
  same order as or slightly larger than PAR.
- `direct_fraction`: fraction of the shortwave forcing treated as direct beam.
  Expected range is `[0, 1]`. Typical values are near `0.8-1.0` for clear-sky
  conditions, around `0.3-0.7` for mixed conditions, and `0.0` for a fully
  diffuse sky.
- `diffuse_fraction`: fraction of the shortwave forcing treated as diffuse sky
  radiation. Expected range is `[0, 1]`. In most use cases, it should satisfy
  `direct_fraction + diffuse_fraction == 1`.

This six-argument constructor stores `ri_sw_f = ri_par_f + ri_nir_f`
automatically.

Examples:

- Nearly zenith, mostly direct forcing:
  `SkyState(180.0, 89.0, 350.0, 250.0, 0.95, 0.05)`
- Lower sun with mixed direct and diffuse light:
  `SkyState(135.0, 35.0, 200.0, 180.0, 0.5, 0.5)`
- Fully diffuse overcast step:
  `SkyState(180.0, 45.0, 120.0, 130.0, 0.0, 1.0)`
"""
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

"""
    TurtleSector

One directional sector of the ARCHIMED turtle, with its direction, weight, and
source (`:sky`, `:sun`, or another source label).
"""
struct TurtleSector
    id::Int
    direction
    weight::Float64
    source::Symbol
end

"""
    TurtleGrid

Collection of [`TurtleSector`](@ref)s used to discretize incoming radiation.
"""
struct TurtleGrid
    sectors::Vector{TurtleSector}
end

"""
    DirectionalFluxes

Per-sector PAR and NIR fluxes aligned with one [`TurtleGrid`](@ref).
"""
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

"""
    FirstOrderResult

Outputs of the first-order interception stage: projected area, incident power,
and hit counts per node.
"""
struct FirstOrderResult
    projected_area_per_node::Dict{Int,Float64}
    incident_power::SpectralNodeValues
    hits_per_node::Dict{Int,Int}
end

"""
    ScatteringResult

Outputs of the multiple-scattering stage: added power per node, iteration
count, and convergence flag.
"""
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

"""
    ScatteringTransferGraph

Compact scene-scale topology used by the scattering solver to move energy
between nodes.
"""
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

"""
    LightBudget

Per-node light budget for one simulation step, storing incident and absorbed
fluxes and energies for PAR, NIR, and optional extra wavebands.
"""
struct LightBudget
    incident_flux::InitialTotalSpectralNodeValues
    incident_energy::InitialTotalSpectralNodeValues
    absorbed_flux::InitialTotalSpectralNodeValues
    absorbed_energy::InitialTotalSpectralNodeValues
    extra_initial_energy_per_band::Dict{String,Dict{Int,Float64}}
    extra_energy_per_band::Dict{String,Dict{Int,Float64}}
end

"""
    LightStepResult

Complete result of one light simulation step, including the sky state, turtle,
directional fluxes, first-order interception, optional scattering, and the
integrated [`LightBudget`](@ref). When requested, the result can also store a
per-node `sky_fraction` map. Results returned by `run_light_step` and
`run_light_series` also carry the render geometry needed by `lightplot`.
"""
struct LightStepResult
    sky::SkyState
    turtle::TurtleGrid
    fluxes::DirectionalFluxes
    first_order::FirstOrderResult
    scattering::Union{Nothing,ScatteringResult}
    budget::LightBudget
    extra_band_irradiance::Dict{String,Float64}
    sky_fraction::Union{Nothing,Dict{Int,Float64}}
    render_geometry::Union{Nothing,LightRenderGeometry}
end

LightStepResult(
    sky::SkyState,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    first_order::FirstOrderResult,
    scattering::Union{Nothing,ScatteringResult},
    budget::LightBudget,
    extra_band_irradiance::Dict{String,Float64},
) = LightStepResult(sky, turtle, fluxes, first_order, scattering, budget, extra_band_irradiance, nothing, nothing)

LightStepResult(
    sky::SkyState,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    first_order::FirstOrderResult,
    scattering::Union{Nothing,ScatteringResult},
    budget::LightBudget,
    extra_band_irradiance::Dict{String,Float64},
    sky_fraction::Union{Nothing,Dict{Int,Float64}},
) = LightStepResult(sky, turtle, fluxes, first_order, scattering, budget, extra_band_irradiance, sky_fraction, nothing)

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
