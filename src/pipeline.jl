function _zero_spectral_node_values(node_ids)
    SpectralNodeValues(
        Dict{Int,Float64}(nid => 0.0 for nid in node_ids),
        Dict{Int,Float64}(nid => 0.0 for nid in node_ids),
    )
end

function _zero_initial_total_spectral_node_values(node_ids)
    InitialTotalSpectralNodeValues(
        _zero_spectral_node_values(node_ids),
        _zero_spectral_node_values(node_ids),
    )
end

function _zero_budget_components(node_ids)
    return (
        incident_flux=_zero_initial_total_spectral_node_values(node_ids),
        incident_energy=_zero_initial_total_spectral_node_values(node_ids),
        absorbed_flux=_zero_initial_total_spectral_node_values(node_ids),
        absorbed_energy=_zero_initial_total_spectral_node_values(node_ids),
    )
end

function _store_node_budget!(
    budget,
    nid::Int,
    area::Float64,
    par0::Float64,
    nir0::Float64,
    par_scat::Float64,
    nir_scat::Float64,
    abs_par::Float64,
    abs_nir::Float64,
    step_duration_seconds::Float64,
)
    ri_par_0_f = par0 / area
    ri_nir_0_f = nir0 / area
    ri_par_f = (par0 + par_scat) / area
    ri_nir_f = (nir0 + nir_scat) / area

    budget.incident_flux.initial.par[nid] = ri_par_0_f
    budget.incident_flux.initial.nir[nid] = ri_nir_0_f
    budget.incident_flux.total.par[nid] = ri_par_f
    budget.incident_flux.total.nir[nid] = ri_nir_f

    budget.incident_energy.initial.par[nid] = par0 * step_duration_seconds
    budget.incident_energy.initial.nir[nid] = nir0 * step_duration_seconds
    budget.incident_energy.total.par[nid] = (par0 + par_scat) * step_duration_seconds
    budget.incident_energy.total.nir[nid] = (nir0 + nir_scat) * step_duration_seconds

    budget.absorbed_flux.initial.par[nid] = ri_par_0_f * abs_par
    budget.absorbed_flux.initial.nir[nid] = ri_nir_0_f * abs_nir
    budget.absorbed_flux.total.par[nid] = ri_par_f * abs_par
    budget.absorbed_flux.total.nir[nid] = ri_nir_f * abs_nir
    budget.absorbed_energy.initial.par[nid] = par0 * abs_par * step_duration_seconds
    budget.absorbed_energy.initial.nir[nid] = nir0 * abs_nir * step_duration_seconds
    budget.absorbed_energy.total.par[nid] = (par0 + par_scat) * abs_par * step_duration_seconds
    budget.absorbed_energy.total.nir[nid] = (nir0 + nir_scat) * abs_nir * step_duration_seconds
    return nothing
end

function _scale_extra_band_energy(extra_q_per_band, step_duration_seconds::Float64)
    Dict{String,Dict{Int,Float64}}(
        band => Dict{Int,Float64}(nid => value * step_duration_seconds for (nid, value) in values)
        for (band, values) in extra_q_per_band
    )
end

"""
    integrate_light(scene, models, first, scat, options; meteo_row=nothing, step_duration_seconds=nothing, check_boundaries=options.check_meteo_boundaries, ...)::LightBudget

Combine first-order interception and scattering into per-node incident and
absorbed light budgets.

The result stores both irradiance-style outputs (`*_f`, W m^-2) and energy
outputs (`*_q`, J component^-1 timestep^-1), plus optional extra-waveband
energies when they were carried through the pipeline.
"""
function integrate_light(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    first::FirstOrderResult,
    scat::Union{Nothing,ScatteringResult},
    options::LightOptions;
    meteo_row=nothing,
    extra_initial_energy_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    extra_energy_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    extra_emitter_escaped_power_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    step_duration_seconds::Union{Nothing,Float64}=nothing,
    component_area_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_par_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_nir_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    dt_seconds =
        if step_duration_seconds !== nothing
            _duration_seconds_strict(
                step_duration_seconds;
                field_name="step_duration_seconds",
            )
        elseif meteo_row !== nothing
            _resolved_meteo_step_or_error(
                meteo_row,
                options;
                check_boundaries=check_boundaries,
            ).duration_seconds
        else
            error("integrate_light requires `step_duration_seconds` or a resolvable `meteo_row`.")
        end
    node_ids = collect(keys(first.projected_area_per_node))
    budget = _zero_budget_components(node_ids)
    area_map = component_area_per_node === nothing ? _interception_area_per_node_local(scene, models, options) : component_area_per_node
    abs_par_map = absorption_par_per_node === nothing ? _node_absorptance_per_band(scene, models, options, "PAR") : absorption_par_per_node
    abs_nir_map = absorption_nir_per_node === nothing ? _node_absorptance_per_band(scene, models, options, "NIR") : absorption_nir_per_node
    default_abs_par = clamp(1.0 - options.scattering_coeff_par, 0.0, 1.0)
    default_abs_nir = clamp(1.0 - options.scattering_coeff_nir, 0.0, 1.0)
    for nid in node_ids
        pa = get(area_map, nid, get(first.projected_area_per_node, nid, 0.0))
        pa = max(pa, eps(Float64))
        p0 = get(first.incident_power.par, nid, 0.0)
        n0 = get(first.incident_power.nir, nid, 0.0)
        ps = scat === nothing ? 0.0 : get(scat.added_power.par, nid, 0.0)
        ns = scat === nothing ? 0.0 : get(scat.added_power.nir, nid, 0.0)
        abs_par =
            clamp(get(abs_par_map, nid, default_abs_par), 0.0, 1.0)
        abs_nir =
            clamp(get(abs_nir_map, nid, default_abs_nir), 0.0, 1.0)

        _store_node_budget!(budget, nid, pa, p0, n0, ps, ns, abs_par, abs_nir, dt_seconds)
    end

    emitter_escaped_power_per_band = Dict{String,Dict{Int,Float64}}(
        "PAR" => Dict{Int,Float64}(first.emitter_escaped_power.par),
        "NIR" => Dict{Int,Float64}(first.emitter_escaped_power.nir),
    )
    for (band, values) in extra_emitter_escaped_power_per_band
        emitter_escaped_power_per_band[uppercase(band)] = values
    end

    return LightBudget(
        budget.incident_flux,
        budget.incident_energy,
        budget.absorbed_flux,
        budget.absorbed_energy,
        _scale_extra_band_energy(extra_initial_energy_per_band, dt_seconds),
        _scale_extra_band_energy(extra_energy_per_band, dt_seconds),
        _scale_extra_band_energy(emitter_escaped_power_per_band, dt_seconds),
    )
end

struct DenseNodeMap{T}
    node_ids::Vector{Int}
    values::Vector{T}
    active_indices::Vector{Int}
end

Base.eltype(::Type{DenseNodeMap{T}}) where {T} = Pair{Int,T}

function _dense_node_map_active_indices(values::Vector{T}) where {T}
    active = Int[]
    sizehint!(active, min(length(values), 64))
    @inbounds for i in eachindex(values)
        iszero(values[i]) && continue
        push!(active, i)
    end
    return active
end

DenseNodeMap(node_ids::Vector{Int}, values::Vector{T}) where {T} =
    DenseNodeMap{T}(node_ids, values, _dense_node_map_active_indices(values))

function Base.iterate(map::DenseNodeMap{T}, state::Int=1) where {T}
    state > length(map.active_indices) && return nothing
    idx = map.active_indices[state]
    @inbounds v = map.values[idx]
    @inbounds nid = map.node_ids[idx]
    return ((nid, v), state + 1)
end

@inline Base.length(map::DenseNodeMap) = length(map.active_indices)
@inline _active_indices(map::DenseNodeMap) = map.active_indices

function _dense_node_map(values::Vector{T}, geometry::InterceptionSceneData) where {T}
    return DenseNodeMap(geometry.node_ids, values)
end

struct SectorResponsesCache
    prepared::PreparedInterceptionData
    projected_area_per_sector::Vector{DenseNodeMap{Float64}}
    hits_per_sector::Vector{DenseNodeMap{Int}}
    node_ids::Vector{Int}
    emitter_incident_power_par::DenseNodeMap{Float64}
    emitter_incident_power_nir::DenseNodeMap{Float64}
    emitter_escaped_power_par::DenseNodeMap{Float64}
    emitter_escaped_power_nir::DenseNodeMap{Float64}
    emitter_transfer::Union{Nothing,EmitterTransferResult}
    scattering_topology::Union{Nothing,ScatteringTopologyCache}
end

function _dense_sector_float(values::Dict{Int,Float64}, geometry::InterceptionSceneData)
    out = zeros(Float64, length(geometry.node_ids))
    for (nid, v) in values
        out[geometry.node_index[nid]] = v
    end
    return _dense_node_map(out, geometry)
end

function _dense_sector_int(values::Dict{Int,Int}, geometry::InterceptionSceneData)
    out = zeros(Int, length(geometry.node_ids))
    for (nid, v) in values
        out[geometry.node_index[nid]] = v
    end
    return _dense_node_map(out, geometry)
end

function _dense_emitter_power(
    transfer::Union{Nothing,EmitterTransferResult},
    prepared::PreparedInterceptionData,
    ;
    emitter_band_par::Union{Nothing,AbstractString}="PAR",
    emitter_band_nir::Union{Nothing,AbstractString}="NIR",
)
    geometry = prepared.geometry
    incident_par = zeros(Float64, length(geometry.node_ids))
    incident_nir = zeros(Float64, length(geometry.node_ids))
    escaped_par = zeros(Float64, length(geometry.node_ids))
    escaped_nir = zeros(Float64, length(geometry.node_ids))
    if transfer !== nothing
        _accumulate_emitter_band_power!(
            incident_par,
            escaped_par,
            transfer,
            _emitter_band_power_by_index(prepared, emitter_band_par),
            geometry,
        )
        _accumulate_emitter_band_power!(
            incident_nir,
            escaped_nir,
            transfer,
            _emitter_band_power_by_index(prepared, emitter_band_nir),
            geometry,
        )
    end
    return (
        _dense_node_map(incident_par, geometry),
        _dense_node_map(incident_nir, geometry),
        _dense_node_map(escaped_par, geometry),
        _dense_node_map(escaped_nir, geometry),
    )
end

function _turtle_cache_key(turtle::TurtleGrid, options::LightOptions)
    h = hash((length(turtle.sectors), options.pixel_size, options.area_ratio))
    for s in turtle.sectors
        d = s.direction
        h = hash(
            (
                round(d[1], digits=8),
                round(d[2], digits=8),
                round(d[3], digits=8),
                round(s.weight, digits=12),
                s.source,
            ),
            h,
        )
    end
    h
end

function _build_sector_responses(
    prepared::PreparedInterceptionData,
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    options::LightOptions,
)
    n = length(turtle.sectors)
    geometry = prepared.geometry
    pa_by_sector = Vector{DenseNodeMap{Float64}}(undef, n)
    hits_by_sector = Vector{DenseNodeMap{Int}}(undef, n)
    emitter_sector_fraction = _lambertian_sector_fractions(turtle)
    emitter_received_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Tuple{Int,Int},Float64}()
    emitter_observed_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Tuple{Int,Int},Float64}()
    emitter_escaped_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Int,Float64}()
    scattering_edge_counts = options.scattering ? Dict{UInt64,Int}() : nothing
    scattering_sun_hits = options.scattering ? Dict{Int,Int}() : nothing
    scattering_scratch = options.scattering ? ScatteringStackScratch() : nothing

    for i in 1:n
        sector = turtle.sectors[i]
        projection = _prepared_direction_projection(prepared, sector.direction, options)
        pa_by_sector[i] =
            _dense_node_map(
                _visible_area_from_projection_dense(
                    projection,
                    options,
                    geometry.plotbox,
                    prepared.virtual_node_mask,
                    prepared.node_transparency_by_index,
                    geometry,
                ),
                geometry,
            )
        hits_by_sector[i] = _dense_node_map(copy(_dense_projection_hits(projection, geometry)), geometry)
        if emitter_received_fraction !== nothing && emitter_sector_fraction[i] > 0.0
            emitter_edge_counts = Dict{UInt64,Int}()
            emitter_observed_edge_counts = Dict{UInt64,Int}()
            emitter_total_from = Dict{Int,Int}()
            if projection isa DenseDirectionProjectionResult
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_observed_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_node_mask,
                    prepared.virtual_node_mask,
                    geometry.node_ids;
                    stacks_sorted=true,
                )
            else
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_observed_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_nodes;
                    virtual_nodes=prepared.virtual_nodes,
                    stacks_sorted=true,
                )
            end
            _merge_emitter_direction_transfer!(
                emitter_received_fraction,
                emitter_observed_fraction,
                emitter_escaped_fraction,
                prepared.emitter_nodes,
                emitter_sector_fraction[i],
                emitter_edge_counts,
                emitter_observed_edge_counts,
                emitter_total_from,
            )
        end
        if scattering_edge_counts !== nothing
            if projection isa DenseDirectionProjectionResult
                _accumulate_scattering_counts!(
                    scattering_edge_counts,
                    scattering_sun_hits,
                    sector,
                    projection,
                    prepared.virtual_node_mask,
                    geometry.pavement_node_mask,
                    scattering_scratch;
                    node_ids=geometry.node_ids,
                    stacks_sorted=true,
                )
            else
                _accumulate_scattering_counts!(
                    scattering_edge_counts,
                    scattering_sun_hits,
                    sector,
                    projection,
                    prepared.virtual_nodes,
                    geometry.node_group,
                    scattering_scratch;
                    node_ids=geometry.node_ids,
                    stacks_sorted=true,
                )
            end
        end
    end
    emitter_transfer =
        emitter_received_fraction === nothing ? nothing :
        _finish_emitter_transfer(
            emitter_received_fraction,
            emitter_observed_fraction,
            emitter_escaped_fraction,
            prepared.emitter_nodes,
            emitter_sector_fraction,
        )
    (
        emitter_incident_power_par,
        emitter_incident_power_nir,
        emitter_escaped_power_par,
        emitter_escaped_power_nir,
    ) = _dense_emitter_power(emitter_transfer, prepared)
    scattering_topology =
        if scattering_edge_counts !== nothing
            _build_scattering_topology_cache(
                scene,
                models,
                prepared,
                _edge_counts_from_packed(scattering_edge_counts),
                scattering_sun_hits,
            )
        else
            nothing
        end
    return SectorResponsesCache(
        prepared,
        pa_by_sector,
        hits_by_sector,
        geometry.node_ids,
        emitter_incident_power_par,
        emitter_incident_power_nir,
        emitter_escaped_power_par,
        emitter_escaped_power_nir,
        emitter_transfer,
        scattering_topology,
    )
end

function _build_sector_responses(scene::PlantGeom.SceneGeometry, models::LightModels, turtle::TurtleGrid, options::LightOptions)
    prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true)
    return _build_sector_responses(prepared, scene, models, turtle, options)
end

function _sky_fraction_from_sector_responses(
    scene::PlantGeom.SceneGeometry,
    responses::SectorResponsesCache,
    turtle::TurtleGrid,
)
    sky_count = 0
    visible_sum = zeros(Float64, length(responses.node_ids))
    for (i, sector) in enumerate(turtle.sectors)
        sector.source == :sun && continue
        sky_count += 1
        sector_area = responses.projected_area_per_sector[i].values
        @inbounds for j in _active_indices(responses.projected_area_per_sector[i])
            visible_sum[j] += sector_area[j]
        end
    end

    out = Dict{Int,Float64}()
    @inbounds for i in eachindex(responses.node_ids)
        nid = responses.node_ids[i]
        area = _scene_area(scene, nid, 0.0)
        out[nid] = area <= 0.0 || sky_count == 0 ? 0.0 : visible_sum[i] / sky_count / area
    end
    return out
end

function _compute_sky_fraction(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    options::LightOptions;
    prepared::Union{Nothing,PreparedInterceptionData}=nothing,
    responses_cache::Union{Nothing,SectorResponsesCache}=nothing,
)
    responses =
        if responses_cache !== nothing
            responses_cache
        elseif prepared !== nothing
            _build_sector_responses(prepared, scene, models, turtle, options)
        else
            _build_sector_responses(scene, models, turtle, options)
        end
    return _sky_fraction_from_sector_responses(scene, responses, turtle)
end

function _combine_sector_responses(
    responses::SectorResponsesCache,
    fluxes::DirectionalFluxes,
    ;
    emitter_band_par::Union{Nothing,AbstractString}="PAR",
    emitter_band_nir::Union{Nothing,AbstractString}="NIR",
)
    node_ids = responses.node_ids
    projected_area_per_node = zeros(Float64, length(node_ids))
    incident_power_par = zeros(Float64, length(node_ids))
    incident_power_nir = zeros(Float64, length(node_ids))
    hits_per_node = zeros(Int, length(node_ids))

    for i in eachindex(responses.projected_area_per_sector)
        pf = fluxes.par[i]
        nf = fluxes.nir[i]
        active_flux = (pf != 0.0) || (nf != 0.0)

        sector_area = responses.projected_area_per_sector[i].values
        if active_flux
            @inbounds for j in _active_indices(responses.projected_area_per_sector[i])
                pa = sector_area[j]
                projected_area_per_node[j] += pa
                pf != 0.0 && (incident_power_par[j] += pf * pa)
                nf != 0.0 && (incident_power_nir[j] += nf * pa)
            end
        end
        sector_hits = responses.hits_per_sector[i].values
        @inbounds for j in _active_indices(responses.hits_per_sector[i])
            hits_per_node[j] += sector_hits[j]
        end
    end

    (
        emitter_incident_power_par,
        emitter_incident_power_nir,
        emitter_escaped_power_par,
        emitter_escaped_power_nir,
    ) =
        emitter_band_par == "PAR" && emitter_band_nir == "NIR" ?
        (
            responses.emitter_incident_power_par,
            responses.emitter_incident_power_nir,
            responses.emitter_escaped_power_par,
            responses.emitter_escaped_power_nir,
        ) :
        _dense_emitter_power(
            responses.emitter_transfer,
            responses.prepared;
            emitter_band_par=emitter_band_par,
            emitter_band_nir=emitter_band_nir,
        )

    emitter_par = emitter_incident_power_par.values
    @inbounds for idx in _active_indices(emitter_incident_power_par)
        incident_power_par[idx] += emitter_par[idx]
    end
    emitter_nir = emitter_incident_power_nir.values
    @inbounds for idx in _active_indices(emitter_incident_power_nir)
        incident_power_nir[idx] += emitter_nir[idx]
    end

    return FirstOrderResult(
        _all_dense_float_node_map(node_ids, projected_area_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(node_ids, incident_power_par),
            _all_dense_float_node_map(node_ids, incident_power_nir),
        ),
        _all_dense_int_node_map(node_ids, hits_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(node_ids, emitter_escaped_power_par.values),
            _all_dense_float_node_map(node_ids, emitter_escaped_power_nir.values),
        ),
    )
end

function _stream_first_order_with_scattering_topology(
    prepared::PreparedInterceptionData,
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    geometry = prepared.geometry
    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    emitter_escaped_power_par = zeros(Float64, length(geometry.node_ids))
    emitter_escaped_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))
    scattering_edge_counts = Dict{UInt64,Int}()
    scattering_sun_hits = Dict{Int,Int}()
    scattering_scratch = ScatteringStackScratch()
    emitter_sector_fraction = _lambertian_sector_fractions(turtle)
    emitter_received_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Tuple{Int,Int},Float64}()
    emitter_observed_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Tuple{Int,Int},Float64}()
    emitter_escaped_fraction =
        isempty(prepared.emitter_nodes) ? nothing : Dict{Int,Float64}()

    for (k, sector) in enumerate(turtle.sectors)
        projection = _prepared_direction_projection(prepared, sector.direction, options)
        visible_area =
            _visible_area_from_projection_dense(
                projection,
                options,
                geometry.plotbox,
                prepared.virtual_node_mask,
                prepared.node_transparency_by_index,
                geometry,
            )

        _accumulate_projection_hits!(hits_per_node, projection, geometry)

        par_flux = fluxes.par[k]
        nir_flux = options.nir_interception ? fluxes.nir[k] : 0.0
        if par_flux != 0.0 || nir_flux != 0.0
            @inbounds for idx in eachindex(visible_area)
                pa = visible_area[idx]
                pa <= 0.0 && continue
                projected_area_per_node[idx] += pa
                par_flux != 0.0 && (incident_power_par[idx] += par_flux * pa)
                nir_flux != 0.0 && (incident_power_nir[idx] += nir_flux * pa)
            end
        end

        if emitter_received_fraction !== nothing && emitter_sector_fraction[k] > 0.0
            emitter_edge_counts = Dict{UInt64,Int}()
            emitter_observed_edge_counts = Dict{UInt64,Int}()
            emitter_total_from = Dict{Int,Int}()
            if projection isa DenseDirectionProjectionResult
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_observed_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_node_mask,
                    prepared.virtual_node_mask,
                    geometry.node_ids;
                    stacks_sorted=true,
                )
            else
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_observed_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_nodes;
                    virtual_nodes=prepared.virtual_nodes,
                    stacks_sorted=true,
                )
            end
            _merge_emitter_direction_transfer!(
                emitter_received_fraction,
                emitter_observed_fraction,
                emitter_escaped_fraction,
                prepared.emitter_nodes,
                emitter_sector_fraction[k],
                emitter_edge_counts,
                emitter_observed_edge_counts,
                emitter_total_from,
            )
        end
        if projection isa DenseDirectionProjectionResult
            _accumulate_scattering_counts!(
                scattering_edge_counts,
                scattering_sun_hits,
                sector,
                projection,
                prepared.virtual_node_mask,
                geometry.pavement_node_mask,
                scattering_scratch;
                node_ids=geometry.node_ids,
                stacks_sorted=true,
            )
        else
            _accumulate_scattering_counts!(
                scattering_edge_counts,
                scattering_sun_hits,
                sector,
                projection,
                prepared.virtual_nodes,
                geometry.node_group,
                scattering_scratch;
                node_ids=geometry.node_ids,
                stacks_sorted=true,
            )
        end
    end

    if emitter_received_fraction !== nothing
        emitter_transfer = _finish_emitter_transfer(
            emitter_received_fraction,
            emitter_observed_fraction,
            emitter_escaped_fraction,
            prepared.emitter_nodes,
            emitter_sector_fraction,
        )
        (
            emitter_incident_power_par,
            emitter_incident_power_nir,
            escaped_power_par,
            escaped_power_nir,
        ) = _dense_emitter_power(
            emitter_transfer,
            prepared;
            emitter_band_nir=options.nir_interception ? "NIR" : nothing,
        )
        @inbounds for idx in _active_indices(emitter_incident_power_par)
            incident_power_par[idx] += emitter_incident_power_par.values[idx]
        end
        @inbounds for idx in _active_indices(emitter_incident_power_nir)
            incident_power_nir[idx] += emitter_incident_power_nir.values[idx]
        end
        emitter_escaped_power_par .= escaped_power_par.values
        emitter_escaped_power_nir .= escaped_power_nir.values
    end

    first = FirstOrderResult(
        _all_dense_float_node_map(geometry.node_ids, projected_area_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(geometry.node_ids, incident_power_par),
            _all_dense_float_node_map(geometry.node_ids, incident_power_nir),
        ),
        _all_dense_int_node_map(geometry.node_ids, hits_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(geometry.node_ids, emitter_escaped_power_par),
            _all_dense_float_node_map(geometry.node_ids, emitter_escaped_power_nir),
        ),
    )
    topology = _build_scattering_topology_cache(
        scene,
        models,
        prepared,
        _edge_counts_from_packed(scattering_edge_counts),
        scattering_sun_hits,
    )
    return first, topology
end

function _row_number_local(row, name::Symbol, default::Float64=0.0)
    if name in _row_propertynames(row)
        v = getproperty(row, name)
        if v isa Number
            return Float64(v)
        elseif v !== missing
            try
                return parse(Float64, string(v))
            catch
            end
        end
    end
    return default
end

function _parse_time_or_default_local(v)
    v isa Dates.Time && return v
    v isa Dates.DateTime && return Dates.Time(v)
    s = strip(string(v))
    if isempty(s)
        return Dates.Time(0)
    end
    try
        Dates.Time(s, Dates.DateFormat("HH:MM:SS"))
    catch
        try
            Dates.Time(s, Dates.DateFormat("HH:MM"))
        catch
            Dates.Time(0)
        end
    end
end

function _step_duration_seconds_local(row)
    names = _row_propertynames(row)
    if (:step_duration in names)
        return _duration_seconds_strict(getproperty(row, :step_duration); field_name="step_duration")
    end
    if (:duration in names)
        return _duration_seconds_strict(getproperty(row, :duration); field_name="duration")
    end
    if (:hour_start in names) && (:hour_end in names)
        t0 = _parse_time_or_default_local(getproperty(row, :hour_start))
        t1 = _parse_time_or_default_local(getproperty(row, :hour_end))
        dt0 = Dates.DateTime(Dates.Date(2000, 1, 1), t0)
        dt1 = Dates.DateTime(Dates.Date(2000, 1, 1), t1)
        dt1 < dt0 && (dt1 += Dates.Day(1))
        dt_seconds = Dates.value(dt1 - dt0) / 1000.0
        dt_seconds > 0.0 || error("Invalid meteo timestep: non-positive duration from hour_start/hour_end.")
        return dt_seconds
    end
    return 1.0
end

function _row_datetime_interval_local(row; index::Int=0)
    errors = String[]
    sources = Dict{Symbol,Symbol}()
    temporal = _meteo_resolve_temporal!(errors, row, index, sources)
    isempty(errors) || error(
        "invalid datetime interval at meteo row $(index): $(join(errors, "; "))",
    )
    temporal.date !== nothing || error(
        "invalid datetime interval at meteo row $(index): a finite date is required for datetime meteo_range selection",
    )
    temporal.start_hour !== nothing && temporal.end_hour !== nothing || error(
        "invalid datetime interval at meteo row $(index): a start time and positive duration are required",
    )

    base = Dates.DateTime(temporal.date)
    start = base + Dates.Millisecond(round(Int, temporal.start_hour * 3_600_000))
    stop = base + Dates.Millisecond(round(Int, temporal.end_hour * 3_600_000))
    return start, stop
end

function _cfg_meteo_range_spec_local(options::LightOptions)
    v = options.meteo_range
    v === nothing && return nothing
    s = strip(string(v))
    isempty(s) ? nothing : s
end

function _parse_bool_strict_local(v, field_name::String)
    if v isa Bool
        return v
    elseif v isa Integer
        (v == 0 || v == 1) && return v == 1
    elseif v isa Number
        isfinite(v) && (abs(v) < eps(Float64) || abs(v - 1.0) < eps(Float64)) && return v > 0
    end

    s = lowercase(strip(string(v)))
    s in ("true", "t", "1", "yes", "y") && return true
    s in ("false", "f", "0", "no", "n") && return false
    error("invalid $(field_name) value: $(repr(v))")
end

function _cfg_bool_override_local(options::LightOptions, key::String)
    key == "nir_scattering" && return options.nir_scattering
    key == "nir_interception" && return options.nir_interception
    return nothing
end

"""
    _nir_scattering_enabled_local(options)::Bool

NIR scattering activation with optional explicit override (`nir_scattering`).
Default behavior remains enabled whenever scattering is enabled.
"""
function _nir_scattering_enabled_local(options::LightOptions)
    options.scattering || return false
    override = _cfg_bool_override_local(options, "nir_scattering")
    override === nothing ? true : override
end

"""
    _nir_interception_enabled_local(options)::Bool

NIR interception activation with optional explicit override (`nir_interception`).
Default behavior remains enabled.
"""
function _nir_interception_enabled_local(options::LightOptions)
    override = _cfg_bool_override_local(options, "nir_interception")
    override === nothing ? true : override
end

function _disable_nir_sky_local(sky::SkyState)
    SkyState(
        sky.sun_azimuth_deg,
        sky.sun_elevation_deg,
        sky.ri_sw_f,
        sky.ri_par_f,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
    )
end

function _disable_nir_first_order_local(first::FirstOrderResult)
    zero_incident = Dict{Int,Float64}(nid => 0.0 for nid in keys(first.incident_power.nir))
    zero_escaped =
        Dict{Int,Float64}(nid => 0.0 for nid in keys(first.emitter_escaped_power.nir))
    return FirstOrderResult(
        first.projected_area_per_node,
        SpectralNodeValues(first.incident_power.par, zero_incident),
        first.hits_per_node,
        SpectralNodeValues(first.emitter_escaped_power.par, zero_escaped),
    )
end

function _parse_range_datetime_token_local(s::AbstractString)
    ss = strip(String(s))
    for fmt in (
        Dates.DateFormat("yyyy/mm/dd-HH:MM:SS"),
        Dates.DateFormat("yyyy/mm/dd-HH:MM"),
        Dates.DateFormat("yyyy-mm-dd-HH:MM:SS"),
        Dates.DateFormat("yyyy-mm-dd-HH:MM"),
    )
        try
            return Dates.DateTime(ss, fmt)
        catch
        end
    end
    error("invalid meteo_range datetime token: $(repr(s))")
end

function _meteo_range_indices_local(
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions,
)
    spec = _cfg_meteo_range_spec_local(options)
    spec === nothing && return collect(1:length(meteo))

    parts = split(spec, ","; limit=2)
    length(parts) == 2 || error("invalid meteo_range format: $(repr(spec))")
    a = strip(parts[1])
    b = strip(parts[2])
    isempty(a) && error("invalid meteo_range format: $(repr(spec))")
    isempty(b) && error("invalid meteo_range format: $(repr(spec))")

    ia = tryparse(Int, a)
    ib = tryparse(Int, b)
    if ia !== nothing && ib !== nothing
        n = length(meteo)
        ia >= 1 || error("invalid meteo_range: start step must be >= 1")
        ib >= ia || error("invalid meteo_range: end step is before start step")
        ib <= n || error("invalid meteo_range: end step exceeds meteo size")
        return collect(ia:ib)
    end

    t0 = _parse_range_datetime_token_local(a)
    t1 = _parse_range_datetime_token_local(b)
    t1 >= t0 || error("invalid meteo_range: end datetime is before start datetime")

    selected = Int[]
    for (i, row) in enumerate(meteo)
        s, e = _row_datetime_interval_local(row; index=i)
        # Java uses closed-interval overlap semantics.
        s <= t1 && e >= t0 && push!(selected, i)
    end
    isempty(selected) && error("invalid meteo_range: selection is empty")
    return selected
end

function _apply_meteo_range_local(
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions,
)
    selected = _meteo_range_indices_local(meteo, options)
    length(selected) == length(meteo) && selected == collect(1:length(meteo)) &&
        return meteo
    return meteo[selected]
end

function _meteo_active_indices_local(
    meteo::PlantMeteo.TimeStepTable;
    source_indices::AbstractVector{<:Integer}=collect(1:length(meteo)),
)
    isempty(meteo) && return Int[]
    length(source_indices) == length(meteo) ||
        error("Internal meteo row-index mismatch during active selection.")
    names = _row_propertynames(first(meteo))
    :active in names || return collect(1:length(meteo))

    selected = Int[]
    for (i, row) in enumerate(meteo)
        flag = try
            _parse_bool_strict_local(getproperty(row, :active), "active")
        catch err
            error(
                "invalid active value at meteo row $(source_indices[i]): $(sprint(showerror, err))",
            )
        end
        flag && push!(selected, i)
    end
    isempty(selected) && error("invalid meteo: no active meteo step")
    return selected
end

function _apply_meteo_active_filter_local(meteo::PlantMeteo.TimeStepTable)
    selected = _meteo_active_indices_local(meteo)
    length(selected) == length(meteo) && selected == collect(1:length(meteo)) &&
        return meteo
    return meteo[selected]
end

function _select_meteo_with_source_indices_local(
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions,
)
    meteo_for_selection = _with_inferred_datetime_durations(meteo)
    range_indices = _meteo_range_indices_local(meteo_for_selection, options)
    ranged =
        length(range_indices) == length(meteo_for_selection) &&
        range_indices == collect(1:length(meteo_for_selection)) ?
        meteo_for_selection : meteo_for_selection[range_indices]
    active_indices = _meteo_active_indices_local(
        ranged;
        source_indices=range_indices,
    )
    selected =
        length(active_indices) == length(ranged) &&
        active_indices == collect(1:length(ranged)) ? ranged : ranged[active_indices]
    selected = _with_inferred_datetime_durations(selected)
    return selected, range_indices[active_indices]
end

function _resolved_datetime_interval_local(step::ResolvedMeteoStep)
    step.start_hour !== nothing && step.end_hour !== nothing || return nothing
    base = Dates.DateTime(something(step.date, Dates.Date(2000, 1, 1)))
    start = base + Dates.Millisecond(round(Int, step.start_hour * 3_600_000))
    stop = base + Dates.Millisecond(round(Int, step.end_hour * 3_600_000))
    return start, stop
end

function _validate_meteo_sequence_local(steps::Vector{ResolvedMeteoStep})
    previous_stop = nothing
    for step in steps
        interval = _resolved_datetime_interval_local(step)
        interval === nothing && continue
        start, stop = interval
        if previous_stop !== nothing && start < previous_stop
            error("invalid overlapping meteo steps at row $(step.row_index)")
        end
        previous_stop = stop
    end
    return nothing
end

_meteo_rows(meteo::PlantMeteo.TimeStepTable) = collect(meteo)

function _prepare_meteo_for_series(
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions;
    check_boundaries::Bool=options.check_meteo_boundaries,
    return_resolved::Bool=false,
)
    isempty(meteo) && error("invalid meteo: no meteo steps")
    source_indices = Int[]
    meteo, source_indices = _select_meteo_with_source_indices_local(meteo, options)
    errors = String[]
    warnings = String[]
    infos = String[]
    resolved_steps = ResolvedMeteoStep[]
    for (i, row) in enumerate(meteo)
        source_index = source_indices[i]
        result = _resolve_meteo_step(
            row,
            options;
            row_index=source_index,
            check_boundaries=check_boundaries,
        )
        append!(errors, result.errors)
        append!(warnings, result.warnings)
        append!(infos, result.infos)
        result.step === nothing || push!(resolved_steps, result.step)
    end
    report = ValidationReport(unique(errors), unique(warnings), unique(infos))
    if !isempty(report.errors)
        error(_validation_error_message(
            "Invalid meteo input",
            report;
            next="Run `check_meteo(meteo; options=options, check_boundaries=$(check_boundaries))` to inspect the failing rows and variables.",
        ))
    end
    length(resolved_steps) == length(meteo) || error("Meteo resolution did not produce every selected step.")
    options.allow_overlapping_meteo_steps ||
        _validate_meteo_sequence_local(resolved_steps)
    return return_resolved ? (meteo, resolved_steps) : meteo
end

"""
    prepare_meteo(meteo, options; check_boundaries=options.check_meteo_boundaries)::PlantMeteo.TimeStepTable

Return the effective meteo table after Java-like meteo controls are applied:
sequence validation, optional `meteo_range`, and optional `active` filtering.

Arguments:

- `meteo`: a `PlantMeteo.TimeStepTable` or Tables.jl-compatible table of meteo
  rows.
- `options`: [`LightOptions`](@ref) controlling overlap validation,
  `meteo_range`, and active-row filtering.

Keywords:

- `check_boundaries`: check physical ranges for used meteo values. Missing,
  nonfinite, derivability, and conflicting-input checks always remain enabled.
"""
function prepare_meteo(
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions;
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    _prepare_meteo_for_series(meteo, options; check_boundaries=check_boundaries)
end

function prepare_meteo(
    meteo,
    options::LightOptions;
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected PlantMeteo.TimeStepTable or a Tables.jl-compatible table.")
    meteo isa PlantMeteo.TimeStepTable &&
        return prepare_meteo(meteo, options; check_boundaries=check_boundaries)
    prepare_meteo(
        _as_plantmeteo_table(meteo),
        options;
        check_boundaries=check_boundaries,
    )
end

function _default_scattering_factor_local(options::LightOptions, band::String)
    b = uppercase(band)
    b == "NIR" && return options.scattering_coeff_nir
    return options.scattering_coeff_par
end

function _compute_scattering_with_flags(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    options::LightOptions;
    mode::Symbol=:raycast,
    backend::Union{Nothing,ScatteringBackend}=nothing,
    nir_scattering::Bool=true,
    responses_cache::Union{Nothing,SectorResponsesCache}=nothing,
    scattering_topology::Union{Nothing,ScatteringTopologyCache}=nothing,
)
    options.scattering || return nothing
    if nir_scattering
        scattering_topology !== nothing && return compute_scattering(scattering_topology, first, options; mode=mode, backend=backend)
        if responses_cache !== nothing && responses_cache.scattering_topology !== nothing
            return compute_scattering(responses_cache.scattering_topology, first, options; mode=mode, backend=backend)
        end
        return compute_scattering(scene, models, turtle, first, options; mode=mode, backend=backend)
    end

    par_only =
        if scattering_topology !== nothing
            compute_scattering_band(
                scattering_topology,
                first,
                options;
                mode=mode,
                backend=backend,
                band="PAR",
                initial_power_per_node=first.incident_power.par,
                default_coeff=options.scattering_coeff_par,
            )
        elseif responses_cache !== nothing && responses_cache.scattering_topology !== nothing
            compute_scattering_band(
                responses_cache.scattering_topology,
                first,
                options;
                mode=mode,
                backend=backend,
                band="PAR",
                initial_power_per_node=first.incident_power.par,
                default_coeff=options.scattering_coeff_par,
            )
        else
            compute_scattering_band(
                scene,
                models,
                turtle,
                first,
                options;
                mode=mode,
                backend=backend,
                band="PAR",
                initial_power_per_node=first.incident_power.par,
                default_coeff=options.scattering_coeff_par,
            )
        end
    return ScatteringResult(SpectralNodeValues(par_only.added_power_per_node, Dict{Int,Float64}()), par_only.iterations, par_only.converged)
end

function _interception_area_per_node_local(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    return _interception_area_per_node_from_geometry(geometry)
end

function _node_absorptance_per_band(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions, band::String)
    geometry = _scene_geometry_for_interception(scene, models, options)
    virtual_nodes = _virtual_sensor_node_ids(geometry.node_group, geometry.node_type, models)
    return _node_absorptance_per_band_from_geometry(scene, models, options, geometry, virtual_nodes, band)
end

function _extra_band_irradiance(meteo_row)
    extras = Dict{String,Float64}()
    for p in _row_propertynames(meteo_row)
        s = String(p)
        su = uppercase(s)
        startswith(su, "RI_") || continue
        endswith(su, "_F") || continue
        # TIR is thermal forcing for the energy-balance model in Java ARCHIMED,
        # not a shortwave interception/scattering band. ArchimedLight is
        # deliberately light-only, so keep it out of the extra-band pipeline.
        su in ("RI_PAR_F", "RI_NIR_F", "RI_SW_F", "RI_TIR_F") && continue
        band = su[4:(end - 2)]
        isempty(band) && continue
        v = _row_number_local(meteo_row, p, NaN)
        isfinite(v) || continue
        extras[band] = max(v, 0.0)
    end
    extras
end

function _extra_band_irradiance(meteo_row, options::LightOptions)
    scale = _radiation_input_effective_scale(meteo_row, options)
    return Dict{String,Float64}(
        band => irradiance * scale
        for (band, irradiance) in _extra_band_irradiance(meteo_row)
    )
end

function _extra_band_irradiance(
    resolved::ResolvedMeteoStep,
    meteo_row,
    options::LightOptions,
)
    scale = _radiation_input_effective_scale(meteo_row, options, resolved)
    return Dict{String,Float64}(
        band => max(irradiance, 0.0) * scale
        for (band, irradiance) in resolved.extra_irradiance
    )
end

function _include_emitter_extra_bands!(
    extras_irradiance::Dict{String,Float64},
    emitter_power_per_band::Dict{String,Dict{Int,Float64}},
)
    for band in keys(emitter_power_per_band)
        band in ("PAR", "NIR", "TIR") && continue
        get!(extras_irradiance, band, 0.0)
    end
    return extras_irradiance
end

function _single_band_flux(
    total_irradiance::Float64,
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    options::LightOptions;
    resolved_step::Union{Nothing,ResolvedMeteoStep}=nothing,
)
    total_irradiance <= 0.0 && return zeros(Float64, length(turtle.sectors))

    # Extra shortwave bands follow the same directional sky distribution as
    # the effective PAR/NIR forcing. Reusing that distribution is important
    # for clearness-only rows: re-running `_radiation_substeps` with a temporary
    # custom-band SkyState would otherwise replace the supplied custom-band
    # magnitude with a newly derived total shortwave irradiance.
    reference = compute_directional_fluxes(
        meteo_row,
        sky,
        turtle,
        options;
        resolved_step=resolved_step,
    )
    reference_sw = reference.par .+ reference.nir
    reference_total = sum(reference_sw)
    if reference_total > 0.0
        return reference_sw .* (total_irradiance / reference_total)
    end

    # Defensive fallback for an interval with custom forcing but no usable
    # PAR/NIR reference (for example, an interval-mean nighttime input).
    tmp = SkyState(
        sky.sun_azimuth_deg,
        sky.sun_elevation_deg,
        total_irradiance,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
    )
    compute_directional_fluxes(tmp, turtle, options).par
end

function _compute_extra_band_light(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    options::LightOptions;
    interception_backend::InterceptionBackend=RasterCPUBackend(),
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    responses_cache::Union{Nothing,SectorResponsesCache}=nothing,
    resolved_step::Union{Nothing,ResolvedMeteoStep}=nothing,
)
    extras_irr =
        meteo_row === nothing ? Dict{String,Float64}() :
        resolved_step === nothing ? _extra_band_irradiance(meteo_row, options) :
        _extra_band_irradiance(resolved_step, meteo_row, options)
    emitter_power_per_band =
        responses_cache === nothing ?
        _emitter_power_per_band_per_node(scene, models) :
        responses_cache.prepared.emitter_power_per_band_per_node
    _include_emitter_extra_bands!(extras_irr, emitter_power_per_band)
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()
    extra_emitter_escaped_power = Dict{String,Dict{Int,Float64}}()

    isempty(extras_irr) &&
        return extra_0_q, extra_q, extras_irr, extra_emitter_escaped_power

    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(
            total_irr,
            meteo_row,
            sky,
            turtle,
            options;
            resolved_step=resolved_step,
        )
        first_band =
            if responses_cache === nothing
                compute_first_order(
                    scene,
                    models,
                    turtle,
                    DirectionalFluxes(ids, flux_band, zeros(Float64, n)),
                    options;
                    backend=interception_backend,
                    emitter_band_par=band,
                    emitter_band_nir=nothing,
                )
            else
                _combine_sector_responses(
                    responses_cache,
                    DirectionalFluxes(ids, flux_band, zeros(Float64, n));
                    emitter_band_par=band,
                    emitter_band_nir=nothing,
                )
            end
        order0 = Dict{Int,Float64}(nid => v for (nid, v) in first_band.incident_power.par)

        added =
            if options.scattering
                if responses_cache !== nothing && responses_cache.scattering_topology !== nothing
                    compute_scattering_band(
                        responses_cache.scattering_topology,
                        first_band,
                        options;
                        mode=scattering_mode,
                        backend=scattering_backend,
                        band=band,
                    ).added_power_per_node
                else
                    compute_scattering_band(
                        scene,
                        models,
                        turtle,
                        first_band,
                        options;
                        mode=scattering_mode,
                        backend=scattering_backend,
                        band=band,
                    ).added_power_per_node
                end
            else
                Dict{Int,Float64}()
            end

        total = Dict{Int,Float64}()
        for nid in union(keys(order0), keys(added))
            total[nid] = get(order0, nid, 0.0) + get(added, nid, 0.0)
        end
        extra_0_q[band] = order0
        extra_q[band] = total
        extra_emitter_escaped_power[band] =
            Dict{Int,Float64}(first_band.emitter_escaped_power.par)
    end
    return extra_0_q, extra_q, extras_irr, extra_emitter_escaped_power
end

function _can_use_series_radiation_cache(::RasterCPUBackend)
    true
end

function _can_use_series_radiation_cache(::InterceptionBackend)
    false
end

mutable struct TurtleLightCacheEntry
    key::UInt64
    turtle::TurtleGrid
    responses_cache::SectorResponsesCache
    par_added_per_sector::Vector{Union{Nothing,Vector{Float64}}}
    nir_added_per_sector::Vector{Union{Nothing,Vector{Float64}}}
    extra_added_per_sector::Dict{String,Vector{Union{Nothing,Vector{Float64}}}}
    par_iterations_per_sector::Vector{Int}
    nir_iterations_per_sector::Vector{Int}
    par_converged_per_sector::Vector{Bool}
    nir_converged_per_sector::Vector{Bool}
    extra_iterations_per_sector::Dict{String,Vector{Int}}
    extra_converged_per_sector::Dict{String,Vector{Bool}}
    resident_bytes::Int
    last_used_tick::Int
end

mutable struct LightSimulationCache
    scene::PlantGeom.SceneGeometry
    models::LightModels
    options::LightOptions
    interception_backend::Any
    resolved_interception_backend::InterceptionBackend
    scattering_mode::Symbol
    scattering_backend::Union{Nothing,ScatteringBackend}
    prepared::Union{Nothing,PreparedInterceptionData}
    render_geometry::LightRenderGeometry
    node_metadata::Union{Nothing,LightNodeMetadata}
    mode::Symbol
    estimated_entry_bytes::Int
    memory_limit_bytes::Int
    resident_bytes::Int
    tick::Int
    entries::Dict{UInt64,TurtleLightCacheEntry}
end

function _node_metadata_narrow(values::Vector{Any})
    isempty(values) && return Any[]
    value_type = foldl((acc, value) -> typejoin(acc, typeof(value)), values; init=Union{})
    value_type === Union{} && return Any[]
    value_type[values...]
end

function _node_metadata_attribute_table(names::Tuple, columns::AbstractVector)
    narrowed = Tuple(_node_metadata_narrow(column) for column in columns)
    NamedTuple{names}(narrowed)
end

function _validate_lightweight_node_metadata_value(value, attribute::Symbol, nid::Int)
    lightweight = value === nothing || value === missing || value isa Number ||
                  value isa AbstractString || value isa Symbol || value isa Char || isbits(value)
    lightweight && return value
    throw(ArgumentError(
        "node_metadata_attributes requested `$attribute`, but node $nid stores " *
        "a non-scalar $(typeof(value)). Only lightweight scalar identifiers can be retained.",
    ))
end

function _build_light_node_metadata(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions,
    geometry_node_ids,
)
    node_ids = sort!(collect(Int, geometry_node_ids))
    output_keys = _interception_output_keys_for_node_ids(scene, node_ids)
    source_topology_ids = Int[]
    object_ids = Int[]
    item_ids = Int[]
    component_ids = Int[]
    groups = String[]
    types = String[]
    symbols = Symbol[]
    scales = Int[]
    attribute_names = options.node_metadata_attributes
    local_columns = [Any[] for _ in attribute_names]
    inherited_columns = [Any[] for _ in attribute_names]

    scene.mtg === nothing && throw(
        ArgumentError("store_node_metadata=true requires an MTG-backed scene."),
    )
    available_attributes = Set(Symbol.(MultiScaleTreeGraph.get_attributes(scene.mtg)))
    missing_attributes = Symbol[name for name in attribute_names if !(name in available_attributes)]
    isempty(missing_attributes) || throw(ArgumentError(
        "Unknown node_metadata_attributes: $(join(string.(missing_attributes), ", ")).",
    ))

    for nid in node_ids
        node = _scene_mtg_node(scene, nid)
        node === nothing && throw(ArgumentError("Geometric node id $nid is absent from the scene MTG."))
        source_topology_id = _scene_source_topology_id(scene, nid, nid)
        object_id = _scene_object_id(scene, nid, -1)
        item_id, component_id = get(output_keys, nid, (object_id, source_topology_id))
        push!(source_topology_ids, source_topology_id)
        push!(object_ids, object_id)
        push!(item_ids, item_id)
        push!(component_ids, component_id)
        push!(groups, _scene_group(scene, nid, ""))
        push!(types, _scene_display_type(scene, models, nid))
        push!(symbols, Symbol(MultiScaleTreeGraph.symbol(node)))
        push!(scales, Int(MultiScaleTreeGraph.scale(node)))

        for (i, attribute) in enumerate(attribute_names)
            local_value = MultiScaleTreeGraph.attribute(node, attribute; default=nothing)
            inherited_value = _inherited_attr(scene, nid, (attribute,), nothing)
            push!(
                local_columns[i],
                _validate_lightweight_node_metadata_value(local_value, attribute, nid),
            )
            push!(
                inherited_columns[i],
                _validate_lightweight_node_metadata_value(inherited_value, attribute, nid),
            )
        end
    end

    LightNodeMetadata(
        node_ids,
        source_topology_ids,
        object_ids,
        item_ids,
        component_ids,
        groups,
        types,
        symbols,
        scales,
        _node_metadata_attribute_table(attribute_names, local_columns),
        _node_metadata_attribute_table(attribute_names, inherited_columns),
    )
end

function _default_light_cache_memory_limit()
    fallback = 512 * 1024 * 1024
    try
        total = Sys.total_memory()
        total > 0 || return fallback
        return max(128 * 1024 * 1024, min(fallback, Int(fld(total, 8))))
    catch
        return fallback
    end
end

function _estimate_light_cache_entry_bytes(prepared::PreparedInterceptionData, options::LightOptions)
    n_nodes = length(prepared.geometry.node_ids)
    n_sectors = max(options.turtle_sectors + (options.all_in_turtle ? 0 : 1), 1)
    base_sector_bytes = n_nodes * (sizeof(Float64) + sizeof(Int))
    scatter_sector_bytes = options.scattering ? n_nodes * 2 * sizeof(Float64) : 0
    base = Int128(n_sectors) * (base_sector_bytes + scatter_sector_bytes)
    if !isempty(prepared.emitter_nodes)
        # A transfer can contain receiver/source and observer/source pairs.
        # Use a conservative pair-dictionary upper bound for mode selection,
        # then replace it with measured storage in each built entry.
        pair_count = Int128(length(prepared.emitter_nodes)) * n_nodes
        edge_bytes = 2 * pair_count * 64
        dense_emitter_bytes = Int128(4) * n_nodes * sizeof(Float64)
        base += edge_bytes + dense_emitter_bytes
    end
    return Int(min(base, Int128(typemax(Int))))
end

function _turtle_cache_entry_retained_bytes(entry::TurtleLightCacheEntry)
    responses = entry.responses_cache
    # Measure every allocation retained only because this entry is resident.
    # `responses.prepared` is deliberately represented by one pointer-sized
    # scalar instead of traversed: the prepared geometry is owned by the parent
    # LightSimulationCache even when no turtle response is cached. Shared node
    # vectors reached through the response maps are counted, making this a
    # conservative retained-size measurement.
    retained = (
        entry.key,
        entry.turtle,
        UInt(0),
        responses.projected_area_per_sector,
        responses.hits_per_sector,
        responses.node_ids,
        responses.emitter_incident_power_par,
        responses.emitter_incident_power_nir,
        responses.emitter_escaped_power_par,
        responses.emitter_escaped_power_nir,
        responses.emitter_transfer,
        responses.scattering_topology,
        entry.par_added_per_sector,
        entry.nir_added_per_sector,
        entry.extra_added_per_sector,
        entry.par_iterations_per_sector,
        entry.nir_iterations_per_sector,
        entry.par_converged_per_sector,
        entry.nir_converged_per_sector,
        entry.extra_iterations_per_sector,
        entry.extra_converged_per_sector,
        entry.resident_bytes,
        entry.last_used_tick,
    )
    return Base.summarysize(retained)
end

@inline function _dense_vector_from_node_values(values::Dict{Int,Float64}, geometry::InterceptionSceneData)
    out = zeros(Float64, length(geometry.node_ids))
    for (nid, v) in values
        out[geometry.node_index[nid]] = v
    end
    return out
end

function _float_dict_from_dense(node_ids::Vector{Int}, values::Vector{Float64})
    out = Dict{Int,Float64}()
    @inbounds for i in eachindex(node_ids)
        v = values[i]
        iszero(v) && continue
        out[node_ids[i]] = v
    end
    return out
end

@inline function _accumulate_scaled!(dst::Vector{Float64}, src::Vector{Float64}, scale::Float64)
    iszero(scale) && return dst
    @inbounds for i in eachindex(dst)
        dst[i] += src[i] * scale
    end
    return dst
end

@inline function _unit_directional_fluxes(turtle::TurtleGrid, sector_idx::Int; par::Float64=0.0, nir::Float64=0.0)
    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    par_flux = zeros(Float64, n)
    nir_flux = zeros(Float64, n)
    par != 0.0 && (par_flux[sector_idx] = par)
    nir != 0.0 && (nir_flux[sector_idx] = nir)
    return DirectionalFluxes(ids, par_flux, nir_flux)
end

function _cache_supports_full_response(
    prepared::Union{Nothing,PreparedInterceptionData},
    ib::InterceptionBackend,
)
    prepared !== nothing || return false
    ib isa RasterCPUBackend || return false
    return true
end

function _cache_mode_for(
    prepared::Union{Nothing,PreparedInterceptionData},
    options::LightOptions,
    ib::InterceptionBackend,
    memory_limit_bytes::Int,
)
    _cache_supports_full_response(prepared, ib) || return :topology_fallback
    estimate = _estimate_light_cache_entry_bytes(prepared, options)
    estimate <= memory_limit_bytes || return :topology_fallback
    return options.all_in_turtle ? :full : :partial
end

function prepare_light_cache(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    memory_limit_bytes=nothing,
)
    ib = _resolve_interception_backend(interception_backend)
    prepared = ib isa RasterCPUBackend ? _prepare_interception_data(scene, models, options; include_budget_maps=true) : nothing
    interception_geometry =
        prepared === nothing ? _scene_geometry_for_interception(scene, models, options) : prepared.geometry
    render_geometry = _light_render_geometry(interception_geometry)
    node_metadata =
        options.store_node_metadata && scene.mtg !== nothing ?
        _build_light_node_metadata(scene, models, options, interception_geometry.node_ids) : nothing
    limit = memory_limit_bytes === nothing ? _default_light_cache_memory_limit() : Int(memory_limit_bytes)
    limit = max(limit, 1)
    estimated = prepared === nothing ? 0 : _estimate_light_cache_entry_bytes(prepared, options)
    mode = _cache_mode_for(prepared, options, ib, limit)
    return LightSimulationCache(
        scene,
        models,
        options,
        interception_backend,
        ib,
        scattering_mode,
        scattering_backend,
        prepared,
        render_geometry,
        node_metadata,
        mode,
        estimated,
        limit,
        0,
        0,
        Dict{UInt64,TurtleLightCacheEntry}(),
    )
end

function cache_summary(cache::LightSimulationCache)
    cached_sector_count = sum(length(entry.turtle.sectors) for entry in values(cache.entries); init=0)
    cached_full_response_sector_count =
        sum(
            count(x -> !isnothing(x), entry.par_added_per_sector) +
            count(x -> !isnothing(x), entry.nir_added_per_sector) +
            sum(count(x -> !isnothing(x), band_values) for band_values in values(entry.extra_added_per_sector); init=0)
            for entry in values(cache.entries);
            init=0,
        )
    return (
        mode=cache.mode,
        estimated_entry_bytes=cache.estimated_entry_bytes,
        resident_bytes=cache.resident_bytes,
        node_metadata_bytes=cache.node_metadata === nothing ? 0 : Base.summarysize(cache.node_metadata),
        node_metadata_count=cache.node_metadata === nothing ? 0 : length(cache.node_metadata.node_id),
        cached_turtle_count=length(cache.entries),
        cached_sector_count=cached_sector_count,
        cached_full_response_sector_count=cached_full_response_sector_count,
        memory_limit_bytes=cache.memory_limit_bytes,
        full_response_enabled=cache.mode != :topology_fallback,
    )
end

mutable struct LightSimulation
    scene::PlantGeom.SceneGeometry
    models::LightModels
    options::LightOptions
    interception_backend::Any
    scattering_mode::Symbol
    scattering_backend::Union{Nothing,ScatteringBackend}
    memory_limit_bytes::Any
    cache::Union{Nothing,LightSimulationCache}
    validation::Union{Nothing,ValidationReport}
end

"""
    LightSimulation(scene, models; options=LightOptions(), kwargs...)

Create a reusable light simulation. Expensive geometry preparation and radiation
caches are built lazily by [`run_light`](@ref).

Arguments:

- `scene`: prepared `PlantGeom.SceneGeometry` used by the solver.
- `models`: model specification accepted by [`prepare_models`](@ref).

Keywords:

- `options`: [`LightOptions`](@ref) controlling interception, scattering, and
  caching.
- `interception_backend`: interception backend selector or backend instance.
  The default is `:raster_cpu`.
- `scattering_mode`: scattering algorithm selector. The default is `:raycast`.
- `scattering_backend`: optional scattering backend instance.
- `memory_limit_bytes`: optional limit for resident directional-response cache
  data. `nothing` uses the package default.
"""
function LightSimulation(
    scene::PlantGeom.SceneGeometry,
    models;
    options::LightOptions=LightOptions(),
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    memory_limit_bytes=nothing,
)
    LightSimulation(
        scene,
        prepare_models(models),
        options,
        interception_backend,
        scattering_mode,
        scattering_backend,
        memory_limit_bytes,
        nothing,
        nothing,
    )
end

function _drop_light_cache!(sim::LightSimulation)
    sim.cache = nothing
    sim.validation = nothing
    return sim
end

function _ensure_light_cache!(sim::LightSimulation)
    sim.cache !== nothing && return sim.cache
    report = check_simulation(sim)
    if !isempty(report.errors)
        error(_validation_error_message(
            "Invalid light simulation",
            report;
            next="Run `check_simulation(sim)` and `summarize_scene(sim.scene; models=sim.models)` to inspect the scene/model inputs.",
        ))
    end
    sim.validation = report
    sim.cache = prepare_light_cache(
        sim.scene,
        sim.models,
        sim.options;
        interception_backend=sim.interception_backend,
        scattering_mode=sim.scattering_mode,
        scattering_backend=sim.scattering_backend,
        memory_limit_bytes=sim.memory_limit_bytes,
    )
    return sim.cache
end

function _validation_error_message(title::AbstractString, report::ValidationReport; next::AbstractString="")
    parts = String[string(title, ":")]
    append!(parts, ["- " * e for e in report.errors])
    if !isempty(report.warnings)
        push!(parts, "Warnings:")
        append!(parts, ["- " * w for w in report.warnings])
    end
    isempty(next) || push!(parts, next)
    return join(parts, "\n")
end

"""
    update_scene!(sim, new_scene)

Replace the scene and immediately release all prepared data tied to the old
scene. The next `run_light` call prepares the new scene lazily.

Arguments:

- `sim`: [`LightSimulation`](@ref) to update in place.
- `new_scene`: replacement `PlantGeom.SceneGeometry`.
"""
function update_scene!(sim::LightSimulation, new_scene::PlantGeom.SceneGeometry)
    _drop_light_cache!(sim)
    sim.scene = new_scene
    return sim
end

"""
    update_models!(sim, new_models)

Replace the model specification in `sim` and release all prepared cache data.
The next [`run_light`](@ref) call prepares the new models lazily.

Arguments:

- `sim`: [`LightSimulation`](@ref) to update in place.
- `new_models`: replacement model specification accepted by
  [`prepare_models`](@ref).
"""
function update_models!(sim::LightSimulation, new_models)
    _drop_light_cache!(sim)
    sim.models = prepare_models(new_models)
    return sim
end

"""
    update_options!(sim, new_options)

Replace the runtime options in `sim` and release all prepared cache data. The
next [`run_light`](@ref) call prepares data using the new options.

Arguments:

- `sim`: [`LightSimulation`](@ref) to update in place.
- `new_options`: replacement [`LightOptions`](@ref).
"""
function update_options!(sim::LightSimulation, new_options::LightOptions)
    _drop_light_cache!(sim)
    sim.options = new_options
    return sim
end

"""
    cache_summary(sim)

Return a named tuple summarizing the current radiation cache state for a
[`LightSimulation`](@ref).

Arguments:

- `sim`: simulation whose prepared cache should be summarized. If no cache has
  been prepared yet, the summary reports `mode=:unprepared`.
"""
function cache_summary(sim::LightSimulation)
    sim.cache === nothing && return (
        mode=:unprepared,
        estimated_entry_bytes=0,
        resident_bytes=0,
        node_metadata_bytes=0,
        node_metadata_count=0,
        cached_turtle_count=0,
        cached_sector_count=0,
        cached_full_response_sector_count=0,
        memory_limit_bytes=sim.memory_limit_bytes === nothing ? _default_light_cache_memory_limit() : Int(sim.memory_limit_bytes),
        full_response_enabled=false,
    )
    return cache_summary(sim.cache)
end

"""
    check_scene(scene)::ValidationReport

Validate basic scene readiness for light interception.

Arguments:

- `scene`: `PlantGeom.SceneGeometry` to check for geometry nodes, faces, and a
  valid xy domain.
"""
function check_scene(scene::PlantGeom.SceneGeometry)
    errors = String[]
    warnings = String[]
    infos = String[]
    isempty(scene.face2node) && push!(errors, "Scene has no geometry faces.")
    isempty(scene.nodes) && push!(errors, "Scene has no geometric nodes.")
    scene.scene_xy_bounds === nothing && push!(errors, "Scene has no XY domain. Pass `domain=` to `PlantGeom.make_scene` or `scene_xy_bounds=` to `PlantGeom.prepare_scene`.")
    if scene.scene_xy_bounds !== nothing
        xmin, ymin, xmax, ymax = scene.scene_xy_bounds
        (xmax > xmin && ymax > ymin) || push!(errors, "Scene XY domain must satisfy xmax > xmin and ymax > ymin.")
    end
    push!(infos, "$(length(scene.nodes)) geometric node(s) found.")
    return ValidationReport(errors, warnings, infos)
end

function _missing_model_pairs(scene::PlantGeom.SceneGeometry, models::LightModels)
    missing = Set{Tuple{String,String}}()
    ignored = _ignored_group_types(models)
    for nid in unique(scene.face2node)
        _is_ignored_node(nid, scene, ignored) && continue
        group = strip(_scene_group(scene, nid, ""))
        type_name = strip(_scene_type(scene, nid, ""))
        _type_model(models, group, type_name) === nothing && push!(missing, (group, type_name))
    end
    return sort!(collect(missing))
end

function _missing_models_snippet(missing::Vector{Tuple{String,String}})
    isempty(missing) && return ""
    grouped = OrderedDict{String,Vector{String}}()
    for (group, type_name) in missing
        push!(get!(grouped, group, String[]), type_name)
    end
    lines = ["Add model entries for these geometric nodes, for example:", "models = models_for("]
    for (group, types) in grouped
        push!(lines, "    $(repr(group)) => (")
        for type_name in types
            push!(lines, "        $(repr(type_name)) => translucent(par=0.15, nir=0.90),")
        end
        push!(lines, "    ),")
    end
    push!(lines, ")")
    return join(lines, "\n")
end

"""
    check_models(scene, models)::ValidationReport

Validate that `models` cover the geometric group/type pairs present in `scene`.

Arguments:

- `scene`: `PlantGeom.SceneGeometry` whose geometric nodes define the required
  group/type pairs.
- `models`: model specification accepted by [`prepare_models`](@ref).
"""
function check_models(scene::PlantGeom.SceneGeometry, models)
    lm = prepare_models(models)
    missing = _missing_model_pairs(scene, lm)
    errors = String[]
    if !isempty(missing)
        details = join(["group=$(repr(g)), type=$(repr(t))" for (g, t) in missing], "; ")
        push!(errors, "Missing model(s) for geometric node group/type pairs: $details.\n" * _missing_models_snippet(missing))
    end
    infos = ["Models are checked only for nodes with geometry."]
    return ValidationReport(errors, String[], infos)
end

"""
    summarize_scene(scene; models=nothing)

Return a `SceneSummary` describing the prepared scene domain, geometric nodes,
faces, group/type pairs, object ids, and missing model pairs.

Pass `models` to include a model coverage check in the summary.

Arguments:

- `scene`: `PlantGeom.SceneGeometry` to summarize.

Keywords:

- `models`: optional model specification used to report missing group/type
  coverage.
"""
function summarize_scene(scene::PlantGeom.SceneGeometry; models=nothing)
    buckets = Dict{Tuple{String,String},NamedTuple}()
    face_counts = Dict{Int,Int}()
    for nid in scene.face2node
        face_counts[nid] = get(face_counts, nid, 0) + 1
    end
    for (nid, node) in scene.nodes
        group = _scene_group(scene, nid, "")
        type_name = _scene_type(scene, nid, "")
        object_id = _scene_object_id(scene, nid, -1)
        key = (group, type_name)
        current = get(buckets, key, (nodes=0, faces=0, area=0.0, object_ids=Set{Int}()))
        object_ids = current.object_ids
        object_id > 0 && push!(object_ids, object_id)
        buckets[key] = (
            nodes=current.nodes + 1,
            faces=current.faces + get(face_counts, nid, 0),
            area=current.area + Float64(_scene_area(scene, nid, 0.0)),
            object_ids=object_ids,
        )
    end
    group_types = NamedTuple[]
    for ((group, type_name), data) in sort!(collect(buckets); by=x -> x[1])
        push!(
            group_types,
            (
                group=group,
                type=type_name,
                nodes=data.nodes,
                faces=data.faces,
                area=data.area,
                object_ids=sort!(collect(data.object_ids)),
            ),
        )
    end
    missing = models === nothing ? Tuple{String,String}[] : _missing_model_pairs(scene, prepare_models(models))
    warnings = String[]
    scene.scene_xy_bounds === nothing && push!(warnings, "Scene has no XY domain. Pass `domain=` to `PlantGeom.make_scene` or `scene_xy_bounds=` to `PlantGeom.prepare_scene`.")
    isempty(scene.nodes) && push!(warnings, "Scene has no geometric nodes.")
    domain = scene.scene_xy_bounds === nothing ? nothing : Tuple(Float64(v) for v in scene.scene_xy_bounds)
    object_ids = Set{Int}()
    for nid in keys(scene.nodes)
        object_id = _scene_object_id(scene, nid, -1)
        object_id > 0 && push!(object_ids, object_id)
    end
    return SceneSummary(domain, length(scene.nodes), length(scene.face2node), length(object_ids), group_types, missing, warnings)
end

function _meteo_rows_and_indices_for_check(
    meteo,
    options::LightOptions;
    apply_selection::Bool=true,
)
    table =
        meteo isa PlantMeteo.TimeStepTable ? meteo :
        Tables.istable(typeof(meteo)) ? _as_plantmeteo_table(meteo) : nothing
    table === nothing && return Any[meteo], [1]
    if apply_selection && !isempty(table)
        table, source_indices = _select_meteo_with_source_indices_local(table, options)
        return _meteo_rows(table), source_indices
    end
    return _meteo_rows(table), collect(1:length(table))
end

function _meteo_rows_for_check(
    meteo,
    options::LightOptions;
    apply_selection::Bool=true,
)
    rows, _ = _meteo_rows_and_indices_for_check(
        meteo,
        options;
        apply_selection=apply_selection,
    )
    return rows
end

"""
    check_meteo(meteo; options=LightOptions(), check_boundaries=options.check_meteo_boundaries)::ValidationReport

Validate meteo rows for required solar geometry, radiation inputs, and timestep
duration data.

Arguments:

- `meteo`: a meteo row, `PlantMeteo.TimeStepTable`, or Tables.jl-compatible
  table.

Keywords:

- `options`: [`LightOptions`](@ref) used for checks that depend on runtime
  meteo handling.
- `check_boundaries`: check physical ranges for used meteo values. Missing,
  nonfinite, derivability, and conflicting-input checks are always performed.
"""
function check_meteo(
    meteo;
    options::LightOptions=LightOptions(),
    check_boundaries::Bool=options.check_meteo_boundaries,
    _apply_selection::Bool=true,
    _row_indices::Union{Nothing,Vector{Int}}=nothing,
)
    rows = Any[]
    selected_indices = Int[]
    errors = String[]
    warnings = String[]
    infos = String[]
    try
        rows, selected_indices = _meteo_rows_and_indices_for_check(
            meteo,
            options;
            apply_selection=_apply_selection,
        )
        _row_indices === nothing || (selected_indices = _row_indices)
    catch err
        push!(errors, "Could not read meteo input as a table or row: $(sprint(showerror, err))")
        return ValidationReport(errors, warnings, infos)
    end
    isempty(rows) && push!(errors, "Meteo input has no rows.")
    isempty(rows) && return ValidationReport(errors, warnings, infos)

    length(selected_indices) == length(rows) || error("Internal meteo row-index mismatch.")
    for (i, row) in enumerate(rows)
        source_index = selected_indices[i]
        result = _resolve_meteo_step(
            row,
            options;
            row_index=source_index,
            check_boundaries=check_boundaries,
        )
        append!(errors, result.errors)
        append!(warnings, result.warnings)
        append!(infos, result.infos)
        if result.step === nothing && isempty(result.errors)
            push!(errors, "Meteo row $source_index could not be resolved for an unknown reason.")
        end
    end
    push!(infos, "$(length(rows)) meteo row(s) checked.")
    return ValidationReport(unique(errors), unique(warnings), unique(infos))
end

"""
    summarize_meteo(meteo; options=LightOptions(), check_boundaries=options.check_meteo_boundaries)

Return a `MeteoSummary` describing row count, columns, timestep duration,
radiation inputs, and the detected solar-geometry path for a meteo table or row.

Arguments:

- `meteo`: a meteo row, `PlantMeteo.TimeStepTable`, or Tables.jl-compatible
  table.

Keywords:

- `options`: [`LightOptions`](@ref) used for checks that depend on runtime
  meteo handling.
- `check_boundaries`: check physical ranges while resolving the summary.
"""
function summarize_meteo(
    meteo;
    options::LightOptions=LightOptions(),
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    warnings = String[]
    rows, source_indices = try
        _meteo_rows_and_indices_for_check(meteo, options)
    catch err
        push!(warnings, "Could not read meteo input as a table or row: $(sprint(showerror, err))")
        return MeteoSummary(0, Symbol[], nothing, false, String[], "unknown", warnings)
    end
    if isempty(rows)
        push!(warnings, "Meteo input has no rows.")
        return MeteoSummary(0, Symbol[], nothing, false, String[], "unknown", warnings)
    end
    first_row = first(rows)
    columns = Symbol[Symbol(name) for name in _row_propertynames(first_row)]
    first_resolution = _resolve_meteo_step(
        first_row,
        options;
        row_index=first(source_indices),
        check_boundaries=check_boundaries,
    )
    append!(warnings, first_resolution.errors)
    append!(warnings, first_resolution.warnings)
    radiation_inputs = String[]
    solar_geometry = "unknown"
    if first_resolution.step !== nothing && isempty(first_resolution.errors)
        step = first_resolution.step
        for (logical_name, label) in (
            (:ri_sw_f, "RI_SW_f"),
            (:ri_par_f, "RI_PAR_f"),
            (:ri_nir_f, "RI_NIR_f"),
            (:ri_sw_direct, "RI_SW_f_direct"),
            (:ri_sw_diffuse, "RI_SW_f_diffuse"),
            (:direct_fraction, "direct_fraction"),
            (:clearness, "clearness"),
        )
            source = get(step.sources, logical_name, nothing)
            source === nothing && continue
            source in getproperty(_METEO_ALIASES, logical_name) || continue
            push!(radiation_inputs, label)
        end
        solar_geometry =
            step.solar_geometry_source == :explicit ? "explicit sun_azimuth/sun_elevation" :
            step.solar_geometry_source == :zero_forcing_default ? "irrelevant for zero radiation" :
            "reconstructed from date/time and latitude"
    else
        has_sun =
            (:sun_azimuth in columns || :sun_azimut in columns) &&
            (:sun_elevation in columns)
        latitude = _row_value(first_row, [:latitude, :lat], NaN)
        solar_geometry =
            has_sun ? "explicit sun_azimuth/sun_elevation" :
            isfinite(latitude) ? "reconstructed from date/time and latitude" :
            "missing latitude or explicit sun position"
    end
    duration_values = Float64[]
    for (i, row) in enumerate(rows)
        source_index = source_indices[i]
        resolution = _resolve_meteo_step(
            row,
            options;
            row_index=source_index,
            check_boundaries=check_boundaries,
        )
        append!(warnings, resolution.errors)
        append!(warnings, resolution.warnings)
        resolution.step === nothing ||
            push!(duration_values, resolution.step.duration_seconds)
    end
    duration_seconds = isempty(duration_values) ? nothing : first(duration_values)
    variable_duration =
        !isempty(duration_values) && any(v -> !isapprox(v, first(duration_values); rtol=0.0, atol=1e-9), duration_values)
    report = check_meteo(
        meteo;
        options=options,
        check_boundaries=check_boundaries,
    )
    append!(warnings, report.errors)
    append!(warnings, report.warnings)
    return MeteoSummary(length(rows), columns, duration_seconds, variable_duration, radiation_inputs, solar_geometry, unique(warnings))
end

function _inferred_radiation_input_columns(row)
    inputs = String[]
    _has_any_column(row, [:clearness, :Kt]) && push!(inputs, "clearness")
    _has_any_column(row, [:RI_SW_f, :Ri_SW_f, :Rg, :rg, :sw_global, :global]) && push!(inputs, "RI_SW_f")
    _has_any_column(row, [:RI_PAR_f, :Ri_PAR_f, :PAR, :par]) && push!(inputs, "RI_PAR_f")
    _has_any_column(row, [:RI_NIR_f, :Ri_NIR_f, :NIR, :nir]) && push!(inputs, "RI_NIR_f")
    return inputs
end

"""
    check_simulation(sim)::ValidationReport
    check_simulation(scene, meteo; models, options=LightOptions(), check_boundaries=options.check_meteo_boundaries)::ValidationReport

Validate a reusable simulation or the separate inputs needed to build and run
one.

Arguments:

- `sim`: [`LightSimulation`](@ref) to validate.
- `scene`: `PlantGeom.SceneGeometry` to validate when checking separate inputs.
- `meteo`: meteo row or table to validate when checking separate inputs.

Keywords:

- `models`: required model specification when checking separate inputs.
- `options`: [`LightOptions`](@ref) used for meteo validation.
- `check_boundaries`: check physical ranges for used meteo values.
"""
function check_simulation(sim::LightSimulation)
    _merge_reports(check_scene(sim.scene), check_models(sim.scene, sim.models))
end

function check_simulation(
    scene::PlantGeom.SceneGeometry,
    meteo;
    models,
    options::LightOptions=LightOptions(),
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    _merge_reports(
        check_scene(scene),
        check_models(scene, models),
        check_meteo(
            meteo;
            options=options,
            check_boundaries=check_boundaries,
        ),
    )
end

@inline _series_progress_enabled(io::IO=stderr) = io isa Base.TTY

function _format_progress_seconds(seconds::Float64)
    if !isfinite(seconds)
        return "?"
    elseif seconds < 60
        return string(round(seconds; digits=1), "s")
    else
        minutes = floor(Int, seconds ÷ 60)
        rem_seconds = seconds - 60 * minutes
        return string(minutes, "m", lpad(string(round(rem_seconds; digits=1)), 4, '0'), "s")
    end
end

function _report_series_progress!(
    io::IO,
    cache::LightSimulationCache,
    i::Int,
    n::Int,
    started_ns::UInt64,
    last_report_ns::Base.RefValue{UInt64};
    force::Bool=false,
)
    _series_progress_enabled(io) || return nothing
    now_ns = time_ns()
    min_interval_ns = 2_000_000_000
    min_step_delta = max(cld(n, 100), 1)
    !force && i < n && (now_ns - last_report_ns[] < min_interval_ns) && (i % min_step_delta != 0) && return nothing
    elapsed = (now_ns - started_ns) / 1.0e9
    rate = elapsed > 0 ? i / elapsed : 0.0
    remaining = rate > 0 ? (n - i) / rate : Inf
    pct = n > 0 ? 100 * i / n : 100.0
    mode_label =
        cache.mode === :full ? "full-cache" :
        cache.mode === :partial ? "partial-cache" :
        "topology-fallback"
    print(
        io,
        "\r[ArchimedLight] run_light_series ",
        i,
        "/",
        n,
        " (",
        round(pct; digits=1),
        "%, ",
        mode_label,
        ", elapsed=",
        _format_progress_seconds(elapsed),
        ", eta=",
        _format_progress_seconds(remaining),
        ")",
    )
    if i == n
        print(io, "\n")
    end
    flush(io)
    last_report_ns[] = now_ns
    return nothing
end

function _touch_cache_entry!(cache::LightSimulationCache, entry::TurtleLightCacheEntry)
    cache.tick += 1
    entry.last_used_tick = cache.tick
    return entry
end

@inline function _cache_can_store_bytes(cache::LightSimulationCache, required_bytes::Int)
    required_bytes >= 0 || return false
    return Int128(cache.resident_bytes) + required_bytes <= cache.memory_limit_bytes
end

function _evict_cache_entries!(cache::LightSimulationCache, required_bytes::Int; protect_key::Union{Nothing,UInt64}=nothing)
    cache.mode == :partial || return nothing
    while cache.resident_bytes + required_bytes > cache.memory_limit_bytes && !isempty(cache.entries)
        victim_key = nothing
        victim_tick = typemax(Int)
        for (key, entry) in cache.entries
            key == protect_key && continue
            if entry.last_used_tick < victim_tick
                victim_key = key
                victim_tick = entry.last_used_tick
            end
        end
        victim_key === nothing && break
        victim = cache.entries[victim_key]
        cache.resident_bytes = max(cache.resident_bytes - victim.resident_bytes, 0)
        delete!(cache.entries, victim_key)
    end
    return nothing
end

function _build_turtle_cache_entry!(
    cache::LightSimulationCache,
    key::UInt64,
    turtle::TurtleGrid,
)
    prepared = cache.prepared
    prepared === nothing && error("Cannot build full-response cache entry without prepared interception data.")
    responses = _build_sector_responses(prepared, cache.scene, cache.models, turtle, cache.options)
    n = length(turtle.sectors)
    entry = TurtleLightCacheEntry(
        key,
        turtle,
        responses,
        fill(nothing, n),
        fill(nothing, n),
        Dict{String,Vector{Union{Nothing,Vector{Float64}}}}(),
        zeros(Int, n),
        zeros(Int, n),
        fill(true, n),
        fill(true, n),
        Dict{String,Vector{Int}}(),
        Dict{String,Vector{Bool}}(),
        0,
        0,
    )
    base_bytes = _turtle_cache_entry_retained_bytes(entry)
    entry.resident_bytes = base_bytes
    cache.mode == :partial && _evict_cache_entries!(cache, base_bytes)
    if !_cache_can_store_bytes(cache, base_bytes)
        # A built entry can exceed its pre-build estimate because it also owns
        # sparse active-index vectors, scattering topology, and container
        # storage. Finish this step with the transient response, then use the
        # uncached topology path later.
        cache.mode = :topology_fallback
        return _touch_cache_entry!(cache, entry)
    end
    cache.entries[key] = entry
    cache.resident_bytes += base_bytes
    return _touch_cache_entry!(cache, entry)
end

function _get_turtle_cache_entry!(
    cache::LightSimulationCache,
    turtle::TurtleGrid,
)
    key = _turtle_cache_key(turtle, cache.options)
    if haskey(cache.entries, key)
        return _touch_cache_entry!(cache, cache.entries[key])
    end
    return _build_turtle_cache_entry!(cache, key, turtle)
end

function _ensure_sector_band_cache!(
    cache::LightSimulationCache,
    entry::TurtleLightCacheEntry,
    sector_idx::Int,
    band::String,
)
    options = cache.options
    geometry = entry.responses_cache.prepared.geometry
    band_u = uppercase(band)
    is_extra_band = false
    if band_u == "PAR"
        target = entry.par_added_per_sector
        iterations = entry.par_iterations_per_sector
        converged = entry.par_converged_per_sector
    elseif band_u == "NIR"
        target = entry.nir_added_per_sector
        iterations = entry.nir_iterations_per_sector
        converged = entry.nir_converged_per_sector
    else
        is_extra_band = true
        target = get(entry.extra_added_per_sector, band_u, nothing)
        iterations = get(entry.extra_iterations_per_sector, band_u, nothing)
        converged = get(entry.extra_converged_per_sector, band_u, nothing)
    end
    existing = target === nothing ? nothing : target[sector_idx]
    existing !== nothing && return existing, iterations[sector_idx], converged[sector_idx]

    unit_fluxes =
        if band_u == "NIR"
            _unit_directional_fluxes(entry.turtle, sector_idx; nir=1.0)
        else
            _unit_directional_fluxes(entry.turtle, sector_idx; par=1.0)
        end
    # Cached unit responses must stay linear in the external sector flux.
    # Artificial-emitter scenes use a combined step-level scattering solve.
    first = _combine_sector_responses(
        entry.responses_cache,
        unit_fluxes;
        emitter_band_par=nothing,
        emitter_band_nir=nothing,
    )
    result =
        compute_scattering_band(
            entry.responses_cache.scattering_topology,
            first,
            options;
            mode=cache.scattering_mode,
            backend=cache.scattering_backend,
            band=band_u,
        )
    dense = _dense_vector_from_node_values(result.added_power_per_node, geometry)
    required_bytes = Base.summarysize(dense)
    cache.mode == :partial &&
        _evict_cache_entries!(cache, required_bytes; protect_key=entry.key)
    entry_is_resident = get(cache.entries, entry.key, nothing) === entry
    # The set of custom wavebands is unbounded and Julia Dict growth is not a
    # stable byte contract. Reuse cached geometry/topology for those bands, but
    # retain only the fixed PAR/NIR response families under the hard cache cap.
    retain_response =
        !is_extra_band &&
        entry_is_resident &&
        _cache_can_store_bytes(cache, required_bytes)
    if retain_response
        target[sector_idx] = dense
        iterations[sector_idx] = result.iterations
        converged[sector_idx] = result.converged
        entry.resident_bytes += required_bytes
        cache.resident_bytes += required_bytes
    end
    return dense, result.iterations, result.converged
end

function _assemble_cached_scattering(
    cache::LightSimulationCache,
    entry::TurtleLightCacheEntry,
    fluxes::DirectionalFluxes,
    nir_scattering::Bool,
)
    options = cache.options
    options.scattering || return nothing

    # The iterative stopping rule is relative to the combined initial field.
    # With emitters, summing independently truncated sky-sector and emitter
    # solves can therefore differ from an uncached combined solve. Reuse the
    # cached geometry/topology, but solve scattering once on the complete
    # first-order field to preserve cache-on/cache-off parity.
    if entry.responses_cache.emitter_transfer !== nothing
        first = _combine_sector_responses(entry.responses_cache, fluxes)
        return _compute_scattering_with_flags(
            cache.scene,
            cache.models,
            entry.turtle,
            first,
            options;
            mode=cache.scattering_mode,
            backend=cache.scattering_backend,
            nir_scattering=nir_scattering,
            responses_cache=entry.responses_cache,
        )
    end

    node_ids = entry.responses_cache.node_ids
    added_par = zeros(Float64, length(node_ids))
    added_nir = zeros(Float64, length(node_ids))
    iterations = 0
    converged = true
    for i in eachindex(entry.turtle.sectors)
        pf = fluxes.par[i]
        if pf != 0.0
            dense, it, conv = _ensure_sector_band_cache!(cache, entry, i, "PAR")
            _accumulate_scaled!(added_par, dense, pf)
            iterations = max(iterations, it)
            converged &= conv
        end
        if nir_scattering
            nf = fluxes.nir[i]
            if nf != 0.0
                dense, it, conv = _ensure_sector_band_cache!(cache, entry, i, "NIR")
                _accumulate_scaled!(added_nir, dense, nf)
                iterations = max(iterations, it)
                converged &= conv
            end
        end
    end

    return ScatteringResult(
        SpectralNodeValues(_float_dict_from_dense(node_ids, added_par), _float_dict_from_dense(node_ids, added_nir)),
        iterations,
        converged,
    )
end

function _compute_extra_band_light_cached(
    cache::LightSimulationCache,
    entry::TurtleLightCacheEntry,
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    resolved_step::ResolvedMeteoStep,
)
    extras_irr = _extra_band_irradiance(resolved_step, meteo_row, cache.options)
    _include_emitter_extra_bands!(
        extras_irr,
        entry.responses_cache.prepared.emitter_power_per_band_per_node,
    )
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()
    extra_emitter_escaped_power = Dict{String,Dict{Int,Float64}}()
    isempty(extras_irr) &&
        return extra_0_q, extra_q, extras_irr, extra_emitter_escaped_power

    node_ids = entry.responses_cache.node_ids
    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(
            total_irr,
            meteo_row,
            sky,
            turtle,
            cache.options;
            resolved_step=resolved_step,
        )
        first_band = _combine_sector_responses(
            entry.responses_cache,
            DirectionalFluxes(ids, flux_band, zeros(Float64, n));
            emitter_band_par=band,
            emitter_band_nir=nothing,
        )
        order0 = Dict{Int,Float64}(nid => v for (nid, v) in first_band.incident_power.par)
        added =
            if cache.options.scattering
                emitter_power = get(
                    entry.responses_cache.prepared.emitter_power_per_band_per_node,
                    uppercase(band),
                    Dict{Int,Float64}(),
                )
                if any(!iszero, values(emitter_power))
                    compute_scattering_band(
                        entry.responses_cache.scattering_topology,
                        first_band,
                        cache.options;
                        mode=cache.scattering_mode,
                        backend=cache.scattering_backend,
                        band=band,
                    ).added_power_per_node
                else
                    dense = zeros(Float64, length(node_ids))
                    iterations = 0
                    converged = true
                    for i in eachindex(flux_band)
                        f = flux_band[i]
                        f == 0.0 && continue
                        unit_dense, it, conv = _ensure_sector_band_cache!(cache, entry, i, band)
                        _accumulate_scaled!(dense, unit_dense, f)
                        iterations = max(iterations, it)
                        converged &= conv
                    end
                    _float_dict_from_dense(node_ids, dense)
                end
            else
                Dict{Int,Float64}()
            end
        total = Dict{Int,Float64}()
        for nid in union(keys(order0), keys(added))
            total[nid] = get(order0, nid, 0.0) + get(added, nid, 0.0)
        end
        extra_0_q[band] = order0
        extra_q[band] = total
        extra_emitter_escaped_power[band] =
            Dict{Int,Float64}(first_band.emitter_escaped_power.par)
    end
    return extra_0_q, extra_q, extras_irr, extra_emitter_escaped_power
end

function _run_light_step_cached(
    cache::LightSimulationCache,
    meteo_row;
    use_full_response::Bool=(cache.mode != :topology_fallback),
    check_boundaries::Bool=cache.options.check_meteo_boundaries,
    resolved_step::Union{Nothing,ResolvedMeteoStep}=nothing,
)
    scene = cache.scene
    models = cache.models
    options = cache.options
    ib = cache.resolved_interception_backend
    nir_interception = _nir_interception_enabled_local(options)
    nir_scattering = _nir_scattering_enabled_local(options) && nir_interception
    resolved =
        resolved_step === nothing ?
        _resolved_meteo_step_or_error(
            meteo_row,
            options;
            check_boundaries=check_boundaries,
        ) : resolved_step
    sky = compute_sky(
        meteo_row,
        options;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
    nir_interception || (sky = _disable_nir_sky_local(sky))
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(
        meteo_row,
        sky,
        turtle,
        options;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )

    prepared = cache.prepared
    responses_cache = nothing
    scattering_topology = nothing
    extra_irr = _extra_band_irradiance(resolved, meteo_row, options)
    first = nothing
    scat = nothing
    if use_full_response && cache.mode != :topology_fallback
        entry = _get_turtle_cache_entry!(cache, turtle)
        responses_cache = entry.responses_cache
        first = _combine_sector_responses(responses_cache, fluxes)
        scat = _assemble_cached_scattering(cache, entry, fluxes, nir_scattering)
        extra_0_q, extra_q, extra_irr, extra_emitter_escaped_power =
            _compute_extra_band_light_cached(
                cache,
                entry,
                meteo_row,
                sky,
                turtle,
                resolved,
            )
    else
        if ib isa RasterCPUBackend && options.scattering && isempty(extra_irr)
            prepared === nothing && (prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true))
            first, scattering_topology =
                _stream_first_order_with_scattering_topology(prepared, scene, models, turtle, fluxes, options)
        elseif ib isa RasterCPUBackend && (options.scattering || !isempty(extra_irr))
            prepared === nothing && (prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true))
            responses_cache = _build_sector_responses(prepared, scene, models, turtle, options)
            first = _combine_sector_responses(responses_cache, fluxes)
        else
            first = compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
        end
        scat = _compute_scattering_with_flags(
            scene,
            models,
            turtle,
            first,
            options;
            mode=cache.scattering_mode,
            backend=cache.scattering_backend,
            nir_scattering=nir_scattering,
            responses_cache=responses_cache,
            scattering_topology=scattering_topology,
        )
        extra_0_q, extra_q, extra_irr, extra_emitter_escaped_power =
            _compute_extra_band_light(
                scene,
                models,
                meteo_row,
                sky,
                turtle,
                options;
                interception_backend=ib,
                scattering_mode=cache.scattering_mode,
                scattering_backend=cache.scattering_backend,
                responses_cache=responses_cache,
                resolved_step=resolved,
            )
    end
    nir_interception || (first = _disable_nir_first_order_local(first))
    budget = integrate_light(
        scene,
        models,
        first,
        scat,
        options;
        step_duration_seconds=resolved.duration_seconds,
        extra_initial_energy_per_band=extra_0_q,
        extra_energy_per_band=extra_q,
        extra_emitter_escaped_power_per_band=extra_emitter_escaped_power,
        component_area_per_node=prepared === nothing ? nothing : prepared.component_area_per_node,
        absorption_par_per_node=prepared === nothing ? nothing : prepared.absorption_par_per_node,
        absorption_nir_per_node=prepared === nothing ? nothing : prepared.absorption_nir_per_node,
    )
    sky_fraction =
        options.include_sky_fraction ?
        _compute_sky_fraction(
            scene,
            models,
            turtle,
            options;
            prepared=prepared,
            responses_cache=responses_cache,
        ) : nothing
    return LightStepResult(
        sky,
        turtle,
        fluxes,
        first,
        scat,
        budget,
        extra_irr,
        sky_fraction,
        cache.render_geometry,
        cache.node_metadata,
    )
end

function _run_light_sky_cached(
    cache::LightSimulationCache,
    sky::SkyState;
    step_duration_seconds::Real,
    use_full_response::Bool=(cache.mode != :topology_fallback),
)
    scene = cache.scene
    models = cache.models
    options = cache.options
    ib = cache.resolved_interception_backend
    nir_interception = _nir_interception_enabled_local(options)
    nir_scattering = _nir_scattering_enabled_local(options) && nir_interception
    nir_interception || (sky = _disable_nir_sky_local(sky))
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(sky, turtle, options)

    prepared = cache.prepared
    responses_cache = nothing
    scattering_topology = nothing
    first = nothing
    scat = nothing
    extra_irr = Dict{String,Float64}()
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()
    if use_full_response && cache.mode != :topology_fallback
        entry = _get_turtle_cache_entry!(cache, turtle)
        responses_cache = entry.responses_cache
        first = _combine_sector_responses(responses_cache, fluxes)
        scat = _assemble_cached_scattering(cache, entry, fluxes, nir_scattering)
    else
        if ib isa RasterCPUBackend && options.scattering
            prepared === nothing && (prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true))
            first, scattering_topology =
                _stream_first_order_with_scattering_topology(prepared, scene, models, turtle, fluxes, options)
        else
            first = compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
        end
        scat = _compute_scattering_with_flags(
            scene,
            models,
            turtle,
            first,
            options;
            mode=cache.scattering_mode,
            backend=cache.scattering_backend,
            nir_scattering=nir_scattering,
            responses_cache=responses_cache,
            scattering_topology=scattering_topology,
        )
    end
    extra_0_q, extra_q, extra_irr, extra_emitter_escaped_power =
        _compute_extra_band_light(
            scene,
            models,
            nothing,
            sky,
            turtle,
            options;
            interception_backend=ib,
            scattering_mode=cache.scattering_mode,
            scattering_backend=cache.scattering_backend,
            responses_cache=responses_cache,
        )
    nir_interception || (first = _disable_nir_first_order_local(first))
    budget = integrate_light(
        scene,
        models,
        first,
        scat,
        options;
        step_duration_seconds=Float64(step_duration_seconds),
        extra_initial_energy_per_band=extra_0_q,
        extra_energy_per_band=extra_q,
        extra_emitter_escaped_power_per_band=extra_emitter_escaped_power,
        component_area_per_node=prepared === nothing ? nothing : prepared.component_area_per_node,
        absorption_par_per_node=prepared === nothing ? nothing : prepared.absorption_par_per_node,
        absorption_nir_per_node=prepared === nothing ? nothing : prepared.absorption_nir_per_node,
    )
    return LightStepResult(
        sky,
        turtle,
        fluxes,
        first,
        scat,
        budget,
        extra_irr,
        nothing,
        cache.render_geometry,
        cache.node_metadata,
    )
end

function _use_full_response_for_sim(sim::LightSimulation, cache::LightSimulationCache)
    sim.options.cache_radiation && cache.mode != :topology_fallback
end

function _is_meteo_series_input(x)
    x isa PlantMeteo.TimeStepTable && return true
    if x isa AbstractVector
        (isempty(x) && eltype(x) <: NamedTuple) && return true
        (!isempty(x) && all(row -> row isa NamedTuple, x)) && return true
        Tables.istable(typeof(x)) && return true
        return false
    end
    if x isa NamedTuple
        return !isempty(x) && all(column -> column isa AbstractVector, values(x))
    end
    Tables.istable(typeof(x)) && return true
    return false
end

function _run_light_series_resolved(
    cache::LightSimulationCache,
    rows_eff::PlantMeteo.TimeStepTable,
    resolved_steps::Vector{ResolvedMeteoStep};
    use_full_response::Bool=(cache.mode != :topology_fallback),
    check_boundaries::Bool=cache.options.check_meteo_boundaries,
)
    length(rows_eff) == length(resolved_steps) ||
        error("Internal meteo/resolved-step length mismatch.")
    out = Vector{LightStepResult}(undef, length(rows_eff))
    io = stderr
    started_ns = time_ns()
    last_report_ns = Ref(started_ns)
    if !isempty(rows_eff)
        _report_series_progress!(
            io,
            cache,
            0,
            length(rows_eff),
            started_ns,
            last_report_ns;
            force=true,
        )
    end
    for (i, row) in enumerate(rows_eff)
        out[i] = _run_light_step_cached(
            cache,
            row;
            use_full_response=use_full_response,
            check_boundaries=check_boundaries,
            resolved_step=resolved_steps[i],
        )
        _report_series_progress!(
            io,
            cache,
            i,
            length(rows_eff),
            started_ns,
            last_report_ns;
            force=(i == length(rows_eff)),
        )
    end
    return out
end

function _run_light_series_sim(
    sim::LightSimulation,
    meteo;
    check_boundaries::Bool=sim.options.check_meteo_boundaries,
)
    meteo_local = meteo isa PlantMeteo.TimeStepTable ? meteo : _as_plantmeteo_table(meteo)
    rows_eff, resolved_steps = _prepare_meteo_for_series(
        meteo_local,
        sim.options;
        check_boundaries=check_boundaries,
        return_resolved=true,
    )
    cache = _ensure_light_cache!(sim)
    return _run_light_series_resolved(
        cache,
        rows_eff,
        resolved_steps;
        use_full_response=_use_full_response_for_sim(sim, cache),
        check_boundaries=check_boundaries,
    )
end

"""
    run_light(sim, meteo_or_row; check_boundaries=sim.options.check_meteo_boundaries)

Run one light step for a meteo row, or a full series for a meteo table.

Arguments:

- `sim`: [`LightSimulation`](@ref) containing the scene, models, options, and
  lazy cache.
- `meteo_or_row`: either a single meteo row for one step, or a
  `PlantMeteo.TimeStepTable`/Tables.jl-compatible table for a series.

Keywords:

- `check_boundaries`: check physical ranges for used meteo values. Derivability
  and conflicting-input checks are always performed.
"""
function run_light(
    sim::LightSimulation,
    meteo_or_row;
    check_boundaries::Bool=sim.options.check_meteo_boundaries,
)
    if _is_meteo_series_input(meteo_or_row)
        return _run_light_series_sim(
            sim,
            meteo_or_row;
            check_boundaries=check_boundaries,
        )
    end
    resolved = _resolved_meteo_step_or_error(
        meteo_or_row,
        sim.options;
        check_boundaries=check_boundaries,
    )
    cache = _ensure_light_cache!(sim)
    return _run_light_step_cached(
        cache,
        meteo_or_row;
        use_full_response=_use_full_response_for_sim(sim, cache),
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
end

"""
    run_light(sim, sky::SkyState; step_duration_seconds)

Run one light step from an already computed sky state. The step duration is
required because there is no meteo row from which to infer it.

Arguments:

- `sim`: [`LightSimulation`](@ref) containing the scene, models, options, and
  lazy cache.
- `sky`: precomputed [`SkyState`](@ref) used for one light step.

Keywords:

- `step_duration_seconds`: duration of the step in seconds. This is required
  for energy integration.
"""
function run_light(sim::LightSimulation, sky::SkyState; step_duration_seconds=nothing)
    if step_duration_seconds === nothing
        error("run_light(sim, sky::SkyState) requires `step_duration_seconds`, for example `run_light(sim, sky; step_duration_seconds=1800.0)`.")
    end
    cache = _ensure_light_cache!(sim)
    return _run_light_sky_cached(
        cache,
        sky;
        step_duration_seconds=step_duration_seconds,
        use_full_response=_use_full_response_for_sim(sim, cache),
    )
end

function _sky_series_durations(skies::AbstractVector{<:SkyState}, step_duration_seconds)
    if step_duration_seconds === nothing
        error(
            "run_light(sim, skies::AbstractVector{<:SkyState}) requires `step_duration_seconds`; " *
            "pass one duration for every state or a vector with length $(length(skies)).",
        )
    elseif step_duration_seconds isa Real
        duration = _duration_seconds_strict(
            step_duration_seconds;
            field_name="step_duration_seconds",
        )
        return fill(duration, length(skies))
    elseif step_duration_seconds isa AbstractVector
        length(step_duration_seconds) == length(skies) || throw(
            ArgumentError(
                "`step_duration_seconds` has length $(length(step_duration_seconds)), " *
                "but the SkyState vector has length $(length(skies)).",
            ),
        )
        durations = Vector{Float64}(undef, length(step_duration_seconds))
        for (i, duration) in enumerate(step_duration_seconds)
            duration isa Real || throw(
                ArgumentError(
                    "`step_duration_seconds[$i]` must be a real number of seconds, " *
                    "got $(typeof(duration)).",
                ),
            )
            durations[i] = _duration_seconds_strict(
                duration;
                field_name="step_duration_seconds[$i]",
            )
        end
        return durations
    end
    throw(
        ArgumentError(
            "`step_duration_seconds` must be a real number or a vector of real numbers, " *
            "got $(typeof(step_duration_seconds)).",
        ),
    )
end

"""
    run_light(sim, skies::AbstractVector{<:SkyState}; step_duration_seconds)

Run a series of already computed sky states. The step duration is required
because sky states do not contain timing metadata.

`step_duration_seconds` may be either one positive duration shared by every
state or a vector of positive durations with the same length as `skies`. The
results are returned as a `Vector{LightStepResult}` in input order, while the
simulation's prepared scene and radiation cache are reused across the series.

# Examples

```julia
series = run_light(sim, skies; step_duration_seconds=1800.0)
series = run_light(sim, skies; step_duration_seconds=[900.0, 1800.0, 900.0])
```
"""
function run_light(
    sim::LightSimulation,
    skies::AbstractVector{<:SkyState};
    step_duration_seconds=nothing,
)
    durations = _sky_series_durations(skies, step_duration_seconds)
    cache = _ensure_light_cache!(sim)
    use_full_response = _use_full_response_for_sim(sim, cache)
    out = Vector{LightStepResult}(undef, length(skies))
    io = stderr
    started_ns = time_ns()
    last_report_ns = Ref(started_ns)
    if !isempty(skies)
        _report_series_progress!(
            io,
            cache,
            0,
            length(skies),
            started_ns,
            last_report_ns;
            force=true,
        )
    end
    for (i, (sky, duration)) in enumerate(zip(skies, durations))
        out[i] = _run_light_sky_cached(
            cache,
            sky;
            step_duration_seconds=duration,
            use_full_response=use_full_response,
        )
        _report_series_progress!(
            io,
            cache,
            i,
            length(skies),
            started_ns,
            last_report_ns;
            force=(i == length(skies)),
        )
    end
    return out
end

"""
    run_light_step(scene, models, meteo_row, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::LightStepResult

Run a complete light computation for one meteo row:
`compute_sky -> build_turtle -> compute_directional_fluxes -> compute_first_order -> compute_scattering -> integrate_light`.
Set `LightOptions(include_sky_fraction=true)` to store the per-node sky-view
fraction needed by downstream MTG attachment or coupled energy-balance models.
When using `read_config`, this option is enabled by requesting `sky_fraction`
in `component_variables` or `opf_variables`.
"""
function run_light_step(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    meteo_row,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    resolved = _resolved_meteo_step_or_error(
        meteo_row,
        options;
        check_boundaries=check_boundaries,
    )
    cache = prepare_light_cache(
        scene,
        models,
        options;
        interception_backend=interception_backend,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
        memory_limit_bytes=0,
    )
    return _run_light_step_cached(
        cache,
        meteo_row;
        use_full_response=false,
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
end

function run_light_step(
    cache::LightSimulationCache,
    meteo_row;
    check_boundaries::Bool=cache.options.check_meteo_boundaries,
)
    resolved = _resolved_meteo_step_or_error(
        meteo_row,
        cache.options;
        check_boundaries=check_boundaries,
    )
    return _run_light_step_cached(
        cache,
        meteo_row;
        check_boundaries=check_boundaries,
        resolved_step=resolved,
    )
end

"""
    run_light_series(scene, models, meteo, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::Vector{LightStepResult}

Run the complete light pipeline for all rows in a `PlantMeteo.TimeStepTable`, with optional directional
response reuse when `LightOptions(cache_radiation=true)` is enabled.
Set `LightOptions(include_sky_fraction=true)` to store `sky_fraction` in each step.
When using `read_config`, this option is enabled by requesting `sky_fraction`
in `component_variables` or `opf_variables`.
"""
function run_light_series(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    rows_eff, resolved_steps = _prepare_meteo_for_series(
        meteo,
        options;
        check_boundaries=check_boundaries,
        return_resolved=true,
    )
    if _can_use_series_radiation_cache(_resolve_interception_backend(interception_backend)) && options.cache_radiation
        cache = prepare_light_cache(
            scene,
            models,
            options;
            interception_backend=interception_backend,
            scattering_mode=scattering_mode,
            scattering_backend=scattering_backend,
        )
        return _run_light_series_resolved(
            cache,
            rows_eff,
            resolved_steps;
            check_boundaries=check_boundaries,
        )
    end

    cache = prepare_light_cache(
        scene,
        models,
        options;
        interception_backend=interception_backend,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
        memory_limit_bytes=0,
    )
    return _run_light_series_resolved(
        cache,
        rows_eff,
        resolved_steps;
        check_boundaries=check_boundaries,
    )
end

function run_light_series(
    cache::LightSimulationCache,
    meteo::PlantMeteo.TimeStepTable;
    check_boundaries::Bool=cache.options.check_meteo_boundaries,
)
    rows_eff, resolved_steps = _prepare_meteo_for_series(
        meteo,
        cache.options;
        check_boundaries=check_boundaries,
        return_resolved=true,
    )
    return _run_light_series_resolved(
        cache,
        rows_eff,
        resolved_steps;
        check_boundaries=check_boundaries,
    )
end

function run_light_series(
    cache::LightSimulationCache,
    meteo;
    check_boundaries::Bool=cache.options.check_meteo_boundaries,
)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected PlantMeteo.TimeStepTable or a Tables.jl-compatible table.")
    meteo isa PlantMeteo.TimeStepTable &&
        return run_light_series(cache, meteo; check_boundaries=check_boundaries)
    return run_light_series(
        cache,
        _as_plantmeteo_table(meteo);
        check_boundaries=check_boundaries,
    )
end

function run_light_series(
    scene::PlantGeom.SceneGeometry,
    models::LightModels,
    meteo,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected PlantMeteo.TimeStepTable or a Tables.jl-compatible table.")
    meteo isa PlantMeteo.TimeStepTable &&
        return run_light_series(
            scene,
            models,
            meteo,
            options;
            interception_backend=interception_backend,
            scattering_mode=scattering_mode,
            scattering_backend=scattering_backend,
            check_boundaries=check_boundaries,
        )
    run_light_series(
        scene,
        models,
        _as_plantmeteo_table(meteo),
        options;
        interception_backend=interception_backend,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
        check_boundaries=check_boundaries,
    )
end
