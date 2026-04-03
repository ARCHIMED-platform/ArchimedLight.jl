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
    integrate_light(scene, models, first, scat, options; meteo_row=nothing, extra_initial_energy_per_band=..., extra_energy_per_band=..., step_duration_seconds=nothing, component_area_per_node=nothing, absorption_par_per_node=nothing, absorption_nir_per_node=nothing)::LightBudget

Combine first-order interception and scattering into per-node incident and
absorbed light budgets.

The result stores both irradiance-style outputs (`*_f`, W m^-2) and energy
outputs (`*_q`, J component^-1 timestep^-1), plus optional extra-waveband
energies when they were carried through the pipeline.
"""
function integrate_light(
    scene::SceneGeometry,
    models::LightModels,
    first::FirstOrderResult,
    scat::Union{Nothing,ScatteringResult},
    options::LightOptions;
    meteo_row=nothing,
    extra_initial_energy_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    extra_energy_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    step_duration_seconds::Union{Nothing,Float64}=nothing,
    component_area_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_par_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_nir_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
)
    dt_seconds = step_duration_seconds === nothing ? _step_duration_seconds_local(meteo_row) : Float64(step_duration_seconds)
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

    return LightBudget(
        budget.incident_flux,
        budget.incident_energy,
        budget.absorbed_flux,
        budget.absorbed_energy,
        _scale_extra_band_energy(extra_initial_energy_per_band, dt_seconds),
        _scale_extra_band_energy(extra_energy_per_band, dt_seconds),
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

function _dense_emitter_incident_power(
    edge_counts::Union{Nothing,Dict{UInt64,Int}},
    total_from::Union{Nothing,Dict{Int,Int}},
    prepared::PreparedInterceptionData,
)
    geometry = prepared.geometry
    par = zeros(Float64, length(geometry.node_ids))
    nir = zeros(Float64, length(geometry.node_ids))
    if edge_counts !== nothing && total_from !== nothing
        for (edge, count) in edge_counts
            src = _unpack_emitter_from(edge)
            n = get(total_from, src, 0)
            n > 0 || continue
            w = count / n
            to = _unpack_emitter_to(edge)
            idx = geometry.node_index[to]
            src_idx = geometry.node_index[src]
            par[idx] += w * prepared.emitter_par_power_by_index[src_idx]
            nir[idx] += w * prepared.emitter_nir_power_by_index[src_idx]
        end
    end
    return _dense_node_map(par, geometry), _dense_node_map(nir, geometry)
end

function _turtle_cache_key(turtle::TurtleGrid, options::LightOptions)
    h = hash((length(turtle.sectors), options.pixel_size, options.area_ratio))
    for s in turtle.sectors
        d = s.direction
        h = hash((round(d[1], digits=8), round(d[2], digits=8), round(d[3], digits=8), s.source), h)
    end
    h
end

function _build_sector_responses(
    prepared::PreparedInterceptionData,
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    options::LightOptions,
)
    n = length(turtle.sectors)
    geometry = prepared.geometry
    pa_by_sector = Vector{DenseNodeMap{Float64}}(undef, n)
    hits_by_sector = Vector{DenseNodeMap{Int}}(undef, n)
    emitter_edge_counts = isempty(prepared.emitter_nodes) ? nothing : Dict{UInt64,Int}()
    emitter_total_from = isempty(prepared.emitter_nodes) ? nothing : Dict{Int,Int}()
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
                    geometry,
                ),
                geometry,
            )
        hits_by_sector[i] = _dense_node_map(copy(_dense_projection_hits(projection, geometry)), geometry)
        if emitter_edge_counts !== nothing
            if projection isa DenseDirectionProjectionResult
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_node_mask,
                    geometry.node_ids;
                    stacks_sorted=true,
                )
            else
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_nodes;
                    stacks_sorted=true,
                )
            end
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
    emitter_incident_power_par, emitter_incident_power_nir =
        _dense_emitter_incident_power(emitter_edge_counts, emitter_total_from, prepared)
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
        scattering_topology,
    )
end

function _build_sector_responses(scene::SceneGeometry, models::LightModels, turtle::TurtleGrid, options::LightOptions)
    prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true)
    return _build_sector_responses(prepared, scene, models, turtle, options)
end

function _combine_sector_responses(
    responses::SectorResponsesCache,
    fluxes::DirectionalFluxes,
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

    emitter_par = responses.emitter_incident_power_par.values
    @inbounds for idx in _active_indices(responses.emitter_incident_power_par)
        incident_power_par[idx] += emitter_par[idx]
    end
    emitter_nir = responses.emitter_incident_power_nir.values
    @inbounds for idx in _active_indices(responses.emitter_incident_power_nir)
        incident_power_nir[idx] += emitter_nir[idx]
    end

    return FirstOrderResult(
        _all_dense_float_node_map(node_ids, projected_area_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(node_ids, incident_power_par),
            _all_dense_float_node_map(node_ids, incident_power_nir),
        ),
        _all_dense_int_node_map(node_ids, hits_per_node),
    )
end

function _stream_first_order_with_scattering_topology(
    prepared::PreparedInterceptionData,
    scene::SceneGeometry,
    models::LightModels,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    options::LightOptions,
)
    geometry = prepared.geometry
    projected_area_per_node = zeros(Float64, length(geometry.node_ids))
    incident_power_par = zeros(Float64, length(geometry.node_ids))
    incident_power_nir = zeros(Float64, length(geometry.node_ids))
    hits_per_node = zeros(Int, length(geometry.node_ids))
    scattering_edge_counts = Dict{UInt64,Int}()
    scattering_sun_hits = Dict{Int,Int}()
    scattering_scratch = ScatteringStackScratch()
    emitter_edge_counts = isempty(prepared.emitter_nodes) ? nothing : Dict{UInt64,Int}()
    emitter_total_from = isempty(prepared.emitter_nodes) ? nothing : Dict{Int,Int}()

    for (k, sector) in enumerate(turtle.sectors)
        projection = _prepared_direction_projection(prepared, sector.direction, options)
        visible_area =
            _visible_area_from_projection_dense(
                projection,
                options,
                geometry.plotbox,
                prepared.virtual_node_mask,
                geometry,
            )

        _accumulate_projection_hits!(hits_per_node, projection, geometry)

        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if par_flux != 0.0 || nir_flux != 0.0
            @inbounds for idx in eachindex(visible_area)
                pa = visible_area[idx]
                pa <= 0.0 && continue
                projected_area_per_node[idx] += pa
                par_flux != 0.0 && (incident_power_par[idx] += par_flux * pa)
                nir_flux != 0.0 && (incident_power_nir[idx] += nir_flux * pa)
            end
        end

        if emitter_edge_counts !== nothing
            if projection isa DenseDirectionProjectionResult
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_node_mask,
                    geometry.node_ids;
                    stacks_sorted=true,
                )
            else
                _accumulate_emitter_transfer_counts!(
                    emitter_edge_counts,
                    emitter_total_from,
                    projection,
                    prepared.emitter_nodes;
                    stacks_sorted=true,
                )
            end
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

    if emitter_edge_counts !== nothing
        emitter_incident_power_par, emitter_incident_power_nir =
            _dense_emitter_incident_power(emitter_edge_counts, emitter_total_from, prepared)
        @inbounds for idx in _active_indices(emitter_incident_power_par)
            incident_power_par[idx] += emitter_incident_power_par.values[idx]
        end
        @inbounds for idx in _active_indices(emitter_incident_power_nir)
            incident_power_nir[idx] += emitter_incident_power_nir.values[idx]
        end
    end

    first = FirstOrderResult(
        _all_dense_float_node_map(geometry.node_ids, projected_area_per_node),
        SpectralNodeValues(
            _all_dense_float_node_map(geometry.node_ids, incident_power_par),
            _all_dense_float_node_map(geometry.node_ids, incident_power_nir),
        ),
        _all_dense_int_node_map(geometry.node_ids, hits_per_node),
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
    if name in propertynames(row)
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
    names = propertynames(row)
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
    try
        return PlantMeteo.row_datetime_interval(
            row;
            index=index,
            date_cols=(:date,),
            start_cols=(:hour_start, :hour, :date),
            end_cols=(:hour_end,),
            duration_cols=(:step_duration, :duration),
            default_date=Dates.Date(2000, 1, 1),
            default_duration_seconds=1.0,
            allow_end_rollover=false,
        )
    catch err
        msg = sprint(showerror, err)
        occursin("end is before start", msg) && error("end is before start at meteo row $(index)")
        rethrow(err)
    end
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
        sky.ri_sw_f,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
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

function _apply_meteo_range_local(rows::Vector{<:NamedTuple}, options::LightOptions)
    spec = _cfg_meteo_range_spec_local(options)
    spec === nothing && return rows

    parts = split(spec, ","; limit=2)
    length(parts) == 2 || error("invalid meteo_range format: $(repr(spec))")
    a = strip(parts[1])
    b = strip(parts[2])
    isempty(a) && error("invalid meteo_range format: $(repr(spec))")
    isempty(b) && error("invalid meteo_range format: $(repr(spec))")

    ia = tryparse(Int, a)
    ib = tryparse(Int, b)
    if ia !== nothing && ib !== nothing
        n = length(rows)
        ia >= 1 || error("invalid meteo_range: start step must be >= 1")
        ib >= ia || error("invalid meteo_range: end step is before start step")
        ib <= n || error("invalid meteo_range: end step exceeds meteo size")
        return rows[ia:ib]
    end

    t0 = _parse_range_datetime_token_local(a)
    t1 = _parse_range_datetime_token_local(b)
    t1 >= t0 || error("invalid meteo_range: end datetime is before start datetime")

    out = PlantMeteo.select_overlapping_timesteps(
        rows,
        t0,
        t1;
        closed=true, # Java uses closed-interval overlap semantics.
        date_cols=(:date,),
        start_cols=(:hour_start, :hour, :date),
        end_cols=(:hour_end,),
        duration_cols=(:step_duration, :duration),
        default_date=Dates.Date(2000, 1, 1),
        default_duration_seconds=1.0,
        allow_end_rollover=false,
    )
    isempty(out) && error("invalid meteo_range: selection is empty")
    return out
end

function _apply_meteo_active_filter_local(rows::Vector{<:NamedTuple})
    isempty(rows) && return rows
    names = propertynames(first(rows))
    :active in names || return rows

    out = NamedTuple[]
    for row in rows
        flag = _parse_bool_strict_local(getproperty(row, :active), "active")
        flag && push!(out, row)
    end
    isempty(out) && error("invalid meteo: no active meteo step")
    return out
end

function _validate_meteo_sequence_local(rows::Vector{<:NamedTuple})
    try
        PlantMeteo.check_non_overlapping_timesteps(
            rows;
            date_cols=(:date,),
            start_cols=(:hour_start, :hour, :date),
            end_cols=(:hour_end,),
            duration_cols=(:step_duration, :duration),
            default_date=Dates.Date(2000, 1, 1),
            default_duration_seconds=1.0,
            allow_end_rollover=false,
        )
    catch err
        msg = sprint(showerror, err)
        if occursin("overlapping timesteps at row", msg)
            m = match(r"row\s+(\d+)", msg)
            i = m === nothing ? 0 : parse(Int, m.captures[1])
            error("invalid overlapping meteo steps at row $(i)")
        end
        rethrow(err)
    end
end

_meteo_rows(meteo::MeteoTable) = collect(meteo.rows)
function _meteo_rows(meteo::PlantMeteo.TimeStepTable)
    meta = _meta_to_namedtuple(getfield(meteo, :metadata))
    raw_rows = _rows_to_namedtuples(meteo)
    [_namedtuple_with_meta(r, meta) for r in raw_rows]
end

_meteo_metadata(meteo::MeteoTable) = meteo.metadata
_meteo_metadata(meteo::PlantMeteo.TimeStepTable) = _meta_to_namedtuple(getfield(meteo, :metadata))

function _prepare_meteo_rows_for_series(meteo::MeteoTable, options::LightOptions)
    rows = _meteo_rows(meteo)
    isempty(rows) && return rows
    _validate_meteo_sequence_local(rows)
    rows = _apply_meteo_range_local(rows, options)
    rows = _apply_meteo_active_filter_local(rows)
    rows
end

function _prepare_meteo_rows_for_series(meteo::PlantMeteo.TimeStepTable, options::LightOptions)
    rows = _meteo_rows(meteo)
    isempty(rows) && return rows
    _validate_meteo_sequence_local(rows)
    rows = _apply_meteo_range_local(rows, options)
    rows = _apply_meteo_active_filter_local(rows)
    rows
end

"""
    prepare_meteo(meteo, options)::MeteoTable

Return the effective meteo table after Java-like meteo controls are applied:
sequence validation, optional `meteo_range`, and optional `active` filtering.
"""
function prepare_meteo(meteo::MeteoTable, options::LightOptions)
    MeteoTable(_prepare_meteo_rows_for_series(meteo, options), meteo.metadata)
end

function prepare_meteo(meteo::PlantMeteo.TimeStepTable, options::LightOptions)
    PlantMeteo.TimeStepTable(_prepare_meteo_rows_for_series(meteo, options), _meteo_metadata(meteo))
end

function prepare_meteo(meteo, options::LightOptions)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected MeteoTable, PlantMeteo.TimeStepTable, or a Tables.jl-compatible table.")
    meteo isa MeteoTable && return prepare_meteo(meteo, options)
    meteo isa PlantMeteo.TimeStepTable && return prepare_meteo(meteo, options)
    prepare_meteo(_as_plantmeteo_table(meteo), options)
end

function _default_scattering_factor_local(options::LightOptions, band::String)
    b = uppercase(band)
    b == "NIR" && return options.scattering_coeff_nir
    return options.scattering_coeff_par
end

function _compute_scattering_with_flags(
    scene::SceneGeometry,
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

function _interception_area_per_node_local(scene::SceneGeometry, models::LightModels, options::LightOptions)
    geometry = _scene_geometry_for_interception(scene, models, options)
    return _interception_area_per_node_from_geometry(geometry)
end

function _node_absorptance_per_band(scene::SceneGeometry, models::LightModels, options::LightOptions, band::String)
    geometry = _scene_geometry_for_interception(scene, models, options)
    virtual_nodes = _virtual_sensor_node_ids(geometry.node_group, geometry.node_type, models)
    return _node_absorptance_per_band_from_geometry(scene, models, options, geometry, virtual_nodes, band)
end

function _extra_band_irradiance(meteo_row)
    extras = Dict{String,Float64}()
    for p in propertynames(meteo_row)
        s = String(p)
        su = uppercase(s)
        startswith(su, "RI_") || continue
        endswith(su, "_F") || continue
        su in ("RI_PAR_F", "RI_NIR_F", "RI_SW_F") && continue
        band = su[4:(end - 2)]
        isempty(band) && continue
        v = _row_number_local(meteo_row, p, NaN)
        isfinite(v) || continue
        extras[band] = max(v, 0.0)
    end
    extras
end

function _single_band_flux(total_irradiance::Float64, meteo_row, sky::SkyState, turtle::TurtleGrid, options::LightOptions)
    tmp = SkyState(
        sky.sun_azimuth_deg,
        sky.sun_elevation_deg,
        total_irradiance,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
    )
    compute_directional_fluxes(meteo_row, tmp, turtle, options).par
end

function _compute_extra_band_light(
    scene::SceneGeometry,
    models::LightModels,
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    options::LightOptions;
    interception_backend::InterceptionBackend=RasterCPUBackend(),
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    responses_cache::Union{Nothing,SectorResponsesCache}=nothing,
)
    extras_irr = _extra_band_irradiance(meteo_row)
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()

    isempty(extras_irr) && return extra_0_q, extra_q, extras_irr

    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(total_irr, meteo_row, sky, turtle, options)
        first_band =
            if responses_cache === nothing
                compute_first_order(
                    scene,
                    models,
                    turtle,
                    DirectionalFluxes(ids, flux_band, zeros(Float64, n)),
                    options;
                    backend=interception_backend,
                )
            else
                _combine_sector_responses(responses_cache, DirectionalFluxes(ids, flux_band, zeros(Float64, n)))
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
    end
    return extra_0_q, extra_q, extras_irr
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
    scene::SceneGeometry
    models::LightModels
    options::LightOptions
    interception_backend::Any
    resolved_interception_backend::InterceptionBackend
    scattering_mode::Symbol
    scattering_backend::Union{Nothing,ScatteringBackend}
    prepared::Union{Nothing,PreparedInterceptionData}
    render_geometry::LightRenderGeometry
    mode::Symbol
    estimated_entry_bytes::Int
    memory_limit_bytes::Int
    resident_bytes::Int
    tick::Int
    entries::Dict{UInt64,TurtleLightCacheEntry}
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
    return n_sectors * (base_sector_bytes + scatter_sector_bytes)
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
    return isempty(prepared.emitter_nodes)
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
    scene::SceneGeometry,
    models::LightModels,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
    memory_limit_bytes=nothing,
)
    ib = _resolve_interception_backend(interception_backend)
    prepared = ib isa RasterCPUBackend ? _prepare_interception_data(scene, models, options; include_budget_maps=true) : nothing
    render_geometry =
        if prepared === nothing
            _light_render_geometry(_scene_geometry_for_interception(scene, models, options))
        else
            _light_render_geometry(prepared.geometry)
        end
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
        cached_turtle_count=length(cache.entries),
        cached_sector_count=cached_sector_count,
        cached_full_response_sector_count=cached_full_response_sector_count,
        memory_limit_bytes=cache.memory_limit_bytes,
        full_response_enabled=cache.mode != :topology_fallback,
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
    if force || i == n
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
    base_bytes =
        length(prepared.geometry.node_ids) * n * (sizeof(Float64) + sizeof(Int))
    cache.mode == :partial && _evict_cache_entries!(cache, base_bytes)
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
        base_bytes,
        0,
    )
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
    n_nodes = length(geometry.node_ids)
    band_u = uppercase(band)
    if band_u == "PAR"
        target = entry.par_added_per_sector
        iterations = entry.par_iterations_per_sector
        converged = entry.par_converged_per_sector
    elseif band_u == "NIR"
        target = entry.nir_added_per_sector
        iterations = entry.nir_iterations_per_sector
        converged = entry.nir_converged_per_sector
    else
        target = get!(entry.extra_added_per_sector, band_u) do
            fill(nothing, length(entry.turtle.sectors))
        end
        iterations = get!(entry.extra_iterations_per_sector, band_u) do
            zeros(Int, length(entry.turtle.sectors))
        end
        converged = get!(entry.extra_converged_per_sector, band_u) do
            fill(true, length(entry.turtle.sectors))
        end
    end
    existing = target[sector_idx]
    existing !== nothing && return existing, iterations[sector_idx], converged[sector_idx]

    cache.mode == :partial && _evict_cache_entries!(cache, n_nodes * sizeof(Float64); protect_key=entry.key)
    unit_fluxes =
        if band_u == "NIR"
            _unit_directional_fluxes(entry.turtle, sector_idx; nir=1.0)
        else
            _unit_directional_fluxes(entry.turtle, sector_idx; par=1.0)
        end
    first = _combine_sector_responses(entry.responses_cache, unit_fluxes)
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
    target[sector_idx] = dense
    iterations[sector_idx] = result.iterations
    converged[sector_idx] = result.converged
    entry.resident_bytes += n_nodes * sizeof(Float64)
    cache.resident_bytes += n_nodes * sizeof(Float64)
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
)
    extras_irr = _extra_band_irradiance(meteo_row)
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()
    isempty(extras_irr) && return extra_0_q, extra_q, extras_irr

    node_ids = entry.responses_cache.node_ids
    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(total_irr, meteo_row, sky, turtle, cache.options)
        first_band = _combine_sector_responses(entry.responses_cache, DirectionalFluxes(ids, flux_band, zeros(Float64, n)))
        order0 = Dict{Int,Float64}(nid => v for (nid, v) in first_band.incident_power.par)
        added =
            if cache.options.scattering
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
            else
                Dict{Int,Float64}()
            end
        total = Dict{Int,Float64}()
        for nid in union(keys(order0), keys(added))
            total[nid] = get(order0, nid, 0.0) + get(added, nid, 0.0)
        end
        extra_0_q[band] = order0
        extra_q[band] = total
    end
    return extra_0_q, extra_q, extras_irr
end

function _run_light_step_cached(
    cache::LightSimulationCache,
    meteo_row;
    use_full_response::Bool=(cache.mode != :topology_fallback),
)
    scene = cache.scene
    models = cache.models
    options = cache.options
    ib = cache.resolved_interception_backend
    nir_interception = _nir_interception_enabled_local(options)
    nir_scattering = _nir_scattering_enabled_local(options) && nir_interception
    sky = compute_sky(meteo_row, options)
    nir_interception || (sky = _disable_nir_sky_local(sky))
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(meteo_row, sky, turtle, options)

    prepared = cache.prepared
    responses_cache = nothing
    scattering_topology = nothing
    extra_irr = _extra_band_irradiance(meteo_row)
    first = nothing
    scat = nothing
    if use_full_response && cache.mode != :topology_fallback
        entry = _get_turtle_cache_entry!(cache, turtle)
        responses_cache = entry.responses_cache
        first = _combine_sector_responses(responses_cache, fluxes)
        scat = _assemble_cached_scattering(cache, entry, fluxes, nir_scattering)
        extra_0_q, extra_q, extra_irr = _compute_extra_band_light_cached(cache, entry, meteo_row, sky, turtle)
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
        extra_0_q, extra_q, extra_irr = _compute_extra_band_light(
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
        )
    end
    budget = integrate_light(
        scene,
        models,
        first,
        scat,
        options;
        meteo_row=meteo_row,
        extra_initial_energy_per_band=extra_0_q,
        extra_energy_per_band=extra_q,
        component_area_per_node=prepared === nothing ? nothing : prepared.component_area_per_node,
        absorption_par_per_node=prepared === nothing ? nothing : prepared.absorption_par_per_node,
        absorption_nir_per_node=prepared === nothing ? nothing : prepared.absorption_nir_per_node,
    )
    return LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr, cache.render_geometry)
end

"""
    run_light_step(scene, models, meteo_row, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::LightStepResult

Run a complete light computation for one meteo row:
`compute_sky -> build_turtle -> compute_directional_fluxes -> compute_first_order -> compute_scattering -> integrate_light`.
"""
function run_light_step(
    scene::SceneGeometry,
    models::LightModels,
    meteo_row,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
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
    return _run_light_step_cached(cache, meteo_row; use_full_response=false)
end

function run_light_step(
    cache::LightSimulationCache,
    meteo_row,
)
    return _run_light_step_cached(cache, meteo_row)
end

"""
    run_light_series(scene, models, meteo, options; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::Vector{LightStepResult}

Run the complete light pipeline for all rows in a `MeteoTable`, with optional directional
response reuse when `LightOptions(cache_radiation=true)` is enabled.
"""
function run_light_series(
    scene::SceneGeometry,
    models::LightModels,
    meteo::MeteoTable,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
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
        return run_light_series(cache, meteo)
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
    return run_light_series(cache, meteo)
end

function run_light_series(
    cache::LightSimulationCache,
    meteo::MeteoTable,
)
    rows_eff = _prepare_meteo_rows_for_series(meteo, cache.options)
    out = Vector{LightStepResult}(undef, length(rows_eff))
    io = stderr
    started_ns = time_ns()
    last_report_ns = Ref(started_ns)
    if !isempty(rows_eff)
        _report_series_progress!(io, cache, 0, length(rows_eff), started_ns, last_report_ns; force=true)
    end
    for i in eachindex(rows_eff)
        out[i] = _run_light_step_cached(cache, rows_eff[i])
        _report_series_progress!(io, cache, i, length(rows_eff), started_ns, last_report_ns; force=(i == length(rows_eff)))
    end
    return out
end

function run_light_series(
    cache::LightSimulationCache,
    meteo::PlantMeteo.TimeStepTable,
)
    meteo_local = MeteoTable(_meteo_rows(meteo), _meteo_metadata(meteo))
    return run_light_series(cache, meteo_local)
end

function run_light_series(
    cache::LightSimulationCache,
    meteo,
)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected MeteoTable, PlantMeteo.TimeStepTable, or a Tables.jl-compatible table.")
    meteo isa MeteoTable && return run_light_series(cache, meteo)
    meteo isa PlantMeteo.TimeStepTable && return run_light_series(cache, meteo)
    return run_light_series(cache, _as_plantmeteo_table(meteo))
end

function run_light_series(
    scene::SceneGeometry,
    models::LightModels,
    meteo::PlantMeteo.TimeStepTable,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
)
    meteo_local = MeteoTable(_meteo_rows(meteo), _meteo_metadata(meteo))
    run_light_series(
        scene,
        models,
        meteo_local,
        options;
        interception_backend=interception_backend,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
    )
end

function run_light_series(
    scene::SceneGeometry,
    models::LightModels,
    meteo,
    options::LightOptions;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
)
    Tables.istable(typeof(meteo)) || error("Unsupported meteo input: expected MeteoTable, PlantMeteo.TimeStepTable, or a Tables.jl-compatible table.")
    meteo isa MeteoTable && return run_light_series(scene, models, meteo, options; interception_backend=interception_backend, scattering_mode=scattering_mode, scattering_backend=scattering_backend)
    meteo isa PlantMeteo.TimeStepTable && return run_light_series(scene, models, meteo, options; interception_backend=interception_backend, scattering_mode=scattering_mode, scattering_backend=scattering_backend)
    run_light_series(scene, models, _as_plantmeteo_table(meteo), options; interception_backend=interception_backend, scattering_mode=scattering_mode, scattering_backend=scattering_backend)
end
