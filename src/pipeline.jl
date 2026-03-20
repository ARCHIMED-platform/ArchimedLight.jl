"""
    integrate_light(scene, models, first, scat, options; meteo_row=nothing, extra_initial_energy_per_band=..., extra_energy_per_band=..., step_duration_seconds=nothing, component_area_per_node=nothing, absorption_par_per_node=nothing, absorption_nir_per_node=nothing)::LightBudget

Combine first-order and scattering results into per-node irradiance (`*_f`, W m^-2)
and energy (`*_q`, J component^-1 timestep^-1) budgets, including absorbed light.
"""
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

struct SectorResponsesCache
    prepared::PreparedInterceptionData
    projections::Union{Nothing,Vector{DirectionProjectionResult}}
    projected_area_per_sector::Vector{Dict{Int,Float64}}
    hits_per_sector::Vector{Dict{Int,Int}}
    node_ids::Vector{Int}
    emitter_weights::Dict{Tuple{Int,Int},Float64}
    scattering_topology::Union{Nothing,ScatteringTopologyCache}
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
    pa_by_sector = Vector{Dict{Int,Float64}}(undef, n)
    hits_by_sector = Vector{Dict{Int,Int}}(undef, n)
    keep_projections = options.scattering || !isempty(prepared.emitter_nodes)
    projections = keep_projections ? Vector{DirectionProjectionResult}(undef, n) : nothing
    for i in 1:n
        projection = _prepared_direction_projection(prepared, turtle.sectors[i].direction, options)
        pa_by_sector[i] = _visible_area_from_projection(projection, options, geometry.plotbox, prepared.virtual_nodes)
        hits_by_sector[i] = projection.node_hits
        keep_projections && (projections[i] = projection)
    end
    emitter_weights =
        projections === nothing ? Dict{Tuple{Int,Int},Float64}() :
        _emitter_transfer_weights_from_projections(projections, turtle, prepared.emitter_nodes)
    scattering_topology =
        if options.scattering && projections !== nothing
            _build_scattering_topology_cache(scene, models, prepared, turtle, projections)
        else
            nothing
        end
    return SectorResponsesCache(prepared, projections, pa_by_sector, hits_by_sector, geometry.node_ids, emitter_weights, scattering_topology)
end

function _build_sector_responses(scene::SceneGeometry, models::LightModels, turtle::TurtleGrid, options::LightOptions)
    prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true)
    return _build_sector_responses(prepared, scene, models, turtle, options)
end

function _combine_sector_responses(
    responses::SectorResponsesCache,
    fluxes::DirectionalFluxes,
)
    node_ids = Set{Int}()
    for k in responses.node_ids
        push!(node_ids, k)
    end
    for d in responses.projected_area_per_sector
        for k in keys(d)
            push!(node_ids, k)
        end
    end

    projected_area_per_node = Dict{Int,Float64}(id => 0.0 for id in node_ids)
    incident_power = _zero_spectral_node_values(node_ids)
    hits_per_node = Dict{Int,Int}(id => 0 for id in node_ids)

    for i in eachindex(responses.projected_area_per_sector)
        pf = fluxes.par[i]
        nf = fluxes.nir[i]
        active_flux = (pf != 0.0) || (nf != 0.0)

        for (nid, pa) in responses.projected_area_per_sector[i]
            if active_flux
                projected_area_per_node[nid] = get(projected_area_per_node, nid, 0.0) + pa
                if pf != 0.0
                    incident_power.par[nid] = get(incident_power.par, nid, 0.0) + pf * pa
                end
                if nf != 0.0
                    incident_power.nir[nid] = get(incident_power.nir, nid, 0.0) + nf * pa
                end
            end
        end
        for (nid, h) in responses.hits_per_sector[i]
            hits_per_node[nid] = get(hits_per_node, nid, 0) + h
        end
    end

    for ((to, src), w) in responses.emitter_weights
        incident_power.par[to] = get(incident_power.par, to, 0.0) + w * get(responses.prepared.emitter_par_power_per_node, src, 0.0)
        incident_power.nir[to] = get(incident_power.nir, to, 0.0) + w * get(responses.prepared.emitter_nir_power_per_node, src, 0.0)
    end

    FirstOrderResult(projected_area_per_node, incident_power, hits_per_node)
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
        return PlantMeteo.positive_duration_seconds(getproperty(row, :step_duration); field_name="step_duration")
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
            start_cols=(:hour_start, :hour),
            end_cols=(:hour_end,),
            duration_cols=(:step_duration,),
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
        start_cols=(:hour_start, :hour),
        end_cols=(:hour_end,),
        duration_cols=(:step_duration,),
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
            start_cols=(:hour_start, :hour),
            end_cols=(:hour_end,),
            duration_cols=(:step_duration,),
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

function _prepare_meteo_rows_for_series(meteo::MeteoTable, options::LightOptions)
    rows = collect(meteo.rows)
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
)
    options.scattering || return nothing
    if nir_scattering
        if responses_cache !== nothing && responses_cache.scattering_topology !== nothing
            return compute_scattering(responses_cache.scattering_topology, first, options; mode=mode, backend=backend)
        end
        return compute_scattering(scene, models, turtle, first, options; mode=mode, backend=backend)
    end

    par_only =
        if responses_cache !== nothing && responses_cache.scattering_topology !== nothing
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
    ib = _resolve_interception_backend(interception_backend)
    nir_interception = _nir_interception_enabled_local(options)
    nir_scattering = _nir_scattering_enabled_local(options) && nir_interception
    sky = compute_sky(meteo_row, options)
    nir_interception || (sky = _disable_nir_sky_local(sky))
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(meteo_row, sky, turtle, options)
    responses_cache = nothing
    if ib isa RasterCPUBackend && (options.scattering || !isempty(_extra_band_irradiance(meteo_row)))
        prepared = _prepare_interception_data(scene, models, options; include_budget_maps=true)
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
        mode=scattering_mode,
        backend=scattering_backend,
        nir_scattering=nir_scattering,
        responses_cache=responses_cache,
    )
    extra_0_q, extra_q, extra_irr = _compute_extra_band_light(
        scene,
        models,
        meteo_row,
        sky,
        turtle,
        options;
        interception_backend=ib,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
        responses_cache=responses_cache,
    )
    budget = integrate_light(
        scene,
        models,
        first,
        scat,
        options;
        meteo_row=meteo_row,
        extra_initial_energy_per_band=extra_0_q,
        extra_energy_per_band=extra_q,
        component_area_per_node=responses_cache === nothing ? nothing : responses_cache.prepared.component_area_per_node,
        absorption_par_per_node=responses_cache === nothing ? nothing : responses_cache.prepared.absorption_par_per_node,
        absorption_nir_per_node=responses_cache === nothing ? nothing : responses_cache.prepared.absorption_nir_per_node,
    )
    LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
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
    ib = _resolve_interception_backend(interception_backend)
    nir_interception = _nir_interception_enabled_local(options)
    nir_scattering = _nir_scattering_enabled_local(options) && nir_interception
    use_cache = _can_use_series_radiation_cache(ib) && (options.cache_radiation || options.scattering)
    cache = Dict{UInt64,SectorResponsesCache}()
    prepared = ib isa RasterCPUBackend ? _prepare_interception_data(scene, models, options; include_budget_maps=true) : nothing
    rows_eff = _prepare_meteo_rows_for_series(meteo, options)
    out = Vector{LightStepResult}(undef, length(rows_eff))
    for i in eachindex(rows_eff)
        row = rows_eff[i]
        sky = compute_sky(row, options)
        nir_interception || (sky = _disable_nir_sky_local(sky))
        turtle = build_turtle(options, sky)
        fluxes = compute_directional_fluxes(row, sky, turtle, options)

        responses = nothing
        first =
            if use_cache
                key = _turtle_cache_key(turtle, options)
                responses = get!(cache, key) do
                    _build_sector_responses(prepared, scene, models, turtle, options)
                end
                _combine_sector_responses(responses, fluxes)
            else
                prepared === nothing ? compute_first_order(scene, models, turtle, fluxes, options; backend=ib) : _compute_first_order(prepared, turtle, fluxes, options)
            end

        scat = _compute_scattering_with_flags(
            scene,
            models,
            turtle,
            first,
            options;
            mode=scattering_mode,
            backend=scattering_backend,
            nir_scattering=nir_scattering,
            responses_cache=responses,
        )
        extra_0_q, extra_q, extra_irr = _compute_extra_band_light(
            scene,
            models,
            row,
            sky,
            turtle,
            options;
            interception_backend=ib,
            scattering_mode=scattering_mode,
            scattering_backend=scattering_backend,
            responses_cache=responses,
        )
        budget = integrate_light(
            scene,
            models,
            first,
            scat,
            options;
            meteo_row=row,
            extra_initial_energy_per_band=extra_0_q,
            extra_energy_per_band=extra_q,
            component_area_per_node=prepared === nothing ? nothing : prepared.component_area_per_node,
            absorption_par_per_node=prepared === nothing ? nothing : prepared.absorption_par_per_node,
            absorption_nir_per_node=prepared === nothing ? nothing : prepared.absorption_nir_per_node,
        )
        out[i] = LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
    end
    out
end
