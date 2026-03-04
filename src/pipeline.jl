"""
    integrate_light(first, scat, cfg; extra_0_q_per_band=..., extra_q_per_band=..., step_duration_seconds=1.0, component_area_per_node=nothing, absorption_par_per_node=nothing, absorption_nir_per_node=nothing)::LightBudget

Combine first-order and scattering results into per-node irradiance (`*_f`, W m^-2)
and energy (`*_q`, J component^-1 timestep^-1) budgets, including absorbed light.
"""
function integrate_light(
    first::FirstOrderResult,
    scat::Union{Nothing,ScatteringResult},
    cfg::LightConfig;
    extra_0_q_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    extra_q_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    step_duration_seconds::Float64=1.0,
    component_area_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_par_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
    absorption_nir_per_node::Union{Nothing,Dict{Int,Float64}}=nothing,
)
    ri_par_0_f_per_node = Dict{Int,Float64}()
    ri_nir_0_f_per_node = Dict{Int,Float64}()
    ri_par_f_per_node = Dict{Int,Float64}()
    ri_nir_f_per_node = Dict{Int,Float64}()
    ri_par_0_q_per_node = Dict{Int,Float64}()
    ri_nir_0_q_per_node = Dict{Int,Float64}()
    ri_par_q_per_node = Dict{Int,Float64}()
    ri_nir_q_per_node = Dict{Int,Float64}()
    ra_par_0_f_per_node = Dict{Int,Float64}()
    ra_nir_0_f_per_node = Dict{Int,Float64}()
    ra_par_f_per_node = Dict{Int,Float64}()
    ra_nir_f_per_node = Dict{Int,Float64}()
    ra_par_0_q_per_node = Dict{Int,Float64}()
    ra_nir_0_q_per_node = Dict{Int,Float64}()
    ra_par_q_per_node = Dict{Int,Float64}()
    ra_nir_q_per_node = Dict{Int,Float64}()

    node_ids = collect(keys(first.projected_area_per_node))
    default_abs_par = clamp(1.0 - cfg.scattering_coeff_par, 0.0, 1.0)
    default_abs_nir = clamp(1.0 - cfg.scattering_coeff_nir, 0.0, 1.0)
    for nid in node_ids
        pa =
            if component_area_per_node === nothing
                get(first.projected_area_per_node, nid, 0.0)
            else
                get(component_area_per_node, nid, get(first.projected_area_per_node, nid, 0.0))
            end
        pa = max(pa, eps(Float64))
        p0 = get(first.incident_par_power_per_node, nid, 0.0)
        n0 = get(first.incident_nir_power_per_node, nid, 0.0)
        ps = scat === nothing ? 0.0 : get(scat.added_par_power_per_node, nid, 0.0)
        ns = scat === nothing ? 0.0 : get(scat.added_nir_power_per_node, nid, 0.0)
        abs_par =
            if absorption_par_per_node === nothing
                default_abs_par
            else
                clamp(get(absorption_par_per_node, nid, default_abs_par), 0.0, 1.0)
            end
        abs_nir =
            if absorption_nir_per_node === nothing
                default_abs_nir
            else
                clamp(get(absorption_nir_per_node, nid, default_abs_nir), 0.0, 1.0)
            end

        ri_par_0_f = p0 / pa
        ri_nir_0_f = n0 / pa
        ri_par_f = (p0 + ps) / pa
        ri_nir_f = (n0 + ns) / pa

        ri_par_0_f_per_node[nid] = ri_par_0_f
        ri_nir_0_f_per_node[nid] = ri_nir_0_f
        ri_par_f_per_node[nid] = ri_par_f
        ri_nir_f_per_node[nid] = ri_nir_f

        ri_par_0_q_per_node[nid] = p0 * step_duration_seconds
        ri_nir_0_q_per_node[nid] = n0 * step_duration_seconds
        ri_par_q_per_node[nid] = (p0 + ps) * step_duration_seconds
        ri_nir_q_per_node[nid] = (n0 + ns) * step_duration_seconds

        ra_par_0_f_per_node[nid] = ri_par_0_f * abs_par
        ra_nir_0_f_per_node[nid] = ri_nir_0_f * abs_nir
        ra_par_f_per_node[nid] = ri_par_f * abs_par
        ra_nir_f_per_node[nid] = ri_nir_f * abs_nir
        ra_par_0_q_per_node[nid] = p0 * abs_par * step_duration_seconds
        ra_nir_0_q_per_node[nid] = n0 * abs_nir * step_duration_seconds
        ra_par_q_per_node[nid] = (p0 + ps) * abs_par * step_duration_seconds
        ra_nir_q_per_node[nid] = (n0 + ns) * abs_nir * step_duration_seconds
    end

    scaled_extra_0_q = Dict{String,Dict{Int,Float64}}()
    for (band, vals) in extra_0_q_per_band
        scaled_extra_0_q[band] = Dict{Int,Float64}(nid => v * step_duration_seconds for (nid, v) in vals)
    end

    scaled_extra_q = Dict{String,Dict{Int,Float64}}()
    for (band, vals) in extra_q_per_band
        scaled_extra_q[band] = Dict{Int,Float64}(nid => v * step_duration_seconds for (nid, v) in vals)
    end

    LightBudget(
        ri_par_0_f_per_node,
        ri_nir_0_f_per_node,
        ri_par_f_per_node,
        ri_nir_f_per_node,
        ri_par_0_q_per_node,
        ri_nir_0_q_per_node,
        ri_par_q_per_node,
        ri_nir_q_per_node,
        ra_par_0_f_per_node,
        ra_nir_0_f_per_node,
        ra_par_f_per_node,
        ra_nir_f_per_node,
        ra_par_0_q_per_node,
        ra_nir_0_q_per_node,
        ra_par_q_per_node,
        ra_nir_q_per_node,
        scaled_extra_0_q,
        scaled_extra_q,
    )
end

function _turtle_cache_key(turtle::TurtleGrid, cfg::LightConfig)
    h = hash((length(turtle.sectors), cfg.pixel_size, cfg.area_ratio))
    for s in turtle.sectors
        d = s.direction
        h = hash((round(d[1], digits=8), round(d[2], digits=8), round(d[3], digits=8), s.source), h)
    end
    h
end

function _build_sector_responses(scene::SceneGeometry, turtle::TurtleGrid, cfg::LightConfig)
    n = length(turtle.sectors)
    vertices, faces, face2node, node_ids, plotbox, _ = _scene_geometry_for_interception(scene, cfg)
    cache_ctx = _projection_cache_context(vertices, faces, face2node, plotbox, cfg)
    upper_hit = _use_upper_hit_pixel_table(cfg)
    pa_by_sector = Vector{Dict{Int,Float64}}(undef, n)
    hits_by_sector = Vector{Dict{Int,Int}}(undef, n)
    for i in 1:n
        pa_by_sector[i], hits_by_sector[i] = _rasterize_direction_java(
            vertices,
            faces,
            face2node,
            turtle.sectors[i].direction,
            cfg,
            plotbox,
            cache_ctx=cache_ctx,
            upper_hit=upper_hit,
        )
    end
    emit_par, emit_nir = _emitter_power_per_node(scene, cfg)
    emitter_nodes = Set(union(keys(emit_par), keys(emit_nir)))
    emitter_weights =
        isempty(emitter_nodes) ? Dict{Tuple{Int,Int},Float64}() :
        _emitter_transfer_weights(vertices, faces, face2node, turtle, cfg, plotbox, emitter_nodes, cache_ctx)
    return pa_by_sector, hits_by_sector, node_ids, emit_par, emit_nir, emitter_weights
end

function _combine_sector_responses(
    pa_by_sector,
    hits_by_sector,
    fluxes::DirectionalFluxes,
    node_ids_all,
    emit_par::Dict{Int,Float64}=Dict{Int,Float64}(),
    emit_nir::Dict{Int,Float64}=Dict{Int,Float64}(),
    emitter_weights::Dict{Tuple{Int,Int},Float64}=Dict{Tuple{Int,Int},Float64}(),
)
    node_ids = Set{Int}()
    for k in node_ids_all
        push!(node_ids, k)
    end
    for d in pa_by_sector
        for k in keys(d)
            push!(node_ids, k)
        end
    end

    projected_area_per_node = Dict{Int,Float64}(id => 0.0 for id in node_ids)
    incident_par_power_per_node = Dict{Int,Float64}(id => 0.0 for id in node_ids)
    incident_nir_power_per_node = Dict{Int,Float64}(id => 0.0 for id in node_ids)
    hits_per_node = Dict{Int,Int}(id => 0 for id in node_ids)

    for i in eachindex(pa_by_sector)
        pf = fluxes.par[i]
        nf = fluxes.nir[i]
        active_flux = (pf != 0.0) || (nf != 0.0)

        for (nid, pa) in pa_by_sector[i]
            if active_flux
                projected_area_per_node[nid] = get(projected_area_per_node, nid, 0.0) + pa
                if pf != 0.0
                    incident_par_power_per_node[nid] = get(incident_par_power_per_node, nid, 0.0) + pf * pa
                end
                if nf != 0.0
                    incident_nir_power_per_node[nid] = get(incident_nir_power_per_node, nid, 0.0) + nf * pa
                end
            end
        end
        for (nid, h) in hits_by_sector[i]
            hits_per_node[nid] = get(hits_per_node, nid, 0) + h
        end
    end

    for ((to, src), w) in emitter_weights
        incident_par_power_per_node[to] = get(incident_par_power_per_node, to, 0.0) + w * get(emit_par, src, 0.0)
        incident_nir_power_per_node[to] = get(incident_nir_power_per_node, to, 0.0) + w * get(emit_nir, src, 0.0)
    end

    FirstOrderResult(projected_area_per_node, incident_par_power_per_node, incident_nir_power_per_node, hits_per_node)
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

function _cfg_lookup_value_local(cfg::LightConfig, keys::Vector{String})
    function _dict_lookup(d, k::String)
        d isa AbstractDict || return nothing
        haskey(d, k) && return d[k]
        lk = lowercase(k)
        for (dk, dv) in d
            lowercase(string(dk)) == lk && return dv
        end
        return nothing
    end

    props = get(cfg.raw, "prop", nothing)
    for d in (props, cfg.raw)
        for k in keys
            v = _dict_lookup(d, k)
            v !== nothing && return v
        end
    end
    return nothing
end

function _cfg_meteo_range_spec_local(cfg::LightConfig)
    v = _cfg_lookup_value_local(cfg, ["meteo_range", "meteorange", "meteoRange"])
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

function _cfg_bool_override_local(cfg::LightConfig, keys::Vector{String})
    v = _cfg_lookup_value_local(cfg, keys)
    v === nothing && return nothing
    _parse_bool_strict_local(v, keys[1])
end

"""
    _nir_scattering_enabled_local(cfg)::Bool

NIR scattering activation with optional explicit override (`nir_scattering`).
Default behavior remains enabled whenever scattering is enabled.
"""
function _nir_scattering_enabled_local(cfg::LightConfig)
    cfg.scattering || return false
    override = _cfg_bool_override_local(cfg, ["nir_scattering", "nirScattering"])
    override === nothing ? true : override
end

"""
    _nir_interception_enabled_local(cfg)::Bool

NIR interception activation with optional explicit override (`nir_interception`).
Default behavior remains enabled.
"""
function _nir_interception_enabled_local(cfg::LightConfig)
    override = _cfg_bool_override_local(cfg, ["nir_interception", "nirInterception"])
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

function _apply_meteo_range_local(rows::Vector{<:NamedTuple}, cfg::LightConfig)
    spec = _cfg_meteo_range_spec_local(cfg)
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

function _prepare_meteo_rows_for_series(meteo::MeteoTable, cfg::LightConfig)
    rows = collect(meteo.rows)
    isempty(rows) && return rows
    _validate_meteo_sequence_local(rows)
    rows = _apply_meteo_range_local(rows, cfg)
    rows = _apply_meteo_active_filter_local(rows)
    rows
end

"""
    prepare_meteo(meteo, cfg)::MeteoTable

Return the effective meteo table after Java-like meteo controls are applied:
sequence validation, optional `meteo_range`, and optional `active` filtering.
"""
function prepare_meteo(meteo::MeteoTable, cfg::LightConfig)
    MeteoTable(_prepare_meteo_rows_for_series(meteo, cfg), meteo.metadata)
end

function _default_scattering_factor_local(cfg::LightConfig, band::String)
    b = uppercase(band)
    b == "NIR" && return cfg.scattering_coeff_nir
    return cfg.scattering_coeff_par
end

function _compute_scattering_with_flags(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    first::FirstOrderResult,
    cfg::LightConfig;
    mode::Symbol=:raycast,
    backend::Union{Nothing,ScatteringBackend}=nothing,
    nir_scattering::Bool=true,
)
    cfg.scattering || return nothing
    nir_scattering && return compute_scattering(scene, turtle, first, cfg; mode=mode, backend=backend)

    par_only = compute_scattering_band(
        scene,
        turtle,
        first,
        cfg;
        mode=mode,
        backend=backend,
        band="PAR",
        initial_power_per_node=first.incident_par_power_per_node,
        default_coeff=cfg.scattering_coeff_par,
    )
    return ScatteringResult(par_only.added_power_per_node, Dict{Int,Float64}(), par_only.iterations, par_only.converged)
end

function _interception_area_per_node_local(scene::SceneGeometry, cfg::LightConfig)
    vertices, faces, face2node, node_ids, _, _ = _scene_geometry_for_interception(scene, cfg)
    area = Dict{Int,Float64}(nid => 0.0 for nid in node_ids)
    @inbounds for i in eachindex(faces)
        f = faces[i]
        nid = face2node[i]
        a = _triangle_area3d(vertices[f[1]], vertices[f[2]], vertices[f[3]])
        area[nid] = get(area, nid, 0.0) + a
    end
    return area
end

function _node_absorptance_per_band(scene::SceneGeometry, cfg::LightConfig, band::String)
    coeffs_by_group_type = _group_optical_coeffs(cfg)
    virtual_groups = _virtual_sensor_groups(cfg)
    b = uppercase(band)
    sf_default = _default_scattering_factor_local(cfg, b)
    out = Dict{Int,Float64}()
    _, _, _, node_ids, _, node_group = _scene_geometry_for_interception(scene, cfg)
    for nid in node_ids
        group = get(node_group, nid, get(scene.node_group, nid, ""))
        if group in virtual_groups
            out[nid] = 0.0
            continue
        end
        default_type = group == "pavement" ? "Cobblestone" : ""
        typ = get(scene.node_type, nid, default_type)
        coeffs = get(
            coeffs_by_group_type,
            (group, typ),
            get(coeffs_by_group_type, (group, "*"), Dict{String,Float64}()),
        )
        sf = get(coeffs, b, sf_default)
        out[nid] = clamp(1.0 - sf, 0.0, 1.0)
    end
    out
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

function _single_band_flux(total_irradiance::Float64, meteo_row, sky::SkyState, turtle::TurtleGrid, cfg::LightConfig)
    tmp = SkyState(
        sky.sun_azimuth_deg,
        sky.sun_elevation_deg,
        total_irradiance,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
    )
    compute_directional_fluxes(meteo_row, tmp, turtle, cfg).par
end

function _compute_extra_band_light(
    scene::SceneGeometry,
    meteo_row,
    sky::SkyState,
    turtle::TurtleGrid,
    cfg::LightConfig;
    interception_backend::InterceptionBackend=RasterCPUBackend(),
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
)
    extras_irr = _extra_band_irradiance(meteo_row)
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()

    isempty(extras_irr) && return extra_0_q, extra_q, extras_irr

    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(total_irr, meteo_row, sky, turtle, cfg)
        first_band = compute_first_order(
            scene,
            turtle,
            DirectionalFluxes(ids, flux_band, zeros(Float64, n)),
            cfg;
            backend=interception_backend,
        )
        order0 = Dict{Int,Float64}(nid => v for (nid, v) in first_band.incident_par_power_per_node)

        added =
            if cfg.scattering
                compute_scattering_band(
                    scene,
                    turtle,
                    first_band,
                    cfg;
                    mode=scattering_mode,
                    backend=scattering_backend,
                    band=band,
                ).added_power_per_node
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
    run_light_step(scene, meteo_row, cfg; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::LightStepResult

Run a complete light computation for one meteo row:
`compute_sky -> build_turtle -> compute_directional_fluxes -> compute_first_order -> compute_scattering -> integrate_light`.
"""
function run_light_step(
    scene::SceneGeometry,
    meteo_row,
    cfg::LightConfig;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
)
    ib = _resolve_interception_backend(interception_backend)
    nir_interception = _nir_interception_enabled_local(cfg)
    nir_scattering = _nir_scattering_enabled_local(cfg) && nir_interception
    dt_seconds = _step_duration_seconds_local(meteo_row)
    area_per_node = _interception_area_per_node_local(scene, cfg)
    abs_par = _node_absorptance_per_band(scene, cfg, "PAR")
    abs_nir = nir_interception ? _node_absorptance_per_band(scene, cfg, "NIR") : Dict{Int,Float64}()
    sky = compute_sky(meteo_row, cfg)
    nir_interception || (sky = _disable_nir_sky_local(sky))
    turtle = build_turtle(cfg, sky)
    fluxes = compute_directional_fluxes(meteo_row, sky, turtle, cfg)
    first = compute_first_order(scene, turtle, fluxes, cfg; backend=ib)
    scat = _compute_scattering_with_flags(
        scene,
        turtle,
        first,
        cfg;
        mode=scattering_mode,
        backend=scattering_backend,
        nir_scattering=nir_scattering,
    )
    extra_0_q, extra_q, extra_irr = _compute_extra_band_light(
        scene,
        meteo_row,
        sky,
        turtle,
        cfg;
        interception_backend=ib,
        scattering_mode=scattering_mode,
        scattering_backend=scattering_backend,
    )
    budget = integrate_light(
        first,
        scat,
        cfg;
        extra_0_q_per_band=extra_0_q,
        extra_q_per_band=extra_q,
        step_duration_seconds=dt_seconds,
        component_area_per_node=area_per_node,
        absorption_par_per_node=abs_par,
        absorption_nir_per_node=abs_nir,
    )
    LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
end

"""
    run_light_series(scene, meteo, cfg; interception_backend=:raster_cpu, scattering_mode=:raycast, scattering_backend=nothing)::Vector{LightStepResult}

Run the complete light pipeline for all rows in a `MeteoTable`, with optional directional
response reuse when `cfg.cache_radiation` is enabled.
"""
function run_light_series(
    scene::SceneGeometry,
    meteo::MeteoTable,
    cfg::LightConfig;
    interception_backend=:raster_cpu,
    scattering_mode::Symbol=:raycast,
    scattering_backend::Union{Nothing,ScatteringBackend}=nothing,
)
    ib = _resolve_interception_backend(interception_backend)
    nir_interception = _nir_interception_enabled_local(cfg)
    nir_scattering = _nir_scattering_enabled_local(cfg) && nir_interception
    area_per_node = _interception_area_per_node_local(scene, cfg)
    abs_par = _node_absorptance_per_band(scene, cfg, "PAR")
    abs_nir = nir_interception ? _node_absorptance_per_band(scene, cfg, "NIR") : Dict{Int,Float64}()
    use_cache = cfg.cache_radiation && _can_use_series_radiation_cache(ib)
    cache = Dict{UInt64,Tuple{
        Vector{Dict{Int,Float64}},
        Vector{Dict{Int,Int}},
        Vector{Int},
        Dict{Int,Float64},
        Dict{Int,Float64},
        Dict{Tuple{Int,Int},Float64},
    }}()
    rows_eff = _prepare_meteo_rows_for_series(meteo, cfg)
    out = Vector{LightStepResult}(undef, length(rows_eff))
    for i in eachindex(rows_eff)
        row = rows_eff[i]
        dt_seconds = _step_duration_seconds_local(row)
        sky = compute_sky(row, cfg)
        nir_interception || (sky = _disable_nir_sky_local(sky))
        turtle = build_turtle(cfg, sky)
        fluxes = compute_directional_fluxes(row, sky, turtle, cfg)

        first =
            if use_cache
                key = _turtle_cache_key(turtle, cfg)
                pa_by_sector, hits_by_sector, node_ids, emit_par, emit_nir, emitter_weights =
                    get!(cache, key) do
                        _build_sector_responses(scene, turtle, cfg)
                    end
                _combine_sector_responses(pa_by_sector, hits_by_sector, fluxes, node_ids, emit_par, emit_nir, emitter_weights)
            else
                compute_first_order(scene, turtle, fluxes, cfg; backend=ib)
            end

        scat = _compute_scattering_with_flags(
            scene,
            turtle,
            first,
            cfg;
            mode=scattering_mode,
            backend=scattering_backend,
            nir_scattering=nir_scattering,
        )
        extra_0_q, extra_q, extra_irr = _compute_extra_band_light(
            scene,
            row,
            sky,
            turtle,
            cfg;
            interception_backend=ib,
            scattering_mode=scattering_mode,
            scattering_backend=scattering_backend,
        )
        budget = integrate_light(
            first,
            scat,
            cfg;
            extra_0_q_per_band=extra_0_q,
            extra_q_per_band=extra_q,
            step_duration_seconds=dt_seconds,
            component_area_per_node=area_per_node,
            absorption_par_per_node=abs_par,
            absorption_nir_per_node=abs_nir,
        )
        out[i] = LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
    end
    out
end
