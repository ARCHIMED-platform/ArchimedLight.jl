function integrate_light(
    first::FirstOrderResult,
    scat::Union{Nothing,ScatteringResult},
    cfg::LightConfig;
    extra_0_q_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
    extra_q_per_band::Dict{String,Dict{Int,Float64}}=Dict{String,Dict{Int,Float64}}(),
)
    ri_par_f_per_node = Dict{Int,Float64}()
    ri_nir_f_per_node = Dict{Int,Float64}()
    ri_par_0_q_per_node = Dict{Int,Float64}()
    ri_nir_0_q_per_node = Dict{Int,Float64}()
    ri_par_q_per_node = Dict{Int,Float64}()
    ri_nir_q_per_node = Dict{Int,Float64}()

    node_ids = collect(keys(first.projected_area_per_node))
    for nid in node_ids
        pa = max(get(first.projected_area_per_node, nid, 0.0), eps(Float64))
        p0 = get(first.incident_par_power_per_node, nid, 0.0)
        n0 = get(first.incident_nir_power_per_node, nid, 0.0)
        ps = scat === nothing ? 0.0 : get(scat.added_par_power_per_node, nid, 0.0)
        ns = scat === nothing ? 0.0 : get(scat.added_nir_power_per_node, nid, 0.0)

        ri_par_f_per_node[nid] = (p0 + ps) / pa
        ri_nir_f_per_node[nid] = (n0 + ns) / pa

        # Unit duration by default. The full time integration is performed by callers if needed.
        ri_par_0_q_per_node[nid] = p0
        ri_nir_0_q_per_node[nid] = n0
        ri_par_q_per_node[nid] = p0 + ps
        ri_nir_q_per_node[nid] = n0 + ns
    end

    LightBudget(
        ri_par_f_per_node,
        ri_nir_f_per_node,
        ri_par_0_q_per_node,
        ri_nir_0_q_per_node,
        ri_par_q_per_node,
        ri_nir_q_per_node,
        extra_0_q_per_band,
        extra_q_per_band,
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
        )
    end
    return pa_by_sector, hits_by_sector, node_ids
end

function _combine_sector_responses(pa_by_sector, hits_by_sector, fluxes::DirectionalFluxes, node_ids_all)
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

function _single_band_flux(total_irradiance::Float64, sky::SkyState, turtle::TurtleGrid, cfg::LightConfig)
    tmp = SkyState(
        sky.sun_azimuth_deg,
        sky.sun_elevation_deg,
        total_irradiance,
        0.0,
        sky.direct_fraction,
        sky.diffuse_fraction,
    )
    compute_directional_fluxes(tmp, turtle, cfg).par
end

function _compute_extra_band_light(scene::SceneGeometry, meteo_row, sky::SkyState, turtle::TurtleGrid, cfg::LightConfig)
    extras_irr = _extra_band_irradiance(meteo_row)
    extra_0_q = Dict{String,Dict{Int,Float64}}()
    extra_q = Dict{String,Dict{Int,Float64}}()

    isempty(extras_irr) && return extra_0_q, extra_q, extras_irr

    ids = [s.id for s in turtle.sectors]
    n = length(ids)
    for (band, total_irr) in extras_irr
        flux_band = _single_band_flux(total_irr, sky, turtle, cfg)
        first_band = compute_first_order(scene, turtle, DirectionalFluxes(ids, flux_band, zeros(Float64, n)), cfg)
        order0 = Dict{Int,Float64}(nid => v for (nid, v) in first_band.incident_par_power_per_node)

        added =
            if cfg.scattering
                compute_scattering_band(scene, turtle, first_band, cfg; band=band).added_power_per_node
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

function run_light_step(scene::SceneGeometry, meteo_row, cfg::LightConfig)
    sky = compute_sky(meteo_row, cfg)
    turtle = build_turtle(cfg, sky)
    fluxes = compute_directional_fluxes(sky, turtle, cfg)
    first = compute_first_order(scene, turtle, fluxes, cfg)
    scat = cfg.scattering ? compute_scattering(scene, turtle, first, cfg) : nothing
    extra_0_q, extra_q, extra_irr = _compute_extra_band_light(scene, meteo_row, sky, turtle, cfg)
    budget = integrate_light(first, scat, cfg; extra_0_q_per_band=extra_0_q, extra_q_per_band=extra_q)
    LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
end

function run_light_series(scene::SceneGeometry, meteo::MeteoTable, cfg::LightConfig)
    cache = Dict{UInt64,Tuple{Vector{Dict{Int,Float64}},Vector{Dict{Int,Int}},Vector{Int}}}()
    out = Vector{LightStepResult}(undef, length(meteo.rows))
    for i in eachindex(meteo.rows)
        sky = compute_sky(meteo.rows[i], cfg)
        turtle = build_turtle(cfg, sky)
        fluxes = compute_directional_fluxes(sky, turtle, cfg)

        first =
            if cfg.cache_radiation
                key = _turtle_cache_key(turtle, cfg)
                pa_by_sector, hits_by_sector, node_ids =
                    get!(cache, key) do
                        _build_sector_responses(scene, turtle, cfg)
                    end
                _combine_sector_responses(pa_by_sector, hits_by_sector, fluxes, node_ids)
            else
                compute_first_order(scene, turtle, fluxes, cfg)
            end

        scat = cfg.scattering ? compute_scattering(scene, turtle, first, cfg) : nothing
        extra_0_q, extra_q, extra_irr = _compute_extra_band_light(scene, meteo.rows[i], sky, turtle, cfg)
        budget = integrate_light(first, scat, cfg; extra_0_q_per_band=extra_0_q, extra_q_per_band=extra_q)
        out[i] = LightStepResult(sky, turtle, fluxes, first, scat, budget, extra_irr)
    end
    out
end
