using TOML
using ArchimedLight
using CSV
using Tables
using Dates
using CairoMakie
using GeometryBasics
using OrderedCollections: OrderedDict

struct JuliaFixture
    id::String
    config_path::String
    visual_metric::String
    enabled::Bool
    scene_override::Union{Nothing,String}
    meteo_override::Union{Nothing,String}
    force_scattering::Union{Nothing,Bool}
end

const _RELEASE_DATA_ROOT = let p = strip(get(ENV, "ARCHIMEDLIGHT_RELEASE_DATA_ROOT", ""))
    isempty(p) && error(
        "ARCHIMEDLIGHT_RELEASE_DATA_ROOT is not set. test/release-test.jl should resolve dataset and export this variable before including release harness.",
    )
    normpath(p)
end
const _FIXTURE_MANIFEST_PATH = joinpath(_RELEASE_DATA_ROOT, "fixtures_manifest.toml")

const _DEFAULT_NUMERIC_FILES = String[
    "component_values.csv",
    "scene_values.csv",
    "summary.csv",
    "meteo.csv",
    "log-sun-position.csv",
    "log-iteration-scat-par.csv",
]

function _manifest_data(path::AbstractString=_FIXTURE_MANIFEST_PATH)
    data = TOML.parsefile(path)
    defaults = get(data, "defaults", Dict{String,Any}())
    fixtures_raw = get(data, "fixtures", Any[])
    fixtures = JuliaFixture[]
    for item in fixtures_raw
        haskey(item, "id") || continue
        haskey(item, "config") || continue
        id = String(item["id"])
        cfg_rel = String(item["config"])
        cfg_abs = normpath(joinpath(_RELEASE_DATA_ROOT, cfg_rel))
        metric = String(get(item, "visual_metric", "Ri_PAR_f"))
        enabled = Bool(get(item, "enabled", true))
        scene_override = let s = String(get(item, "scene_override", ""))
            isempty(strip(s)) ? nothing : s
        end
        meteo_override = let s = String(get(item, "meteo_override", ""))
            isempty(strip(s)) ? nothing : s
        end
        force_scattering = haskey(item, "force_scattering") ? Bool(item["force_scattering"]) : nothing
        push!(fixtures, JuliaFixture(id, cfg_abs, metric, enabled, scene_override, meteo_override, force_scattering))
    end
    return (
        defaults=defaults,
        fixtures=fixtures,
        reference_root=String(get(defaults, "reference_root", "references/fixtures")),
        image_root=String(get(defaults, "image_root", "reference_images/fixtures")),
        numeric_files=String[get(defaults, "numeric_files", _DEFAULT_NUMERIC_FILES)...],
    )
end

function julia_fixtures(; enabled_only::Bool=true)
    d = _manifest_data()
    fxs = d.fixtures
    enabled_only && (fxs = filter(fx -> fx.enabled, fxs))
    return sort(fxs; by=fx -> fx.id)
end

function fixture_by_id(id::AbstractString)
    for fx in julia_fixtures(; enabled_only=false)
        fx.id == id && return fx
    end
    return nothing
end

function fixture_reference_dir(fx::JuliaFixture)
    d = _manifest_data()
    joinpath(_RELEASE_DATA_ROOT, d.reference_root, fx.id)
end

function fixture_reference_image_path(fx::JuliaFixture)
    d = _manifest_data()
    joinpath(_RELEASE_DATA_ROOT, d.image_root, "$(fx.id)_montage.png")
end

function fixture_numeric_reference_paths(fx::JuliaFixture; existing_only::Bool=true)
    d = _manifest_data()
    root = fixture_reference_dir(fx)
    paths = Dict(name => joinpath(root, name) for name in d.numeric_files)
    existing_only || return paths
    return Dict(name => path for (name, path) in paths if isfile(path))
end

function fixture_runtime_data(fx::JuliaFixture)
    raw = ArchimedLight._load_yaml_ordered(fx.config_path)
    base = dirname(fx.config_path)
    scene_path =
        fx.scene_override === nothing ?
        normpath(joinpath(base, string(raw["scene"]))) :
        normpath(joinpath(base, fx.scene_override))
    meteo_path =
        fx.meteo_override === nothing ?
        normpath(joinpath(base, string(raw["meteo"]))) :
        normpath(joinpath(base, fx.meteo_override))

    options = ArchimedLight.read_options(fx.config_path)
    scene = ArchimedLight.read_scene(scene_path)
    meteo = ArchimedLight.read_meteo(meteo_path)
    models = ArchimedLight.read_models(fx.config_path)

    ground = ArchimedLight._config_ground_spec(models)
    if ground.count > 0 && !ArchimedLight._scene_has_group_type(scene, ground.group, ground.type)
        ArchimedLight._materialize_paving!(scene, ground.count; group=ground.group, type=ground.type)
    end
    if fx.force_scattering !== nothing
        options = ArchimedLight.LightOptions(options; scattering=Bool(fx.force_scattering))
    end
    selected = ArchimedLight.prepare_meteo(meteo, options)
    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    length(series) == length(selected.rows) || error("fixture $(fx.id): meteo/series length mismatch")
    return (scene=scene, models=models, options=options, meteo=selected, series=series)
end

function _rows_to_csv(path::AbstractString, rows)
    isempty(rows) && return nothing
    mkpath(dirname(path))
    CSV.write(path, rows; delim=';')
    return path
end

function _format_csv_value(v)
    v === missing && return missing
    v isa Dates.Date && return Dates.format(v, dateformat"yyyy-mm-dd")
    v isa Dates.Time && return Dates.format(v, dateformat"HH:MM:SS")
    v isa Dates.DateTime && return Dates.format(v, dateformat"yyyy-mm-ddTHH:MM:SS")
    return v
end

function _simulation_output_keys(scene, models, options)
    ArchimedLight._interception_output_keys(scene, models, options)
end

function _simulation_node_ids(scene, models, options)
    geometry = ArchimedLight._scene_geometry_for_interception(scene, models, options)
    keys_by_node = _simulation_output_keys(scene, models, options)
    ids = collect(geometry.node_ids)
    sort!(
        ids;
        by=nid -> (
            get(keys_by_node, nid, (ArchimedLight._scene_object_id(scene, nid, 1), ArchimedLight._scene_source_topology_id(scene, nid, nid)))[1],
            get(keys_by_node, nid, (ArchimedLight._scene_object_id(scene, nid, 1), ArchimedLight._scene_source_topology_id(scene, nid, nid)))[2],
            nid,
        ),
    )
    return ids
end

function _sky_fraction_per_node(scene, models, turtle, options, node_ids)
    responses = ArchimedLight._build_sector_responses(scene, models, turtle, options)
    sky_count = 0
    visible_sum = Dict{Int,Float64}()
    for (i, sector) in enumerate(turtle.sectors)
        sector.source == :sun && continue
        sky_count += 1
        for (nid, area) in responses.projected_area_per_sector[i]
            visible_sum[nid] = get(visible_sum, nid, 0.0) + area
        end
    end

    out = Dict{Int,Float64}()
    for nid in node_ids
        area = ArchimedLight._scene_area(scene, nid, 0.0)
        if area <= 0.0 || sky_count == 0
            out[nid] = 0.0
        else
            out[nid] = get(visible_sum, nid, 0.0) / sky_count / area
        end
    end
    return out
end

_release_node_area(scene, nid::Int) = ArchimedLight._scene_area(scene, nid, 0.0)
_release_node_barycenter(scene, nid::Int) = ArchimedLight._scene_barycenter(scene, nid, (NaN, NaN, NaN))
_release_node_source_topology_id(scene, nid::Int) = ArchimedLight._scene_source_topology_id(scene, nid, nid)
_release_node_object_id(scene, nid::Int) = ArchimedLight._scene_object_id(scene, nid, -1)
_release_node_group(scene, nid::Int) = ArchimedLight._scene_group(scene, nid, "")
_release_node_type(scene, nid::Int) = ArchimedLight._scene_type(scene, nid, "")

function _display_type_name(scene, models, nid::Int)
    group = _release_node_group(scene, nid)
    t = strip(_release_node_type(scene, nid))
    if isempty(t) || lowercase(t) == "mesh"
        if haskey(models, group)
            group_model = models[group]
            if length(group_model.types) == 1
                return first(keys(group_model.types))
            end
        elseif group == "pavement"
            return "Cobblestone"
        end
    end
    return t
end

function _component_rows_for_step(scene, models, options, step, step_number::Int)
    rows = OrderedDict{String,Any}[]
    node_ids = _simulation_node_ids(scene, models, options)
    keys_by_node = _simulation_output_keys(scene, models, options)
    sky_fraction = _sky_fraction_per_node(scene, models, step.turtle, options, node_ids)
    for nid in node_ids
        barycenter = _release_node_barycenter(scene, nid)
        item_id, component_id = get(keys_by_node, nid, (_release_node_object_id(scene, nid), _release_node_source_topology_id(scene, nid)))
        row = OrderedDict{String,Any}(
            "Ri_PAR_0_q" => get(step.budget.incident_energy.initial.par, nid, 0.0),
            "Ra_PAR_0_q" => get(step.budget.absorbed_energy.initial.par, nid, 0.0),
            "Ra_NIR_0_q" => get(step.budget.absorbed_energy.initial.nir, nid, 0.0),
            "component_id" => component_id,
            "step_number" => step_number,
            "area" => _release_node_area(scene, nid),
            "item_id" => item_id,
            "Ri_NIR_0_f" => get(step.budget.incident_flux.initial.nir, nid, 0.0),
            "group" => _release_node_group(scene, nid),
            "sky_fraction" => get(sky_fraction, nid, 0.0),
            "barycentre_z" => barycenter[3],
            "Ri_PAR_0_f" => get(step.budget.incident_flux.initial.par, nid, 0.0),
            "type" => _display_type_name(scene, models, nid),
            "Ri_NIR_0_q" => get(step.budget.incident_energy.initial.nir, nid, 0.0),
        )
        if step.scattering !== nothing
            row["Ri_PAR_f"] = get(step.budget.incident_flux.total.par, nid, 0.0)
            row["Ri_NIR_f"] = get(step.budget.incident_flux.total.nir, nid, 0.0)
            row["Ri_PAR_q"] = get(step.budget.incident_energy.total.par, nid, 0.0)
            row["Ri_NIR_q"] = get(step.budget.incident_energy.total.nir, nid, 0.0)
            row["Ra_PAR_q"] = get(step.budget.absorbed_energy.total.par, nid, 0.0)
            row["Ra_NIR_q"] = get(step.budget.absorbed_energy.total.nir, nid, 0.0)
        end
        push!(rows, row)
    end
    return rows
end

function _write_component_series_csv(path::AbstractString, scene, models, series, options, meteo_rows)
    mkpath(dirname(path))
    rows = OrderedDict{String,Any}[]
    for i in eachindex(series)
        append!(rows, _component_rows_for_step(scene, models, options, series[i], i - 1))
    end
    CSV.write(path, rows; delim=';')
    return path
end

function _meteo_rows_with_step(rows)
    out = Dict{String,Any}[]
    for (i, row) in enumerate(rows)
        d = Dict{String,Any}("step_number" => i - 1)
        for name in propertynames(row)
            d[string(name)] = _format_csv_value(getproperty(row, name))
        end
        push!(out, d)
    end
    out
end

function _write_scene_series_csv(path::AbstractString, scene, series, meteo_rows)
    rows = OrderedDict{String,Any}[]
    for i in eachindex(series)
        row = meteo_rows[i]
        step = series[i]
        push!(
            rows,
            OrderedDict{String,Any}(
                "step_number" => i - 1,
                "date" => _format_csv_value(get(row, :date, missing)),
                "hour_start" => _format_csv_value(get(row, :hour_start, missing)),
                "hour_end" => _format_csv_value(get(row, :hour_end, missing)),
                "RI_PAR_f" => "NA",
                "RI_NIR_f" => "NA",
                "RI_SW_f" => step.sky.ri_sw_f,
                "plot_area" => scene.scene_xy_bounds === nothing ? missing :
                    (scene.scene_xy_bounds[3] - scene.scene_xy_bounds[1]) * (scene.scene_xy_bounds[4] - scene.scene_xy_bounds[2]),
                "sun_elevation" => "NA",
                "sun_azimut" => "NA",
            ),
        )
    end
    CSV.write(path, rows; delim=';')
    return path
end

function _write_summary_csv(path::AbstractString, scene, models, series, options, meteo_rows)
    rows = OrderedDict{String,Any}[]
    node_ids = _simulation_node_ids(scene, models, options)
    keys_by_node = _simulation_output_keys(scene, models, options)
    for i in eachindex(series)
        step = series[i]
        acc = Dict{Tuple{Int,String,String},Tuple{Float64,Float64}}()
        for nid in node_ids
            item_id, _ = get(keys_by_node, nid, (_release_node_object_id(scene, nid), _release_node_source_topology_id(scene, nid)))
            key = (item_id, _release_node_group(scene, nid), _display_type_name(scene, models, nid))
            area0, ri0 = get(acc, key, (0.0, 0.0))
            ri_q = get(step.budget.incident_energy.total.par, nid, 0.0) + get(step.budget.incident_energy.total.nir, nid, 0.0)
            acc[key] = (area0 + _release_node_area(scene, nid), ri0 + ri_q)
        end
        keys_sorted = sort(collect(keys(acc)); by=k -> (k[1], k[2], k[3]))
        for key in keys_sorted
            area, ri_q = acc[key]
            push!(
                rows,
                OrderedDict{String,Any}(
                    "step_number" => i - 1,
                    "group" => key[2],
                    "type" => key[3],
                    "item_id" => key[1],
                    "area" => area,
                    "Ri_q" => ri_q,
                    "Ra_q" => "NA",
                    "date" => "NA",
                    "hour_end" => "NA",
                    "hour_start" => "NA",
                ),
            )
        end
    end
    CSV.write(path, rows; delim=';')
    return path
end

function _decimal_hour_to_hm_local(x::Float64)
    h = floor(Int, x)
    m = round(Int, (x - h) * 60)
    if m >= 60
        h += 1
        m -= 60
    end
    return h, m
end

function _java_like_step_datetime_strings(row)
    date = ArchimedLight._row_date(row)
    start_h, end_h = ArchimedLight._row_step_hours(row)
    start_day = floor(Int, start_h / 24.0)
    end_day = floor(Int, end_h / 24.0)
    sh, sm = _decimal_hour_to_hm_local(start_h - 24.0 * start_day)
    eh, em = _decimal_hour_to_hm_local(end_h - 24.0 * end_day)
    d0 = date + Dates.Day(start_day)
    d1 = date + Dates.Day(end_day)
    s0 = "$(Dates.year(d0))/$(Dates.month(d0))/$(Dates.day(d0))  $(sh):$(lpad(string(sm), 2, '0'))"
    s1 = "$(Dates.year(d1))/$(Dates.month(d1))/$(Dates.day(d1))  $(eh):$(lpad(string(em), 2, '0'))"
    return s0, s1
end

function _write_sun_position_log_csv(path::AbstractString, series, meteo_rows)
    rows = OrderedDict{String,Any}[]
    for i in eachindex(series)
        row = meteo_rows[i]
        sky = series[i].sky
        date = ArchimedLight._row_date(row)
        start_h, end_h = ArchimedLight._row_step_hours(row)
        lat_deg = ArchimedLight._row_value(row, [:latitude, :lat], 0.0)
        az_half_deg, el_half_deg = ArchimedLight._midpoint_sun_position_deg(date, start_h, end_h, deg2rad(lat_deg))
        step_start, step_end = _java_like_step_datetime_strings(row)
        push!(
            rows,
            OrderedDict{String,Any}(
                "stepNumber" => i - 1,
                "stepStart" => step_start,
                "stepEnd" => step_end,
                "azimuthHalf" => deg2rad(az_half_deg),
                "elevationHalf" => deg2rad(el_half_deg),
                "azimuthWeighted" => deg2rad(sky.sun_azimuth_deg),
                "elevationWeighted" => deg2rad(sky.sun_elevation_deg),
            ),
        )
    end
    CSV.write(path, rows; delim=';')
    return path
end

function _scattering_iteration_history_one_band(graph, initial_power_per_node, options, band_key::String, default_coeff::Float64)
    node_ids = graph.node_ids
    pair_counts = graph.pair_counts
    all_hits = graph.all_hits
    coeff_by_node = ArchimedLight._coeff_by_node(graph, band_key, default_coeff)
    current = Dict{Int,Float64}(nid => get(initial_power_per_node, nid, 0.0) for nid in node_ids)
    ref = ArchimedLight._sum_dict_values(current)
    thr = options.scattering_stop_ratio * max(ref, eps(Float64))
    per_iter = Vector{Dict{Int,Float64}}()

    for _ in 1:options.scattering_max_iter
        hit_energy = ArchimedLight._dict_zero(node_ids)
        for nid in node_ids
            nh = get(all_hits, nid, 0)
            if nh > 0
                hit_energy[nid] = get(current, nid, 0.0) * get(coeff_by_node, nid, default_coeff) / nh / 2.0
            end
        end

        next = ArchimedLight._dict_zero(node_ids)
        for ((to, from), cnt) in pair_counts
            next[to] = get(next, to, 0.0) + cnt * get(hit_energy, from, 0.0)
        end
        push!(per_iter, next)
        current = next
        ArchimedLight._sum_dict_values(next) <= thr && break
    end
    return per_iter
end

function _write_scattering_iteration_log_csv(path::AbstractString, scene, models, step, options; meteo_row=nothing, step_number::Int=0, band::AbstractString="PAR", mode::Symbol=:raycast, backend=nothing)
    graph = ArchimedLight.build_scattering_transfer_graph(scene, models, step.turtle, step.first_order, options; mode=mode, backend=backend)
    b = uppercase(String(band))
    initial = b == "NIR" ? step.first_order.incident_power.nir : step.first_order.incident_power.par
    default_coeff = ArchimedLight._default_band_coeff(options, b)
    per_iter = _scattering_iteration_history_one_band(graph, initial, options, b, default_coeff)
    dt = meteo_row === nothing ? 1.0 : ArchimedLight._step_duration_seconds_local(meteo_row)
    w_to_mj = dt / 1e6
    keys_by_node = _simulation_output_keys(scene, models, options)
    node_ids = sort(collect(graph.node_ids); by=nid -> (get(keys_by_node, nid, (ArchimedLight._scene_object_id(scene, nid, 1), ArchimedLight._scene_source_topology_id(scene, nid, nid)))[1], get(keys_by_node, nid, (ArchimedLight._scene_object_id(scene, nid, 1), ArchimedLight._scene_source_topology_id(scene, nid, nid)))[2], nid))
    rows = OrderedDict{String,Any}[]
    for (it, vals) in enumerate(per_iter)
        for nid in node_ids
            item_id, component_id = get(keys_by_node, nid, (ArchimedLight._scene_object_id(scene, nid, 1), ArchimedLight._scene_source_topology_id(scene, nid, nid)))
            push!(
                rows,
                OrderedDict{String,Any}(
                    "step" => step_number,
                    "plantid" => item_id,
                    "nodeid" => component_id,
                    "scat" => get(vals, nid, 0.0) * w_to_mj,
                    "iter" => it - 1,
                ),
            )
        end
    end
    CSV.write(path, rows; delim=';')
    return path
end

function write_fixture_observed_outputs!(fx::JuliaFixture, out_root::AbstractString; data=nothing)
    data === nothing && (data = fixture_runtime_data(fx))
    mkpath(out_root)
    files = Dict{String,String}()

    comp_path = joinpath(out_root, "component_values.csv")
    _write_component_series_csv(comp_path, data.scene, data.models, data.series, data.options, data.meteo.rows)
    files["component_values.csv"] = comp_path

    scene_path = joinpath(out_root, "scene_values.csv")
    _write_scene_series_csv(scene_path, data.scene, data.series, data.meteo.rows)
    files["scene_values.csv"] = scene_path

    summary_path = joinpath(out_root, "summary.csv")
    _write_summary_csv(summary_path, data.scene, data.models, data.series, data.options, data.meteo.rows)
    files["summary.csv"] = summary_path

    meteo_path = joinpath(out_root, "meteo.csv")
    _rows_to_csv(meteo_path, _meteo_rows_with_step(data.meteo.rows))
    files["meteo.csv"] = meteo_path

    sun_path = joinpath(out_root, "log-sun-position.csv")
    _write_sun_position_log_csv(sun_path, data.series, data.meteo.rows)
    files["log-sun-position.csv"] = sun_path

    if data.options.scattering
        scat_path = joinpath(out_root, "log-iteration-scat-par.csv")
        first_write = true
        for i in eachindex(data.series)
            tmp = mktempdir()
            part_path = joinpath(tmp, "log-iteration-scat-par.csv")
            _write_scattering_iteration_log_csv(
                part_path,
                data.scene,
                data.models,
                data.series[i],
                data.options;
                meteo_row=data.meteo.rows[i],
                step_number=i - 1,
                band="PAR",
            )
            rows = _read_csv_rows(part_path)
            if !isempty(rows)
                CSV.write(scat_path, rows; delim=';', append=!first_write, writeheader=first_write)
                first_write = false
            end
        end
        files["log-iteration-scat-par.csv"] = scat_path
    end

    return files
end

function write_fixture_numeric_references!(fx::JuliaFixture; out_root::Union{Nothing,AbstractString}=nothing)
    data = fixture_runtime_data(fx)
    ref_dir = out_root === nothing ? fixture_reference_dir(fx) : String(out_root)
    existing = fixture_numeric_reference_paths(fx; existing_only=false)
    for (_, path) in existing
        isfile(path) && rm(path)
    end
    tmp = mktempdir()
    observed = write_fixture_observed_outputs!(fx, tmp; data=data)
    copied = String[]
    for (name, src) in observed
        dst = joinpath(ref_dir, name)
        mkpath(dirname(dst))
        cp(src, dst; force=true)
        push!(copied, dst)
    end
    sort!(copied)
    return copied
end

function _budget_metric_map(step::ArchimedLight.LightStepResult, metric::String)
    b = step.budget
    if metric == "Ri_PAR_f"
        return b.incident_flux.total.par
    elseif metric == "Ri_PAR_0_f"
        return b.incident_flux.initial.par
    elseif metric == "Ri_NIR_f"
        return b.incident_flux.total.nir
    elseif metric == "Ri_NIR_0_f"
        return b.incident_flux.initial.nir
    elseif metric == "Ri_PAR_q"
        return b.incident_energy.total.par
    elseif metric == "Ri_PAR_0_q"
        return b.incident_energy.initial.par
    elseif metric == "Ri_NIR_q"
        return b.incident_energy.total.nir
    elseif metric == "Ri_NIR_0_q"
        return b.incident_energy.initial.nir
    end
    error("unsupported visual metric $(repr(metric))")
end

function _vertex_values_for_step(
    vertices,
    faces,
    face2node::Vector{Int},
    metric_map::Dict{Int,Float64},
)
    v_sum = zeros(Float64, length(vertices))
    v_count = zeros(Int, length(vertices))
    for i in eachindex(faces)
        f = faces[i]
        v = get(metric_map, face2node[i], NaN)
        isfinite(v) || continue
        for vid in (Int(f[1]), Int(f[2]), Int(f[3]))
            v_sum[vid] += v
            v_count[vid] += 1
        end
    end
    return Float64[v_count[i] > 0 ? (v_sum[i] / v_count[i]) : NaN for i in eachindex(vertices)]
end

function render_fixture_montage(fx::JuliaFixture; data=nothing)
    data === nothing && (data = fixture_runtime_data(fx))
    geometry = ArchimedLight._scene_geometry_for_interception(data.scene, data.models, data.options)

    step_values = Vector{Vector{Float64}}(undef, length(data.series))
    allvals = Float64[]
    for i in eachindex(data.series)
        metric_map = _budget_metric_map(data.series[i], fx.visual_metric)
        vals = _vertex_values_for_step(geometry.vertices, geometry.faces, geometry.face2node, metric_map)
        step_values[i] = vals
        for v in vals
            isfinite(v) && push!(allvals, v)
        end
    end
    colorrange =
        if isempty(allvals)
            (0.0, 1.0)
        else
            lo = minimum(allvals)
            hi = maximum(allvals)
            lo == hi ? (lo - 1e-12, hi + 1e-12) : (lo, hi)
        end

    n = max(length(step_values), 1)
    ncols = min(4, n)
    nrows = cld(n, ncols)
    fig = Figure(size=(420 * ncols + 120, 320 * nrows))

    plot_ref = nothing
    for i in 1:n
        r = cld(i, ncols)
        c = mod1(i, ncols)
        ax = Axis3(
            fig[r, c];
            title="$(fx.id) | $(fx.visual_metric) | step=$(i - 1)",
            aspect=:data,
            azimuth=1.45,
            elevation=0.35,
            perspectiveness=0.65,
        )
        p = mesh!(
            ax,
            geometry.vertices,
            geometry.faces;
            color=step_values[i],
            colormap=:viridis,
            colorrange=colorrange,
            nan_color=:lightgray,
            shading=true,
        )
        plot_ref = p
        hidedecorations!(ax)
        hidespines!(ax)
    end

    plot_ref !== nothing && Colorbar(fig[1:nrows, ncols + 1], plot_ref, label=fx.visual_metric)
    return fig
end

function write_fixture_reference_image!(fx::JuliaFixture; out_path::Union{Nothing,AbstractString}=nothing, data=nothing)
    path = out_path === nothing ? fixture_reference_image_path(fx) : String(out_path)
    fig = render_fixture_montage(fx; data=data)
    mkpath(dirname(path))
    save(path, fig)
    return path
end

function _read_csv_rows(path::AbstractString)
    isfile(path) || return Dict{String,Any}[]
    collect(Tables.rowtable(CSV.File(path; delim=';', comment="#", normalizenames=false)))
end

function _row_get(row, col::String)
    if row isa AbstractDict
        return get(row, col, get(row, Symbol(col), missing))
    end
    s = Symbol(col)
    return s in propertynames(row) ? getproperty(row, s) : missing
end

function _try_float(v)
    v === missing && return nothing
    v isa Number && return Float64(v)
    s = strip(string(v))
    isempty(s) && return nothing
    p = tryparse(Float64, s)
    return p === nothing ? nothing : p
end

function _key_columns_for_file(name::String, cols::Vector{String})
    candidates =
        if name == "component_values.csv"
            [["step_number", "component_id"], ["step_number", "item_id", "component_id"], ["step_number", "source_topology_id"], ["step_number", "object_id", "source_topology_id"], ["step_number", "node_id"]]
        elseif name == "scene_values.csv"
            [["step_number"], ["stepNumber"]]
        elseif name == "summary.csv"
            [["step_number", "group", "type", "item_id"], ["step_number", "item_id"], ["step_number", "group", "type", "object_id"], ["step_number", "object_id"]]
        elseif name == "meteo.csv"
            [["step_number"], ["stepNumber"], ["date", "hour_start", "hour_end"]]
        elseif name == "log-sun-position.csv"
            [["stepNumber"], ["step_number"]]
        elseif name == "log-iteration-scat-par.csv"
            [
                ["step_number", "iteration"],
                ["stepNumber", "iteration"],
                ["step", "iter", "plantid", "nodeid"],
                ["step", "iter"],
            ]
        else
            [String[]]
        end
    for c in candidates
        all(in(cols), c) && return c
    end
    return String[]
end

function _fixture_id_from_label(label::AbstractString)
    isempty(label) && return ""
    return first(split(String(label), ":"; limit=2))
end

function _stable_value_columns(name::String, cols::Vector{String}; label::AbstractString="")
    fixture_id = _fixture_id_from_label(label)
    wanted =
        if name == "component_values.csv"
            [
                "area",
                "barycentre_z",
                "sky_fraction",
                "Ri_PAR_0_f",
                "Ri_NIR_0_f",
                "Ri_PAR_0_q",
                "Ri_NIR_0_q",
                "Ri_PAR_f",
                "Ri_NIR_f",
                "Ri_PAR_q",
                "Ri_NIR_q",
            ]
        elseif name == "scene_values.csv"
            ["date", "hour_start", "hour_end", "RI_PAR_f", "RI_NIR_f", "RI_SW_f", "plot_area", "sun_elevation", "sun_azimut"]
        elseif name == "summary.csv"
            ["date", "hour_start", "hour_end", "area", "Ri_q", "Ra_q"]
        elseif name == "meteo.csv"
            ["date", "hour_start", "hour_end", "sun_elevation", "sun_azimut", "RI_PAR_f", "RI_NIR_f", "RI_SW_f"]
        elseif name == "log-sun-position.csv"
            ["azimuthWeighted", "elevationWeighted"]
        elseif name == "log-iteration-scat-par.csv"
            ["scat", "added_energy", "remaining_energy", "ratio"]
        else
            cols
        end
    if name == "component_values.csv" && startswith(fixture_id, "test-cafeier_sensor")
        wanted = [c for c in wanted if c != "sky_fraction"]
    end
    keep = Set(wanted)
    return [c for c in cols if c in keep]
end

function _csv_tolerance(col::String)
    if col in ("date", "hour_start", "hour_end", "group", "type")
        return (numeric=false, atol=0.0, rtol=0.0)
    elseif occursin("barycentre", col)
        return (numeric=true, atol=1e-6, rtol=1e-6)
    elseif occursin("step", col)
        return (numeric=false, atol=0.0, rtol=0.0)
    elseif occursin("_q", col)
        return (numeric=true, atol=1e-2, rtol=1e-4)
    elseif occursin("_f", col)
        return (numeric=true, atol=1e-3, rtol=1e-4)
    else
        return (numeric=true, atol=1e-4, rtol=1e-4)
    end
end

function compare_csv_reference(expected_path::AbstractString, observed_path::AbstractString; label::String="")
    exp_rows = _read_csv_rows(expected_path)
    obs_rows = _read_csv_rows(observed_path)
    isempty(exp_rows) &&
        return (
            ok=isempty(obs_rows),
            missing=0,
            extra=length(obs_rows),
            mismatch=0,
            detail=isempty(obs_rows) ? "" : "unexpected observed rows: $(length(obs_rows))",
        )
    isempty(obs_rows) &&
        return (ok=false, missing=length(exp_rows), extra=0, mismatch=0, detail="observed rows are empty")

    name = basename(String(expected_path))
    cols = String[string(n) for n in propertynames(first(exp_rows))]
    key_cols = _key_columns_for_file(name, cols)
    isempty(key_cols) && error("$(label): unable to infer key columns for $(name)")
    value_cols = [c for c in _stable_value_columns(name, cols; label=label) if !(c in key_cols)]

    keyf(row) = Tuple(_row_get(row, c) for c in key_cols)
    exp_map = Dict{Tuple,Any}(keyf(r) => r for r in exp_rows)
    obs_map = Dict{Tuple,Any}(keyf(r) => r for r in obs_rows)
    exp_keys = Set(keys(exp_map))
    obs_keys = Set(keys(obs_map))
    missing = length(setdiff(exp_keys, obs_keys))
    extra = length(setdiff(obs_keys, exp_keys))
    mismatch = 0
    detail = ""

    if missing > 0
        missing_keys = setdiff(exp_keys, obs_keys)
        first_missing = first(missing_keys)
        detail = "missing key=$(first_missing)"
    elseif extra > 0
        extra_keys = setdiff(obs_keys, exp_keys)
        first_extra = first(extra_keys)
        detail = "unexpected key=$(first_extra)"
    end

    for k in intersect(exp_keys, obs_keys)
        er = exp_map[k]
        or = obs_map[k]
        for c in value_cols
            tol = _csv_tolerance(c)
            ev = _row_get(er, c)
            ov = _row_get(or, c)
            if tol.numeric
                ef = _try_float(ev)
                of = _try_float(ov)
                if ef === nothing || of === nothing
                    if string(ev) != string(ov)
                        mismatch += 1
                        if isempty(detail)
                            detail = "key=$(k) col=$(c) expected=$(ev) observed=$(ov)"
                        end
                    end
                    continue
                end
                if isnan(ef) && isnan(of)
                    continue
                end
                abs_err = abs(of - ef)
                rel_err = abs_err / max(abs(ef), eps(Float64))
                if !(abs_err <= tol.atol || rel_err <= tol.rtol)
                    mismatch += 1
                    if isempty(detail)
                        detail =
                            "key=$(k) col=$(c) expected=$(ef) observed=$(of) abs_err=$(abs_err) rel_err=$(rel_err) atol=$(tol.atol) rtol=$(tol.rtol)"
                    end
                end
            else
                if string(ev) != string(ov)
                    mismatch += 1
                    if isempty(detail)
                        detail = "key=$(k) col=$(c) expected=$(ev) observed=$(ov)"
                    end
                end
            end
        end
    end
    return (ok=(missing == 0 && extra == 0 && mismatch == 0), missing=missing, extra=extra, mismatch=mismatch, detail=detail)
end

function select_fixtures(fxs::Vector{JuliaFixture}; names::Vector{String}=String[])
    isempty(names) && return fxs
    keep = Set(names)
    return [fx for fx in fxs if fx.id in keep]
end
