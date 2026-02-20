import CSV
import Tables
import Dates

const PARITY_LIMITS = Dict{Symbol,Float64}(
    # IDs/counts/order invariants.
    :exact => 0.0,
    :hitcount_total_rel => 1e-2,
    :hitcount_component_rel => 3e-3,
    :hitcount_component_hi_rel_turtle => 2.2e-3,
    :hitcount_component_hi_rel_raycast => 5e-3,
    :hitcount_hist_abs_turtle => 500.0,
    :hitcount_hist_abs_raycast => 900.0,
    :hitcount_hist_rel_turtle => 1.2e-3,
    :hitcount_hist_rel_raycast => 2e-3,
    # Irradiance / energy parity.
    :irr_component_rel_strict => 1e-6,
    :irr_component_rel_dense => 9e-3,
    :irr_component_abs_sparse => 7.5,
    :scene_riq_rel => 2e-3,
    :scattering_total_rel_loose => 5e-2,
    :scattering_total_rel_strict => 1e-2,
    # Sun / sky metrics.
    :sun_deg_snapshot => 5e-2,
    :sun_deg_az_series => 2e-2,
    :sun_deg_el_series => 5e-2,
)

parity_limit(metric::Symbol) = get(PARITY_LIMITS, metric) do
    error("Missing parity limit for metric=$(metric)")
end

relerr(a, b) = abs(a - b) / max(abs(b), eps(Float64))

"""
    parity_rel_with_limit(actual, expected, metric)

Return `(err, limit)` for a relative parity comparison controlled by `PARITY_LIMITS`.
"""
function parity_rel_with_limit(actual::Real, expected::Real, metric::Symbol)
    err = relerr(Float64(actual), Float64(expected))
    return err, parity_limit(metric)
end

"""
    parity_abs_with_limit(actual, expected, metric)

Return `(err, limit)` for an absolute parity comparison controlled by `PARITY_LIMITS`.
"""
function parity_abs_with_limit(actual::Real, expected::Real, metric::Symbol)
    err = abs(Float64(actual) - Float64(expected))
    return err, parity_limit(metric)
end

struct ParityFixture
    name::String
    config_path::String
    scene_path::String
    meteo_path::String
    expected_dir::Union{Nothing,String}
end

function _repo_root()
    normpath(joinpath(@__DIR__, ".."))
end

function _fixture_root()
    joinpath(_repo_root(), "java_implementation", "archimed-lib-2018", "tests")
end

function _fixture_expected_dir(test_dir::String)
    d = joinpath(test_dir, "expected")
    isdir(d) ? d : nothing
end

function _fixture_join(base::AbstractString, p::AbstractString)
    isabspath(p) ? normpath(p) : normpath(joinpath(base, p))
end

function _fixture_row(
    name::String;
    fixture_name::String=name,
    config_rel::String="config.yml",
    scene_rel::Union{Nothing,String}=nothing,
    meteo_rel::Union{Nothing,String}=nothing,
)
    test_dir = joinpath(_fixture_root(), name)
    cfg = joinpath(test_dir, config_rel)
    scene = ""
    meteo = ""
    if isfile(cfg)
        c = ArchimedLight.read_light_config(cfg)
        scene = scene_rel === nothing ? c.scene : _fixture_join(dirname(cfg), scene_rel)
        meteo = meteo_rel === nothing ? c.meteo : _fixture_join(dirname(cfg), meteo_rel)
    end
    ParityFixture(fixture_name, cfg, scene, meteo, _fixture_expected_dir(test_dir))
end

function light_parity_fixtures()
    entries = [
        ("test-hitcount", "test-hitcount", "config.yml"),
        ("test-hitcount2", "test-hitcount2", "config.yml"),
        ("test-hitcount3", "test-hitcount3", "config.yml"),
        ("test-weighted-sun", "test-weighted-sun", "config.yml"),
        ("test-save_on_disk1", "test-save_on_disk1", "config.yml"),
        ("test-save_on_disk6", "test-save_on_disk6", "config.yml"),
        ("test-skytir", "test-skytir", "config.yml"),
        ("test-area_ratio", "test-area_ratio", "config.yml"),
        ("test-area_ratio2", "test-area_ratio2", "config.yml"),
        ("test-area_ratio3", "test-area_ratio3", "config.yml"),
        ("test-area_ratio4", "test-area_ratio4", "config.yml"),
        ("test-area_ratio5", "test-area_ratio5", "config.yml"),
        ("test-cafeier", "test-cafeier", "config.yml"),
        ("test-cafeier2", "test-cafeier2", "config.yml"),
        ("test-cafeier_sensor", "test-cafeier_sensor", "config.yml"),
        ("test-cafeier_sensor2", "test-cafeier_sensor2", "config.yml"),
        ("test-cafeier_sensor3", "test-cafeier_sensor3", "config.yml"),
        ("test-cafeier_sensor4", "test-cafeier_sensor4", "config.yml"),
        ("test-scattering-one-plate", "test-scattering-one-plate", "config.yml"),
        ("test-scattering-two-plates", "test-scattering-two-plates", "config.yml"),
        ("test-scattering-divergence", "test-scattering-divergence", "config.yml"),
        ("test-links", "test-links", "config.yml"),
        ("test-links2", "test-links2", "config.yml"),
        ("test-links3", "test-links3", "config.yml"),
        ("test-links4", "test-links4", "config.yml"),
        ("test-links-stats", "test-links-stats", "config.yml"),
        ("test-links-pixeltable", "test-links-pixeltable", "config.yml"),
        ("test-links-pixeltable2", "test-links-pixeltable2", "config.yml"),
        ("test-links-sensor-plates", "test-links-sensor-plates", "config.yml"),
        ("test-links-sensor-plates2", "test-links-sensor-plates2", "config.yml"),
        ("test-cached-radiation", "test-cached-radiation", "config.yml"),
        ("test-cached-radiation2", "test-cached-radiation2", "config.yml"),
        ("test-cached-radiation3", "test-cached-radiation3", "config.yml"),
        ("test-cached-radiation4", "test-cached-radiation4", "config.yml"),
        ("test-absorb", "test-absorb", "config.yml"),
        ("test-absorb2", "test-absorb2", "config.yml"),
        ("test-independant-steps", "test-independant-steps", "config.yml", nothing, "meteo1.csv"),
        ("test-independant-steps2", "test-independant-steps2", "config.yml", nothing, "meteo1.csv"),
        ("test-meteo-stepduration", "test-meteo-stepduration", "config.yml", nothing, "meteo1.csv"),
        ("test-timestep1", "test-timestep1-onestep", "onestep/config.yml"),
        ("test-timestep1", "test-timestep1-manysteps", "manysteps/config.yml"),
        ("test-timestep2", "test-timestep2-onestep", "onestep/config.yml"),
        ("test-timestep2", "test-timestep2-manysteps", "manysteps/config.yml"),
        ("test-timestep3", "test-timestep3-onestep", "onestep/config.yml"),
        ("test-timestep3", "test-timestep3-manysteps", "manysteps/config.yml"),
        ("test-save_on_disk2", "test-save_on_disk2", "config.yml"),
        ("test-save_on_disk3", "test-save_on_disk3", "config.yml"),
        ("test-save_on_disk4", "test-save_on_disk4", "config.yml"),
        ("test-save_on_disk5", "test-save_on_disk5", "config.yml"),
        ("test-customband", "test-customband", "config.yml"),
    ]
    [
        if length(e) == 3
            _fixture_row(e[1]; fixture_name=e[2], config_rel=e[3])
        else
            _fixture_row(e[1]; fixture_name=e[2], config_rel=e[3], scene_rel=e[4], meteo_rel=e[5])
        end for e in entries
    ]
end

function read_java_csv(path::AbstractString)
    Tables.rowtable(CSV.File(path; delim=';', comment="#")) |> collect
end

function fixture_path(fx::ParityFixture)
    dirname(fx.config_path)
end

function _with_overrides(cfg::ArchimedLight.LightConfig; all_in_turtle=nothing, cache_radiation=nothing)
    ArchimedLight.LightConfig(
        cfg.scene,
        cfg.meteo,
        all_in_turtle === nothing ? cfg.all_in_turtle : Bool(all_in_turtle),
        cfg.turtle_sectors,
        cfg.pixel_size,
        cfg.area_ratio,
        cfg.scattering,
        cfg.scattering_max_iter,
        cfg.scattering_stop_ratio,
        cfg.scattering_coeff_par,
        cfg.scattering_coeff_nir,
        cache_radiation === nothing ? cfg.cache_radiation : Bool(cache_radiation),
        cfg.raw
    )
end

function run_fixture_once(fx::ParityFixture; all_in_turtle=nothing, cache_radiation=nothing)
    cfg = ArchimedLight.read_light_config(fx.config_path)
    if !isempty(fx.scene_path) || !isempty(fx.meteo_path)
        raw = copy(cfg.raw)
        !isempty(fx.scene_path) && (raw["scene"] = fx.scene_path)
        !isempty(fx.meteo_path) && (raw["meteo"] = fx.meteo_path)
        cfg = ArchimedLight.LightConfig(
            isempty(fx.scene_path) ? cfg.scene : fx.scene_path,
            isempty(fx.meteo_path) ? cfg.meteo : fx.meteo_path,
            cfg.all_in_turtle,
            cfg.turtle_sectors,
            cfg.pixel_size,
            cfg.area_ratio,
            cfg.scattering,
            cfg.scattering_max_iter,
            cfg.scattering_stop_ratio,
            cfg.scattering_coeff_par,
            cfg.scattering_coeff_nir,
            cfg.cache_radiation,
            raw,
        )
    end
    cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle, cache_radiation=cache_radiation)
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    row = first(meteo.rows)
    ArchimedLight.run_light_step(scene, row, cfg)
end

function _first_existing(paths::Vector{String})
    for p in paths
        isfile(p) && return p
    end
    return nothing
end

function _expected_component_values_path(fx::ParityFixture; all_in_turtle=nothing)
    base = fixture_path(fx)
    if fx.name == "test-hitcount2" || fx.name == "test-hitcount3"
        if all_in_turtle === nothing
            return nothing
        end
        folder = all_in_turtle ? "expected-turtletrue-6dirs" : "expected-turtlefalse-6dirs"
        p = joinpath(base, folder, "component_values.csv")
        return isfile(p) ? p : nothing
    end

    if fx.expected_dir === nothing
        return nothing
    end

    candidates = [
        joinpath(fx.expected_dir, "component_values.csv"),
    ]
    p = _first_existing(candidates)
    p !== nothing && return p

    # Some tests (e.g. weighted-sun) store expected outputs under nested timestamped directories.
    nested = collect(filter(isfile, [joinpath(root, "component_values.csv") for (root, _, _) in walkdir(fx.expected_dir)]))
    isempty(nested) ? nothing : first(sort(nested))
end

function _expected_scene_values_path(fx::ParityFixture)
    fx.expected_dir === nothing && return nothing
    p = joinpath(fx.expected_dir, "scene_values.csv")
    isfile(p) && return p
    nested = collect(filter(isfile, [joinpath(root, "scene_values.csv") for (root, _, _) in walkdir(fx.expected_dir)]))
    isempty(nested) ? nothing : first(sort(nested))
end

function _expected_sun_log_path(fx::ParityFixture)
    fx.expected_dir === nothing && return nothing
    p = joinpath(fx.expected_dir, "log-sun-position.csv")
    isfile(p) && return p
    nested = collect(filter(isfile, [joinpath(root, "log-sun-position.csv") for (root, _, _) in walkdir(fx.expected_dir)]))
    isempty(nested) ? nothing : first(sort(nested))
end

function _expected_scattering_log_path(fx::ParityFixture)
    fx.expected_dir === nothing && return nothing
    p = joinpath(fx.expected_dir, "log-iteration-scat-par.csv")
    isfile(p) ? p : nothing
end

function _expected_summary_path(fx::ParityFixture)
    fx.expected_dir === nothing && return nothing
    p = joinpath(fx.expected_dir, "summary.csv")
    isfile(p) && return p
    nested = collect(filter(isfile, [joinpath(root, "summary.csv") for (root, _, _) in walkdir(fx.expected_dir)]))
    isempty(nested) ? nothing : first(sort(nested))
end

function _sum_float_col(rows, col::Symbol)
    vals = Float64[]
    for r in rows
        if col in propertynames(r)
            v = getproperty(r, col)
            if v isa Number
                push!(vals, Float64(v))
            elseif v !== missing
                try
                    push!(vals, parse(Float64, string(v)))
                catch
                end
            end
        end
    end
    sum(vals)
end

function _mean_float_col(rows, col::Symbol)
    vals = Float64[]
    for r in rows
        if col in propertynames(r)
            v = getproperty(r, col)
            if v isa Number
                push!(vals, Float64(v))
            elseif v !== missing
                try
                    push!(vals, parse(Float64, string(v)))
                catch
                end
            end
        end
    end
    isempty(vals) ? NaN : sum(vals) / length(vals)
end

function _to_degrees_if_radians(x)
    isnan(x) && return x
    abs(x) <= 2pi + 1e-9 ? rad2deg(x) : x
end

function _parse_time_or_default(v)
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

function _step_duration_seconds(row)
    names = propertynames(row)
    if (:hour_start in names) && (:hour_end in names)
        t0 = _parse_time_or_default(getproperty(row, :hour_start))
        t1 = _parse_time_or_default(getproperty(row, :hour_end))
        dt0 = Dates.DateTime(Dates.Date(2000, 1, 1), t0)
        dt1 = Dates.DateTime(Dates.Date(2000, 1, 1), t1)
        dt1 < dt0 && (dt1 += Dates.Day(1))
        return Dates.value(dt1 - dt0) / 1000.0
    end
    return 1.0
end

function _hitcount_total_from_component_values(rows)
    total = 0
    for r in rows
        names = propertynames(r)
        (:area in names && :surface_hits in names) || continue
        area = Float64(getproperty(r, :area))
        hits = Float64(getproperty(r, :surface_hits))
        total += Int(floor(area * hits))
    end
    total
end

function parity_snapshot(step::ArchimedLight.LightStepResult)
    total_hits = sum(values(step.first_order.hits_per_node))
    total_proj_area = sum(values(step.first_order.projected_area_per_node))
    total_par0 = sum(values(step.first_order.incident_par_power_per_node))
    total_nir0 = sum(values(step.first_order.incident_nir_power_per_node))
    total_par = sum(values(step.budget.ri_par_q_per_node))
    total_nir = sum(values(step.budget.ri_nir_q_per_node))
    (
        total_hits=total_hits,
        total_projected_area=total_proj_area,
        total_par0=total_par0,
        total_nir0=total_nir0,
        total_par=total_par,
        total_nir=total_nir,
        sun_azimuth=step.sky.sun_azimuth_deg,
        sun_elevation=step.sky.sun_elevation_deg,
    )
end

function fixture_parity_report(fx::ParityFixture; all_in_turtle=nothing, cache_radiation=nothing)
    step = run_fixture_once(fx; all_in_turtle=all_in_turtle, cache_radiation=cache_radiation)
    snap = parity_snapshot(step)
    cfg = ArchimedLight.read_light_config(fx.config_path)
    cfg = _with_overrides(cfg; all_in_turtle=all_in_turtle, cache_radiation=cache_radiation)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    row = first(meteo.rows)
    dt_seconds = _step_duration_seconds(row)
    julia_scattering_total =
        if step.scattering === nothing
            0.0
        else
            sum(values(step.scattering.added_par_power_per_node)) * dt_seconds / 1e6
        end

    component_expected_path = _expected_component_values_path(fx; all_in_turtle=all_in_turtle)
    scene_expected_path = _expected_scene_values_path(fx)
    sun_log_path = _expected_sun_log_path(fx)
    scat_path = _expected_scattering_log_path(fx)

    expected_hitcount_total =
        if component_expected_path === nothing
            nothing
        else
            rows = read_java_csv(component_expected_path)
            _hitcount_total_from_component_values(rows)
        end

    expected_scene_ri_sw =
        if scene_expected_path === nothing
            nothing
        else
            rows = read_java_csv(scene_expected_path)
            :RI_SW_f in propertynames(first(rows)) ? Float64(first(rows).RI_SW_f) : nothing
        end

    expected_sun_az =
        if sun_log_path === nothing
            nothing
        else
            rows = read_java_csv(sun_log_path)
            if isempty(rows) || !(:azimuthWeighted in propertynames(first(rows)))
                nothing
            else
                _to_degrees_if_radians(Float64(first(rows).azimuthWeighted))
            end
        end

    expected_sun_el =
        if sun_log_path === nothing
            nothing
        else
            rows = read_java_csv(sun_log_path)
            if isempty(rows) || !(:elevationWeighted in propertynames(first(rows)))
                nothing
            else
                _to_degrees_if_radians(Float64(first(rows).elevationWeighted))
            end
        end

    expected_scat_total =
        if scat_path === nothing
            nothing
        else
            rows = read_java_csv(scat_path)
            _sum_float_col(rows, :scat)
        end

    return (
        fixture=fx.name,
        all_in_turtle=all_in_turtle,
        snapshot=snap,
        expected_hitcount_total=expected_hitcount_total,
        expected_scene_ri_sw=expected_scene_ri_sw,
        expected_sun_azimuth=expected_sun_az,
        expected_sun_elevation=expected_sun_el,
        expected_scattering_total=expected_scat_total,
        julia_scattering_total=julia_scattering_total,
        paths=(
            component_values=component_expected_path,
            scene_values=scene_expected_path,
            sun_log=sun_log_path,
            scat_log=scat_path,
        ),
    )
end

_row_get_by_col(row::AbstractDict{String,Any}, col::String) = get(row, col, missing)
function _row_get_by_col(row::NamedTuple, col::String)
    s = Symbol(col)
    s in propertynames(row) ? getproperty(row, s) : missing
end

function _canonical_cell_value(v)
    if v === missing || v === nothing
        return missing
    elseif v isa Dates.Date
        return v
    elseif v isa Dates.DateTime
        return v
    elseif v isa Dates.Time
        return Dates.format(v, Dates.DateFormat("HH:MM:SS"))
    elseif v isa AbstractString
        s = strip(String(v))
        for fmt in (Dates.DateFormat("yyyy/mm/dd"), Dates.DateFormat("yyyy-mm-dd"))
            try
                return Dates.Date(s, fmt)
            catch
            end
        end
        for fmt in (
            Dates.DateFormat("HH:MM:SS"),
            Dates.DateFormat("H:MM:SS"),
            Dates.DateFormat("HH:MM"),
            Dates.DateFormat("H:MM"),
        )
            try
                return Dates.format(Dates.Time(s, fmt), Dates.DateFormat("HH:MM:SS"))
            catch
            end
        end
        return s
    end
    return v
end

function _row_key_from_dict(row, key_cols::Vector{String})
    Tuple(_row_get_by_col(row, c) for c in key_cols)
end

function _row_float_value(row, col::String)
    v = _row_get_by_col(row, col)
    v === missing && return nothing
    v isa Number && return Float64(v)
    p = tryparse(Float64, strip(string(v)))
    p === nothing ? nothing : p
end

"""
    compare_rows_by_key(expected_rows, observed_rows; key_cols, value_cols, atol=1e-6, rtol=1e-3, top_n=5)

Compare two row sets keyed by `key_cols`. Returns:
- `ok::Bool`
- `missing_keys::Vector`
- `extra_keys::Vector`
- `mismatches::Vector{NamedTuple}`

Each mismatch entry contains `(key, column, expected, observed, abs_err, rel_err)`.
"""
function compare_rows_by_key(
    expected_rows::Vector,
    observed_rows::Vector;
    key_cols::Vector{String},
    value_cols::Vector{String},
    atol::Float64=1e-6,
    rtol::Float64=1e-3,
    top_n::Int=5,
)
    exp_map = Dict{Tuple,Any}(_row_key_from_dict(r, key_cols) => r for r in expected_rows)
    obs_map = Dict{Tuple,Any}(_row_key_from_dict(r, key_cols) => r for r in observed_rows)

    missing_keys = collect(setdiff(Set(keys(exp_map)), Set(keys(obs_map))))
    extra_keys = collect(setdiff(Set(keys(obs_map)), Set(keys(exp_map))))

    mismatches = NamedTuple[]
    for k in intersect(Set(keys(exp_map)), Set(keys(obs_map)))
        e = exp_map[k]
        o = obs_map[k]
        for c in value_cols
            evf = _row_float_value(e, c)
            ovf = _row_float_value(o, c)
            if evf !== nothing && ovf !== nothing
                if isnan(evf) && isnan(ovf)
                    continue
                end
                abs_err = abs(ovf - evf)
                rel_err = abs_err / max(abs(evf), eps(Float64))
                if !(abs_err <= atol || rel_err <= rtol)
                    push!(mismatches, (key=k, column=c, expected=evf, observed=ovf, abs_err=abs_err, rel_err=rel_err))
                end
            else
                ev = _canonical_cell_value(_row_get_by_col(e, c))
                ov = _canonical_cell_value(_row_get_by_col(o, c))
                ev == ov || push!(mismatches, (key=k, column=c, expected=ev, observed=ov, abs_err=NaN, rel_err=NaN))
            end
        end
    end

    # Keep diagnostics concise.
    sort!(mismatches; by=x -> (isnan(x.rel_err) ? Inf : -x.rel_err))
    if length(mismatches) > top_n
        resize!(mismatches, top_n)
    end

    ok = isempty(missing_keys) && isempty(extra_keys) && isempty(mismatches)
    return (ok=ok, missing_keys=missing_keys, extra_keys=extra_keys, mismatches=mismatches)
end
