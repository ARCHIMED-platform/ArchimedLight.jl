using CSV
using Dates
using Statistics
using Tables
using OrderedCollections: OrderedDict

const _REGRESSION_FILTER = Set(
    filter(!isempty, strip.(split(get(ENV, "ARCHIMEDLIGHT_REGRESSION_FILTER", ""), ","))),
)

_regression_update_enabled() = lowercase(strip(get(ENV, "ARCHIMEDLIGHT_REGRESSION_UPDATE", ""))) in ("1", "true", "yes", "on")

function _regression_case_enabled(id::AbstractString)
    isempty(_REGRESSION_FILTER) && return true
    return String(id) in _REGRESSION_FILTER
end

function _regression_report_root()
    raw = strip(get(ENV, "ARCHIMEDLIGHT_REGRESSION_REPORT_DIR", ""))
    if !isempty(raw)
        return normpath(raw)
    end
    return joinpath(@__DIR__, "reports", "latest")
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

function _ordered_nt_to_dict(nt::NamedTuple)
    out = OrderedDict{String,Any}()
    for name in propertynames(nt)
        out[string(name)] = getproperty(nt, name)
    end
    return out
end

function _stable_value_string(v)
    if v isa Bool
        return v ? "true" : "false"
    elseif v isa Integer
        return string(v)
    elseif v isa AbstractFloat
        if isfinite(v)
            return replace(string(round(v; digits=6)), "." => "p", "-" => "m")
        end
        return lowercase(string(v))
    elseif v isa Symbol
        return lowercase(String(v))
    end
    s = lowercase(strip(string(v)))
    s = replace(s, r"[^a-z0-9]+" => "_")
    s = replace(s, r"_+" => "_")
    return strip(s, '_')
end

function _case_id(prefix::AbstractString, options::OrderedDict{String,Any})
    parts = String[prefix]
    for (k, v) in options
        push!(parts, "$(k)-$(_stable_value_string(v))")
    end
    return join(parts, "__")
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
    return s in ArchimedLight._row_propertynames(row) ? getproperty(row, s) : missing
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
            [
                ["step_number", "item_id", "component_id"],
                ["step_number", "object_id", "source_topology_id"],
                ["step_number", "source_topology_id"],
                ["step_number", "component_id"],
                ["step_number", "node_id"],
            ]
        elseif name == "scene_values.csv"
            [["step_number"], ["stepNumber"]]
        elseif name == "summary.csv"
            [["step_number", "group", "type", "item_id"], ["step_number", "item_id"], ["step_number", "group", "type"], ["step_number", "group", "type", "object_id"], ["step_number", "object_id"]]
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
        elseif name == "sky_summary.csv"
            [["step_number"]]
        elseif name == "sector_flux.csv"
            [["step_number", "sector_id"]]
        else
            [String[]]
        end
    for c in candidates
        all(in(cols), c) && return c
    end
    return String[]
end

function _stable_value_columns(name::String, cols::Vector{String})
    wanted =
        if name == "component_values.csv"
            [
                "area",
                "Ri_PAR_0_f",
                "Ri_NIR_0_f",
                "Ri_PAR_0_q",
                "Ri_NIR_0_q",
                "Ra_PAR_0_q",
                "Ra_NIR_0_q",
                "Ri_PAR_f",
                "Ri_NIR_f",
                "Ri_PAR_q",
                "Ri_NIR_q",
                "Ra_PAR_q",
                "Ra_NIR_q",
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
        elseif name == "sky_summary.csv"
            ["sky_mode", "turtle_sectors", "all_in_turtle", "sun_sector_count", "ri_par_f", "ri_nir_f", "sum_par", "sum_nir"]
        elseif name == "sector_flux.csv"
            ["source", "par", "nir", "weight", "dir_x", "dir_y", "dir_z"]
        else
            cols
        end
    keep = Set(wanted)
    return [c for c in cols if c in keep]
end

function _csv_tolerance(col::String)
    if col in ("date", "hour_start", "hour_end", "group", "type", "sky_mode", "source")
        return (numeric=false, atol=0.0, rtol=0.0)
    elseif occursin("barycentre", col) || startswith(col, "dir_")
        return (numeric=true, atol=1e-6, rtol=1e-6)
    elseif occursin("step", col)
        return (numeric=false, atol=0.0, rtol=0.0)
    elseif occursin("_q", col)
        return (numeric=true, atol=1e-2, rtol=1e-4)
    elseif occursin("_f", col) || col in ("par", "nir", "weight", "sum_par", "sum_nir")
        return (numeric=true, atol=1e-3, rtol=1e-4)
    else
        return (numeric=true, atol=1e-4, rtol=1e-4)
    end
end

function compare_stable_csv_paths(expected_path::AbstractString, observed_path::AbstractString; label::String="")
    exp_rows = _read_csv_rows(expected_path)
    obs_rows = _read_csv_rows(observed_path)
    isempty(exp_rows) &&
        return (
            ok=isempty(obs_rows),
            missing=0,
            extra=length(obs_rows),
            mismatch=0,
            detail=isempty(obs_rows) ? "" : "unexpected observed rows: $(length(obs_rows))",
            max_abs_error=0.0,
            max_rel_error=0.0,
        )
    isempty(obs_rows) &&
        return (
            ok=false,
            missing=length(exp_rows),
            extra=0,
            mismatch=0,
            detail="observed rows are empty",
            max_abs_error=Inf,
            max_rel_error=Inf,
        )

    name = basename(String(expected_path))
    cols = String[string(n) for n in propertynames(first(exp_rows))]
    key_cols = _key_columns_for_file(name, cols)
    isempty(key_cols) && error("$(label): unable to infer key columns for $(name)")
    value_cols = [c for c in _stable_value_columns(name, cols) if !(c in key_cols)]

    keyf(row) = Tuple(_row_get(row, c) for c in key_cols)
    exp_map = Dict{Tuple,Any}(keyf(r) => r for r in exp_rows)
    obs_map = Dict{Tuple,Any}(keyf(r) => r for r in obs_rows)
    exp_keys = Set(keys(exp_map))
    obs_keys = Set(keys(obs_map))
    missing = length(setdiff(exp_keys, obs_keys))
    extra = length(setdiff(obs_keys, exp_keys))
    mismatch = 0
    detail = ""
    max_abs_error = 0.0
    max_rel_error = 0.0

    if missing > 0
        missing_keys = setdiff(exp_keys, obs_keys)
        detail = "missing key=$(first(missing_keys))"
    elseif extra > 0
        extra_keys = setdiff(obs_keys, exp_keys)
        detail = "unexpected key=$(first(extra_keys))"
    end

    for k in intersect(exp_keys, obs_keys)
        er = exp_map[k]
        orow = obs_map[k]
        for c in value_cols
            tol = _csv_tolerance(c)
            ev = _row_get(er, c)
            ov = _row_get(orow, c)
            if tol.numeric
                ef = _try_float(ev)
                of = _try_float(ov)
                if ef === nothing || of === nothing
                    if string(ev) != string(ov)
                        mismatch += 1
                        isempty(detail) && (detail = "key=$(k) col=$(c) expected=$(ev) observed=$(ov)")
                    end
                    continue
                end
                if isnan(ef) && isnan(of)
                    continue
                end
                abs_err = abs(of - ef)
                rel_err = abs_err / max(abs(ef), eps(Float64))
                max_abs_error = max(max_abs_error, abs_err)
                max_rel_error = max(max_rel_error, rel_err)
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
                    isempty(detail) && (detail = "key=$(k) col=$(c) expected=$(ev) observed=$(ov)")
                end
            end
        end
    end

    return (
        ok=(missing == 0 && extra == 0 && mismatch == 0),
        missing=missing,
        extra=extra,
        mismatch=mismatch,
        detail=detail,
        max_abs_error=max_abs_error,
        max_rel_error=max_rel_error,
    )
end

function _write_component_series_csv(path::AbstractString, scene, models, series, options, meteo_rows)
    sim = ArchimedLight.LightSimulation(scene, models; options=options)
    return ArchimedLight.write_component_values(path, sim, series; step_index_base=0)
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

function _write_sky_summary_csv(path::AbstractString, sky, turtle, fluxes, options; step_number::Int=0, sky_mode::AbstractString)
    row = OrderedDict{String,Any}(
        "step_number" => step_number,
        "sky_mode" => String(sky_mode),
        "turtle_sectors" => options.turtle_sectors,
        "all_in_turtle" => options.all_in_turtle,
        "sun_sector_count" => count(s -> s.source == :sun, turtle.sectors),
        "ri_par_f" => sky.ri_par_f,
        "ri_nir_f" => sky.ri_nir_f,
        "sum_par" => sum(fluxes.par),
        "sum_nir" => sum(fluxes.nir),
    )
    return _rows_to_csv(path, [row])
end

function _write_sector_flux_csv(path::AbstractString, turtle, fluxes; step_number::Int=0)
    rows = OrderedDict{String,Any}[]
    for i in eachindex(turtle.sectors)
        sec = turtle.sectors[i]
        push!(
            rows,
            OrderedDict{String,Any}(
                "step_number" => step_number,
                "sector_id" => sec.id,
                "source" => String(sec.source),
                "par" => fluxes.par[i],
                "nir" => fluxes.nir[i],
                "weight" => sec.weight,
                "dir_x" => Float64(sec.direction[1]),
                "dir_y" => Float64(sec.direction[2]),
                "dir_z" => Float64(sec.direction[3]),
            ),
        )
    end
    return _rows_to_csv(path, rows)
end

function _render_ri_par_f_figure(scene, models, options, step; title::String)
    geometry = ArchimedLight._scene_geometry_for_interception(scene, models, options)
    metric = step.budget.incident_flux.total.par

    v_sum = zeros(Float64, length(geometry.vertices))
    v_count = zeros(Int, length(geometry.vertices))
    for i in eachindex(geometry.faces)
        f = geometry.faces[i]
        v = get(metric, geometry.face2node[i], NaN)
        isfinite(v) || continue
        for vid in (Int(f[1]), Int(f[2]), Int(f[3]))
            v_sum[vid] += v
            v_count[vid] += 1
        end
    end
    vertex_values = Float64[v_count[i] > 0 ? (v_sum[i] / v_count[i]) : NaN for i in eachindex(geometry.vertices)]

    colorrange = (0.0, max(step.sky.ri_par_f, eps(Float64)))
    fig = Figure(size=(960, 720))
    ax = Axis3(
        fig[1, 1];
        title=title,
        aspect=:data,
        azimuth=1.45,
        elevation=0.35,
        perspectiveness=0.65,
    )
    p = mesh!(
        ax,
        geometry.vertices,
        geometry.faces;
        color=vertex_values,
        colormap=:viridis,
        colorrange=colorrange,
        nan_color=:lightgray,
        shading=false,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    Colorbar(fig[1, 2], p, label="Ri_PAR_f")
    return fig
end

function _save_figure_png(path::AbstractString, fig)
    mkpath(dirname(path))
    CairoMakie.save(path, fig)
    return path
end

function _psnr_from_reference(reference_path::AbstractString, actual)
    reference_file = ReferenceTests.query_extended(reference_path)
    actual_img = ReferenceTests._convert(ReferenceTests.DataFormat{:PNG}, actual)
    reference_img = ReferenceTests.loadfile(typeof(actual_img), reference_file)
    return ReferenceTests._psnr(reference_img, actual_img)
end

function _compare_png_reference(reference_path::AbstractString, actual; threshold::Real=35)
    psnr = _psnr_from_reference(reference_path, actual)
    return (ok=psnr >= threshold, psnr=psnr, threshold=Float64(threshold))
end

function _timed_median(f::Function; samples::Int=3)
    times = Float64[]
    bytes = Int[]
    value = nothing
    for _ in 1:samples
        GC.gc()
        t = @timed f()
        value = t.value
        push!(times, Float64(t.time))
        push!(bytes, Int(t.bytes))
    end
    return (
        value=value,
        runtime_median_seconds=median(times),
        alloc_median_bytes=median(bytes),
        runtime_min_seconds=minimum(times),
        alloc_min_bytes=minimum(bytes),
    )
end

function _copy_tree!(src::AbstractString, dst::AbstractString)
    mkpath(dst)
    for name in readdir(src)
        src_path = joinpath(src, name)
        dst_path = joinpath(dst, name)
        if isdir(src_path)
            _copy_tree!(src_path, dst_path)
        else
            mkpath(dirname(dst_path))
            cp(src_path, dst_path; force=true)
        end
    end
    return dst
end

function _update_baseline_dir!(baseline_dir::AbstractString, observed_dir::AbstractString)
    isdir(baseline_dir) && rm(baseline_dir; recursive=true, force=true)
    return _copy_tree!(observed_dir, baseline_dir)
end

function _write_report(path::AbstractString, rows)
    mkpath(dirname(path))
    CSV.write(path, rows; delim=';')
    return path
end
