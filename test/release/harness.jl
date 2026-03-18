using TOML
using CSV
using Tables
using Dates
using CairoMakie
using GeometryBasics

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
        "ARCHIMEDLIGHT_RELEASE_DATA_ROOT is not set. test/test_release.jl should resolve dataset and export this variable before including release harness.",
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
    cfg0 = ArchimedLight.read_light_config(fx.config_path)
    cfg = deepcopy(cfg0)
    if fx.scene_override !== nothing
        scene_path = normpath(joinpath(dirname(fx.config_path), fx.scene_override))
        cfg.source_files.scene = scene_path
    end
    if fx.meteo_override !== nothing
        meteo_path = normpath(joinpath(dirname(fx.config_path), fx.meteo_override))
        cfg.source_files.meteo = meteo_path
    end
    scattering = fx.force_scattering === nothing ? get(cfg.general, "scattering", true) : Bool(fx.force_scattering)
    cfg.general["scattering"] = scattering
    ArchimedLight.refresh_light_config!(cfg)

    scene = ArchimedLight.read_scene(cfg.source_files.scene)
    meteo = ArchimedLight.read_meteo(cfg.source_files.meteo)
    selected = ArchimedLight.prepare_meteo(meteo, cfg)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)
    length(series) == length(selected.rows) || error("fixture $(fx.id): meteo/series length mismatch")
    return (cfg=cfg, scene=scene, meteo=selected, series=series)
end

function _rows_to_csv(path::AbstractString, rows)
    isempty(rows) && return nothing
    mkpath(dirname(path))
    CSV.write(path, rows; delim=';')
    return path
end

function _write_component_series_csv(path::AbstractString, scene, series, cfg, meteo_rows)
    mkpath(dirname(path))
    first_write = true
    cols = [
        "step_number",
        "node_id",
        "source_topology_id",
        "object_id",
        "group",
        "type",
        "area",
        "barycentre_z",
        "sky_fraction",
        "Ri_PAR_0_f",
        "Ri_NIR_0_f",
        "Ri_PAR_0_q",
        "Ri_NIR_0_q",
        "Ra_PAR_0_q",
        "Ra_NIR_0_q",
    ]
    if get(cfg.general, "scattering", true)
        append!(cols, ["Ri_PAR_f", "Ri_NIR_f", "Ri_PAR_q", "Ri_NIR_q", "Ra_PAR_q", "Ra_NIR_q"])
    end
    for i in eachindex(series)
        t = ArchimedLight.component_values_table(
            scene,
            series[i],
            cfg;
            meteo_row=meteo_rows[i],
            step_number=i - 1,
            columns=cols,
            strict=false,
        )
        CSV.write(path, t.rows; delim=';', append=!first_write, writeheader=first_write)
        first_write = false
    end
    return path
end

function _meteo_rows_with_step(rows)
    out = Dict{String,Any}[]
    for (i, row) in enumerate(rows)
        d = Dict{String,Any}("step_number" => i - 1)
        for name in propertynames(row)
            d[string(name)] = getproperty(row, name)
        end
        push!(out, d)
    end
    out
end

function write_fixture_observed_outputs!(fx::JuliaFixture, out_root::AbstractString; data=nothing)
    data === nothing && (data = fixture_runtime_data(fx))
    mkpath(out_root)
    files = Dict{String,String}()

    comp_path = joinpath(out_root, "component_values.csv")
    _write_component_series_csv(comp_path, data.scene, data.series, data.cfg, data.meteo.rows)
    files["component_values.csv"] = comp_path

    scene_path = joinpath(out_root, "scene_values.csv")
    ArchimedLight.write_scene_values_csv(
        scene_path,
        data.scene,
        data.series,
        data.cfg;
        meteo_rows=data.meteo.rows,
        start_step_number=0,
        columns=[
            "step_number",
            "date",
            "hour_start",
            "hour_end",
            "RI_PAR_f",
            "RI_NIR_f",
            "RI_SW_f",
            "plot_area",
            "sun_elevation",
            "sun_azimut",
        ],
        strict=false,
    )
    files["scene_values.csv"] = scene_path

    summary_path = joinpath(out_root, "summary.csv")
    try
        ArchimedLight.write_summary_csv(
            summary_path,
            data.scene,
            data.series,
            data.cfg;
            meteo_rows=data.meteo.rows,
            start_step_number=0,
            columns=[
                "step_number",
                "date",
                "hour_start",
                "hour_end",
                "object_id",
                "group",
                "type",
                "area",
                "Ri_q",
                "Ra_q",
            ],
        )
        files["summary.csv"] = summary_path
    catch
        # Some fixtures intentionally exclude summary fields; keep regression runnable without summary.
    end

    sun_path = joinpath(out_root, "log-sun-position.csv")
    ArchimedLight.write_sun_position_log_csv(sun_path, data.series, data.meteo.rows; start_step_number=0)
    files["log-sun-position.csv"] = sun_path

    if get(data.cfg.general, "scattering", true)
        scat_path = joinpath(out_root, "log-iteration-scat-par.csv")
        first_write = true
        for i in eachindex(data.series)
            t = ArchimedLight.scattering_iteration_log_table(
                data.scene,
                data.series[i],
                data.cfg;
                meteo_row=data.meteo.rows[i],
                step_number=i - 1,
                band="PAR",
            )
            CSV.write(scat_path, t.rows; delim=';', append=!first_write, writeheader=first_write)
            first_write = false
        end
        files["log-iteration-scat-par.csv"] = scat_path
    end

    meteo_path = joinpath(out_root, "meteo.csv")
    _rows_to_csv(meteo_path, _meteo_rows_with_step(data.meteo.rows))
    files["meteo.csv"] = meteo_path

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
        return b.incident.par.flux_per_node
    elseif metric == "Ri_PAR_0_f"
        return b.incident.par.initial_flux_per_node
    elseif metric == "Ri_NIR_f"
        return b.incident.nir.flux_per_node
    elseif metric == "Ri_NIR_0_f"
        return b.incident.nir.initial_flux_per_node
    elseif metric == "Ri_PAR_q"
        return b.incident.par.energy_per_node
    elseif metric == "Ri_PAR_0_q"
        return b.incident.par.initial_energy_per_node
    elseif metric == "Ri_NIR_q"
        return b.incident.nir.energy_per_node
    elseif metric == "Ri_NIR_0_q"
        return b.incident.nir.initial_energy_per_node
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
    vertices, faces, face2node, _, _, _ = ArchimedLight._scene_geometry_for_interception(data.scene, data.cfg)

    step_values = Vector{Vector{Float64}}(undef, length(data.series))
    allvals = Float64[]
    for i in eachindex(data.series)
        metric_map = _budget_metric_map(data.series[i], fx.visual_metric)
        vals = _vertex_values_for_step(vertices, faces, face2node, metric_map)
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
            vertices,
            faces;
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
            [["step_number", "node_id"]]
        elseif name == "scene_values.csv"
            [["step_number"], ["stepNumber"]]
        elseif name == "summary.csv"
            [["step_number", "object_id", "group", "type"], ["step_number", "group", "type"], ["step_number", "group"]]
        elseif name == "meteo.csv"
            [["step_number"], ["stepNumber"], ["date", "hour_start", "hour_end"]]
        elseif name == "log-sun-position.csv"
            [["stepNumber"], ["step_number"]]
        elseif name == "log-iteration-scat-par.csv"
            [
                ["step_number", "iteration"],
                ["stepNumber", "iteration"],
                ["step", "iter", "node_id"],
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

function _stable_value_columns(name::String, cols::Vector{String})
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

function fixture_filter_from_env()
    raw = strip(get(ENV, "ARCHIMEDLIGHT_FIXTURE_FILTER", ""))
    isempty(raw) && return String[]
    String[s for s in split(raw, ',') if !isempty(strip(s))]
end

function select_fixtures(fxs::Vector{JuliaFixture}; names::Vector{String}=fixture_filter_from_env())
    isempty(names) && return fxs
    keep = Set(names)
    return [fx for fx in fxs if fx.id in keep]
end
