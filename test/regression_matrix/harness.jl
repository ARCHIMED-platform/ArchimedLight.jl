using Test
using Artifacts: artifact_hash, artifact_path
using Pkg.Artifacts: ensure_artifact_installed
using ArchimedLight
using CairoMakie
using CSV
using Dates
using GeometryBasics
using ReferenceTests
using SHA
using Tables
import PlantGeom

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "synthetic_support.jl"))
include(joinpath(@__DIR__, "..", "support.jl"))

struct RegressionScenario
    id::String
    source_kind::Symbol
    source_id::String
end

struct RegressionCase
    id::String
    scenario::RegressionScenario
    options::OrderedDict{String,Any}
    strict::Bool
    visual::Bool
end

const _FAST_FIXTURE_CACHE = Dict{String,Any}()
const _KNOWN_OPTION_COLUMNS = [
    "sky_mode",
    "toricity",
    "area_ratio",
    "scattering",
    "pixel_size_m",
    "cache_radiation",
    "cache_pixel_table",
    "nir_interception",
    "nir_scattering",
    "radiation_timestep",
    "java_logged_turtle_dirs",
    "scattering_mode",
]

function _regression_profile()
    lowercase(strip(get(ENV, "ARCHIMEDLIGHT_REGRESSION_PROFILE", "fast")))
end

function _baseline_root()
    joinpath(@__DIR__, "baselines")
end

function _default_case_options()
    OrderedDict{String,Any}(
        "sky_mode" => "46_direct",
        "toricity" => false,
        "area_ratio" => true,
        "scattering" => false,
        "pixel_size_m" => 0.05,
        "cache_radiation" => false,
        "cache_pixel_table" => false,
        "nir_interception" => true,
        "nir_scattering" => true,
        "radiation_timestep" => 15,
        "java_logged_turtle_dirs" => false,
        "scattering_mode" => :raycast,
    )
end

function _merge_case_options(overrides::AbstractDict{String,<:Any})
    out = _default_case_options()
    for (k, v) in overrides
        out[k] = v
    end
    return out
end

function _push_case!(cases::Vector{RegressionCase}, prefix::AbstractString, scenario::RegressionScenario, overrides::AbstractDict{String,<:Any}; strict::Bool=false, visual::Bool=false)
    opts = _merge_case_options(overrides)
    push!(cases, RegressionCase(_case_id(prefix, opts), scenario, opts, strict, visual))
    return cases
end

function _sky_mode_params(mode::AbstractString)
    if mode == "1_direct"
        return (sectors=1, all_in_turtle=false)
    elseif mode == "16_all_in_turtle"
        return (sectors=16, all_in_turtle=true)
    elseif mode == "46_direct"
        return (sectors=46, all_in_turtle=false)
    end
    error("Unsupported sky_mode $(repr(mode))")
end

function _apply_case_options(options0::ArchimedLight.LightOptions, options::OrderedDict{String,Any})
    sky = _sky_mode_params(String(options["sky_mode"]))
    ArchimedLight.LightOptions(
        options0;
        turtle_sectors=sky.sectors,
        all_in_turtle=sky.all_in_turtle,
        scattering=Bool(options["scattering"]),
        pixel_size=Float64(options["pixel_size_m"]),
        toricity=Bool(options["toricity"]),
        area_ratio=Bool(options["area_ratio"]),
        cache_radiation=Bool(options["cache_radiation"]),
        cache_pixel_table=Bool(options["cache_pixel_table"]),
        radiation_timestep_minutes=Float64(options["radiation_timestep"]),
        nir_interception=Bool(options["nir_interception"]),
        nir_scattering=Bool(options["nir_scattering"]),
        java_logged_turtle_dirs=Bool(options["java_logged_turtle_dirs"]),
    )
end

function _scattering_mode(options::OrderedDict{String,Any})
    value = options["scattering_mode"]
    return value isa Symbol ? value : Symbol(value)
end

function _fast_fixture_source(source_id::String)
    get!(_FAST_FIXTURE_CACHE, source_id) do
        case_root, input_path =
            if source_id == "coffee_example_2"
                root = normpath(joinpath(@__DIR__, "..", "..", "example_2"))
                (root, joinpath(root, "config.yml"))
            else
                root = joinpath(@__DIR__, "..", "fast_fixtures", source_id)
                (root, joinpath(root, "input"))
            end
        fixture = load_fixture_inputs(input_path)
        merge((case_root=case_root,), fixture)
    end
end

function _core_matrix_option_sets()
    rows = OrderedDict{String,Any}[]
    for sky_mode in ("1_direct", "16_all_in_turtle", "46_direct"),
        toricity in (false, true),
        area_ratio in (false, true),
        scattering in (false, true),
        pixel_size_m in (0.01, 0.05, 0.40)
        push!(
            rows,
            OrderedDict{String,Any}(
                "sky_mode" => sky_mode,
                "toricity" => toricity,
                "area_ratio" => area_ratio,
                "scattering" => scattering,
                "pixel_size_m" => pixel_size_m,
            ),
        )
    end
    return rows
end

function _targeted_option_sets()
    return OrderedDict{String,Any}[
        OrderedDict("cache_radiation" => true),
        OrderedDict("cache_pixel_table" => true),
        OrderedDict("nir_interception" => false),
        OrderedDict("scattering" => true, "nir_scattering" => false),
        OrderedDict("radiation_timestep" => 5),
        OrderedDict("radiation_timestep" => 30),
        OrderedDict("java_logged_turtle_dirs" => true, "sky_mode" => "16_all_in_turtle"),
        OrderedDict("scattering" => true, "scattering_mode" => :raycast),
        OrderedDict("sky_mode" => "1_direct", "cache_pixel_table" => true, "toricity" => true),
    ]
end

function regression_cases(; profile::String=_regression_profile())
    cases = RegressionCase[]

    single_plate = RegressionScenario("single_plate_direct", :synthetic, "single_plate_direct")
    partial_overlap = RegressionScenario("partial_overlap_direct", :synthetic, "partial_overlap_direct")
    toricity_wrap = RegressionScenario("toricity_wraparound", :synthetic, "toricity_wraparound")
    stacked_scattering = RegressionScenario("stacked_scattering", :synthetic, "stacked_scattering")
    cached_series = RegressionScenario("cached_series_parity", :synthetic, "cached_series_parity")
    simpleplant = RegressionScenario("simpleplant", :fast_fixture, "simpleplant_16_notoric")
    sky_fixture = RegressionScenario("sky_fixture", :fast_fixture, "sky_46_direct")
    coffee_dense = RegressionScenario("coffee_dense_default", :fast_fixture, "coffee_example_2")

    _push_case!(cases, "strict_single_plate", single_plate, OrderedDict("sky_mode" => "1_direct", "pixel_size_m" => 0.01); strict=true)
    _push_case!(cases, "strict_partial_overlap", partial_overlap, OrderedDict("sky_mode" => "1_direct", "pixel_size_m" => 0.01); strict=true)
    _push_case!(cases, "strict_toricity_off", toricity_wrap, OrderedDict("sky_mode" => "1_direct", "pixel_size_m" => 0.01, "toricity" => false); strict=true)
    _push_case!(cases, "strict_toricity_on", toricity_wrap, OrderedDict("sky_mode" => "1_direct", "pixel_size_m" => 0.01, "toricity" => true); strict=true)
    _push_case!(cases, "strict_stacked_scattering_off", stacked_scattering, OrderedDict("sky_mode" => "1_direct", "scattering" => false, "pixel_size_m" => 0.01); strict=true)
    _push_case!(cases, "strict_stacked_scattering_on", stacked_scattering, OrderedDict("sky_mode" => "1_direct", "scattering" => true, "pixel_size_m" => 0.01); strict=true)
    _push_case!(cases, "strict_cached_series", cached_series, OrderedDict("sky_mode" => "46_direct", "pixel_size_m" => 0.01); strict=true)
    _push_case!(cases, "strict_simpleplant", simpleplant, OrderedDict("sky_mode" => "46_direct", "pixel_size_m" => 0.40, "area_ratio" => true, "toricity" => false); strict=true, visual=true)
    _push_case!(cases, "strict_sky_fixture", sky_fixture, OrderedDict("sky_mode" => "46_direct"); strict=true)
    _push_case!(
        cases,
        "strict_coffee_dense",
        coffee_dense,
        OrderedDict(
            "sky_mode" => "16_all_in_turtle",
            "toricity" => true,
            "area_ratio" => true,
            "scattering" => true,
            "pixel_size_m" => 0.003,
            "cache_radiation" => false,
            "cache_pixel_table" => false,
            "nir_interception" => true,
            "nir_scattering" => true,
            "radiation_timestep" => 5,
            "java_logged_turtle_dirs" => false,
            "scattering_mode" => :raycast,
        );
        strict=true,
    )

    for opts in _core_matrix_option_sets()
        _push_case!(cases, "matrix_partial_overlap", partial_overlap, opts)
    end

    for opts in _targeted_option_sets()
        _push_case!(cases, "matrix_simpleplant", simpleplant, opts)
        _push_case!(cases, "matrix_sky_fixture", sky_fixture, opts)
    end

    if profile == "release"
        append!(cases, _release_regression_cases())
    end

    seen = Set{String}()
    unique_cases = RegressionCase[]
    for case in cases
        case.id in seen && continue
        push!(seen, case.id)
        _regression_case_enabled(case.id) || continue
        push!(unique_cases, case)
    end
    return unique_cases
end

function _release_dataset_root()
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    artifacts_toml = joinpath(repo_root, "Artifacts.toml")
    artifact_name = get(ENV, "ARCHIMEDLIGHT_RELEASE_ARTIFACT_NAME", "archimedlight-release-fixtures")
    dataset_root = strip(get(ENV, "ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR", ""))
    if isempty(dataset_root) && isfile(artifacts_toml)
        hash = artifact_hash(artifact_name, artifacts_toml)
        if hash !== nothing
            ensure_artifact_installed(artifact_name, artifacts_toml)
            dataset_root = artifact_path(hash)
        end
    end
    isempty(dataset_root) && error(
        "Release dataset not found. Set ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR or register artifact $(repr(artifact_name)) in Artifacts.toml.",
    )
    isdir(dataset_root) || error("Release dataset directory does not exist: $(dataset_root)")
    manifest_path = joinpath(dataset_root, "fixtures_manifest.toml")
    isfile(manifest_path) || error("Missing fixtures manifest in release dataset: $(manifest_path)")
    return normpath(dataset_root)
end

function _ensure_release_harness_loaded!()
    if !isdefined(@__MODULE__, :JuliaFixture)
        old_root = get(ENV, "ARCHIMEDLIGHT_RELEASE_DATA_ROOT", nothing)
        ENV["ARCHIMEDLIGHT_RELEASE_DATA_ROOT"] = _release_dataset_root()
        try
            Base.include(@__MODULE__, joinpath(@__DIR__, "..", "release", "harness.jl"))
        finally
            if old_root === nothing
                delete!(ENV, "ARCHIMEDLIGHT_RELEASE_DATA_ROOT")
            else
                ENV["ARCHIMEDLIGHT_RELEASE_DATA_ROOT"] = old_root
            end
        end
    end
    return nothing
end

_release_func(name::Symbol) = Base.invokelatest(getfield, @__MODULE__, name)
_release_call(name::Symbol, args...; kwargs...) = Base.invokelatest(_release_func(name), args...; kwargs...)
_release_cache_pair_fixture_ids() = Set(["test-cached-radiation", "test-cached-radiation3"])

function _release_cached_pair_fixture(fx)
    cfg2 = joinpath(dirname(fx.config_path), "config2.yml")
    isfile(cfg2) || error("$(fx.id): missing paired cached-radiation config2.yml")
    return _release_call(
        :JuliaFixture,
        fx.id * "-config2",
        cfg2,
        fx.visual_metric,
        true,
        fx.scene_override,
        fx.meteo_override,
        fx.force_scattering,
    )
end

function _release_regression_cases()
    _ensure_release_harness_loaded!()
    cases = RegressionCase[]
    fixtures = _release_call(:select_fixtures, _release_call(:julia_fixtures))
    for fx in fixtures
        scenario = RegressionScenario("release_$(fx.id)", :release_fixture, fx.id)
        opts = _default_case_options()
        push!(cases, RegressionCase("release_fixture__$(fx.id)", scenario, opts, true, true))
    end
    return cases
end

function _synthetic_scene_for_source(source_id::String)
    with_ground(specs) = vcat(
        specs,
        [
            (x0=0.0, x1=1.0 / 3.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=0.0, x1=1.0 / 3.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=0.0, x1=1.0 / 3.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
        ],
    )

    if source_id == "single_plate_direct"
        return _synthetic_horizontal_scene(with_ground([
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
        ]))
    elseif source_id == "partial_overlap_direct"
        return _synthetic_horizontal_scene(with_ground([
            (x0=0.0, x1=0.5, y0=0.0, y1=1.0, z=1.0, group="upper_half", type="plate", object_id=1),
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower_full", type="plate", object_id=2),
        ]))
    elseif source_id == "toricity_wraparound"
        scene = _synthetic_horizontal_scene(with_ground([
            (x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1),
        ]))
        return PlantGeom.SceneGeometry(
            scene.mtg,
            scene.merged_mesh,
            scene.face2node,
            scene.nodes,
            scene.source_path,
            (0.0, 0.0, 1.0, 1.0),
        )
    elseif source_id == "stacked_scattering"
        return _synthetic_horizontal_scene(with_ground([
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
        ]))
    elseif source_id == "cached_series_parity"
        return _synthetic_horizontal_scene(with_ground([
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
        ]))
    end
    error("Unsupported synthetic source $(repr(source_id))")
end

function _synthetic_meteo_for_source(source_id::String)
    if source_id == "cached_series_parity"
        rows = [
            _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
            _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
            _synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
        ]
        return ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_cached_series"))
    elseif source_id == "toricity_wraparound"
        return _synthetic_meteo_row(; duration_seconds=1.0, ri_par_f=100.0, ri_nir_f=0.0, direct_fraction=1.0, sun_azimut=270.0, sun_elevation=45.0)
    end
    return _synthetic_meteo_row(; duration_seconds=1.0, ri_par_f=100.0, ri_nir_f=0.0, direct_fraction=1.0)
end

function _compute_synthetic_case(case::RegressionCase)
    scene = _synthetic_scene_for_source(case.scenario.source_id)
    models = _default_synthetic_models()
    options = _apply_case_options(ArchimedLight.LightOptions(), case.options)
    scattering_mode = _scattering_mode(case.options)

    if case.scenario.source_id == "cached_series_parity"
        meteo = _synthetic_meteo_for_source(case.scenario.source_id)
        options_uncached = ArchimedLight.LightOptions(options; cache_radiation=false)
        options_cached = ArchimedLight.LightOptions(options; cache_radiation=true)
        series_uncached = ArchimedLight.run_light_series(scene, models, meteo, options_uncached; scattering_mode=scattering_mode)
        series_cached = ArchimedLight.run_light_series(scene, models, meteo, options_cached; scattering_mode=scattering_mode)
        diffs = OrderedDict{String,Float64}()
        for i in eachindex(series_uncached)
            diffs["step$(i)_ri_par_f"] = _max_abs_float_dict_diff(series_uncached[i].budget.incident_flux.total.par, series_cached[i].budget.incident_flux.total.par)
            diffs["step$(i)_ri_nir_f"] = _max_abs_float_dict_diff(series_uncached[i].budget.incident_flux.total.nir, series_cached[i].budget.incident_flux.total.nir)
            diffs["step$(i)_ri_par_q"] = _max_abs_float_dict_diff(series_uncached[i].budget.incident_energy.total.par, series_cached[i].budget.incident_energy.total.par)
            diffs["step$(i)_ri_nir_q"] = _max_abs_float_dict_diff(series_uncached[i].budget.incident_energy.total.nir, series_cached[i].budget.incident_energy.total.nir)
        end
        strict_result = _synthetic_exact_check(case.scenario.source_id, (diffs=diffs,))
        return (
            kind=:scene_outputs,
            scene=scene,
            models=models,
            options=options_cached,
            series=series_cached,
            meteo_rows=collect(meteo),
            figure=nothing,
            strict_result=strict_result,
            meta=OrderedDict{String,Any}("cache_radiation" => true),
        )
    end

    meteo_row = _synthetic_meteo_for_source(case.scenario.source_id)
    step = ArchimedLight.run_light_step(scene, models, meteo_row, options; scattering_mode=scattering_mode)
    strict_result = _synthetic_exact_check(
        case.scenario.source_id,
        (step=step, meta=OrderedDict{String,Any}("toricity" => case.options["toricity"], "scattering" => case.options["scattering"])),
    )
    return (
        kind=:scene_outputs,
        scene=scene,
        models=models,
        options=options,
        series=[step],
        meteo_rows=[meteo_row],
        figure=nothing,
        strict_result=strict_result,
        meta=OrderedDict{String,Any}("source" => case.scenario.source_id),
    )
end

function _compute_fast_fixture_case(case::RegressionCase)
    src = _fast_fixture_source(case.scenario.source_id)
    options = _apply_case_options(src.options, case.options)
    scattering_mode = _scattering_mode(case.options)
    if startswith(case.scenario.source_id, "sky_")
        row = first(ArchimedLight.prepare_meteo(src.meteo, options))
        sky = ArchimedLight.compute_sky(row, options)
        turtle = ArchimedLight.build_turtle(options, sky)
        fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
        strict_ok =
            isapprox(sum(fluxes.par), sky.ri_par_f; atol=1e-6, rtol=1e-6) &&
            isapprox(sum(fluxes.nir), sky.ri_nir_f; atol=1e-6, rtol=1e-6)
        return (
            kind=:sky_outputs,
            options=options,
            sky=sky,
            turtle=turtle,
            fluxes=fluxes,
            sky_mode=String(case.options["sky_mode"]),
            figure=nothing,
            strict_result=(ok=strict_ok, detail=strict_ok ? "" : "flux sum mismatch"),
            meta=OrderedDict{String,Any}("fixture" => case.scenario.source_id),
        )
    end

    selected = ArchimedLight.prepare_meteo(src.meteo, options)
    series = ArchimedLight.run_light_series(src.scene, src.models, src.meteo, options; scattering_mode=scattering_mode)
    figure = case.visual ? _render_ri_par_f_figure(src.scene, src.models, options, first(series); title="$(case.scenario.source_id) | $(case.id)") : nothing
    strict_result = (ok=true, detail="")
    return (
        kind=:scene_outputs,
        scene=src.scene,
        models=src.models,
        options=options,
        series=series,
        meteo_rows=collect(selected),
        figure=figure,
        strict_result=strict_result,
        meta=OrderedDict{String,Any}("fixture" => case.scenario.source_id),
    )
end

function _compute_release_case(case::RegressionCase)
    _ensure_release_harness_loaded!()
    fx = _release_call(:fixture_by_id, case.scenario.source_id)
    fx === nothing && error("Unknown release fixture $(repr(case.scenario.source_id))")
    data = _release_call(:fixture_runtime_data, fx)
    if fx.id in _release_cache_pair_fixture_ids()
        pair_fx = _release_cached_pair_fixture(fx)
        pair_data = _release_call(:fixture_runtime_data, pair_fx)
        fig = _release_call(:render_fixture_montage, fx; data=data)
        return (
            kind=:release_pair_outputs,
            fx=fx,
            data=data,
            pair_fx=pair_fx,
            pair_data=pair_data,
            figure=fig,
            strict_result=(ok=true, detail=""),
            meta=OrderedDict{String,Any}("fixture" => fx.id, "pair_fixture" => pair_fx.id),
        )
    end
    fig = _release_call(:render_fixture_montage, fx; data=data)
    return (
        kind=:release_outputs,
        fx=fx,
        data=data,
        figure=fig,
        strict_result=(ok=true, detail=""),
        meta=OrderedDict{String,Any}("fixture" => fx.id),
    )
end

function _compute_case(case::RegressionCase)
    if case.scenario.source_kind == :synthetic
        return _compute_synthetic_case(case)
    elseif case.scenario.source_kind == :fast_fixture
        return _compute_fast_fixture_case(case)
    elseif case.scenario.source_kind == :release_fixture
        return _compute_release_case(case)
    end
    error("Unsupported source kind $(repr(case.scenario.source_kind))")
end

function _write_observed_outputs!(case::RegressionCase, observed_dir::AbstractString, data)
    files = Dict{String,String}()
    image_path = nothing

    if data.kind == :scene_outputs
        if case.scenario.source_id != "coffee_example_2"
            comp_path = joinpath(observed_dir, "component_values.csv")
            _write_component_series_csv(comp_path, data.scene, data.models, data.series, data.options, data.meteo_rows)
            files["component_values.csv"] = comp_path
        end

        scene_path = joinpath(observed_dir, "scene_values.csv")
        _write_scene_series_csv(scene_path, data.scene, data.series, data.meteo_rows)
        files["scene_values.csv"] = scene_path

        if data.figure !== nothing
            image_path = _save_figure_png(joinpath(observed_dir, "ri_par_f.png"), data.figure)
        end
    elseif data.kind == :sky_outputs
        sky_summary_path = joinpath(observed_dir, "sky_summary.csv")
        _write_sky_summary_csv(sky_summary_path, data.sky, data.turtle, data.fluxes, data.options; step_number=0, sky_mode=data.sky_mode)
        files["sky_summary.csv"] = sky_summary_path
        sector_path = joinpath(observed_dir, "sector_flux.csv")
        _write_sector_flux_csv(sector_path, data.turtle, data.fluxes; step_number=0)
        files["sector_flux.csv"] = sector_path
    elseif data.kind == :release_outputs
        files = _release_call(:write_fixture_observed_outputs!, data.fx, observed_dir; data=data.data)
        image_path = _save_figure_png(joinpath(observed_dir, "$(data.fx.id)_montage.png"), data.figure)
    elseif data.kind == :release_pair_outputs
        out1 = joinpath(observed_dir, "config")
        out2 = joinpath(observed_dir, "config2")
        files1 = _release_call(:write_fixture_observed_outputs!, data.fx, out1; data=data.data)
        files2 = _release_call(:write_fixture_observed_outputs!, data.pair_fx, out2; data=data.pair_data)
        image_path = _save_figure_png(
            joinpath(observed_dir, "$(data.fx.id)_montage.png"),
            data.figure,
        )
        return (files=files1, pair_files=files2, image_path=image_path)
    else
        error("Unsupported observed data kind $(repr(data.kind))")
    end

    return (files=files, image_path=image_path)
end

function _case_baseline_dir(case::RegressionCase)
    if case.scenario.source_kind == :release_fixture
        _ensure_release_harness_loaded!()
        fx = _release_call(:fixture_by_id, case.scenario.source_id)
        fx === nothing && error("Unknown release fixture $(repr(case.scenario.source_id))")
        return _release_call(:fixture_reference_dir, fx)
    end
    return joinpath(_baseline_root(), _case_path_token(case.id))
end

function _case_report_dir(case::RegressionCase)
    return joinpath(_regression_report_root(), "cases", _case_path_token(case.id))
end

function _case_path_token(id::AbstractString)
    digest = bytes2hex(sha1(codeunits(String(id))))[1:12]
    prefix = first(split(String(id), "__"))
    short_prefix = first(prefix, min(length(prefix), 24))
    return "$(short_prefix)__$(digest)"
end

function _baseline_image_path(case::RegressionCase, observed)
    if case.scenario.source_kind == :release_fixture
        _ensure_release_harness_loaded!()
        fx = _release_call(:fixture_by_id, case.scenario.source_id)
        fx === nothing && return nothing
        return _release_call(:fixture_reference_image_path, fx)
    elseif observed.image_path !== nothing
        return joinpath(_case_baseline_dir(case), basename(observed.image_path))
    end
    return nothing
end

function _compare_case_against_baseline(case::RegressionCase, data, observed, observed_dir::AbstractString)
    pair_cmp = nothing
    if data.kind == :release_pair_outputs
        name = "component_values.csv"
        obs_path = get(observed.files, name, "")
        pair_path = get(observed.pair_files, name, "")
        pair_cmp =
            if isempty(obs_path) || !isfile(obs_path)
                (ok=false, missing=0, extra=0, mismatch=0, detail="observed file missing: $(name)", max_abs_error=Inf, max_rel_error=Inf)
            elseif isempty(pair_path) || !isfile(pair_path)
                (ok=false, missing=0, extra=0, mismatch=0, detail="paired observed file missing: $(name)", max_abs_error=Inf, max_rel_error=Inf)
            else
                compare_stable_csv_paths(obs_path, pair_path; label="$(case.id):$(name)")
            end
    end

    baseline_dir = _case_baseline_dir(case)
    update = _regression_update_enabled() && case.scenario.source_kind != :release_fixture
    if update
        _update_baseline_dir!(baseline_dir, observed_dir)
    elseif !isdir(baseline_dir)
        return (
            ok=!case.strict,
            status=case.strict ? "missing_baseline" : "report_missing_baseline",
            detail="baseline directory missing",
            total_missing=0,
            total_extra=0,
            total_mismatch=0,
            max_abs_error=0.0,
            max_rel_error=0.0,
            psnr=missing,
        )
    end

    total_missing = 0
    total_extra = 0
    total_mismatch = 0
    max_abs_error = 0.0
    max_rel_error = 0.0
    details = String[]
    ok = true

    for (name, obs_path) in sort(collect(observed.files); by=first)
        exp_path = joinpath(baseline_dir, name)
        if !isfile(exp_path)
            ok &= !case.strict
            push!(details, "missing baseline file $(name)")
            continue
        end
        cmp = compare_stable_csv_paths(exp_path, obs_path; label="$(case.id):$(name)")
        ok &= cmp.ok || !case.strict
        total_missing += cmp.missing
        total_extra += cmp.extra
        total_mismatch += cmp.mismatch
        max_abs_error = max(max_abs_error, cmp.max_abs_error)
        max_rel_error = max(max_rel_error, cmp.max_rel_error)
        !isempty(cmp.detail) && push!(details, "$(name): $(cmp.detail)")
    end

    psnr = missing
    if observed.image_path !== nothing && data.figure !== nothing
        exp_img = _baseline_image_path(case, observed)
        if exp_img !== nothing && isfile(exp_img)
            img_cmp = _compare_png_reference(exp_img, data.figure)
            psnr = img_cmp.psnr
            ok &= img_cmp.ok || !case.strict
            img_cmp.ok || push!(details, "image PSNR $(img_cmp.psnr) < $(img_cmp.threshold)")
        elseif case.strict
            ok = false
            push!(details, "missing baseline image")
        end
    end

    if pair_cmp !== nothing
        ok &= pair_cmp.ok || !case.strict
        total_missing += pair_cmp.missing
        total_extra += pair_cmp.extra
        total_mismatch += pair_cmp.mismatch
        max_abs_error = max(max_abs_error, pair_cmp.max_abs_error)
        max_rel_error = max(max_rel_error, pair_cmp.max_rel_error)
        pair_cmp.ok || push!(details, "cached radiation pair: $(pair_cmp.detail)")
    end

    status =
        if update
            "updated"
        elseif case.strict
            ok ? "strict_pass" : "strict_fail"
        else
            ok ? "report_ok" : "report_drift"
        end
    return (
        ok=ok,
        status=status,
        detail=join(details, " | "),
        total_missing=total_missing,
        total_extra=total_extra,
        total_mismatch=total_mismatch,
        max_abs_error=max_abs_error,
        max_rel_error=max_rel_error,
        psnr=psnr,
    )
end

function _timing_samples_for_case(case::RegressionCase)
    if case.scenario.source_kind == :release_fixture
        return 1
    elseif case.strict
        return 2
    end
    return 2
end

function _report_row(case::RegressionCase, timing, cmp, strict_detail::String)
    row = OrderedDict{String,Any}()
    row["case_id"] = case.id
    row["scenario_id"] = case.scenario.id
    row["source_kind"] = String(case.scenario.source_kind)
    row["source_id"] = case.scenario.source_id
    row["strict"] = case.strict
    row["status"] = cmp.status
    row["runtime_median_seconds"] = timing.runtime_median_seconds
    row["alloc_median_bytes"] = timing.alloc_median_bytes
    row["runtime_min_seconds"] = timing.runtime_min_seconds
    row["alloc_min_bytes"] = timing.alloc_min_bytes
    row["total_missing"] = cmp.total_missing
    row["total_extra"] = cmp.total_extra
    row["total_mismatch"] = cmp.total_mismatch
    row["max_abs_error"] = cmp.max_abs_error
    row["max_rel_error"] = cmp.max_rel_error
    row["psnr"] = cmp.psnr
    row["detail"] = cmp.detail
    row["strict_detail"] = strict_detail
    for key in _KNOWN_OPTION_COLUMNS
        row[key] = get(case.options, key, missing)
    end
    return row
end

function run_regression_case!(case::RegressionCase)
    report_dir = _case_report_dir(case)
    observed_dir = joinpath(report_dir, "observed")
    mkpath(observed_dir)

    timing = _timed_median(() -> _compute_case(case); samples=_timing_samples_for_case(case))
    data = timing.value
    observed = _write_observed_outputs!(case, observed_dir, data)
    cmp = _compare_case_against_baseline(case, data, observed, observed_dir)

    strict_ok = getproperty(data.strict_result, :ok)
    strict_detail = getproperty(data.strict_result, :detail)
    cmp_ok = cmp.ok

    if case.strict
        @test strict_ok
        @test cmp_ok
    else
        @test true
    end

    return _report_row(case, timing, cmp, strict_detail)
end

function run_regression_matrix!()
    profile = _regression_profile()
    cases = regression_cases(; profile=profile)
    @test !isempty(cases)
    rows = OrderedDict{String,Any}[]
    @testset "Regression matrix ($(profile))" begin
        for case in cases
            @testset "$(case.id)" begin
                push!(rows, run_regression_case!(case))
            end
        end
    end
    report_path = joinpath(_regression_report_root(), "regression_report.csv")
    _write_report(report_path, rows)
    @info "Regression matrix report written" path=report_path cases=length(rows)
    return rows
end
