#!/usr/bin/env julia

using ArchimedLight
using CairoMakie
using Dates
using GeometryBasics
using MultiScaleTreeGraph
using PlantGeom
using OrderedCollections: OrderedDict
using StaticArrays: SVector

include(joinpath(dirname(@__DIR__), "test", "regression_matrix", "synthetic_support.jl"))

function _render_metric!(slot, scene, models, options, step; title::String, colorrange)
    geometry = ArchimedLight._scene_geometry_for_interception(scene, models, options)
    metric = step.budget.incident_flux.total.par

    v_sum = zeros(Float64, length(geometry.vertices))
    v_count = zeros(Int, length(geometry.vertices))
    for i in eachindex(geometry.faces)
        face = geometry.faces[i]
        value = get(metric, geometry.face2node[i], NaN)
        isfinite(value) || continue
        for vid in (Int(face[1]), Int(face[2]), Int(face[3]))
            v_sum[vid] += value
            v_count[vid] += 1
        end
    end
    vertex_values = Float64[v_count[i] > 0 ? v_sum[i] / v_count[i] : NaN for i in eachindex(geometry.vertices)]

    ax = Axis3(
        slot;
        title=title,
        aspect=:data,
        azimuth=1.45,
        elevation=0.35,
        perspectiveness=0.65,
    )
    plot = mesh!(
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
    return plot
end

function _scene_specs_for_source(source_id::String)
    if source_id == "single_plate_direct"
        return [
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
        ], nothing
    elseif source_id == "partial_overlap_direct"
        return [
            (x0=0.0, x1=0.5, y0=0.0, y1=1.0, z=1.0, group="upper_half", type="plate", object_id=1),
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower_full", type="plate", object_id=2),
        ], nothing
    elseif source_id == "toricity_wraparound"
        return [
            (x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1),
        ], (0.0, 0.0, 1.0, 1.0)
    elseif source_id == "stacked_scattering"
        return [
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
        ], nothing
    elseif source_id == "cached_series_parity"
        return [
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
        ], nothing
    end
    error("Unsupported synthetic source $(repr(source_id))")
end

function _synthetic_scene_for_source(source_id::String; include_ground::Bool=false)
    specs, bounds = _scene_specs_for_source(source_id)
    all_specs = collect(specs)
    if include_ground
        append!(all_specs, [
            (x0=0.0, x1=1.0 / 3.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=0.0, y1=1.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=0.0, x1=1.0 / 3.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=1.0 / 3.0, y1=2.0 / 3.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=0.0, x1=1.0 / 3.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=1.0 / 3.0, x1=2.0 / 3.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
            (x0=2.0 / 3.0, x1=1.0, y0=2.0 / 3.0, y1=1.0, z=0.0, group="pavement", type="Cobblestone", object_id=-1),
        ])
    end
    scene = _synthetic_horizontal_scene(all_specs)
    if bounds !== nothing
        return PlantGeom.SceneGeometry(
            scene.mtg,
            scene.merged_mesh,
            scene.face2node,
            scene.nodes,
            scene.source_path,
            bounds,
        )
    end
    return scene
end

function _synthetic_meteo_for_source(source_id::String)
    if source_id == "cached_series_parity"
        rows = [
            _synthetic_meteo_row(; date=Date(2020, 6, 21), start_time=Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
            _synthetic_meteo_row(; date=Date(2020, 6, 21), start_time=Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
            _synthetic_meteo_row(; date=Date(2020, 6, 22), start_time=Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
        ]
        return ArchimedLight.MeteoTable(rows, (; source="synthetic_cached_series"))
    elseif source_id == "toricity_wraparound"
        return _synthetic_meteo_row(; duration_seconds=1.0, ri_par_f=100.0, ri_nir_f=0.0, direct_fraction=1.0, sun_azimut=270.0, sun_elevation=45.0)
    end
    return _synthetic_meteo_row(; duration_seconds=1.0, ri_par_f=100.0, ri_nir_f=0.0, direct_fraction=1.0)
end

function _default_options(; toricity::Bool=false, scattering::Bool=false, sky_mode::String="1_direct")
    sectors = sky_mode == "46_direct" ? 46 : 1
    ArchimedLight.LightOptions(
        turtle_sectors=sectors,
        toricity=toricity,
        scattering=scattering,
        area_ratio=true,
        pixel_size=0.01,
        cache_radiation=false,
        cache_pixel_table=false,
        nir_interception=true,
        nir_scattering=true,
        radiation_timestep_minutes=15.0,
        java_logged_turtle_dirs=false,
    )
end

function _render_pair(name::String, source_id::String; toricity::Bool=false, scattering::Bool=false, sky_mode::String="1_direct")
    models = _default_synthetic_models()
    options = _default_options(; toricity=toricity, scattering=scattering, sky_mode=sky_mode)

    current_scene = _synthetic_scene_for_source(source_id; include_ground=false)
    baseline_scene = _synthetic_scene_for_source(source_id; include_ground=true)

    if source_id == "cached_series_parity"
        meteo = _synthetic_meteo_for_source(source_id)
        current_step = first(ArchimedLight.run_light_series(current_scene, models, meteo, options))
        baseline_step = first(ArchimedLight.run_light_series(baseline_scene, models, meteo, options))
    else
        row = _synthetic_meteo_for_source(source_id)
        current_step = ArchimedLight.run_light_step(current_scene, models, row, options)
        baseline_step = ArchimedLight.run_light_step(baseline_scene, models, row, options)
    end

    max_flux = maximum((
        maximum(values(current_step.budget.incident_flux.total.par); init=0.0),
        maximum(values(baseline_step.budget.incident_flux.total.par); init=0.0),
        eps(Float64),
    ))

    fig = Figure(size=(1700, 760))
    plot1 = _render_metric!(
        fig[1, 1],
        baseline_scene,
        models,
        options,
        baseline_step;
        title="Baseline-like scene (with 3x3 ground)",
        colorrange=(0.0, max_flux),
    )
    _render_metric!(
        fig[1, 2],
        current_scene,
        models,
        options,
        current_step;
        title="Current scene (no ground)",
        colorrange=(0.0, max_flux),
    )
    Label(
        fig[0, 1:2],
        "$(name): old baseline topology vs current topology",
        fontsize=24,
    )
    Colorbar(fig[1, 3], plot1, label="Ri_PAR_f")

    outdir = joinpath(dirname(@__DIR__), "test", "regression_matrix", "reports", "latest", "visual_diffs")
    mkpath(outdir)
    outpath = joinpath(outdir, "$(name).png")
    save(outpath, fig)
    println(outpath)
    return outpath
end

function main()
    _render_pair("strict_single_plate_topology_diff", "single_plate_direct")
    _render_pair("strict_partial_overlap_topology_diff", "partial_overlap_direct")
    _render_pair("strict_toricity_off_topology_diff", "toricity_wraparound"; toricity=false)
    _render_pair("strict_toricity_on_topology_diff", "toricity_wraparound"; toricity=true)
    _render_pair("strict_stacked_scattering_off_topology_diff", "stacked_scattering"; scattering=false)
    _render_pair("strict_stacked_scattering_on_topology_diff", "stacked_scattering"; scattering=true)
    _render_pair("strict_cached_series_topology_diff", "cached_series_parity"; sky_mode="46_direct")
end

main()
