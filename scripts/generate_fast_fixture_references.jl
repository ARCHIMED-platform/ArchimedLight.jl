#!/usr/bin/env julia

using ArchimedLight
using CSV
using Tables
using CairoMakie

function render_ri_par_f_figure(scene, step, cfg; title::String)
    vertices, faces, face2node, _, _, _ = ArchimedLight._scene_geometry_for_interception(scene, cfg)
    metric = step.budget.ri_par_f_per_node

    v_sum = zeros(Float64, length(vertices))
    v_count = zeros(Int, length(vertices))
    for i in eachindex(faces)
        f = faces[i]
        v = get(metric, face2node[i], NaN)
        isfinite(v) || continue
        for vid in (Int(f[1]), Int(f[2]), Int(f[3]))
            v_sum[vid] += v
            v_count[vid] += 1
        end
    end
    vertex_values = Float64[v_count[i] > 0 ? (v_sum[i] / v_count[i]) : NaN for i in eachindex(vertices)]

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
        vertices,
        faces;
        color=vertex_values,
        colormap=:viridis,
        colorrange=(0.0, max(step.sky.ri_par_f, eps(Float64))),
        nan_color=:lightgray,
        shading=false,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    Colorbar(fig[1, 2], p, label="Ri_PAR_f")
    return fig
end

function write_case(case_name::String)
    case_root = joinpath(dirname(@__DIR__), "test", "fast_fixtures", case_name)
    cfg = ArchimedLight.read_light_config(joinpath(case_root, "input", "config.yml"))
    scene = ArchimedLight.read_scene(cfg.source_files.scene)
    meteo = ArchimedLight.read_meteo(cfg.source_files.meteo)
    selected = ArchimedLight.prepare_meteo(meteo, cfg)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)

    step = first(series)
    meteo_row = first(selected.rows)
    expected_root = joinpath(case_root, "expected")
    mkpath(expected_root)

    csv_path = joinpath(expected_root, "component_values.csv")
    ArchimedLight.write_component_values_csv(
        csv_path,
        scene,
        step,
        cfg;
        meteo_row=meteo_row,
        step_number=0,
        columns=["step_number", "node_id", "area", "Ri_PAR_0_q"],
        strict=false,
    )

    png_path = joinpath(expected_root, "ri_par_f_step0.png")
    fig = render_ri_par_f_figure(scene, step, cfg; title="$(case_name) | Ri_PAR_f")
    save(png_path, fig)

    println("updated ", case_name)
    println("  ", csv_path)
    println("  ", png_path)
end

for case_name in ("simpleplant_16_notoric", "simpleplant_16_toric")
    write_case(case_name)
end
