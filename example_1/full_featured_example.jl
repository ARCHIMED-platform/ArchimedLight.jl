const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight
using CSV
using PlantGeom
using CairoMakie
using MultiScaleTreeGraph
using Dates
using LinearAlgebra
using Statistics

function max_abs_float_dict_diff(a::Dict{Int,Float64}, b::Dict{Int,Float64})
    ids = union(keys(a), keys(b))
    isempty(ids) && return 0.0
    maximum(abs(get(a, id, 0.0) - get(b, id, 0.0)) for id in ids)
end

function max_abs_int_dict_diff(a::Dict{Int,Int}, b::Dict{Int,Int})
    ids = union(keys(a), keys(b))
    isempty(ids) && return 0
    maximum(abs(get(a, id, 0) - get(b, id, 0)) for id in ids)
end

to_float(x) = x isa Number ? Float64(x) : parse(Float64, string(x))
to_int(x) = x isa Integer ? Int(x) : Int(round(to_float(x)))

function triangle_area3d(p1, p2, p3)
    v1 = p2 - p1
    v2 = p3 - p1
    0.5 * norm(cross(v1, v2))
end

function node_area_from_interception_geometry(scene, models, options)
    geometry = ArchimedLight._scene_geometry_for_interception(scene, models, options)
    area = Dict{Int,Float64}()
    for i in eachindex(geometry.faces)
        f = geometry.faces[i]
        a = triangle_area3d(geometry.vertices[f[1]], geometry.vertices[f[2]], geometry.vertices[f[3]])
        nid = geometry.face2node[i]
        area[nid] = get(area, nid, 0.0) + a
    end
    area
end

function read_java_component_values(path::String)
    rows = collect(CSV.File(path; delim=';'))
    d = Dict{Tuple{Int,Int},NamedTuple}()
    for r in rows
        object_id = :object_id in propertynames(r) ? to_int(getproperty(r, :object_id)) : to_int(getproperty(r, :item_id))
        source_topology_id =
            :source_topology_id in propertynames(r) ? to_int(getproperty(r, :source_topology_id)) : to_int(getproperty(r, :component_id))
        key = (object_id, source_topology_id)
        names = propertynames(r)
        d[key] = (
            area=(:area in names ? to_float(getproperty(r, :area)) : NaN),
            Ri_PAR_0_q=(:Ri_PAR_0_q in names ? to_float(getproperty(r, :Ri_PAR_0_q)) : NaN),
            Ri_PAR_q=(:Ri_PAR_q in names ? to_float(getproperty(r, :Ri_PAR_q)) : NaN),
            Ri_custom_0_q=(:Ri_custom_0_q in names ? to_float(getproperty(r, :Ri_custom_0_q)) : NaN),
            Ri_custom_q=(:Ri_custom_q in names ? to_float(getproperty(r, :Ri_custom_q)) : NaN),
        )
    end
    return d
end

function build_julia_component_values(scene, models, options, step)
    key_by_node = ArchimedLight._interception_output_keys(scene, models, options)
    area_by_node = node_area_from_interception_geometry(scene, models, options)
    custom0 = get(step.budget.extra_initial_energy_per_band, "CUSTOM", Dict{Int,Float64}())
    customn = get(step.budget.extra_energy_per_band, "CUSTOM", Dict{Int,Float64}())
    d = Dict{Tuple{Int,Int},NamedTuple}()
    for (nid, key) in key_by_node
        d[key] = (
            area=get(area_by_node, nid, NaN),
            Ri_PAR_0_q=get(step.budget.incident_energy.initial.par, nid, NaN),
            Ri_PAR_q=get(step.budget.incident_energy.total.par, nid, NaN),
            Ri_custom_0_q=get(custom0, nid, NaN),
            Ri_custom_q=get(customn, nid, NaN),
        )
    end
    d
end

function write_component_comparison(out_csv::String, java_vals, julia_vals)
    keys_common = sort(collect(intersect(Set(keys(java_vals)), Set(keys(julia_vals)))))
    missing_java = setdiff(Set(keys(julia_vals)), Set(keys(java_vals)))
    missing_julia = setdiff(Set(keys(java_vals)), Set(keys(julia_vals)))

    cols = (:area, :Ri_PAR_0_q, :Ri_PAR_q, :Ri_custom_0_q, :Ri_custom_q)
    stats = Dict{Symbol,NamedTuple{(:max_abs,:max_rel,:mean_abs),Tuple{Float64,Float64,Float64}}}()

    mkpath(dirname(out_csv))
    open(out_csv, "w") do io
        println(io, "object_id;source_topology_id;java_area;julia_area;java_Ri_PAR_0_q;julia_Ri_PAR_0_q;java_Ri_PAR_q;julia_Ri_PAR_q;java_Ri_custom_0_q;julia_Ri_custom_0_q;java_Ri_custom_q;julia_Ri_custom_q;absdiff_Ri_PAR_q;reldiff_Ri_PAR_q")
        for key in keys_common
            jv = java_vals[key]
            lv = julia_vals[key]
            abs_par = abs(lv.Ri_PAR_q - jv.Ri_PAR_q)
            rel_par = abs_par / max(abs(jv.Ri_PAR_q), eps(Float64))
            println(
                io,
                string(
                    key[1], ";", key[2], ";",
                    jv.area, ";", lv.area, ";",
                    jv.Ri_PAR_0_q, ";", lv.Ri_PAR_0_q, ";",
                    jv.Ri_PAR_q, ";", lv.Ri_PAR_q, ";",
                    jv.Ri_custom_0_q, ";", lv.Ri_custom_0_q, ";",
                    jv.Ri_custom_q, ";", lv.Ri_custom_q, ";",
                    abs_par, ";", rel_par,
                ),
            )
        end
    end

    for c in cols
        absd = Float64[]
        reld = Float64[]
        for key in keys_common
            jv = getproperty(java_vals[key], c)
            lv = getproperty(julia_vals[key], c)
            if isfinite(jv) && isfinite(lv)
                d = abs(lv - jv)
                push!(absd, d)
                push!(reld, d / max(abs(jv), eps(Float64)))
            end
        end
        if !isempty(absd)
            stats[c] = (max_abs=maximum(absd), max_rel=maximum(reld), mean_abs=mean(absd))
        end
    end

    return (keys_common=keys_common, missing_java=missing_java, missing_julia=missing_julia, stats=stats)
end

function step_duration_seconds(row)
    if (:step_duration in propertynames(row))
        return to_float(getproperty(row, :step_duration))
    end
    if (:hour_start in propertynames(row)) && (:hour_end in propertynames(row))
        t0 = Time(string(getproperty(row, :hour_start)))
        t1 = Time(string(getproperty(row, :hour_end)))
        dt0 = DateTime(Date(2000, 1, 1), t0)
        dt1 = DateTime(Date(2000, 1, 1), t1)
        dt1 < dt0 && (dt1 += Day(1))
        return Dates.value(dt1 - dt0) / 1000.0
    end
    return 1.0
end

function main()
    base = @__DIR__

    config_path = joinpath(base, "config.yml")
    options, scene, meteo, models = read_config(config_path)
    row = first(prepare_meteo(meteo, options).rows)

    # Stage-by-stage execution with explicit backend objects.
    sky = compute_sky(row, options)
    turtle = build_turtle(options, sky)
    fluxes = compute_directional_fluxes(row, sky, turtle, options)
    first_order = compute_first_order(scene, models, turtle, fluxes, options; backend=RasterCPUBackend())
    scat = compute_scattering(scene, models, turtle, first_order, options; backend=RaycastScatteringBackend())
    budget = integrate_light(scene, models, first_order, scat, options; meteo_row=row)

    # Single-call pipeline with backend kwargs.
    step = run_light_step(
        scene,
        models,
        row,
        options;
        interception_backend=RasterCPUBackend(),
        scattering_backend=RaycastScatteringBackend(),
    )

    # Sanity checks: staged and pipeline outputs should match exactly on CPU.
    @assert max_abs_float_dict_diff(first_order.projected_area_per_node, step.first_order.projected_area_per_node) == 0.0
    @assert max_abs_int_dict_diff(first_order.hits_per_node, step.first_order.hits_per_node) == 0
    staged_pipeline_par_diff = max_abs_float_dict_diff(budget.incident_energy.total.par, step.budget.incident_energy.total.par)
    staged_pipeline_nir_diff = max_abs_float_dict_diff(budget.incident_energy.total.nir, step.budget.incident_energy.total.nir)

    # The fixture contains an extra band (RI_custom_f).
    @assert haskey(step.extra_band_irradiance, "CUSTOM")

    # Series execution (same inputs) with explicit scattering backend selection.
    series = run_light_series(
        scene,
        models,
        meteo,
        options;
        interception_backend=RasterCPUBackend(),
        scattering_backend=RaycastScatteringBackend(),
    )

    par_vals = collect(values(step.budget.incident_energy.total.par))
    nir_vals = collect(values(step.budget.incident_energy.total.nir))

    println("ArchimedLight full example completed")
    println("- meteo rows: ", length(meteo.rows))
    println("- turtle sectors: ", length(step.turtle.sectors))
    println("- nodes in budget: ", length(step.budget.incident_energy.total.par))
    println("- custom bands: ", join(sort(collect(keys(step.extra_band_irradiance))), ", "))
    println("- mean Ri_PAR_q: ", mean(par_vals))
    println("- mean Ri_NIR_q: ", mean(nir_vals))
    println("- series steps: ", length(series))
    println("- staged vs pipeline max abs Ri_PAR_q diff: ", staged_pipeline_par_diff)
    println("- staged vs pipeline max abs Ri_NIR_q diff: ", staged_pipeline_nir_diff)

    # Map intercepted PAR (per-step quantity) back to MTG nodes for visualization.
    par_q_by_source_topology = Dict{Int,Float64}()
    for (nid, q) in step.budget.incident_energy.total.par
        node = PlantGeom.scene_node(scene, nid)
        node === nothing && continue
        par_q_by_source_topology[node.source_topology_id] = Float64(q)
    end
    traverse!(scene.mtg) do node
        haskey(node, :source_topology_id) || return true
        source_topology_id = node[:source_topology_id]
        source_topology_id === nothing && return true
        node[:Ri_PAR_q] = get(par_q_by_source_topology, to_int(source_topology_id), nothing)
        return true
    end

    # 3D visualization using PlantGeom + CairoMakie, colored by intercepted PAR.
    CairoMakie.activate!()
    fig, _, p = plantviz(scene.mtg, color=:Ri_PAR_q)
    PlantGeom.colorbar(fig[1, 2], p)
    out_dir = joinpath(base, "output")
    mkpath(out_dir)
    out_png = joinpath(out_dir, "scene_3d_par_intercepted.png")
    CairoMakie.save(out_png, fig)
    println("- 3D visualization (PAR intercepted): ", out_png)

    # Compare Java vs Julia outputs if Java component_values.csv is available.
    java_component_path = joinpath(base, "output", "000001", "component_values.csv")
    if isfile(java_component_path)
        java_vals = read_java_component_values(java_component_path)
        julia_vals = build_julia_component_values(scene, models, options, step)
        cmp_path = joinpath(out_dir, "component_values_java_vs_julia.csv")
        cmp = write_component_comparison(cmp_path, java_vals, julia_vals)
        println("- Java vs Julia comparison CSV: ", cmp_path)
        println("- Components matched: ", length(cmp.keys_common))
        println("- Components missing in Java file: ", length(cmp.missing_java))
        println("- Components missing in Julia output: ", length(cmp.missing_julia))
        for c in (:area, :Ri_PAR_0_q, :Ri_PAR_q, :Ri_custom_0_q, :Ri_custom_q)
            haskey(cmp.stats, c) || continue
            s = cmp.stats[c]
            println("  * ", c, " | max_abs=", s.max_abs, " | max_rel=", s.max_rel, " | mean_abs=", s.mean_abs)
        end

        # Map Java Ri_PAR_q on MTG nodes and render a second comparison figure.
        java_par_by_source_topology = Dict{Int,Float64}()
        for ((object_id, source_topology_id), vals) in java_vals
            object_id == 1 || continue
            java_par_by_source_topology[source_topology_id] = vals.Ri_PAR_q
        end
        traverse!(scene.mtg) do node
            haskey(node, :source_topology_id) || return true
            source_topology_id = node[:source_topology_id]
            source_topology_id === nothing && return true
            node[:Ri_PAR_q_java] = get(java_par_by_source_topology, to_int(source_topology_id), nothing)
            return true
        end
        fig_java, _, p_java = plantviz(scene.mtg, color=:Ri_PAR_q_java)
        PlantGeom.colorbar(fig_java[1, 2], p_java)
        out_java_png = joinpath(out_dir, "scene_3d_par_java.png")
        CairoMakie.save(out_java_png, fig_java)
        println("- 3D visualization (Java PAR): ", out_java_png)
    else
        println("- Java component_values.csv not found at: ", java_component_path)
        println("  Run the reference implementation once and place outputs under example_1/output/000001 to enable comparison.")
    end
end

main()
