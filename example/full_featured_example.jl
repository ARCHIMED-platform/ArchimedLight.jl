using ArchimedLight
using PlantGeom
using CairoMakie
using MultiScaleTreeGraph
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

function main()
    base = @__DIR__

    cfg = read_light_config(joinpath(base, "config.yml"))
    scene = read_scene(cfg.scene)
    meteo = read_meteo(cfg.meteo)
    row = first(meteo.rows)

    # Stage-by-stage execution with explicit backend objects.
    sky = compute_sky(row, cfg)
    turtle = build_turtle(cfg, sky)
    fluxes = compute_directional_fluxes(sky, turtle, cfg)
    first_order = compute_first_order(scene, turtle, fluxes, cfg; backend=RasterCPUBackend())
    graph = build_scattering_transfer_graph(scene, turtle, first_order, cfg; backend=RaycastScatteringBackend())
    scat = compute_scattering(graph, first_order, cfg; backend=RaycastScatteringBackend())
    budget = integrate_light(first_order, scat, cfg)

    # Single-call pipeline with backend kwargs.
    step = run_light_step(
        scene,
        row,
        cfg;
        interception_backend=RasterCPUBackend(),
        scattering_backend=RaycastScatteringBackend(),
    )

    # Sanity checks: staged and pipeline outputs should match exactly on CPU.
    @assert max_abs_float_dict_diff(first_order.projected_area_per_node, step.first_order.projected_area_per_node) == 0.0
    @assert max_abs_int_dict_diff(first_order.hits_per_node, step.first_order.hits_per_node) == 0
    @assert max_abs_float_dict_diff(budget.ri_par_q_per_node, step.budget.ri_par_q_per_node) == 0.0
    @assert max_abs_float_dict_diff(budget.ri_nir_q_per_node, step.budget.ri_nir_q_per_node) == 0.0

    # The fixture contains an extra band (RI_custom_f).
    @assert haskey(step.extra_band_irradiance, "CUSTOM")

    # Series execution (same inputs) with links backend selection.
    series = run_light_series(
        scene,
        meteo,
        cfg;
        interception_backend=RasterCPUBackend(),
        scattering_mode=:links,
        scattering_backend=LinksScatteringBackend(),
    )

    par_vals = collect(values(step.budget.ri_par_q_per_node))
    nir_vals = collect(values(step.budget.ri_nir_q_per_node))

    println("ArchimedLight full example completed")
    println("- meteo rows: ", length(meteo.rows))
    println("- turtle sectors: ", length(step.turtle.sectors))
    println("- nodes in budget: ", length(step.budget.ri_par_q_per_node))
    println("- custom bands: ", join(sort(collect(keys(step.extra_band_irradiance))), ", "))
    println("- mean Ri_PAR_q: ", mean(par_vals))
    println("- mean Ri_NIR_q: ", mean(nir_vals))
    println("- series steps: ", length(series))

    # Map intercepted PAR back to MTG nodes for visualization.
    par_q_by_component = Dict{Int,Float64}()
    for (nid, q) in step.budget.ri_par_q_per_node
        comp_id = get(scene.java_component_id_per_node, nid, nothing)
        comp_id === nothing && continue
        par_q_by_component[Int(comp_id)] = Float64(q)
    end
    traverse!(scene.mtg) do node
        comp_id = Int(node_id(node)) + 1
        node[:Ri_PAR_q] = get(par_q_by_component, comp_id, nothing)
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
end

main()
