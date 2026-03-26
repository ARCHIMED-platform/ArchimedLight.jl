using ArchimedLight

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

function _source_topology_area_q(scene, models, options, step)
    out = Dict{Tuple{Int,Int},NamedTuple{(:area, :Ri_PAR_0_q),Tuple{Float64,Float64}}}()
    keys_by_node = ArchimedLight._interception_output_keys(scene, models, options)
    for nid in keys(keys_by_node)
        node = ArchimedLight.scene_node(scene, nid)
        node === nothing && continue
        object_id, source_topology_id = keys_by_node[nid]
        out[(source_topology_id, object_id)] = (
            area=node.area,
            Ri_PAR_0_q=get(step.budget.incident_energy.initial.par, nid, 0.0),
        )
    end
    out
end
