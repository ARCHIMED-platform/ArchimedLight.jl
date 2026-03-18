const _DEFAULT_BUDGET_ATTRS = Dict{Symbol,Symbol}(
    :ri_par_0_f_per_node => :Ri_PAR_0_f,
    :ri_nir_0_f_per_node => :Ri_NIR_0_f,
    :ri_par_f_per_node => :Ri_PAR_f,
    :ri_nir_f_per_node => :Ri_NIR_f,
    :ri_par_0_q_per_node => :Ri_PAR_0_q,
    :ri_nir_0_q_per_node => :Ri_NIR_0_q,
    :ri_par_q_per_node => :Ri_PAR_q,
    :ri_nir_q_per_node => :Ri_NIR_q,
    :ra_par_0_f_per_node => :Ra_PAR_0_f,
    :ra_nir_0_f_per_node => :Ra_NIR_0_f,
    :ra_par_f_per_node => :Ra_PAR_f,
    :ra_nir_f_per_node => :Ra_NIR_f,
    :ra_par_0_q_per_node => :Ra_PAR_0_q,
    :ra_nir_0_q_per_node => :Ra_NIR_0_q,
    :ra_par_q_per_node => :Ra_PAR_q,
    :ra_nir_q_per_node => :Ra_NIR_q,
)

function _attach_node_values!(
    mtg,
    node_ids,
    attr::Symbol,
    values::AbstractDict{<:Integer};
    fill_value=nothing,
)
    MultiScaleTreeGraph.traverse!(mtg) do node
        nid = Int(MultiScaleTreeGraph.node_id(node))
        if nid in node_ids
            node[attr] = get(values, nid, fill_value)
        end
        return true
    end
    return mtg
end

function _geometry_node_ids(scene::SceneGeometry)
    Set{Int}(keys(scene.total_area_per_node))
end

function _budget_node_field(step::LightStepResult, field::Symbol)
    field in propertynames(step.budget) || error("Unknown LightBudget field: $field")
    values = getproperty(step.budget, field)
    values isa AbstractDict || error("LightBudget field `$field` is not a per-node dictionary.")
    return values
end

function _budget_attr_name(field::Symbol, names::AbstractDict{Symbol,Symbol})
    haskey(names, field) && return names[field]
    haskey(_DEFAULT_BUDGET_ATTRS, field) && return _DEFAULT_BUDGET_ATTRS[field]
    error("No default MTG attribute name for `$field`. Pass it explicitly in `names=Dict($field => :YourAttr)`.")
end

function _resolved_step_fields(step::LightStepResult, fields, names)
    Dict{Symbol,AbstractDict{Int,Float64}}(
        _budget_attr_name(field, names) => _budget_node_field(step, field) for field in fields
    )
end

function _resolved_series_fields(steps::AbstractVector{<:LightStepResult}, fields, names, fill_value::Float64)
    out = Dict{Symbol,Dict{Int,Vector{Float64}}}()
    for field in fields
        attr = _budget_attr_name(field, names)
        all_node_ids = Set{Int}()
        for step in steps
            union!(all_node_ids, keys(_budget_node_field(step, field)))
        end
        by_node = Dict{Int,Vector{Float64}}(nid => Float64[] for nid in all_node_ids)
        for step in steps
            vals = _budget_node_field(step, field)
            for nid in all_node_ids
                push!(by_node[nid], Float64(get(vals, nid, fill_value)))
            end
        end
        out[attr] = by_node
    end
    return out
end

function _paving_tile_mesh(points, faces_for_tile)
    used = sort!(unique(vcat([[Int(f[1]), Int(f[2]), Int(f[3])] for f in faces_for_tile]...)))
    remap = Dict{Int,Int}(old => new for (new, old) in enumerate(used))
    tile_points = GeometryBasics.Point3f[GeometryBasics.Point3f(points[idx]...) for idx in used]
    tile_faces = PlantGeom.Face3[
        PlantGeom.Face3(remap[Int(f[1])], remap[Int(f[2])], remap[Int(f[3])]) for f in faces_for_tile
    ]
    return GeometryBasics.Mesh(tile_points, tile_faces), tile_points
end

function _paving_nodes_data(scene::SceneGeometry, cfg::LightConfig; xy_bounds=nothing)
    plot_paving = _cfg_plot_paving(cfg)
    plot_paving > 0 || return NamedTuple[]

    raw_vertices = GeometryBasics.decompose(GeometryBasics.Point3, scene.merged_mesh)
    vertices = [StaticArrays.SVector{3,Float64}(v[1], v[2], v[3]) for v in raw_vertices]
    plotbox = _plotbox(scene, vertices, cfg.general.pixel_size)
    first_node = isempty(scene.total_area_per_node) ? 1 : (maximum(keys(scene.total_area_per_node)) + 1)
    paving_vertices, paving_faces, paving_face2node, _ = _paving_mesh(plotbox, plot_paving, first_node)

    faces_by_node = Dict{Int,Vector{PlantGeom.Face3}}()
    for (face, nid) in zip(paving_faces, paving_face2node)
        push!(get!(faces_by_node, nid, PlantGeom.Face3[]), face)
    end

    out = NamedTuple[]
    for nid in sort!(collect(keys(faces_by_node)))
        tile_mesh, tile_points = _paving_tile_mesh(paving_vertices, faces_by_node[nid])
        xs = Float64[p[1] for p in tile_points]
        ys = Float64[p[2] for p in tile_points]
        zs = Float64[p[3] for p in tile_points]
        bary = ((minimum(xs) + maximum(xs)) / 2, (minimum(ys) + maximum(ys)) / 2, sum(zs) / length(zs))
        if xy_bounds !== nothing
            xmin, ymin, xmax, ymax = xy_bounds
            xmin <= bary[1] <= xmax || continue
            ymin <= bary[2] <= ymax || continue
        end
        push!(out, (node_id=nid, mesh=tile_mesh, barycenter=bary))
    end
    return out
end

function _register_ref_mesh!(mtg, ref_mesh)
    ref_meshes = haskey(mtg, :ref_meshes) ? mtg[:ref_meshes] : nothing
    ref_meshes isa AbstractDict || return nothing
    next_id = isempty(ref_meshes) ? 0 : (maximum(Int.(keys(ref_meshes))) + 1)
    ref_meshes[next_id] = ref_mesh
    return next_id
end

function _append_paving_nodes!(
    mtg,
    scene::SceneGeometry,
    cfg::LightConfig,
    field_values::AbstractDict{Symbol,<:AbstractDict};
    fill_value=nothing,
    xy_bounds=nothing,
)
    tiles = _paving_nodes_data(scene, cfg; xy_bounds=xy_bounds)
    isempty(tiles) && return mtg

    root_scale = MultiScaleTreeGraph.scale(mtg)
    for (idx, tile) in enumerate(tiles)
        ref_mesh = PlantGeom.RefMesh("Cobblestone_$(tile.node_id)", tile.mesh)
        _register_ref_mesh!(mtg, ref_mesh)
        attrs = Dict{Symbol,Any}(
            :geometry => PlantGeom.Geometry(ref_mesh=ref_mesh),
            :functional_group => "pavement",
            :type => "Cobblestone",
            :item_id => -1,
            :source_topology_id => tile.node_id,
        )
        for (attr, values) in field_values
            attrs[attr] = get(values, tile.node_id, fill_value)
        end
        MultiScaleTreeGraph.Node(
            tile.node_id,
            mtg,
            MultiScaleTreeGraph.MutableNodeMTG(:+, :Cobblestone, idx, root_scale + 1),
            attrs,
        )
    end
    return mtg
end

function _drop_scene_surface!(mtg)
    haskey(mtg, :scene_dimensions) || return mtg
    haskey(mtg, :geometry) || return mtg
    mtg[:geometry] = nothing
    return mtg
end

"""
    attach_node_values!(scene, attr, values; fill_value=nothing)

Attach node-indexed values back onto the geometry nodes of `scene.mtg`.

`values` must be keyed by the internal scene node ids used by `SceneGeometry`.
The stored values may be scalars, vectors, or any other Julia object accepted as an MTG attribute.
This makes it possible to attach either a single simulated output or a time series per node.
"""
function attach_node_values!(
    scene::SceneGeometry,
    attr::Symbol,
    values::AbstractDict{<:Integer};
    fill_value=nothing,
)
    _attach_node_values!(scene.mtg, _geometry_node_ids(scene), attr, values; fill_value=fill_value)
end

"""
    attach_light_step!(scene, step; fields=[:ri_par_f_per_node], names=Dict(), fill_value=nothing)

Attach one or more per-node `LightBudget` fields from a single `LightStepResult`
onto `scene.mtg`, using Java-style attribute names by default.

Example:

```julia
attach_light_step!(scene, step; fields=[:ri_par_f_per_node])
plantviz(scene.mtg, color=:Ri_PAR_f)
```
"""
function attach_light_step!(
    scene::SceneGeometry,
    step::LightStepResult;
    fields::AbstractVector{Symbol}=[:ri_par_f_per_node],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value=nothing,
)
    for field in fields
        attach_node_values!(
            scene,
            _budget_attr_name(field, names),
            _budget_node_field(step, field);
            fill_value=fill_value,
        )
    end
    return scene.mtg
end

"""
    visual_scene_mtg(scene, cfg; include_paving=true, keep_scene_surface=false, xy_bounds=nothing)
    visual_scene_mtg(scene, cfg, step; fields=[:ri_par_f_per_node], names=Dict(), include_paving=true, keep_scene_surface=false, xy_bounds=nothing, fill_value=nothing)
    visual_scene_mtg(scene, cfg, steps; fields=[:ri_par_f_per_node], names=Dict(), include_paving=true, keep_scene_surface=false, xy_bounds=nothing, fill_value=NaN)

Return a copied MTG prepared for visualization with `plantviz`.

The copy keeps the original plant geometry, optionally materializes ARCHIMED paving as
`Cobblestone` geometry nodes, and can attach one or several simulated light outputs directly
onto the nodes so `plantviz(..., color=:Ri_PAR_f)` works on the full simulated scene.

`xy_bounds=(xmin, ymin, xmax, ymax)` can be used to keep only a local paving neighborhood in
the visualization copy while leaving the plant untouched.
"""
function visual_scene_mtg(
    scene::SceneGeometry,
    cfg::LightConfig;
    include_paving::Bool=true,
    keep_scene_surface::Bool=false,
    xy_bounds=nothing,
)
    mtg = deepcopy(scene.mtg)
    keep_scene_surface || _drop_scene_surface!(mtg)
    include_paving && _append_paving_nodes!(mtg, scene, cfg, Dict{Symbol,Dict{Int,Float64}}(); xy_bounds=xy_bounds)
    return mtg
end

function visual_scene_mtg(
    scene::SceneGeometry,
    cfg::LightConfig,
    step::LightStepResult;
    fields::AbstractVector{Symbol}=[:ri_par_f_per_node],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    include_paving::Bool=true,
    keep_scene_surface::Bool=false,
    xy_bounds=nothing,
    fill_value=nothing,
)
    mtg = deepcopy(scene.mtg)
    field_values = _resolved_step_fields(step, fields, names)
    for (attr, values) in field_values
        _attach_node_values!(mtg, _geometry_node_ids(scene), attr, values; fill_value=fill_value)
    end
    keep_scene_surface || _drop_scene_surface!(mtg)
    include_paving && _append_paving_nodes!(mtg, scene, cfg, field_values; fill_value=fill_value, xy_bounds=xy_bounds)
    return mtg
end

function visual_scene_mtg(
    scene::SceneGeometry,
    cfg::LightConfig,
    steps::AbstractVector{<:LightStepResult};
    fields::AbstractVector{Symbol}=[:ri_par_f_per_node],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    include_paving::Bool=true,
    keep_scene_surface::Bool=false,
    xy_bounds=nothing,
    fill_value::Float64=NaN,
)
    mtg = deepcopy(scene.mtg)
    field_values = _resolved_series_fields(steps, fields, names, fill_value)
    for (attr, values) in field_values
        _attach_node_values!(mtg, _geometry_node_ids(scene), attr, values; fill_value=fill(fill_value, length(steps)))
    end
    keep_scene_surface || _drop_scene_surface!(mtg)
    include_paving && _append_paving_nodes!(mtg, scene, cfg, field_values; fill_value=fill(fill_value, length(steps)), xy_bounds=xy_bounds)
    return mtg
end

"""
    attach_light_series!(scene, steps; fields=[:ri_par_f_per_node], names=Dict(), fill_value=NaN)

Attach one or more `LightBudget` fields from several time steps onto `scene.mtg`.
Each attached attribute becomes a vector ordered like `steps`, which can then be
visualized with `plantviz(..., color=:Attr, index=timestep)`.
"""
function attach_light_series!(
    scene::SceneGeometry,
    steps::AbstractVector{<:LightStepResult};
    fields::AbstractVector{Symbol}=[:ri_par_f_per_node],
    names::AbstractDict{Symbol,Symbol}=Dict{Symbol,Symbol}(),
    fill_value::Float64=NaN,
)
    node_ids = sort!(collect(keys(scene.total_area_per_node)))
    for field in fields
        by_node = Dict{Int,Vector{Float64}}(nid => Float64[] for nid in node_ids)
        for step in steps
            values = _budget_node_field(step, field)
            for nid in node_ids
                push!(by_node[nid], Float64(get(values, nid, fill_value)))
            end
        end
        attach_node_values!(
            scene,
            _budget_attr_name(field, names),
            by_node;
            fill_value=fill(fill_value, length(steps)),
        )
    end
    return scene.mtg
end
