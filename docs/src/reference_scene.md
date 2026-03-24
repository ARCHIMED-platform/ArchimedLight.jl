# Scene Files And Semantics

The scene definition is spread across three related concepts:

- the plot, described in the `.ops` terrain line or passed as `scene_xy_bounds`
- the objects placed in that plot
- the per-component geometry and metadata stored in `.opf` or `.gwa`

## The Three Entry Formats

`read_scene(path)` accepts:

- `.ops`: open plant scene, usually the main ARCHIMED entry point
- `.opf`: plant topology plus geometry
- `.gwa`: geometry-only objects

In practice, `.ops` is the most common because it lets you place several items in one plot and assign them to functional groups.

## `.ops`: Plot, Placement, Functional Groups

Example:

```text
# T xmin ymin zmin xmax ymax flat
T 0 0 0 2 2 flat

#[Archimed] coffee

#sceneId objectId FilePath x y z scale inclinationAzimut inclinationAngle rotation
1	1	opf/coffee.opf	0	0	0	0	0	0	0
```

The essential parts are:

- terrain line `T ...`: defines the plot bounds used as the projection area
- `#[Archimed] group_name`: ARCHIMED-specific functional-group tag
- object rows: path, position, and object-level transform

### What The Terrain Line Means

Historically, the `.ops` terrain line serves two roles at once:

- it defines the usable plot
- it defines the support over which projected pixel tables are built

That second role is why this line matters so much for light calculations. The plot is not just a visualization box. It determines the projection domain.

### Functional Groups In `.ops`

ARCHIMED introduces a convention on top of the plain `.ops` format:

```text
#[Archimed] coffee
```

Every following object line inherits that functional group until another `#[Archimed] ...` line appears.
Those group names are how scene items are matched to model YAML files.

## `.opf`: Plant Topology Plus Geometry

An `.opf` stores:

- reference meshes
- materials for visualization
- attributes attached to topology nodes
- plant topology, usually as an MTG-like hierarchy

![Simple OPF illustration](assets/simple_opf.png)

For `ArchimedLight.jl`, the most important consequence is that an `.opf` can preserve organ-scale identity. That means the solver can keep outputs at the component level and attach them back to the corresponding nodes.

### What `ArchimedLight.jl` Extracts From The Scene

When you call `read_scene` or `prepare_scene`, the package turns the original scene into a `SceneGeometry` that contains:

- one merged mesh for fast geometric processing
- `face2node`, the map from triangles back to component ids
- node areas
- node barycentres
- functional group and type per node
- `source_topology_id` and `object_id`
- the scene XY bounds

This is the core runtime representation used by first-order interception and scattering.

## `.gwa`: Geometry Without Plant Topology

A `.gwa` is useful for:

- soil or paving geometry
- virtual structures such as walls, panels, or support frames
- simple geometric sensors or emitters

Unlike `.opf`, it does not carry a plant topology. That means it is lighter, but it also means fewer topological semantics are available.

## Scene Metadata That Matter For Light

### Plot Bounds

The solver needs an XY domain for rasterization. With file-based scenes, that comes from the `.ops` terrain line. With in-memory scenes, provide it explicitly:

```julia
scene = prepare_scene(mtg; scene_xy_bounds=(0.0, 0.0, 2.0, 2.0))
```

### `group`

This is the functional group used to resolve optical models.

### `type`

This is the component type name inside the functional group.
By default, `prepare_scene` derives it from node attributes if present, otherwise from the node symbol.

### `object_id`

Higher-level entity id used to regroup component outputs.

### `source_topology_id`

Stable component id imported from the original source when available.

## In-Memory Scenes

The same semantics apply if the scene never touches disk:

```julia
scene = prepare_scene(
    mtg;
    source_path="interactive.opf",
    scene_xy_bounds=(-1.0, -1.0, 1.0, 1.0),
)
```

The critical rule is that the MTG nodes with geometry must expose enough metadata for model lookup:

- a functional group, typically through `:functional_group`
- a type, either explicit or implied by the node symbol

### Interactive Equivalent

There is no mandatory scene file in an interactive workflow. The equivalent of
the `.ops` + `.opf` information is:

- the MTG topology and geometry you build in Julia
- the per-node metadata you attach to that MTG
- the explicit `scene_xy_bounds` you pass to `prepare_scene`

Typical pattern:

```julia
mtg = Node(MutableNodeMTG(:/, :Plant, 1, 1))
leaf = Node(mtg, MutableNodeMTG(:+, :Leaf, 1, 2))
leaf[:geometry] = some_geometry
leaf[:functional_group] = "coffee"
leaf[:object_id] = 1

scene = prepare_scene(
    mtg;
    source_path="interactive.opf",
    scene_xy_bounds=(-1.0, -1.0, 1.0, 1.0),
)
```

The correspondences are:

- `.ops` terrain line -> `scene_xy_bounds=(xmin, ymin, xmax, ymax)`
- `#[Archimed] coffee` -> `node[:functional_group] = "coffee"`
- object placement in `.ops` -> the geometry coordinates or transformations you build in Julia
- organ identity from `.opf` -> the MTG nodes you keep in memory

So the interactive question is not “which file format should I write?”, but
“does my MTG expose geometry, group, type, and plot bounds clearly enough for
`prepare_scene`?”

## Ground And Paving

Two mechanisms exist:

- `read_config` can materialize paving automatically when a model file contains `plot_paving`
- `add_ground!` can add explicit tiles to an already-loaded scene

Explicit ground geometry is especially useful when:

- you want inspectable soil outputs
- you want realistic scattering exchanges near the bottom of the canopy
- you want the ground exported with the final scene

## Practical Checks Before Running

- make sure the plot bounds are representative if `toricity=true`
- make sure each simulated node resolves to a valid `(group, type)` model pair
- use `.opf` when component identity matters and `.gwa` when geometry alone is enough
- keep the scene footprint consistent with the real repeated pattern when using toricity
