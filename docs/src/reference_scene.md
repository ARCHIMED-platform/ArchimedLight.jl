# Scene Reference

The light solver needs a scene description with four kinds of information:

- an XY projection domain
- geometry placed in that domain
- per-component metadata used for model matching
- stable node ids that let outputs be traced back to the original scene

Historically that information often comes from `.ops`, `.opf`, and `.gwa` files. In an interactive workflow, the same information can be created directly in Julia and passed to `prepare_scene`.

## What Matters In A Scene

Whatever the source format, the prepared [`SceneGeometry`](@ref) used by the
solver must be able to answer these questions:

- what horizontal domain should be rasterized?
- where is each geometric object located?
- which functional group does each geometric component belong to?
- which component type does each geometric component belong to?
- which plant or object does each component belong to?
- which original topology id did this component come from?

## Recommended Interactive Builder

For new interactive workflows, start with `light_scene`:

```julia
using FileIO, MeshIO

sensor_mesh = load("sensor.obj")

scene = light_scene(domain=(0.0, 0.0, 2.0, 2.0)) do s
    add_plant!(s, "plant.opf"; group="coffee", id=1, at=(0.0, 0.0, 0.0), rotate=(z=15.0,), deg=true)
    add_plant!(s, plant_mtg; group="banana", id=2, at=(1.0, 0.0, 0.0), scale=0.8)
    add_object!(s, sensor_mesh; group="sensor", type="panel", id=10, at=(0.5, 0.0, 1.2), scale=0.05)
    add_ground!(s; group="soil", type="ground", nx=20, ny=20)
end
```

This creates a prepared `SceneGeometry`. The important user-facing concepts are:

- `domain`: the XY plot footprint used by projection and toricity
- `add_plant!`: imports or places one plant, and assigns its functional group and object id
- `add_object!`: imports or places a non-plant MTG, `GeometryBasics` mesh, `.opf`, or `.gwa` object
- `group`: the high-level model key, such as `"coffee"` or `"soil"`
- `type`: the component type, usually read from the OPF/GWA/MTG node symbol or `:type` attribute
- `id`: the plant or object id used for grouping outputs

Both `add_plant!` and `add_object!` accept placement keywords:

```julia
add_plant!(s, plant_mtg; group="coffee", id=1, at=(1.0, 0.0, 0.0), scale=0.8)
add_object!(s, mesh; group="sensor", type="panel", id=10, rotate=(z=90.0,), deg=true)
```

Use `at` for translation, `scale` for uniform or axis-wise scaling, and
`rotate` for local rotations. Tuple rotations use fixed X, then Y, then Z order.
Named-tuple rotations preserve the field order:

```julia
add_object!(s, mesh; group="sensor", type="panel", id=10, rotate=(x=10, y=20, z=30), deg=true)
add_object!(s, mesh; group="sensor", type="panel", id=11, rotate=(y=20, z=30, x=10), deg=true)
```

For OPS-compatible placement, use scalar `scale` with `rotation`,
`inclination_azimut`, and `inclination_angle`.

For mesh files, let MeshIO/FileIO do the format-specific IO and pass the loaded
mesh to `add_object!`:

```julia
using FileIO, MeshIO
sensor_mesh = load("sensor.obj")
add_object!(s, sensor_mesh; group="sensor", type="panel", id=10)
```

If you already have a complete MTG, you can still use `prepare_scene` directly.

## 1. Plot Domain

The raster interception model needs an XY domain. The domain defines the minimum and maximum scene dimensions, which may be different from the bounds defined by the objects in the scene. For example, when we activate the toricity feature, it is common to define the scene as a voronoi tile with a single plant that represents the whole plot. The domain of the scene is then the tile footprint, which may be larger than the actual plant geometry, and is usually equal to the inter and intra-row spacing between plants.

The domain is used internally to build the projected pixel tables.

### From Files

In an `.ops`, the plot domain is usually defined by the terrain line:

```text
# T xmin ymin zmin xmax ymax flat
T 0 0 0 2 2 flat
```

This is what `read_scene` uses to recover the scene XY bounds.

### Dynamically In Julia

In a dynamic workflow, the equivalent is `scene_xy_bounds=` in `prepare_scene`:

```julia
scene = prepare_scene(
    mtg;
    scene_xy_bounds=(0.0, 0.0, 2.0, 2.0),
)
```

This is the direct equivalent of the `.ops` terrain rectangle.

If you want to add explicit soil or paving later on, the same bounds can also be
reused by `add_ground!`:

```julia
add_ground!(
    scene;
    nx=20, # Number of tiles in X direction
    ny=20, # Number of tiles in Y direction
    xy_bounds=(0.0, 0.0, 2.0, 2.0), # Define the domain here if not already defined in prepare_scene
    group="pavement",
    type="Cobblestone",
)
```

## 2. Object Placement And Geometry

The light model needs actual geometry in 3D space. It does not care whether
that geometry came from an imported file or from code, as long as the final MTG
contains geometric nodes with valid transformations.

### From Files

In an `.ops`, each object row places one imported object in the scene:

```text
#sceneId objectId FilePath x y z scale inclinationAzimut inclinationAngle rotation
1	1	opf/coffee.opf	0	0	0	0	0	0	0
1	2	opf/coffee.opf	1	0	0	0	0	0	0
```

Those rows provide:

- the object file to import
- its translation, scaling and orientation parameters
- an object-level id

The detailed organ geometry usually comes from the referenced `.opf` or `.gwa`.

### Dynamically In Julia

In a dynamic workflow, placement is simply the geometry and transformation you construct on the MTG nodes.

For example, a node can carry geometry directly:

```julia
leaf = Node(mtg, MutableNodeMTG(:+, :Leaf, 1, 2))
leaf[:geometry] = some_geometry
```

or a transformed geometry:

```julia
leaf[:geometry] = PlantGeom.Geometry(
    ref_mesh=leaf_ref_mesh,
    transformation=PlantGeom.Translation(0.5, 0.0, 0.3),
)
```

If you build the plant procedurally with PlantGeom's growth API, the same
placement information is implicit in the internode and leaf parameters you pass
to the constructors.

## 3. Functional Groups

The functional group is the high-level key used to resolve optical models. It is
one half of the `(group, type)` lookup pair.

### From Files

In `.ops`, ARCHIMED adds the convention:

```text
#[Archimed] coffee
```

All following object rows inherit that functional group until another
`#[Archimed] ...` line appears.

### Dynamically In Julia

In an MTG built in memory, the equivalent is to set the `:functional_group`
attribute on the relevant nodes:

```julia
plant[:functional_group] = "coffee"
```

This can be attached at the plant level or directly on geometric nodes,
depending on how your MTG is organized. What matters is that `prepare_scene`
can recover the right group for each geometric component.

Example with two plants carrying different groups:

```julia
coffee = Node(scene_root, MutableNodeMTG(:/, :Plant, 1, 1))
coffee[:functional_group] = "coffee"

banana = Node(scene_root, MutableNodeMTG(:/, :Plant, 2, 1))
banana[:functional_group] = "banana"
```

## 4. Component Types

The type is the lower-level key used to distinguish leaves, stems, paving,
sensors, and other components inside one functional group.

### From Files

In `.opf` or `.gwa`, the type typically comes from the topology node type or
from explicit type metadata stored on the geometric nodes.

### Dynamically In Julia

In a dynamic workflow, `prepare_scene` derives the type from:

1. explicit type-like attributes such as `:type`
2. otherwise the node symbol

So these are equivalent:

```julia
leaf = Node(axis, MutableNodeMTG(:+, :Leaf, 1, 2))
```

and

```julia
leaf = Node(axis, MutableNodeMTG(:+, :Organ, 1, 2))
leaf[:type] = "Leaf"
```

This is what lets one group such as `"coffee"` resolve different optical models
for `"Leaf"`, `"Internode"`, `"Fruit"`, or `"Sensor"`.

## 5. Object Identity

`object_id` is used to regroup several geometric components under the same plant
or object.

### From Files

In `.ops`, the object rows already carry an object id column.

### Dynamically In Julia

In an interactive workflow, you can assign it directly:

```julia
plant[:object_id] = 1
```

If you create several plants in one scene, giving them different `object_id`
values makes it easier to aggregate outputs by individual later on.

## 6. Stable Source Ids

`source_topology_id` is the stable id of the original source component when that
information is available.

### From Files

An imported `.opf` often already carries stable topology ids.

### Dynamically In Julia

You can also set them manually when your own simulator already has stable organ
ids:

```julia
leaf[:source_topology_id] = 42
```

If you do not provide them, `prepare_scene` and `add_ground!` create consistent
fallback ids automatically.

## 7. The Runtime Representation

Once the scene metadata and geometry are in place, `prepare_scene` converts the
original MTG into the dense representation used by the solver:

```julia
scene = prepare_scene(
    mtg;
    source_path="interactive.opf",
    scene_xy_bounds=(-1.0, -1.0, 1.0, 1.0),
)
```

The resulting [`SceneGeometry`](@ref) stores:

- one merged mesh for efficient geometric processing
- a `face2node` map linking triangles back to scene node ids
- node areas and barycentres
- per-node `group`, `type`, `object_id`, and `source_topology_id`
- the scene XY bounds

This is the object that the light pipeline actually consumes.

## Scene Containers

The file formats are still useful, but they are just ways to deliver the scene
semantics described above.

### `.ops`

Useful when you want one file to define:

- the plot domain
- the placement of several imported objects
- ARCHIMED functional-group tags

### `.opf`

Useful when you want:

- plant topology
- organ-scale identity
- geometry with reusable reference meshes
- outputs that can be attached back to plant nodes

### `.gwa`

Useful when you want:

- geometry without plant topology
- simple objects such as paving, walls, frames, or sensors
- lightweight geometric containers

## Ground And Paving

Ground is not a special scene concept. It is just geometry with a group and a
type, often something like `(group="pavement", type="Cobblestone")`.

You can obtain it in two ways:

- from files, for example with automatic paving materialized by `read_simulation`
- dynamically with `add_ground!` on an existing prepared scene

The dynamic route is often preferable when you want explicit inspectable ground
geometry in the final MTG.

## Practical Checks

Before running a simulation, check that:

- `scene_xy_bounds` really match the intended plot footprint
- every simulated component resolves to a valid `(group, type)` pair
- `object_id` is meaningful when several plants share the same scene
- `toricity=true` is used only with a footprint that really represents the
  repeated tile

If those semantics are correct, the solver does not particularly care whether
they came from `.ops` / `.opf` / `.gwa` files or from Julia code.
