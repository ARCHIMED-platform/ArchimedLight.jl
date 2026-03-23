# Tutorial: Interactive Workflow

This tutorial shows the in-memory workflow: build or modify a plant directly in Julia, convert it to a `SceneGeometry`, then run the same light pipeline as in the file-based workflow.

The key bridge is:

- `PlantGeom` builds or edits an MTG with geometry
- `prepare_scene` turns that MTG into the dense scene representation used by `ArchimedLight.jl`

## When To Use This Workflow

- you generate plants procedurally
- you grow or prune plants in a simulation loop
- you want to test light behavior on synthetic scenes without writing files first
- you want to couple light with a geometry-building API such as the `PlantGeom` growth API

## Minimal Example With The Growth API

```julia
using ArchimedLight
using GeometryBasics
using MultiScaleTreeGraph
using OrderedCollections: OrderedDict
using PlantGeom

tri = GeometryBasics.TriangleFace{Int}

stem_mesh = GeometryBasics.mesh(
    GeometryBasics.Cylinder(
        Point(0.0, 0.0, 0.0),
        Point(1.0, 0.0, 0.0),
        0.08,
    ),
)

leaf_mesh = GeometryBasics.Mesh(
    [
        Point(0.0, -0.08, 0.0),
        Point(0.0,  0.08, 0.0),
        Point(0.35, 0.0,  0.0),
    ],
    [tri(1, 2, 3)],
)

prototypes = Dict(
    :Internode => RefMesh("stem", stem_mesh),
    :Leaf => RefMesh("leaf", leaf_mesh),
)

mtg = Node(NodeMTG(:/, :Plant, 1, 1))

axis = emit_internode!(
    mtg;
    link=:/,
    index=1,
    length=0.25,
    width=0.02,
    y_euler=90.0,
    attributes=(functional_group="coffee", object_id=1),
    bump_scene=false,
)

emit_leaf!(
    axis;
    index=1,
    length=0.14,
    width=0.05,
    y_insertion_angle=55.0,
    attributes=(functional_group="coffee", object_id=1),
    bump_scene=false,
)

emit_leaf!(
    axis;
    index=2,
    length=0.12,
    width=0.045,
    y_insertion_angle=-50.0,
    phyllotaxy=160.0,
    attributes=(functional_group="coffee", object_id=1),
    bump_scene=false,
)

rebuild_geometry!(mtg, prototypes; bump_scene=false)

scene = prepare_scene(
    mtg;
    source_path="interactive_growth.opf",
    scene_xy_bounds=(-0.5, -0.5, 0.5, 0.5),
)
```

A few details matter here:

- `functional_group="coffee"` is what links the geometry to the ARCHIMED models
- the node symbols `:Internode` and `:Leaf` become the component type names unless you override them
- `prepare_scene` computes the merged mesh, node areas, barycentres, and node-to-face mapping needed by the light solver

## Add Models In Memory

You do not need model files if you already know the optical parameters:

```julia
models = prepare_models([
    GroupModel(
        "coffee";
        types=OrderedDict(
            "Internode" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.90),
                ),
            ),
            "Leaf" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.90),
                ),
            ),
        ),
    ),
])
```

Wildcard models are also supported:

```julia
models = prepare_models([
    GroupModel(
        "*";
        types=OrderedDict(
            "*" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.15, 0.30),
                ),
            ),
        ),
    ),
])
```

That pattern is convenient for synthetic scenes or tests.

## Run Light Without A Meteo File

If you already know the sky state you want to test, you can bypass file I/O entirely:

```julia
options = LightOptions(
    turtle_sectors=16,
    pixel_size=0.01,
    toricity=false,
    scattering=true,
)

sky = SkyState(
    180.0,  # sun azimuth in degrees
    55.0,   # sun elevation in degrees
    350.0,  # PAR irradiance on horizontal ground, W m^-2
    250.0,  # NIR irradiance on horizontal ground, W m^-2
    0.65,   # direct fraction
    0.35,   # diffuse fraction
)

turtle = build_turtle(options, sky)
fluxes = compute_directional_fluxes(sky, turtle, options)
first = compute_first_order(scene, models, turtle, fluxes, options)
scat = compute_scattering(scene, models, turtle, first, options)
budget = integrate_light(scene, models, first, scat, options; step_duration_seconds=1800.0)
```

This is the most direct way to couple a geometry generator and the light model.

## Modify Geometry And Re-run

Once the geometry stays in memory, you can grow or edit it and reuse the same pattern:

```julia
grow_length!(axis; delta=0.05, bump_scene=false)
set_growth_attributes!(axis; Width=0.025, Thickness=0.025, bump_scene=false)
rebuild_geometry!(mtg, prototypes; bump_scene=false)

scene = prepare_scene(
    mtg;
    source_path="interactive_growth.opf",
    scene_xy_bounds=(-0.5, -0.5, 0.5, 0.5),
)

step = run_light_step(
    scene,
    models,
    (
        date=Dates.Date(2020, 6, 21),
        hour_start="12:00:00",
        hour_end="12:30:00",
        latitude=15.0,
        relativeHumidity=60.0,
        RI_PAR_f=350.0,
        RI_NIR_f=250.0,
        use="relativeHumidity RI_PAR_f RI_NIR_f",
    ),
    options,
)
```

## Practical Advice

- pass `scene_xy_bounds=` explicitly when your scene was not read from an `.ops`
- add ground explicitly with `add_ground!` when you want inspectable paving or better scattering realism near the soil
- prefer wildcard models for quick synthetic tests and explicit group/type models for production scenes
- use `attach_light_step!` once you want to visualize the results through `PlantGeom`
