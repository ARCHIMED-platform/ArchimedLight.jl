# Tutorial: Interactive Workflow

This tutorial shows the in-memory workflow: build or modify a plant directly in Julia, convert it to a `PlantGeom.SceneGeometry`, then run the same light pipeline as in the file-based workflow.

The key bridge is:

- `PlantGeom` builds or edits an MTG with geometry
- `PlantGeom.prepare_scene` turns that MTG into the dense scene representation used by `ArchimedLight.jl`

```@setup interactive_workflow
using CairoMakie
CairoMakie.activate!(type = "png")
```

## When To Use This Workflow

- you generate plants procedurally
- you grow or prune plants in a simulation loop
- you want to test light behavior on synthetic scenes without writing files first
- you want to couple light with a geometry-building API such as the `PlantGeom` growth API

## Creating a Plant With PlantGeom's Growth API

First we import the necessary packages:

```@example interactive_workflow
using ArchimedLight
using CairoMakie # For plotting
using Colors # For color definitions
using GeometryBasics # For geometry
using MultiScaleTreeGraph # For the MTG data structure
using OrderedCollections: OrderedDict
using PlantGeom # For the growth and visualization API
```

Then we build a simple plant with two successive internode–leaf pairs. The internodes are cylinders and the leaves are laminae with a midrib.

First, we define mesh prototypes for the internode and leaf shapes:

```@example interactive_workflow
stem_ref = RefMesh(
    "stem",
    GeometryBasics.mesh(
        GeometryBasics.Cylinder(
            Point(0.0, 0.0, 0.0),
            Point(1.0, 0.0, 0.0),
            0.5,
        ),
    ),
    RGB(0.48, 0.36, 0.25),
)

leaf_ref = lamina_refmesh(
    "leaf";
    length=1.0,
    max_width=1.0,
    n_long=36,
    n_half=7,
    material=RGB(0.19, 0.61, 0.29),
)

prototypes = Dict(
    :Internode => RefMeshPrototype(stem_ref, true),
    :Leaf => PointMapPrototype(
        leaf_ref;
        defaults=(base_angle_deg=42.0, bend=0.30, tip_drop=0.08),
        intrinsic_shape=params -> LaminaMidribMap(
            base_angle_deg=params.base_angle_deg,
            bend=params.bend,
            tip_drop=params.tip_drop,
        ),
    ),
)
```

Then, we add two internode–leaf pairs to the MTG. This low-level PlantGeom helper emits component nodes; it does not create botanical phytomere nodes:

```@example interactive_workflow
scene_mtg = Node(NodeMTG(:/, :Scene, 1, 1))
small_plant = Node(scene_mtg, NodeMTG(:/, :Plant, 1, 1))
small_plant[:functional_group] = "example_plant"
small_plant[:object_id] = 1

first_pair = emit_internode_leaf!(
    small_plant;
    internode=(link=:/, index=1, length=0.20, width=0.022),
    leaf=(index=1, offset=0.15, length=0.22, width=0.05, thickness=0.02, y_insertion_angle=52.0),
)

second_pair = emit_internode_leaf!(
    first_pair.internode;
    internode=(index=2, length=0.18, width=0.020),
    leaf=(index=2, offset=0.14, length=0.24, width=0.055, thickness=0.02, phyllotaxy=180.0, y_insertion_angle=54.0),
)

scene_mtg
```

By default, the growth API from `PlantGeom` uses `:Internode` and `:Leaf` symbols for the component types. The component type names are important because they link the geometry to the optical models later on.

This example is taken from PlantGeom's documentation, and annotated with ARCHIMED metadata (`functional_group`, `object_id`) so it can be passed to `ArchimedLight.jl`.

The only addition we had to make was to declare the `functional_group` of the plant as an attribute, because `ArchimedLight.jl` uses the functional group to link the geometry to the models. The `object_id` is not mandatory, it is used to identify the plants in the scene. If not specified, the plants are numbered by their order as scene children.

We can now visualize the MTG with `plantviz`:

```@example interactive_workflow
rebuild_geometry!(scene_mtg, prototypes)
plantviz(scene_mtg, figure=(size=(820, 560),))
```

## Convert To SceneGeometry

The MTG is not yet in the format needed by `ArchimedLight.jl`. We need to convert it to a `PlantGeom.SceneGeometry` with the `PlantGeom.prepare_scene` function:

```@example interactive_workflow
scene = PlantGeom.prepare_scene(scene_mtg;scene_xy_bounds=(-0.8, -0.8, 0.8, 0.8))
```

`PlantGeom.prepare_scene` computes the merged mesh, node areas, barycentres, and node-to-face mapping needed by the light solver. The `scene_xy_bounds` are used for the rasterization stage, so they should be chosen to fit the plant geometry.

## Add Models

ArchimedLight.jl needs optical models to compute light interception. We can define them directly in Julia as `GroupModel` and `TypeModel` objects, then prepare them with `prepare_models`. The example below defines two group models: one for the plant and one for the pavement. Each group model has one type model with `Translucent` interception and specific optical properties.

The pavement model is included here because we will add an explicit ground later on, and we want it to have different optical properties from the plant. Optical properties are defined as `(scattering_coeff_par, scattering_coeff_nir)`, where the absorbed fraction is `1 - scattering_coeff`.

```@example interactive_workflow
models = prepare_models([
    GroupModel(
        "example_plant";
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
    GroupModel(
        "pavement";
        types=OrderedDict(
            "Cobblestone" => TypeModel(
                interception=InterceptionModel(
                    model="Translucent",
                    transparency=0.0,
                    optical_properties=OpticalProperties(0.12, 0.60),
                ),
            ),
        ),
    ),
])

keys(models.groups) |> collect
```

That pattern is convenient for synthetic scenes, tests, and gradual model setup.

## Add Ground Explicitly

The ground can be added as a special case of scene geometry with `PlantGeom.add_ground!`. This is useful when:

- you want to include the ground in the directional visibility and scattering calculations for better realism near the soil
- you want to visualize the soil interception with an explicit ground patch

```@example interactive_workflow
add_ground!(
    scene;
    nx=60,
    ny=60,
    xy_bounds=(-0.8, -0.8, 0.8, 0.8),
    group="pavement",
    type="Cobblestone",
)
```

We can now visualize the plant and the ground together:

```@example interactive_workflow
plantviz(scene.mtg, figure=(size=(820, 560),))
```

## Run One Simulation And Visualize It

The example below runs a single light simulation on that plant, and visualizes the result directly in 3D.

First, we define the sky conditions with a `SkyState`:

```@example interactive_workflow
sky = SkyState(
    135.0,  # sun azimuth in degrees
    60.0,   # sun elevation in degrees
    350.0,  # PAR irradiance on horizontal ground, W m^-2
    250.0,  # NIR irradiance on horizontal ground, W m^-2
    0.95,   # direct fraction
    0.05,   # diffuse fraction
)
```

This can be replaced with a row from a `PlantMeteo.TimeStepTable`, and ArchimedLight.jl will compute the `SkyState` from the solar geometry and irradiance values in the meteo file.

Then we define the `LightOptions`, which defines simulations parameters such as the number of turtle sectors, pixel size, toricity, and scattering behavior:

```@example interactive_workflow
options = LightOptions(
    turtle_sectors=16,
    pixel_size=0.01,
    toricity=true,
    scattering=false,
)
```

Finally, we create a `LightSimulation` and run one light step:

```@example interactive_workflow
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

The same `LightSimulation` can also be reused with meteo rows or meteo tables.
Use the staged pipeline only when you need to inspect or customize internal
solver stages.

To visualize the spatial distribution of light, we need to attach the results back to the MTG with `attach_light_step!`:

```@example interactive_workflow
attach_light_step!(scene, step; fields=[:incident_par_flux])
```

Then we can plot the scene colored by `Ri_PAR_f`:

```@example interactive_workflow
par_max = maximum(values(step.budget.incident_flux.total.par))

fig, ax, p = plantviz(
    scene.mtg;
    color=:Ri_PAR_f,
    colormap=:thermal,
    colorrange=(0.0, par_max),
    color_missing=:gray85,
    figure=(size=(900, 700),),
)

ax.show_axis[] = false
PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")
fig
```

## Practical Advice

- pass `scene_xy_bounds=` explicitly when your scene was not read from an `.ops`
- add ground explicitly with `PlantGeom.add_ground!` when you want inspectable paving or better scattering realism near the soil
- use `attach_light_step!` once you want to visualize the results through `PlantGeom`
