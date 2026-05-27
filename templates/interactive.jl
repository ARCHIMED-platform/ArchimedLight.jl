using ArchimedLight
using DataFrames
using PlantGeom

bounds = (-1.0, -1.0, 1.0, 1.0)

scene = make_scene(domain=bounds) do s
    add_plant!(s, "scene/plant.opf"; group="coffee", id=1, at=(0.0, 0.0, 0.0))
    add_ground!(s; group="soil", type="ground", nx=20, ny=20)
end

models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
    "soil" => (
        "ground" => translucent(par=0.10, nir=0.40),
    ),
)

meteo = DataFrame(
    date=["2020/06/21"],
    hour_start=["12:00:00"],
    hour_end=["13:00:00"],
    latitude=[15.0],
    RI_PAR_f=[350.0],
    RI_NIR_f=[250.0],
)

options = LightOptions(pixel_size=0.01, scattering=true)
sim = LightSimulation(scene, models; options=options)

summarize_scene(scene; models=models)
summarize_meteo(meteo; options=options)

step = run_light(sim, first(eachrow(meteo)))
