using ArchimedLight

sim = LightSimulation(scene, models; options=LightOptions(cache_radiation=true))

for row in meteo_rows
    light = run_light(sim, row)
    # pass `light` to the host model
end

# When the host model changes geometry:
update_scene!(sim, new_scene)

for row in next_meteo_rows
    light = run_light(sim, row)
    # pass `light` to the host model
end
