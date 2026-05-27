using ArchimedLight

sky = SkyState(
    135.0,  # sun azimuth, degrees
    60.0,   # sun elevation, degrees
    350.0,  # PAR irradiance, W m^-2
    250.0,  # NIR irradiance, W m^-2
    0.60,   # direct fraction
    0.40,   # diffuse fraction
)

sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
