function load_fixture_inputs(path::AbstractString)
    config_path = isdir(path) ? joinpath(path, "config.yml") : String(path)
    raw = ArchimedLight._load_yaml_ordered(config_path)
    base = dirname(config_path)
    scene_path = normpath(joinpath(base, string(raw["scene"])))
    meteo_path = normpath(joinpath(base, string(raw["meteo"])))
    model_paths =
        if haskey(raw, "models") && raw["models"] isa AbstractVector
            [normpath(joinpath(base, string(p))) for p in raw["models"]]
        else
            String[]
        end

    options, scene, meteo, models = ArchimedLight.read_config(config_path)

    return (
        config_path=config_path,
        scene_path=scene_path,
        meteo_path=meteo_path,
        model_paths=model_paths,
        scene=scene,
        models=models,
        options=options,
        meteo=meteo,
    )
end
