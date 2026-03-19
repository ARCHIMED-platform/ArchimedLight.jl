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

    scene = ArchimedLight.read_scene(scene_path)
    models = ArchimedLight.read_models(model_paths)
    options = ArchimedLight.read_options(config_path)
    meteo = ArchimedLight.read_meteo(meteo_path)

    paving_count = 0
    for group_model in values(models.groups)
        for type_model in values(group_model.types)
            if haskey(type_model.extras, "plot_paving")
                paving_count = max(paving_count, Int(type_model.extras["plot_paving"]))
            end
        end
    end
    if paving_count > 0
        n = max(1, round(Int, sqrt(paving_count)))
        # Legacy fixtures used the Java paving mesh, which sits 0.5 cm above z=0.
        ArchimedLight.add_ground!(scene; nx=n, ny=n, z=0.005)
    end

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
