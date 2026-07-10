using Pkg
Pkg.activate(joinpath(dirname(@__DIR__), "benchmark"))

using PlantGeom # For the growth and visualization API
using GLMakie
using ArchimedLight
using PlantMeteo
using Dates
using Artifacts: artifact_hash, artifact_path
using Pkg.Artifacts: ensure_artifact_installed

const REPO_ROOT = dirname(@__DIR__)
const ARTIFACT_NAME = "archimedlight-benchmark-scenes"
const ARTIFACTS_TOML = joinpath(REPO_ROOT, "Artifacts.toml")
const BENCH_STEP_SECONDS = 1800.0
const BENCH_STEPS = 24

ensure_artifact_installed(ARTIFACT_NAME, ARTIFACTS_TOML)
path_scenes = joinpath(artifact_path(artifact_hash(ARTIFACT_NAME, ARTIFACTS_TOML)), "benchmark_scenes")
scenes = readdir(path_scenes, join=true)

function meteo_generator(n::Int)
    n > 0 || error("Step count must be positive, got $n.")
    # Generate a synthetic meteo series with PAR that varies sinusoidally over the day, peaking at noon.
    return [
        (
            date=DateTime(2023, 6, 1, 0, 0) + Dates.Hour(i - 1),
            latitude=48.8566,
            Ri_PAR_f=350.0 * max(0, sin(π * (i - 1) / (n - 1))),
        )
        for i in 1:n
    ]
end

function models_()
    vegetation = ArchimedLight.translucent(par=0.15, nir=0.90)
    panel = ArchimedLight.translucent(par=0.0, nir=0.0)
    return ArchimedLight.models_for(
        "panel" => (
            "*" => panel,
            "Panel" => panel,
        ),
        "*" => (
            "*" => vegetation,
            "Panel" => panel,
        ),
    )
end

function generate_video(scene, output_video_path::AbstractString; toricity=true, nx=1, ny=1)
    scene = read_scene(scene)
    add_ground!(scene; nx=20, ny=20, group="pavement", type="Cobblestone")
    meteo_gen = meteo_generator(BENCH_STEPS)

    options = LightOptions(
        turtle_sectors=46,
        pixel_size=0.01,
        toricity=toricity,
        scattering=true,
        cache_radiation=true,
        all_in_turtle=true,
        nir_interception=false,
    )


    models = models_()
    sim = LightSimulation(scene, models; options=options)
    series = run_light(sim, meteo_gen)

    # Make a video of the light simulation:
    fig = Figure(size=(900, 700))
    ax = Axis3(
        fig[1, 1],
        aspect=:data,
        title="Absorbed PAR (W m⁻²) $(meteo_gen[1].date)",
        xlabel="x (m)",
        ylabel="y (m)",
        zlabel="z (m)",
    )

    tiled_scene = tile_light_geometry(scene, series; nx=nx, ny=ny, centered=true, xperiod=nothing, yperiod=nothing)

    p = ArchimedLight.lightplot!(ax, tiled_scene, series; color=:incident_par_flux, colormap=:thermal, timestep=1)
    Colorbar(fig[1, 2], p, label="aPAR (W m⁻²)")

    record(fig, output_video_path, 2:length(series), framerate=1) do frame
        Makie.update!(p, timestep=frame)
        ax.title = "Absorbed PAR (W m⁻²) $(meteo_gen[frame].date)"
    end

    Makie.update!(p, timestep=12)
    ax.title = "Absorbed PAR (W m⁻²) $(meteo_gen[12].date)"
    save(joinpath(dirname(output_video_path), string(first(split(basename(output_video_path), ".")), ".png")), fig)
end

# Sole wheat plant:
ASSETS_DIR = joinpath(REPO_ROOT, "docs", "src", "assets")
@time generate_video(scenes[7], joinpath(ASSETS_DIR, "archimedlight_day_cycle_1.mp4"); toricity=false)
@time generate_video(scenes[1], joinpath(ASSETS_DIR, "archimedlight_day_cycle_2.mp4"); toricity=true, nx=10, ny=3)
