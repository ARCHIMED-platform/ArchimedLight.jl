#!/usr/bin/env julia

# Generate the docs landing-page animation from the benchmark-scene artifact.
#
# Default:
#   julia --project=docs scripts/generate_home_video.jl
#
# Select one or several scenes:
#   julia --project=docs scripts/generate_home_video.jl coffee wheat elaeis
#
# Useful overrides:
#   ARCHIMEDLIGHT_HOME_VIDEO_FRAME_MINUTES=60
#   ARCHIMEDLIGHT_HOME_VIDEO_GROUND_TILES=3600
#   ARCHIMEDLIGHT_HOME_VIDEO_PIXELS_PER_SQRT_AREA=180
#   ARCHIMEDLIGHT_HOME_VIDEO_FPS=2
#   ARCHIMEDLIGHT_HOME_VIDEO_OUTPUT=/tmp/day_cycle.mp4

using ArchimedLight
using Artifacts: artifact_hash, artifact_path
using CairoMakie
using Pkg.Artifacts: ensure_artifact_installed
using PlantGeom
using Printf

const REPO_ROOT = dirname(@__DIR__)
const ARTIFACT_NAME = "archimedlight-benchmark-scenes"
const ARTIFACTS_TOML = joinpath(REPO_ROOT, "Artifacts.toml")
const DEFAULT_SCENES = ["coffee", "agrivoltaics_wheat"]
const DEFAULT_OUTPUT = joinpath(
    REPO_ROOT,
    "docs",
    "src",
    "assets",
    "archimedlight_day_cycle.mp4",
)
const AVAILABLE_SCENES = [
    "simple_plant",
    "coffee",
    "wheat",
    "agrivoltaics_wheat",
    "elaeis",
    "elaeis_two_plants",
]
const SCENE_TITLES = Dict(
    "simple_plant" => "Simple plant",
    "coffee" => "Coffee canopy",
    "wheat" => "Wheat",
    "agrivoltaics_wheat" => "Agrivoltaic wheat",
    "elaeis" => "Oil palm",
    "elaeis_two_plants" => "Two oil palms",
)
const BACKGROUND = RGBf(0.982, 0.978, 0.969)
const COLOR_RANGE = (0.0, 500.0)

_env_int(name::AbstractString, default::Integer) = parse(Int, get(ENV, name, string(default)))
_env_float(name::AbstractString, default::Real) = parse(Float64, get(ENV, name, string(default)))

function _scene_names(args)
    names = if !isempty(args)
        String.(args)
    else
        raw = strip(get(ENV, "ARCHIMEDLIGHT_HOME_VIDEO_SCENES", ""))
        isempty(raw) ? copy(DEFAULT_SCENES) : [strip(name) for name in split(raw, ',')]
    end
    names = [name for name in names if !isempty(name)]
    isempty(names) && error("Select at least one scene.")
    unknown = setdiff(names, AVAILABLE_SCENES)
    isempty(unknown) || error(
        "Unknown scene(s): $(join(unknown, ", ")). Available scenes: " *
        join(AVAILABLE_SCENES, ", "),
    )
    length(unique(names)) == length(names) || error("Each scene may be selected only once.")
    return names
end

function _normalize_scene_root(root::AbstractString)
    if isdir(root) && any(endswith(name, ".ops") for name in readdir(root))
        return normpath(root)
    end
    nested = joinpath(root, "benchmark_scenes")
    if isdir(nested) && any(endswith(name, ".ops") for name in readdir(nested))
        return normpath(nested)
    end
    error("Benchmark scene root does not contain OPS files: $root")
end

function benchmark_scene_root()
    override = strip(get(ENV, "ARCHIMEDLIGHT_BENCHMARK_SCENES_ROOT", ""))
    isempty(override) || return _normalize_scene_root(override)

    hash = artifact_hash(ARTIFACT_NAME, ARTIFACTS_TOML)
    hash === nothing && error("Artifact $(repr(ARTIFACT_NAME)) is not registered in $ARTIFACTS_TOML.")
    isdir(artifact_path(hash)) || ensure_artifact_installed(ARTIFACT_NAME, ARTIFACTS_TOML)
    return _normalize_scene_root(artifact_path(hash))
end

function _models()
    vegetation = translucent(par=0.15, nir=0.90)
    pavement = translucent(par=0.10, nir=0.40)
    panel = translucent(par=0.0, nir=0.0)
    return models_for(
        "pavement" => (
            "Cobblestone" => pavement,
            "*" => pavement,
        ),
        "panel" => (
            "Panel" => panel,
            "*" => panel,
        ),
        "*" => (
            "Panel" => panel,
            "*" => vegetation,
        ),
    )
end

function _ground_grid(bounds, target_tiles::Int)
    target_tiles > 0 || error("ARCHIMEDLIGHT_HOME_VIDEO_GROUND_TILES must be positive.")
    xmin, ymin, xmax, ymax = bounds
    width = Float64(xmax - xmin)
    height = Float64(ymax - ymin)
    width > 0.0 && height > 0.0 || error("Scene bounds must have positive width and height.")
    tile_edge = sqrt(width * height / target_tiles)
    return max(1, round(Int, width / tile_edge)), max(1, round(Int, height / tile_edge))
end

function _add_dense_ground!(scene, target_tiles::Int)
    bounds = scene.scene_xy_bounds
    bounds === nothing && error("The scene has no XY plot bounds for ground paving.")
    nx, ny = _ground_grid(bounds, target_tiles)
    PlantGeom.add_ground!(
        scene;
        z=0.005,
        nx=nx,
        ny=ny,
        xy_bounds=bounds,
        group="pavement",
        type="Cobblestone",
    )
    return nx, ny
end

function _pixel_size(bounds, pixels_per_sqrt_area::Float64)
    pixels_per_sqrt_area > 0.0 || error(
        "ARCHIMEDLIGHT_HOME_VIDEO_PIXELS_PER_SQRT_AREA must be positive.",
    )
    xmin, ymin, xmax, ymax = bounds
    return sqrt(Float64((xmax - xmin) * (ymax - ymin))) / pixels_per_sqrt_area
end

function _day_minutes(frame_minutes::Int)
    frame_minutes > 0 || error("ARCHIMEDLIGHT_HOME_VIDEO_FRAME_MINUTES must be positive.")
    start_minute = 7 * 60
    end_minute = 17 * 60
    minutes = collect(start_minute:frame_minutes:end_minute)
    last(minutes) == end_minute || push!(minutes, end_minute)
    return minutes
end

function _clear_sky_day(minutes)
    start_minute, end_minute = first(minutes), last(minutes)
    return map(minutes) do minute
        phase = (minute - start_minute) / (end_minute - start_minute)
        daylight = sinpi(phase)
        SkyState(
            90.0 + 180.0 * phase,
            5.0 + 65.0 * daylight,
            80.0 + 420.0 * daylight,
            0.0,
            0.90,
            0.10,
        )
    end
end

function _simulate_scene(
    scene_root::AbstractString,
    name::AbstractString,
    skies;
    ground_tiles::Int,
    pixels_per_sqrt_area::Float64,
    step_seconds::Float64,
)
    scene_path = joinpath(scene_root, string(name, ".ops"))
    isfile(scene_path) || error("Artifact scene is missing: $scene_path")
    scene = read_scene(scene_path)
    nx, ny = _add_dense_ground!(scene, ground_tiles)
    pixel_size = _pixel_size(scene.scene_xy_bounds, pixels_per_sqrt_area)
    options = LightOptions(
        all_in_turtle=false,
        turtle_sectors=16,
        pixel_size=pixel_size,
        scattering=false,
        cache_radiation=false,
        toricity=true,
        nir_interception=false,
        nir_scattering=false,
    )
    sim = LightSimulation(
        scene,
        _models();
        options=options,
        interception_backend=:raster_cpu,
        scattering_mode=:raycast,
    )
    @info "Simulating landing-page scene" scene=name frames=length(skies) ground="$(nx)x$(ny)" pixel_size
    steps = [run_light(sim, sky; step_duration_seconds=step_seconds) for sky in skies]
    return (name=String(name), scene=scene, steps=steps, geometry=light_render_geometry(steps))
end

function _format_time(minute::Integer)
    @sprintf("%02d:%02d", minute ÷ 60, minute % 60)
end

function _fit_camera!(axis, geometry)
    xs = getindex.(geometry.vertices, 1)
    ys = getindex.(geometry.vertices, 2)
    zs = getindex.(geometry.vertices, 3)
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    zmin, zmax = extrema(zs)
    center = Vec3f((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
    span = max(xmax - xmin, ymax - ymin, 1.5 * (zmax - zmin), eps(Float64))
    eye = Vec3f(center[1] + 0.85 * span, center[2] - 1.05 * span, center[3] + 1.15 * span)
    lookat = Vec3f(center[1], center[2], zmin + 0.25 * (zmax - zmin))
    update_cam!(axis.scene, eye, lookat, Vec3f(0, 0, 1))
    return axis
end

function _build_figure(results, minutes)
    count = length(results)
    ncols = min(count, 3)
    nrows = cld(count, ncols)
    fig = Figure(
        size=(520 * ncols, 500 * nrows + 135),
        backgroundcolor=BACKGROUND,
        fontsize=22,
    )
    time_label = Observable("Representative clear-sky day  ·  $(_format_time(first(minutes)))")
    Label(fig[1, 1:ncols], time_label, fontsize=30, font=:bold, tellwidth=false)

    plots = Any[]
    for (i, result) in enumerate(results)
        row = cld(i, ncols)
        col = mod1(i, ncols)
        title_row = 2 * row
        scene_row = title_row + 1
        Label(fig[title_row, col], SCENE_TITLES[result.name], fontsize=24, font=:bold)
        axis = LScene(
            fig[scene_row, col];
            show_axis=false,
            scenekw=(backgroundcolor=BACKGROUND, clear=true),
        )
        plot = lightplot!(
            axis,
            result.geometry,
            result.steps;
            color=:Ri_PAR_f,
            timestep=1,
            colormap=:thermal,
            colorrange=COLOR_RANGE,
            lowclip=:black,
        )
        _fit_camera!(axis, result.geometry)
        push!(plots, plot)
    end

    colorbar_row = 2 * nrows + 2
    Colorbar(
        fig[colorbar_row, 1:ncols];
        limits=COLOR_RANGE,
        colormap=:thermal,
        vertical=false,
        label="Incident PAR irradiance (W m⁻²)",
        flipaxis=false,
    )
    return fig, plots, time_label
end

function generate_home_video(
    scene_names::AbstractVector{<:AbstractString}=_scene_names(String[]);
    out_path::AbstractString=get(ENV, "ARCHIMEDLIGHT_HOME_VIDEO_OUTPUT", DEFAULT_OUTPUT),
    frame_minutes::Int=_env_int("ARCHIMEDLIGHT_HOME_VIDEO_FRAME_MINUTES", 30),
    ground_tiles::Int=_env_int("ARCHIMEDLIGHT_HOME_VIDEO_GROUND_TILES", 2500),
    pixels_per_sqrt_area::Float64=_env_float(
        "ARCHIMEDLIGHT_HOME_VIDEO_PIXELS_PER_SQRT_AREA",
        150.0,
    ),
    fps::Int=_env_int("ARCHIMEDLIGHT_HOME_VIDEO_FPS", 2),
)
    fps > 0 || error("ARCHIMEDLIGHT_HOME_VIDEO_FPS must be positive.")
    minutes = _day_minutes(frame_minutes)
    skies = _clear_sky_day(minutes)
    scene_root = benchmark_scene_root()
    results = [
        _simulate_scene(
            scene_root,
            name,
            skies;
            ground_tiles=ground_tiles,
            pixels_per_sqrt_area=pixels_per_sqrt_area,
            step_seconds=60.0 * frame_minutes,
        ) for name in scene_names
    ]

    CairoMakie.activate!()
    fig, plots, time_label = _build_figure(results, minutes)
    mkpath(dirname(out_path))
    record(fig, out_path, eachindex(minutes); framerate=fps) do frame
        time_label[] = "Representative clear-sky day  ·  $(_format_time(minutes[frame]))"
        for plot in plots
            plot[:timestep][] = frame
        end
    end
    @info "Wrote landing-page video" path=out_path scenes=join(scene_names, ",") frames=length(minutes)
    return out_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_home_video(_scene_names(ARGS))
end
