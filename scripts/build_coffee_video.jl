const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "benchmark"))

using ArchimedLight
using KernelAbstractions
using Metal
using Dates
using GLMakie

metal_backend = KernelAbstractions.get_backend(MtlArray(zeros(Float32, 1)))

options = LightOptions(
    pixel_size=0.003,          # meters = 10 cm, keeps elaeis fast
    radiation_timestep_minutes=5,
    turtle_sectors=16,
    scattering=true,
    cache_radiation=true,
    all_in_turtle=true,
    toricity=true,
    allow_overlapping_meteo_steps=true
)

tile_face_capacity = options.turtle_sectors < 136 ? 128 : 256
max_hits_per_pixel = options.turtle_sectors < 136 ? 128 : 512
tile_face_capacity = options.turtle_sectors < 136 ? 128 : 512
interception = ArchimedLight.RasterGPUBackend(backend=metal_backend, tile_size=1, tile_face_capacity=tile_face_capacity, max_hits_per_pixel=max_hits_per_pixel, edge_accumulation=:auto)
scattering = ArchimedLight.RasterGPUScatteringBackend(interception)

config_path = joinpath("example_2", "config.yml")
sim, _ = read_simulation(
    config_path;
    interception_backend=interception,
    scattering_backend=scattering,
)
sim.options = options

date_start = DateTime(2020, 6, 21, 8, 0, 0)
increment = Minute(30)
rows = [
    (
        date=date_start+increment*(i-1),
        duration=increment,
        step_duration=increment,
        latitude=15.0,
        Ri_PAR_f=320.0 + 80.0 * sin(pi * i / 24),
        Ri_NIR_f=180.0 + 50.0 * sin(pi * i / 24),
    )
    for i in 1:24
]

step = run_light(sim, rows)
out = ArchimedLight.tile_light_geometry(sim.scene, step; nx=3, ny=3)
fig = Figure(size=(900, 700))
ax = Axis3(
    fig[1, 1],
    aspect=:data,
    title="Absorbed PAR (W m⁻²) $(rows[1].date)",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    elevation=π/4,
    azimuth=(-π/2)
)
p = ArchimedLight.lightplot!(ax, out, step; color=:incident_par_flux, colormap=:thermal, timestep=1, colorrange=(0.0, 450.0))
Colorbar(fig[1, 2], p, label="aPAR (W m⁻²)")
record(fig, "light_rasterizer_new.mp4", 2:length(step), framerate=1) do frame
    Makie.update!(p, timestep=frame)
    ax.title = "Absorbed PAR (W m⁻²) $(rows[frame].date)"
end