using ArchimedLight

# Run with:
#   julia --project=. examples/voxel_scattering.jl
#
# These optical values are pedagogical and are not a species calibration.

grid = VoxelGrid(
    (0.0, 0.0, 0.0),
    (2.0, 2.0, 3.0),
    fill(0.4, 2, 2, 3),
)

backend = VoxelCPUBackend(
    rays_per_voxel=4,
    boundary=:periodic,
    traversal=:dda,
)

row = (
    step_duration=3600.0,
    Ri_PAR_f=100.0,
    Ri_NIR_f=50.0,
    sun_azimuth=180.0,
    sun_elevation=35.0,
    direct_fraction=0.5,
)

base_options = LightOptions(
    turtle_sectors=6,
    all_in_turtle=true,
    scattering=false,
    radiation_timestep_minutes=5.0,
)

optics = VoxelOpticalProperties(
    grid;
    par_reflectance=0.08,
    par_transmittance=0.07,
    nir_reflectance=0.45,
    nir_transmittance=0.40,
)
ground = VoxelGroundOptics(grid; par_reflectance=0.10, nir_reflectance=0.30)

first_only = run_voxel_light_step(grid, row, base_options; backend=backend)
one_order = run_voxel_light_step(
    grid,
    row,
    LightOptions(
        base_options;
        scattering=true,
        scattering_max_iter=1,
        scattering_stop_ratio=0.0,
    );
    backend=backend,
    optics=optics,
    ground=ground,
)
converged = run_voxel_light_step(
    grid,
    row,
    LightOptions(
        base_options;
        scattering=true,
        scattering_max_iter=30,
        scattering_stop_ratio=1e-8,
    );
    backend=backend,
    optics=optics,
    ground=ground,
)

println("First-order intercepted energy (J)")
println("  PAR: ", sum(first_only.first_order.par_energy))
println("  NIR: ", sum(first_only.first_order.nir_energy))
println()
println("Absorbed energy after one order and after convergence (J)")
println("  PAR: ", sum(one_order.scattering.par.absorbed_energy),
        " -> ", sum(converged.scattering.par.absorbed_energy))
println("  NIR: ", sum(one_order.scattering.nir.absorbed_energy),
        " -> ", sum(converged.scattering.nir.absorbed_energy))
println()

function print_budget(name, band, external)
    accounted = sum(band.absorbed_energy) + band.ground_absorbed_energy +
                band.escaped_top_energy + band.escaped_side_energy +
                band.escaped_bottom_energy + band.unresolved_energy
    println(name)
    println("  foliage absorbed: ", sum(band.absorbed_energy))
    println("  ground absorbed:  ", band.ground_absorbed_energy)
    println("  top escape:       ", band.escaped_top_energy)
    println("  side escape:      ", band.escaped_side_energy)
    println("  bottom escape:    ", band.escaped_bottom_energy)
    println("  unresolved:       ", band.unresolved_energy)
    println("  external:         ", external)
    println("  accounted:        ", accounted)
    println("  iterations:       ", band.iterations)
    println("  converged:        ", band.converged)
    println("  relative residual:", band.relative_balance_residual)
    @assert isapprox(accounted, external; atol=1e-7, rtol=1e-12)
end

print_budget(
    "PAR energy budget",
    converged.scattering.par,
    converged.first_order.par_incoming_energy + converged.first_order.par_injected_energy,
)
println()
print_budget(
    "NIR energy budget",
    converged.scattering.nir,
    converged.first_order.nir_incoming_energy + converged.first_order.nir_injected_energy,
)
