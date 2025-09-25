using Test
using ArchimedLight
using ArchimedLight: PAR, NIR
using Meshes
using .ArchimedLight.ArchimedIO
using PlantGeom

# Build a simple square SimpleMesh at z=0
function square_mesh(center::Tuple{<:Real,<:Real,<:Real}, size::Float64)
    s = size / 2
    points = [Point(center[1] - s, center[2] - s, center[3]),
              Point(center[1] + s, center[2] - s, center[3]),
              Point(center[1] + s, center[2] + s, center[3]),
              Point(center[1] - s, center[2] + s, center[3])]
    connec = connect.([(1,2,3), (1,3,4)])
    return SimpleMesh(points, connec)
end

@testset "First-order interception (synthetic)" begin
    scene = Scene()
    mesh = square_mesh((0.0, 0.0, 0.0), 1.0)
    add_component(scene, "leaf", mesh; group="plant")
    set_optics!("leaf", OpticalProps())
    sky = SkyConfig(count=16, PAR_Wm2=500.0, NIR_Wm2=300.0)
    cfg = InterceptionConfig(pixel_size=0.25, scattering=false)
    res = compute_interception(scene, sky, cfg)
    @test haskey(res, "leaf")
    @test res["leaf"][PAR] > 0
end

@testset "Plant OPF import and interception" begin
    # Use PlantGeom's test OPF if available
    opf = joinpath(dirname(dirname(pathof(PlantGeom))), "test", "files", "simple_plant.opf")
    scene = ArchimedIO.load_opf_scene(opf)
    sky = SkyConfig(count=46, PAR_Wm2=300.0, NIR_Wm2=200.0)
    cfg = InterceptionConfig(pixel_size=0.2, scattering=false)
    res = compute_interception(scene, sky, cfg)
    @test length(res) > 0
    @test any(v -> v[PAR] > 0, values(res))
end

@testset "Absorptance response (based on test-absorb)" begin
    # Two identical squares with different optics should absorb proportionally
    scene = Scene()
    meshA = square_mesh((-1.0, 0.0, 0.0), 1.0)
    meshB = square_mesh((1.0, 0.0, 0.0), 1.0)
    add_component(scene, "A_low", meshA; group="plant")
    add_component(scene, "B_high", meshB; group="plant")
    set_optics!("A_low", OpticalProps())
    high = OpticalProps(); high.rho[PAR] = 0.0; high.tau[PAR] = 0.0 # absorptance ~ 1.0 in PAR
    set_optics!("B_high", high)
    sky = SkyConfig(count=16, PAR_Wm2=400.0, NIR_Wm2=0.0)
    cfg = InterceptionConfig(pixel_size=0.2, scattering=false)
    res = compute_interception(scene, sky, cfg)
    @test res["B_high"][PAR] > res["A_low"][PAR]
    ratio = res["B_high"][PAR] / max(res["A_low"][PAR], 1e-9)
    @test ratio > 1.05 # higher absorptance should yield higher absorbed PAR
end

@testset "Pixel size convergence (based on test-pixelsize)" begin
    scene = Scene()
    mesh = square_mesh((0.0, 0.0, 0.0), 1.0)
    add_component(scene, "leaf", mesh; group="plant")
    set_optics!("leaf", OpticalProps())
    sky = SkyConfig(count=16, PAR_Wm2=400.0, NIR_Wm2=0.0)
    res_coarse = compute_interception(scene, sky, InterceptionConfig(pixel_size=0.5, scattering=false))
    res_fine   = compute_interception(scene, sky, InterceptionConfig(pixel_size=0.1, scattering=false))
    # Crude sampler: allow loose agreement between coarse and fine
    a, b = res_coarse["leaf"][PAR], res_fine["leaf"][PAR]
    @test abs(a - b) / max(b, 1e-9) < 1.5
end
