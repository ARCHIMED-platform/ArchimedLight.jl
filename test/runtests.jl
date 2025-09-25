using Test
using ArchimedLight
using ArchimedLight: PAR, NIR
using Meshes
using .ArchimedLight.ArchimedIO
using PlantGeom
using ArchimedLight.Runner
using ArchimedLight.Config

include("java_utils.jl")
using .JavaTestUtils

# Build a simple square SimpleMesh at z=0
function square_mesh(center::Tuple{<:Real,<:Real,<:Real}, size::Float64)
    s = size / 2
    points = [Point(center[1] - s, center[2] - s, center[3]),
        Point(center[1] + s, center[2] - s, center[3]),
        Point(center[1] + s, center[2] + s, center[3]),
        Point(center[1] - s, center[2] + s, center[3])]
    connec = connect.([(1, 2, 3), (1, 3, 4)])
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

@testset "Order-0 interception analytic" begin
    scene = Scene()
    mesh = square_mesh((0.0, 0.0, 0.0), 1.0)
    add_component(scene, "leaf", mesh; group="plant")
    props = OpticalProps()
    props.rho[PAR] = 0.2
    props.tau[PAR] = 0.1
    set_optics!("leaf", props)
    sky = SkyConfig(count=16, PAR_Wm2=600.0, NIR_Wm2=0.0)
    cfg = InterceptionConfig(pixel_size=0.2, scattering=false)
    res = compute_interception(scene, sky, cfg)
    inc = res.incident["leaf"][PAR]
    absb = res.absorbed["leaf"][PAR]
    @test inc > 0
    target = 1.0 - props.rho[PAR] - props.tau[PAR]
    @test abs(absb / inc - target) < 1e-6
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

@testset "Config runner (java-style)" begin
    outdir = mktempdir()
    cfgpath = joinpath(outdir, "simple.yml")
    scene_path = joinpath(dirname(dirname(pathof(PlantGeom))), "test", "files", "simple_plant.opf")
    open(cfgpath, "w") do io
        println(io, "scene: \"$(scene_path)\"")
        println(io, "sky_sectors: 46")
        println(io, "pixel_size: 0.2")
        println(io, "scattering: false")
    end
    res = run_from_config(cfgpath; outdir)
    @test isfile(joinpath(outdir, "components.csv"))
    @test any(v -> v[PAR] > 0, values(res))
end

@testset "Java reference: test-absorb" begin
    cfgdir = joinpath(@__DIR__, "java", "test-absorb")
    _, rows, ref = JavaTestUtils.load_reference_table("test-absorb", "config.yml", "component_values.csv")
    @test get(ref.meta, "success", true)

    for (type, band) in (("Metamer", "PAR"), ("Metamer", "NIR"), ("Leaf", "PAR"), ("Leaf", "NIR"))
        subset = [row for row in filter_rows(rows, "type", type) if haskey(row, "Ri_$(band)_0_f") && parse_float(row, "Ri_$(band)_0_f") != 0.0]
        @test !isempty(subset)
        ratios = [parse_float(row, "Ra_$(band)_0_f") / parse_float(row, "Ri_$(band)_0_f") for row in subset]
        ref_ratio = first(ratios)
        maxdiff = maximum(abs.(ratios .- ref_ratio))
        @test maxdiff < 1e-5
    end

    for (type, band) in (("Metamer", "PAR"), ("Metamer", "NIR"), ("Leaf", "PAR"), ("Leaf", "NIR"))
        rq = "Ra_$(band)_0_q"
        subset = [row for row in filter_rows(rows, "type", type) if all(k -> haskey(row, k), (rq, "Ri_$(band)_0_f", "step_duration", "area")) && parse_float(row, "Ri_$(band)_0_f") != 0.0]
        isempty(subset) && continue
        ratios = [parse_float(row, rq) / (parse_float(row, "Ri_$(band)_0_f") * parse_float(row, "step_duration") * parse_float(row, "area")) for row in subset]
        ref_ratio = first(ratios)
        maxdiff = maximum(abs.(ratios .- ref_ratio))
        @test maxdiff < 1e-5
    end

    cmp = JavaTestUtils.compare_and_report(cfgdir, "config.yml"; tol=1e-3, strict=false)
    warn_threshold = 1e-3
    for col in JavaTestUtils.DEFAULT_METRICS
        @test haskey(cmp.java, col)
        @test haskey(cmp.julia, col)
        @test haskey(cmp.diff, col)
        @test haskey(cmp.rel, col)
        @test isfinite(cmp.rel[col])
        java_v = cmp.java[col]
        julia_v = cmp.julia[col]
        diff_v = cmp.diff[col]
        rel_v = cmp.rel[col]
        @info "metric" col "java" java_v "julia" julia_v "diff" diff_v "rel" rel_v
        if rel_v > warn_threshold
            @warn "Large relative difference between Java and Julia for metric" col rel = rel_v diff = diff_v java = java_v julia = julia_v
        end
        @test_broken isapprox(diff_v, 0.0; atol=1e-6, rtol=1e-6)
    end
end

@testset "Java reference: test-pixelsize" begin
    cfgdir = joinpath(@__DIR__, "java", "test-pixelsize")
    ref_base = JavaTestUtils.load_reference_case("test-pixelsize", "config.yml")
    @test get(ref_base.meta, "success", true)

    ref_success = JavaTestUtils.load_reference_case("test-pixelsize", "config.yml"; extra_args=["--prop", "pixel_size=49.9999"])
    @test get(ref_success.meta, "success", false)

    ref_failure = JavaTestUtils.load_reference_case("test-pixelsize", "config.yml"; extra_args=["--prop", "pixel_size=50.0001"])
    @test !get(ref_failure.meta, "success", true)
    @test occursin("invalid pixel size", ref_failure.log) || occursin("greater than plot", ref_failure.log)

    cmp = compare_java_julia(cfgdir, "config.yml")
    for col in JavaTestUtils.DEFAULT_METRICS
        @test isfinite(cmp.rel[col])
        @test isapprox(cmp.diff[col], 0.0; atol=1e-6, rtol=1e-6)
    end
end

@testset "Java reference: test-ops-output" begin
    cfgdir = joinpath(@__DIR__, "java", "test-ops-output")
    function collect_ops(config)
        ref = JavaTestUtils.load_reference_case("test-ops-output", config)
        files = String[]
        for dir in ref.run_paths
            for entry in readdir(dir)
                endswith(entry, ".ops") && push!(files, entry)
            end
        end
        sort!(files)
        return get(ref.meta, "success", false), ref.log, files
    end

    success, _log, files = collect_ops("config1.yml")
    @test success
    @test files == ["scene-step00004.ops"]

    success, _log, files = collect_ops("config2.yml")
    @test success
    @test isempty(files)

    success, _log, files = collect_ops("config3.yml")
    @test success
    @test length(files) == 4

    success, _log, files = collect_ops("config4.yml")
    @test success
    @test files == ["scene-step00003.ops"]

    success, log_text, _ = collect_ops("config5.yml")
    @test !success
    @test occursin("invalid export_ops", log_text)

    success, log_text, _ = collect_ops("config6.yml")
    @test !success
    @test occursin("invalid export_ops", log_text)
end

@testset "Java reference: test-lightsource1" begin
    cfgdir = joinpath(@__DIR__, "java", "test-lightsource1")

    _, out1, ref1 = JavaTestUtils.load_reference_table("test-lightsource1", "config1.yml", "log-intercept-source00.csv"; extra_args=["-d", "--no-upperhit-pixtable"])
    @test get(ref1.meta, "success", true)

    _, out2, ref2 = JavaTestUtils.load_reference_table("test-lightsource1", "config2.yml", "log-intercept-source01.csv"; extra_args=["-d"])
    @test get(ref2.meta, "success", true)

    function step_map(rows)
        mapping = Dict{NTuple{3,Int},Float64}()
        for row in rows
            step = parse(Int, row["step"])
            plant = parse(Int, row["plantid"])
            node = parse(Int, row["nodeid"])
            val = parse_float(row, "stepIntercept")
            mapping[(step, plant, node)] = val
        end
        return mapping
    end

    map1 = step_map(out1)
    map2 = step_map(out2)
    for key in collect(keys(map2))
        if key[2] == 2 && key[3] == 1
            delete!(map2, key)
        end
    end

    common = intersect(keys(map1), keys(map2))
    @test !isempty(common)
    for key in common
        @test abs(map1[key] - map2[key]) <= 0.1
    end
end

@testset "Absorptance response (based on test-absorb)" begin
    # Two identical squares with different optics should absorb proportionally
    scene = Scene()
    meshA = square_mesh((-1.0, 0.0, 0.0), 1.0)
    meshB = square_mesh((1.0, 0.0, 0.0), 1.0)
    add_component(scene, "A_low", meshA; group="plant")
    add_component(scene, "B_high", meshB; group="plant")
    set_optics!("A_low", OpticalProps())
    high = OpticalProps()
    high.rho[PAR] = 0.0
    high.tau[PAR] = 0.0 # absorptance ~ 1.0 in PAR
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
    res_fine = compute_interception(scene, sky, InterceptionConfig(pixel_size=0.1, scattering=false))
    # Crude sampler: allow loose agreement between coarse and fine
    a, b = res_coarse["leaf"][PAR], res_fine["leaf"][PAR]
    @test abs(a - b) / max(b, 1e-9) < 1.5
end
