if lowercase(get(ENV, "ARCHIMEDLIGHT_DEVICE_BACKEND", "cpu")) == "metal" && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using GeometryBasics
using KernelAbstractions
using LinearAlgebra
using MultiScaleTreeGraph
using PlantGeom
using StaticArrays

function _load_optional_package(name::Symbol)
    Base.find_package(String(name)) === nothing &&
        error("Package $name is not available in the active environment.")
    return Base.require(Main, name)
end

function _backend_from_env()
    name = lowercase(get(ENV, "ARCHIMEDLIGHT_DEVICE_BACKEND", "cpu"))
    name == "cpu" && return KernelAbstractions.CPU()

    if name == "cuda"
        mod = _load_optional_package(:CUDA)
        functional = getproperty(mod, :functional)
        Base.invokelatest(functional) || error("CUDA is available but CUDA.functional() is false.")
        array_type = getproperty(mod, :CuArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif name == "metal"
        mod = _load_optional_package(:Metal)
        array_type = getproperty(mod, :MtlArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif name == "oneapi"
        mod = _load_optional_package(:oneAPI)
        array_type = getproperty(mod, :oneArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    elseif name == "amdgpu"
        mod = _load_optional_package(:AMDGPU)
        array_type = getproperty(mod, :ROCArray)
        return Base.invokelatest(KernelAbstractions.get_backend, Base.invokelatest(array_type, zeros(Float32, 1)))
    end

    error(
        "Unsupported ARCHIMEDLIGHT_DEVICE_BACKEND=$name. " *
        "Use cpu, cuda, metal, oneapi, or amdgpu.",
    )
end

function _smoke_scene()
    specs = (
        (z=1.0, group="upper", object_id=1, source_topology_id=1),
        (z=0.1, group="lower", object_id=2, source_topology_id=2),
    )
    points = GeometryBasics.Point{3,Float32}[]
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    nodes = Dict{Int,PlantGeom.SceneNodeData{Float64}}()
    mtg = MultiScaleTreeGraph.Node(
        MultiScaleTreeGraph.MutableNodeMTG(:/, :Scene, 0, 0),
        Dict{Symbol,Any}(),
    )

    for (i, spec) in pairs(specs)
        p1 = (0.0, 0.0, spec.z)
        p2 = (1.0, 0.0, spec.z)
        p3 = (1.0, 1.0, spec.z)
        p4 = (0.0, 1.0, spec.z)
        base = length(points)
        append!(
            points,
            GeometryBasics.Point{3,Float32}[
                GeometryBasics.Point{3,Float32}(Float32(p1[1]), Float32(p1[2]), Float32(p1[3])),
                GeometryBasics.Point{3,Float32}(Float32(p2[1]), Float32(p2[2]), Float32(p2[3])),
                GeometryBasics.Point{3,Float32}(Float32(p3[1]), Float32(p3[2]), Float32(p3[3])),
                GeometryBasics.Point{3,Float32}(Float32(p4[1]), Float32(p4[2]), Float32(p4[3])),
            ],
        )
        append!(faces, GeometryBasics.TriangleFace{Int}[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, (i, i))

        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        barycenter = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            spec.z,
        )
        nodes[i] = PlantGeom.SceneNodeData(area1 + area2, barycenter, spec.source_topology_id)
        mesh = _plate_mesh(Float32(spec.z))
        MultiScaleTreeGraph.Node(
            i,
            mtg,
            MultiScaleTreeGraph.MutableNodeMTG(:+, :plate, i, 1),
            Dict{Symbol,Any}(
                :geometry => PlantGeom.Geometry(ref_mesh=PlantGeom.RefMesh("raycore_device_smoke_$i", mesh)),
                :group => spec.group,
                :functional_group => spec.group,
                :type => "plate",
                :object_id => spec.object_id,
                :source_topology_id => spec.source_topology_id,
            ),
        )
    end

    return PlantGeom.SceneGeometry(
        mtg,
        GeometryBasics.Mesh(points, faces),
        face2node,
        nodes,
        "raycore_device_smoke",
        (0.0, 0.0, 1.0, 1.0),
    )
end

function _plate_mesh(z::Float32)
    return GeometryBasics.Mesh(
        GeometryBasics.Point3f[
            GeometryBasics.Point3f(0, 0, z),
            GeometryBasics.Point3f(1, 0, z),
            GeometryBasics.Point3f(1, 1, z),
            GeometryBasics.Point3f(0, 1, z),
        ],
        GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)],
    )
end

function _smoke_models()
    return ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=ArchimedLight.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
end

function _dict_close(a::Dict{Int,Float64}, b::Dict{Int,Float64}; atol::Float64, rtol::Float64)
    keys_union = union(keys(a), keys(b))
    return all(keys_union) do key
        isapprox(get(a, key, 0.0), get(b, key, 0.0); atol=atol, rtol=rtol)
    end
end

function _assert_close(label::AbstractString, actual::Real, expected::Real; atol::Float64, rtol::Float64)
    isapprox(actual, expected; atol=atol, rtol=rtol) ||
        error("$label $actual does not match CPU reference $expected.")
    return nothing
end

function _assert_graph_close(raycore_graph, cpu_graph)
    raycore_graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes ||
        error("Raycore scattering graph to_nodes do not match CPU raycast.")
    raycore_graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes ||
        error("Raycore scattering graph from_nodes do not match CPU raycast.")
    raycore_graph.pair_counts.counts == cpu_graph.pair_counts.counts ||
        error("Raycore scattering graph counts do not match CPU raycast.")
    raycore_graph.all_hits == cpu_graph.all_hits ||
        error("Raycore scattering all_hits do not match CPU raycast.")
    return nothing
end

function run_raycore_device_smoke(;
    backend=_backend_from_env(),
    first_order_area_atol::Float64=1e-10,
    first_order_area_rtol::Float64=1e-10,
    first_order_power_atol::Float64=1e-8,
    first_order_power_rtol::Float64=1e-8,
    scattering_atol::Float64=1e-5,
    scattering_rtol::Float64=1e-5,
)
    scene = _smoke_scene()
    models = _smoke_models()
    options = ArchimedLight.LightOptions(
        turtle_sectors=1,
        all_in_turtle=false,
        scattering=false,
        pixel_size=0.01,
        toricity=false,
    )
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    raster = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    ib = ArchimedLight.RaycoreInterceptionBackend(backend=backend)
    raycore_data = ArchimedLight._prepare_raycore_interception_data(scene, models, options, ib)
    reduction_capabilities = ArchimedLight._raycore_device_reduction_capabilities(raycore_data)
    raycore = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        fluxes,
        options;
        backend=ib,
    )

    raster_area = sum(values(raster.projected_area_per_node); init=0.0)
    raycore_area = sum(values(raycore.projected_area_per_node); init=0.0)
    raster_par = sum(values(raster.incident_power.par); init=0.0)
    raycore_par = sum(values(raycore.incident_power.par); init=0.0)

    raycore_area > 0.0 || error("Raycore smoke scene did not intercept any projected area.")
    _assert_close(
        "Raycore projected area",
        raycore_area,
        raster_area;
        atol=first_order_area_atol,
        rtol=first_order_area_rtol,
    )
    _assert_close(
        "Raycore PAR",
        raycore_par,
        raster_par;
        atol=first_order_power_atol,
        rtol=first_order_power_rtol,
    )

    scatter_options = ArchimedLight.LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        scattering=true,
        pixel_size=0.01,
        toricity=false,
        scattering_max_iter=10,
    )
    scatter_sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    scatter_turtle = ArchimedLight.build_turtle(scatter_options, scatter_sky)
    scatter_fluxes = ArchimedLight.compute_directional_fluxes(scatter_sky, scatter_turtle, scatter_options)
    scatter_first = ArchimedLight.compute_first_order(scene, models, scatter_turtle, scatter_fluxes, scatter_options; backend=ib)
    cpu_graph = ArchimedLight.build_scattering_transfer_graph(
        scene,
        models,
        scatter_turtle,
        scatter_first,
        scatter_options;
        backend=ArchimedLight.RaycastScatteringBackend(),
    )
    raycore_scattering_backend = ArchimedLight.RaycoreScatteringBackend(ib)
    raycore_graph = ArchimedLight.build_scattering_transfer_graph(
        scene,
        models,
        scatter_turtle,
        scatter_first,
        scatter_options;
        backend=raycore_scattering_backend,
    )
    _assert_graph_close(raycore_graph, cpu_graph)
    cpu_scattering = ArchimedLight.compute_scattering(cpu_graph, scatter_first, scatter_options; backend=ArchimedLight.RaycastScatteringBackend())
    raycore_scattering = ArchimedLight.compute_scattering(raycore_graph, scatter_first, scatter_options; backend=raycore_scattering_backend)

    _dict_close(raycore_scattering.added_power.par, cpu_scattering.added_power.par; atol=scattering_atol, rtol=scattering_rtol) ||
        error("Raycore scattering PAR does not match CPU raycast scattering: raycore=$(raycore_scattering.added_power.par), cpu=$(cpu_scattering.added_power.par).")
    _dict_close(raycore_scattering.added_power.nir, cpu_scattering.added_power.nir; atol=scattering_atol, rtol=scattering_rtol) ||
        error("Raycore scattering NIR does not match CPU raycast scattering: raycore=$(raycore_scattering.added_power.nir), cpu=$(cpu_scattering.added_power.nir).")

    @info "Raycore device smoke passed" backend_type=typeof(backend) raycore_area raycore_par scattering_eltype=ib.config.scattering_eltype reduction_capabilities
    return (
        backend=backend,
        first_order=(area=raycore_area, par=raycore_par),
        scattering=raycore_scattering,
        scattering_eltype=ib.config.scattering_eltype,
        reduction_capabilities=reduction_capabilities,
    )
end

function main()
    run_raycore_device_smoke()
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
