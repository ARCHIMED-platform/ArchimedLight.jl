using OrderedCollections: OrderedDict

const _SYNTHETIC_CASE_FILTERS = Set(filter(!isempty, strip.(lowercase.(split(get(ENV, "ARCHIMEDLIGHT_SYNTHETIC_CASE", ""), ",")))))

function _synthetic_case_enabled(name::String)
    isempty(_SYNTHETIC_CASE_FILTERS) && return true
    lowercase(name) in _SYNTHETIC_CASE_FILTERS
end

function _synthetic_cfg(
    cfg::ArchimedLight.LightConfig;
    sectors::Int=1,
    all_in_turtle::Bool=false,
    scattering::Bool=false,
    pixel_size::Float64=0.01,
    toricity::Bool=true,
    cache_radiation::Bool=false,
    models::Vector{String}=String[],
)
    out = deepcopy(cfg)
    out.general["all_in_turtle"] = all_in_turtle
    out.general["sky_sectors"] = sectors
    out.general["scattering"] = scattering
    out.general["pixel_size"] = pixel_size
    out.general["toricity"] = toricity
    out.general["cache_radiation"] = cache_radiation
    base = out.source_files.base_dir
    out.source_files.models = OrderedDict{String,String}()
    if !isempty(models)
        for rel in models
            abs_path = isabspath(rel) ? rel : normpath(joinpath(base, rel))
            group = ArchimedLight._model_group_from_file(abs_path)
            out.source_files.models[group] = abs_path
        end
    end
    out.models = OrderedDict{String,OrderedDict{String,Any}}()
    ArchimedLight.refresh_light_config!(out; reload_models=!isempty(models))
    return out
end

function _synthetic_horizontal_scene(specs::AbstractVector{<:NamedTuple})
    quad_specs = map(specs) do spec
        (
            p1=(spec.x0, spec.y0, spec.z),
            p2=(spec.x1, spec.y0, spec.z),
            p3=(spec.x1, spec.y1, spec.z),
            p4=(spec.x0, spec.y1, spec.z),
            group=get(spec, :group, "plate"),
            type=get(spec, :type, "plate"),
        )
    end
    _synthetic_quad_scene(quad_specs)
end

function _synthetic_quad_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float32}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    total_area_per_node = Dict{Int,Float64}()
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    source_topology_id_per_node = Dict{Int,Int}()
    object_id_per_node = Dict{Int,Int}()
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()

    xs = Float64[]
    ys = Float64[]
    for (i, spec) in enumerate(specs)
        p1 = ntuple(j -> Float64(spec.p1[j]), 3)
        p2 = ntuple(j -> Float64(spec.p2[j]), 3)
        p3 = ntuple(j -> Float64(spec.p3[j]), 3)
        p4 = ntuple(j -> Float64(spec.p4[j]), 3)
        append!(xs, (p1[1], p2[1], p3[1], p4[1]))
        append!(ys, (p1[2], p2[2], p3[2], p4[2]))

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
        append!(faces, PlantGeom.Face3[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, [i, i])

        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        total_area_per_node[i] = area1 + area2
        barycenter_per_node[i] = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            (p1[3] + p2[3] + p3[3] + p4[3]) / 4,
        )
        source_topology_id_per_node[i] = Int(get(spec, :source_topology_id, i))
        object_id_per_node[i] = Int(get(spec, :object_id, source_topology_id_per_node[i]))
        node_group[i] = String(get(spec, :group, "plate"))
        node_type[i] = String(get(spec, :type, "plate"))
    end

    ArchimedLight.SceneGeometry(
        nothing,
        GeometryBasics.Mesh(points, faces),
        face2node,
        total_area_per_node,
        barycenter_per_node,
        source_topology_id_per_node,
        object_id_per_node,
        node_group,
        node_type,
        "synthetic_scene_cases",
        (minimum(xs), minimum(ys), maximum(xs), maximum(ys)),
    )
end

function _synthetic_meteo_row(;
    date::Dates.Date=Dates.Date(2020, 6, 21),
    start_time::Dates.Time=Dates.Time(12),
    duration_seconds::Float64=1.0,
    latitude::Float64=0.0,
    relative_humidity::Float64=60.0,
    ri_par_f::Float64=100.0,
    ri_nir_f::Float64=0.0,
    direct_fraction::Float64=1.0,
    sun_azimut::Float64=180.0,
    sun_elevation::Float64=90.0,
    use::String="relativeHumidity RI_PAR_f",
)
    start_dt = Dates.DateTime(date, start_time)
    end_dt = start_dt + Dates.Millisecond(round(Int, duration_seconds * 1000))
    (
        date=date,
        hour_start=Dates.format(Dates.Time(start_dt), Dates.DateFormat("HH:MM:SS")),
        hour_end=Dates.format(Dates.Time(end_dt), Dates.DateFormat("HH:MM:SS")),
        step_duration=duration_seconds,
        latitude=latitude,
        relativeHumidity=relative_humidity,
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
        use=use,
    )
end

function _max_abs_float_dict_diff(a::Dict{Int,Float64}, b::Dict{Int,Float64})
    maximum(abs(get(a, id, 0.0) - get(b, id, 0.0)) for id in union(keys(a), keys(b)); init=0.0)
end

@testset "Synthetic scene cases" begin
    cfg_ref = ArchimedLight.read_light_config(
        joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml"),
    )

        if _synthetic_case_enabled("single_plate_direct")
            @testset "Scenario: single 1 m² plate under zenith direct PAR" begin
                inputs = (
                    scene=[(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0),
                    cfg=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01),
                )
                expected = (
                    projected_area=1.0,
                    incident_par_q=100.0,
                    ri_par_0_f=100.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                cfg = _synthetic_cfg(cfg_ref; inputs.cfg...)
                turtle = ArchimedLight.build_turtle(cfg, inputs.sky)
                fluxes = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle, cfg)
                first = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
                budget = ArchimedLight.integrate_light(first, nothing, cfg; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(first.projected_area_per_node, 1, 0.0), expected.projected_area; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), expected.incident_par_q; atol=1e-10, rtol=1e-10)
                @test isapprox(get(budget.ri_par_0_f_per_node, 1, 0.0), expected.ri_par_0_f; atol=1e-10, rtol=1e-10)
            end
        end

        if _synthetic_case_enabled("stacked_scattering")
            @testset "Scenario: two aligned plates, lower receives only scattered PAR" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0),
                    cfg_no_scat=(; sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01),
                    cfg_scat=(; sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01),
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                cfg_no_scat = _synthetic_cfg(cfg_ref; inputs.cfg_no_scat...)
                turtle_no_scat = ArchimedLight.build_turtle(cfg_no_scat, inputs.sky)
                flux_no_scat = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_no_scat, cfg_no_scat)
                first_no_scat = ArchimedLight.compute_first_order(scene, turtle_no_scat, flux_no_scat, cfg_no_scat)
                budget_no_scat = ArchimedLight.integrate_light(first_no_scat, nothing, cfg_no_scat; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(first_no_scat.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
                @test isapprox(get(budget_no_scat.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget_no_scat.ri_par_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

                cfg_scat = _synthetic_cfg(cfg_ref; inputs.cfg_scat...)
                turtle_scat = ArchimedLight.build_turtle(cfg_scat, inputs.sky)
                flux_scat = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_scat, cfg_scat)
                first_scat = ArchimedLight.compute_first_order(scene, turtle_scat, flux_scat, cfg_scat)
                scat = ArchimedLight.compute_scattering(scene, turtle_scat, first_scat, cfg_scat)
                budget_scat = ArchimedLight.integrate_light(first_scat, scat, cfg_scat; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test get(scat.added_par_power_per_node, 2, 0.0) > 0.0
                @test isapprox(get(budget_scat.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
                @test get(budget_scat.ri_par_q_per_node, 2, 0.0) > 0.0
                @test get(budget_scat.ri_par_q_per_node, 2, 0.0) ≈ get(scat.added_par_power_per_node, 2, 0.0) atol = 1e-10 rtol = 1e-10
            end
        end

        if _synthetic_case_enabled("partial_overlap_direct")
            @testset "Scenario: upper half-plate shadows half of lower plate" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=0.5, y0=0.0, y1=1.0, z=1.0, group="upper_half", type="plate", object_id=1),
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower_full", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0),
                    cfg=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01),
                )
                expected = (
                    upper_projected=0.5,
                    lower_projected=0.5,
                    upper_ri_par_0_q=50.0,
                    lower_ri_par_0_q=50.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                cfg = _synthetic_cfg(cfg_ref; inputs.cfg...)
                turtle = ArchimedLight.build_turtle(cfg, inputs.sky)
                fluxes = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle, cfg)
                first = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
                budget = ArchimedLight.integrate_light(first, nothing, cfg; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(first.projected_area_per_node, 1, 0.0), expected.upper_projected; atol=1e-10, rtol=1e-10)
                @test isapprox(get(first.projected_area_per_node, 2, 0.0), expected.lower_projected; atol=1e-10, rtol=1e-10)
                @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), expected.upper_ri_par_0_q; atol=1e-9, rtol=1e-9)
                @test isapprox(get(budget.ri_par_0_q_per_node, 2, 0.0), expected.lower_ri_par_0_q; atol=1e-9, rtol=1e-9)
            end
        end

        if _synthetic_case_enabled("tilted_plate_projection")
            @testset "Scenario: 1 m² plate tilted by 60 degrees" begin
                tilt_deg = 60.0
                c = cosd(tilt_deg)
                s = sind(tilt_deg)
                inputs = (
                    scene=[(
                        p1=(0.0, 0.0, 0.0),
                        p2=(1.0, 0.0, 0.0),
                        p3=(1.0, c, s),
                        p4=(0.0, c, s),
                        group="tilted",
                        type="plate",
                        object_id=1,
                    )],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0),
                    cfg=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01),
                )
                expected = (area=1.0, projected_area=c, ri_par_0_q=100.0 * c)

                scene = _synthetic_quad_scene(inputs.scene)
                cfg = _synthetic_cfg(cfg_ref; inputs.cfg...)
                turtle = ArchimedLight.build_turtle(cfg, inputs.sky)
                fluxes = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle, cfg)
                first = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
                budget = ArchimedLight.integrate_light(first, nothing, cfg; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(scene.total_area_per_node[1], expected.area; atol=1e-10, rtol=1e-10)
                @test isapprox(get(first.projected_area_per_node, 1, 0.0), expected.projected_area; atol=1e-10, rtol=1e-10)
                @test isapprox(get(budget.ri_par_0_q_per_node, 1, 0.0), expected.ri_par_0_q; atol=1e-9, rtol=1e-9)
            end
        end

        if _synthetic_case_enabled("oblique_shadow")
            @testset "Scenario: oblique sun shadow with lateral offset plates" begin
                inputs = (
                    scene=[
                        (x0=1.0, x1=2.0, y0=0.0, y1=1.0, z=1.0, group="upper_offset", type="plate", object_id=1),
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0, group="lower_target", type="plate", object_id=2),
                    ],
                    sky_shadow=ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0),
                    sky_clear=ArchimedLight.SkyState(0.0, 45.0, 100.0, 0.0, 1.0, 0.0),
                    cfg=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01),
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                cfg = _synthetic_cfg(cfg_ref; inputs.cfg...)

                turtle_shadow = ArchimedLight.build_turtle(cfg, inputs.sky_shadow)
                flux_shadow = ArchimedLight.compute_directional_fluxes(inputs.sky_shadow, turtle_shadow, cfg)
                first_shadow = ArchimedLight.compute_first_order(scene, turtle_shadow, flux_shadow, cfg)

                turtle_clear = ArchimedLight.build_turtle(cfg, inputs.sky_clear)
                flux_clear = ArchimedLight.compute_directional_fluxes(inputs.sky_clear, turtle_clear, cfg)
                first_clear = ArchimedLight.compute_first_order(scene, turtle_clear, flux_clear, cfg)

                @test isapprox(get(first_shadow.incident_par_power_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
                @test isapprox(get(first_clear.incident_par_power_per_node, 2, 0.0), 100.0; atol=1e-5, rtol=1e-7)
                @test get(first_shadow.incident_par_power_per_node, 2, 0.0) < 0.02 * get(first_clear.incident_par_power_per_node, 2, 0.0)
            end
        end

        if _synthetic_case_enabled("toricity_wraparound")
            @testset "Scenario: toricity wraps edge-crossing plate across plot border" begin
                inputs = (
                    scene=[
                        (x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1),
                    ],
                    sky=ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0),
                    cfg_notoric=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false),
                    cfg_toric=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true),
                )
                expected = (
                    no_toric_projected_area=0.0,
                    no_toric_ri_par_0_q=0.0,
                    toric_projected_area=0.4,
                    toric_ri_par_0_q=40.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)

                cfg_notoric = _synthetic_cfg(cfg_ref; inputs.cfg_notoric...)
                turtle_notoric = ArchimedLight.build_turtle(cfg_notoric, inputs.sky)
                flux_notoric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_notoric, cfg_notoric)
                first_notoric = ArchimedLight.compute_first_order(scene, turtle_notoric, flux_notoric, cfg_notoric)
                budget_notoric = ArchimedLight.integrate_light(first_notoric, nothing, cfg_notoric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                cfg_toric = _synthetic_cfg(cfg_ref; inputs.cfg_toric...)
                turtle_toric = ArchimedLight.build_turtle(cfg_toric, inputs.sky)
                flux_toric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_toric, cfg_toric)
                first_toric = ArchimedLight.compute_first_order(scene, turtle_toric, flux_toric, cfg_toric)
                budget_toric = ArchimedLight.integrate_light(first_toric, nothing, cfg_toric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(first_notoric.projected_area_per_node, 1, 0.0), expected.no_toric_projected_area; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget_notoric.ri_par_0_q_per_node, 1, 0.0), expected.no_toric_ri_par_0_q; atol=1e-12, rtol=1e-12)
                @test isapprox(get(first_toric.projected_area_per_node, 1, 0.0), expected.toric_projected_area; atol=1e-6, rtol=1e-9)
                @test isapprox(get(budget_toric.ri_par_0_q_per_node, 1, 0.0), expected.toric_ri_par_0_q; atol=1e-5, rtol=1e-9)
                @test get(first_toric.projected_area_per_node, 1, 0.0) > get(first_notoric.projected_area_per_node, 1, 0.0)
            end
        end

        if _synthetic_case_enabled("toricity_cross_border_shadow")
            @testset "Scenario: toricity lets a border-adjacent upper plate shadow across the plot edge" begin
                inputs = (
                    scene=[
                        (x0=0.85, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper_edge", type="plate", object_id=1),
                        (x0=0.0, x1=0.9, y0=0.0, y1=1.0, z=0.0, group="lower_target", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0),
                    cfg_notoric=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false),
                    cfg_toric=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true),
                )
                expected = (
                    no_toric_upper_ri_par_0_q=0.0,
                    no_toric_lower_ri_par_0_q=90.0,
                    toric_upper_projected_area=0.15,
                    toric_upper_ri_par_0_q=15.0,
                    toric_lower_projected_area=0.85,
                    toric_lower_ri_par_0_q=85.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)

                cfg_notoric = _synthetic_cfg(cfg_ref; inputs.cfg_notoric...)
                turtle_notoric = ArchimedLight.build_turtle(cfg_notoric, inputs.sky)
                flux_notoric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_notoric, cfg_notoric)
                first_notoric = ArchimedLight.compute_first_order(scene, turtle_notoric, flux_notoric, cfg_notoric)
                budget_notoric = ArchimedLight.integrate_light(first_notoric, nothing, cfg_notoric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                cfg_toric = _synthetic_cfg(cfg_ref; inputs.cfg_toric...)
                turtle_toric = ArchimedLight.build_turtle(cfg_toric, inputs.sky)
                flux_toric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_toric, cfg_toric)
                first_toric = ArchimedLight.compute_first_order(scene, turtle_toric, flux_toric, cfg_toric)
                budget_toric = ArchimedLight.integrate_light(first_toric, nothing, cfg_toric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(budget_notoric.ri_par_0_q_per_node, 1, 0.0), expected.no_toric_upper_ri_par_0_q; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget_notoric.ri_par_0_q_per_node, 2, 0.0), expected.no_toric_lower_ri_par_0_q; atol=1e-5, rtol=1e-9)

                @test get(first_toric.projected_area_per_node, 1, 0.0) > 0.0
                @test get(first_toric.projected_area_per_node, 2, 0.0) < get(first_notoric.projected_area_per_node, 2, 0.0)
                @test isapprox(
                    get(first_toric.projected_area_per_node, 1, 0.0) + get(first_toric.projected_area_per_node, 2, 0.0),
                    1.0;
                    atol=2e-2,
                    rtol=0.0,
                )
                @test get(budget_toric.ri_par_0_q_per_node, 1, 0.0) > get(budget_notoric.ri_par_0_q_per_node, 1, 0.0)
                @test get(budget_toric.ri_par_0_q_per_node, 2, 0.0) < get(budget_notoric.ri_par_0_q_per_node, 2, 0.0)
                @test isapprox(
                    get(budget_toric.ri_par_0_q_per_node, 1, 0.0) + get(budget_toric.ri_par_0_q_per_node, 2, 0.0),
                    100.0;
                    atol=2.0,
                    rtol=0.0,
                )
            end
        end

        if _synthetic_case_enabled("toricity_diffuse_cross_border_shadow")
            @testset "Scenario: toricity changes diffuse sky interception across the plot edge" begin
                inputs = (
                    scene=[
                        (x0=0.85, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper_edge", type="plate", object_id=1),
                        (x0=0.0, x1=0.9, y0=0.0, y1=1.0, z=0.0, group="lower_target", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 0.0, 1.0),
                    cfg_notoric=(; sectors=16, all_in_turtle=true, scattering=false, pixel_size=0.01, toricity=false),
                    cfg_toric=(; sectors=16, all_in_turtle=true, scattering=false, pixel_size=0.01, toricity=true),
                )
                expected = (
                    no_toric_upper_ri_par_0_q=4.678040747015135,
                    no_toric_lower_ri_par_0_q=87.3445321040145,
                    toric_upper_ri_par_0_q=14.999998069798712,
                    toric_lower_ri_par_0_q=79.13276575225333,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)

                cfg_notoric = _synthetic_cfg(cfg_ref; inputs.cfg_notoric...)
                turtle_notoric = ArchimedLight.build_turtle(cfg_notoric, inputs.sky)
                flux_notoric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_notoric, cfg_notoric)
                first_notoric = ArchimedLight.compute_first_order(scene, turtle_notoric, flux_notoric, cfg_notoric)
                budget_notoric = ArchimedLight.integrate_light(first_notoric, nothing, cfg_notoric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                cfg_toric = _synthetic_cfg(cfg_ref; inputs.cfg_toric...)
                turtle_toric = ArchimedLight.build_turtle(cfg_toric, inputs.sky)
                flux_toric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_toric, cfg_toric)
                first_toric = ArchimedLight.compute_first_order(scene, turtle_toric, flux_toric, cfg_toric)
                budget_toric = ArchimedLight.integrate_light(first_toric, nothing, cfg_toric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test get(budget_notoric.ri_par_0_q_per_node, 1, 0.0) > 0.0
                @test get(budget_notoric.ri_par_0_q_per_node, 2, 0.0) > get(budget_notoric.ri_par_0_q_per_node, 1, 0.0)
                @test get(budget_toric.ri_par_0_q_per_node, 1, 0.0) > get(budget_notoric.ri_par_0_q_per_node, 1, 0.0)
                @test get(budget_toric.ri_par_0_q_per_node, 2, 0.0) < get(budget_notoric.ri_par_0_q_per_node, 2, 0.0)
                @test get(budget_toric.ri_par_0_q_per_node, 1, 0.0) > 10.0
                @test get(budget_toric.ri_par_0_q_per_node, 2, 0.0) > 70.0
            end
        end

        if _synthetic_case_enabled("toricity_scattering_cross_border")
            @testset "Scenario: toricity creates a wrapped scattering source across the plot edge" begin
                inputs = (
                    scene=[
                        (x0=0.8, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper_edge", type="plate", object_id=1),
                        (x0=0.0, x1=0.8, y0=0.0, y1=1.0, z=0.0, group="lower_target", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0),
                    cfg_notoric=(; sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false),
                    cfg_toric=(; sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=true),
                )
                expected = (
                    no_toric_upper_ri_par_0_q=0.0,
                    no_toric_upper_ri_par_q=0.0,
                    no_toric_lower_ri_par_0_q=80.00000119209906,
                    no_toric_lower_ri_par_q=80.00000119209906,
                    toric_upper_ri_par_0_q=20.00000476837205,
                    toric_upper_ri_par_q=20.678303998693803,
                    toric_lower_ri_par_0_q=80.00000119209906,
                    toric_lower_ri_par_q=80.6975728110636,
                    toric_lower_scattered_par_q=0.6975716189645285,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)

                cfg_notoric = _synthetic_cfg(cfg_ref; inputs.cfg_notoric...)
                turtle_notoric = ArchimedLight.build_turtle(cfg_notoric, inputs.sky)
                flux_notoric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_notoric, cfg_notoric)
                first_notoric = ArchimedLight.compute_first_order(scene, turtle_notoric, flux_notoric, cfg_notoric)
                scat_notoric = ArchimedLight.compute_scattering(scene, turtle_notoric, first_notoric, cfg_notoric)
                budget_notoric = ArchimedLight.integrate_light(first_notoric, scat_notoric, cfg_notoric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                cfg_toric = _synthetic_cfg(cfg_ref; inputs.cfg_toric...)
                turtle_toric = ArchimedLight.build_turtle(cfg_toric, inputs.sky)
                flux_toric = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_toric, cfg_toric)
                first_toric = ArchimedLight.compute_first_order(scene, turtle_toric, flux_toric, cfg_toric)
                scat_toric = ArchimedLight.compute_scattering(scene, turtle_toric, first_toric, cfg_toric)
                budget_toric = ArchimedLight.integrate_light(first_toric, scat_toric, cfg_toric; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                @test isapprox(get(budget_notoric.ri_par_0_q_per_node, 1, 0.0), expected.no_toric_upper_ri_par_0_q; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget_notoric.ri_par_q_per_node, 1, 0.0), expected.no_toric_upper_ri_par_q; atol=1e-12, rtol=1e-12)
                @test isapprox(get(budget_notoric.ri_par_0_q_per_node, 2, 0.0), expected.no_toric_lower_ri_par_0_q; atol=1e-5, rtol=1e-9)
                @test isapprox(get(budget_notoric.ri_par_q_per_node, 2, 0.0), expected.no_toric_lower_ri_par_q; atol=1e-5, rtol=1e-9)

                @test get(budget_toric.ri_par_0_q_per_node, 1, 0.0) > 0.0
                @test get(budget_toric.ri_par_q_per_node, 1, 0.0) > get(budget_toric.ri_par_0_q_per_node, 1, 0.0)
                @test get(budget_toric.ri_par_0_q_per_node, 2, 0.0) >= get(budget_notoric.ri_par_0_q_per_node, 2, 0.0)
                @test get(budget_toric.ri_par_q_per_node, 2, 0.0) > get(budget_toric.ri_par_0_q_per_node, 2, 0.0)
                @test get(scat_toric.added_par_power_per_node, 2, 0.0) > 0.0
                @test isapprox(
                    get(budget_toric.ri_par_q_per_node, 2, 0.0) - get(budget_toric.ri_par_0_q_per_node, 2, 0.0),
                    get(scat_toric.added_par_power_per_node, 2, 0.0);
                    atol=1e-10,
                    rtol=1e-10,
                )
                @test get(budget_toric.ri_par_q_per_node, 1, 0.0) > get(budget_notoric.ri_par_q_per_node, 1, 0.0)
                @test get(budget_toric.ri_par_q_per_node, 2, 0.0) > get(budget_notoric.ri_par_q_per_node, 2, 0.0)
            end
        end

        if _synthetic_case_enabled("virtual_sensor_transparency")
            @testset "Scenario: virtual sensor receives light without blocking lower plate" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensors", type="pave", object_id=1),
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
                    ],
                    sky=ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0),
                    cfg=(; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01),
                    sensor_model="""
                    ---
                    Group: sensors
                    Type:
                      pave:
                          Interception:
                              model: VirtualSensor
                    ...
                    """,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                cfg_opaque = _synthetic_cfg(cfg_ref; inputs.cfg...)
                turtle_opaque = ArchimedLight.build_turtle(cfg_opaque, inputs.sky)
                flux_opaque = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_opaque, cfg_opaque)
                first_opaque = ArchimedLight.compute_first_order(scene, turtle_opaque, flux_opaque, cfg_opaque)
                budget_opaque = ArchimedLight.integrate_light(first_opaque, nothing, cfg_opaque; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                mktempdir() do tmp
                    sensor_model = joinpath(tmp, "model_sensor.yml")
                    write(sensor_model, inputs.sensor_model)
                    cfg_sensor = _synthetic_cfg(cfg_ref; models=[sensor_model], inputs.cfg...)
                    turtle_sensor = ArchimedLight.build_turtle(cfg_sensor, inputs.sky)
                    flux_sensor = ArchimedLight.compute_directional_fluxes(inputs.sky, turtle_sensor, cfg_sensor)
                    first_sensor = ArchimedLight.compute_first_order(scene, turtle_sensor, flux_sensor, cfg_sensor)
                    budget_sensor = ArchimedLight.integrate_light(first_sensor, nothing, cfg_sensor; step_duration_seconds=1.0, component_area_per_node=scene.total_area_per_node)

                    @test isapprox(get(budget_opaque.ri_par_0_q_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
                    @test isapprox(get(budget_sensor.ri_par_0_q_per_node, 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(budget_sensor.ri_par_0_q_per_node, 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
                end
            end
        end

        if _synthetic_case_enabled("single_plate_absorptance")
            @testset "Scenario: single 1 m² plate with explicit absorptance" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="absorber", type="plate", object_id=1),
                    ],
                    meteo=_synthetic_meteo_row(; duration_seconds=2.0, ri_par_f=100.0, ri_nir_f=80.0, direct_fraction=1.0),
                    model_yaml="""
                    ---
                    Group: absorber
                    Type:
                      plate:
                          Interception:
                              model: Translucent
                              transparency: 0
                              optical_properties:
                                  PAR: 0.2
                                  NIR: 0.7
                    ...
                    """,
                )
                expected = (
                    ri_par_0_q=200.0,
                    ri_nir_0_q=160.0,
                    ra_par_0_q=160.0,
                    ra_nir_0_q=48.0,
                    ra_par_0_f=80.0,
                    ra_nir_0_f=24.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                mktempdir() do tmp
                    model_path = joinpath(tmp, "model_absorber.yml")
                    write(model_path, inputs.model_yaml)
                    cfg = _synthetic_cfg(cfg_ref; models=[model_path])
                    step = ArchimedLight.run_light_step(scene, inputs.meteo, cfg)

                    @test isapprox(get(step.budget.ri_par_0_q_per_node, 1, 0.0), expected.ri_par_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(step.budget.ri_nir_0_q_per_node, 1, 0.0), expected.ri_nir_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(step.budget.ra_par_0_q_per_node, 1, 0.0), expected.ra_par_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(step.budget.ra_nir_0_q_per_node, 1, 0.0), expected.ra_nir_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(step.budget.ra_par_0_f_per_node, 1, 0.0), expected.ra_par_0_f; atol=1e-9, rtol=1e-9)
                    @test isapprox(get(step.budget.ra_nir_0_f_per_node, 1, 0.0), expected.ra_nir_0_f; atol=1e-9, rtol=1e-9)
                end
            end
        end

        if _synthetic_case_enabled("two_planes_shadow_absorptance")
            @testset "Scenario: two 1 m² plates, upper blocks direct light" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
                    ],
                    meteo=_synthetic_meteo_row(; duration_seconds=2.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
                    model_yaml="""
                    ---
                    Group: upper
                    Type:
                      plate:
                          Interception:
                              model: Translucent
                              transparency: 0
                              optical_properties:
                                  PAR: 0.25
                                  NIR: 0.5
                    ---
                    Group: lower
                    Type:
                      plate:
                          Interception:
                              model: Translucent
                              transparency: 0
                              optical_properties:
                                  PAR: 0.4
                                  NIR: 0.2
                    ...
                    """,
                )
                expected = (
                    upper_ri_par_0_q=240.0,
                    upper_ri_nir_0_q=160.0,
                    upper_ra_par_0_q=180.0,
                    upper_ra_nir_0_q=80.0,
                    lower_ri_par_0_q=0.0,
                    lower_ri_nir_0_q=0.0,
                    lower_ra_par_0_q=0.0,
                    lower_ra_nir_0_q=0.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                mktempdir() do tmp
                    model_path = joinpath(tmp, "model_two_plates.yml")
                    write(model_path, inputs.model_yaml)
                    cfg = _synthetic_cfg(cfg_ref; sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, models=[model_path])
                    step = ArchimedLight.run_light_step(scene, inputs.meteo, cfg)

                    rows = ArchimedLight.component_values_table(
                        scene,
                        step,
                        cfg;
                        meteo_row=inputs.meteo,
                        step_number=1,
                        step_duration_seconds=2.0,
                        columns=["node_id", "Ri_PAR_0_q", "Ri_NIR_0_q", "Ra_PAR_0_q", "Ra_NIR_0_q"],
                    ).rows
                    row_by_node = Dict(Int(r["node_id"]) => r for r in rows)

                    @test isapprox(Float64(row_by_node[1]["Ri_PAR_0_q"]), expected.upper_ri_par_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(Float64(row_by_node[1]["Ri_NIR_0_q"]), expected.upper_ri_nir_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(Float64(row_by_node[1]["Ra_PAR_0_q"]), expected.upper_ra_par_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(Float64(row_by_node[1]["Ra_NIR_0_q"]), expected.upper_ra_nir_0_q; atol=1e-9, rtol=1e-9)
                    @test isapprox(Float64(row_by_node[2]["Ri_PAR_0_q"]), expected.lower_ri_par_0_q; atol=1e-12, rtol=1e-12)
                    @test isapprox(Float64(row_by_node[2]["Ri_NIR_0_q"]), expected.lower_ri_nir_0_q; atol=1e-12, rtol=1e-12)
                    @test isapprox(Float64(row_by_node[2]["Ra_PAR_0_q"]), expected.lower_ra_par_0_q; atol=1e-12, rtol=1e-12)
                    @test isapprox(Float64(row_by_node[2]["Ra_NIR_0_q"]), expected.lower_ra_nir_0_q; atol=1e-12, rtol=1e-12)
                end
            end
        end

        if _synthetic_case_enabled("cached_series_parity")
            @testset "Scenario: repeated simple meteo rows with and without radiation cache" begin
                inputs = (
                    scene=[
                        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
                    ],
                    meteo_rows=[
                        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
                        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
                        _synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0, direct_fraction=1.0),
                    ],
                )
                expected = (
                    ri_par_f=120.0,
                    ri_nir_f=80.0,
                    step1_ri_par_q=120.0 * 600.0,
                    step1_ri_nir_q=80.0 * 600.0,
                    step2_scale=3.0,
                )

                scene = _synthetic_horizontal_scene(inputs.scene)
                meteo = ArchimedLight.MeteoTable(inputs.meteo_rows, (; file="synthetic-cached-series"))
                cfg_uncached = _synthetic_cfg(cfg_ref; cache_radiation=false)
                cfg_cached = _synthetic_cfg(cfg_ref; cache_radiation=true)
                series_uncached = ArchimedLight.run_light_series(scene, meteo, cfg_uncached)
                series_cached = ArchimedLight.run_light_series(scene, meteo, cfg_cached)

                @test length(series_uncached) == 3
                @test length(series_cached) == 3
                for i in eachindex(series_uncached)
                    @test _max_abs_float_dict_diff(series_uncached[i].budget.ri_par_f_per_node, series_cached[i].budget.ri_par_f_per_node) == 0.0
                    @test _max_abs_float_dict_diff(series_uncached[i].budget.ri_nir_f_per_node, series_cached[i].budget.ri_nir_f_per_node) == 0.0
                    @test _max_abs_float_dict_diff(series_uncached[i].budget.ri_par_q_per_node, series_cached[i].budget.ri_par_q_per_node) == 0.0
                    @test _max_abs_float_dict_diff(series_uncached[i].budget.ri_nir_q_per_node, series_cached[i].budget.ri_nir_q_per_node) == 0.0
                end

                @test isapprox(get(series_uncached[1].budget.ri_par_f_per_node, 1, 0.0), expected.ri_par_f; atol=1e-9, rtol=1e-9)
                @test isapprox(get(series_uncached[1].budget.ri_nir_f_per_node, 1, 0.0), expected.ri_nir_f; atol=1e-9, rtol=1e-9)
                @test isapprox(get(series_uncached[1].budget.ri_par_q_per_node, 1, 0.0), expected.step1_ri_par_q; atol=1e-6, rtol=1e-12)
                @test isapprox(get(series_uncached[1].budget.ri_nir_q_per_node, 1, 0.0), expected.step1_ri_nir_q; atol=1e-6, rtol=1e-12)
                @test isapprox(
                    get(series_uncached[2].budget.ri_par_q_per_node, 1, 0.0),
                    expected.step2_scale * expected.step1_ri_par_q;
                    atol=1e-6,
                    rtol=1e-12,
                )
                @test isapprox(
                    get(series_uncached[2].budget.ri_nir_q_per_node, 1, 0.0),
                    expected.step2_scale * expected.step1_ri_nir_q;
                    atol=1e-6,
                    rtol=1e-12,
                )
            end
        end
end
