using Dates
using LinearAlgebra: norm, cross
using StaticArrays: SVector
using GeometryBasics
import PlantGeom

const _SYNTHETIC_TEST_PROFILE = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_PROFILE", "all"))
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
    cache_radiation::Bool=false,
    models::Vector{String}=String[],
)
    raw = copy(cfg.raw)
    raw["all_in_turtle"] = all_in_turtle
    raw["sky_sectors"] = sectors
    raw["scattering"] = scattering
    raw["pixel_size"] = pixel_size * 100.0
    raw["cache_radiation"] = cache_radiation
    raw["models"] = models
    ArchimedLight.LightConfig(
        cfg.scene,
        cfg.meteo,
        all_in_turtle,
        sectors,
        pixel_size,
        cfg.area_ratio,
        scattering,
        cfg.scattering_max_iter,
        cfg.scattering_stop_ratio,
        cfg.scattering_coeff_par,
        cfg.scattering_coeff_nir,
        cache_radiation,
        raw,
    )
end

function _synthetic_horizontal_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float32}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    total_area_per_node = Dict{Int,Float64}()
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    java_item_id_per_node = Dict{Int,Int}()
    java_component_id_per_node = Dict{Int,Int}()

    xs = Float64[]
    ys = Float64[]
    for (i, spec) in enumerate(specs)
        p1 = (spec.x0, spec.y0, spec.z)
        p2 = (spec.x1, spec.y0, spec.z)
        p3 = (spec.x1, spec.y1, spec.z)
        p4 = (spec.x0, spec.y1, spec.z)
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
        node_group[i] = String(get(spec, :group, "plate"))
        node_type[i] = String(get(spec, :type, "plate"))
        java_item_id_per_node[i] = Int(get(spec, :item_id, i))
        java_component_id_per_node[i] = Int(get(spec, :component_id, 1))
    end

    ArchimedLight.SceneGeometry(
        nothing,
        GeometryBasics.Mesh(points, faces),
        face2node,
        total_area_per_node,
        barycenter_per_node,
        node_group,
        node_type,
        java_item_id_per_node,
        java_component_id_per_node,
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

if _SYNTHETIC_TEST_PROFILE == "synthetic"
@testset "Synthetic scene cases" begin
    cfg_ref = ArchimedLight.read_light_config(
        joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests", "test-links-stats", "config.yml"),
    )

    if _synthetic_case_enabled("single_plate_absorptance")
    @testset "Scenario: single 1 m² plate with explicit absorptance" begin
        inputs = (
            scene=[
                (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="absorber", type="plate", item_id=1),
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
                (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", item_id=1),
                (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", item_id=2),
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
                columns=["item_id", "Ri_PAR_0_q", "Ri_NIR_0_q", "Ra_PAR_0_q", "Ra_NIR_0_q"],
            ).rows
            row_by_item = Dict(Int(r["item_id"]) => r for r in rows)

            @test isapprox(Float64(row_by_item[1]["Ri_PAR_0_q"]), expected.upper_ri_par_0_q; atol=1e-9, rtol=1e-9)
            @test isapprox(Float64(row_by_item[1]["Ri_NIR_0_q"]), expected.upper_ri_nir_0_q; atol=1e-9, rtol=1e-9)
            @test isapprox(Float64(row_by_item[1]["Ra_PAR_0_q"]), expected.upper_ra_par_0_q; atol=1e-9, rtol=1e-9)
            @test isapprox(Float64(row_by_item[1]["Ra_NIR_0_q"]), expected.upper_ra_nir_0_q; atol=1e-9, rtol=1e-9)
            @test isapprox(Float64(row_by_item[2]["Ri_PAR_0_q"]), expected.lower_ri_par_0_q; atol=1e-12, rtol=1e-12)
            @test isapprox(Float64(row_by_item[2]["Ri_NIR_0_q"]), expected.lower_ri_nir_0_q; atol=1e-12, rtol=1e-12)
            @test isapprox(Float64(row_by_item[2]["Ra_PAR_0_q"]), expected.lower_ra_par_0_q; atol=1e-12, rtol=1e-12)
            @test isapprox(Float64(row_by_item[2]["Ra_NIR_0_q"]), expected.lower_ra_nir_0_q; atol=1e-12, rtol=1e-12)
        end
    end
    end

    if _synthetic_case_enabled("cached_series_parity")
    @testset "Scenario: repeated simple meteo rows with and without radiation cache" begin
        inputs = (
            scene=[
                (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", item_id=1),
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
end
