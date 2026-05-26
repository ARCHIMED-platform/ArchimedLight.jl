using Dates
using GeometryBasics
using LinearAlgebra: cross, norm
using MultiScaleTreeGraph
using OrderedCollections: OrderedDict
using PlantGeom
using StaticArrays: SVector

function _synthetic_options(;
    sectors::Int=1,
    all_in_turtle::Bool=false,
    scattering::Bool=false,
    pixel_size_m::Float64=0.01,
    toricity::Bool=false,
    area_ratio::Bool=true,
    cache_radiation::Bool=false,
    cache_pixel_table::Bool=false,
    radiation_timestep::Int=15,
    nir_interception::Bool=true,
    nir_scattering::Bool=true,
    java_logged_turtle_dirs::Bool=false,
)
    ArchimedLight.LightOptions(
        turtle_sectors=sectors,
        all_in_turtle=all_in_turtle,
        scattering=scattering,
        pixel_size=pixel_size_m,
        toricity=toricity,
        area_ratio=area_ratio,
        cache_radiation=cache_radiation,
        cache_pixel_table=cache_pixel_table,
        radiation_timestep_minutes=radiation_timestep,
        nir_interception=nir_interception,
        nir_scattering=nir_scattering,
        java_logged_turtle_dirs=java_logged_turtle_dirs,
    )
end

function _default_synthetic_models()
    ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=OrderedDict(
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

function _synthetic_horizontal_scene(specs::AbstractVector{<:NamedTuple})
    quad_specs = map(specs) do spec
        (
            p1=(spec.x0, spec.y0, spec.z),
            p2=(spec.x1, spec.y0, spec.z),
            p3=(spec.x1, spec.y1, spec.z),
            p4=(spec.x0, spec.y1, spec.z),
            source_topology_id=get(spec, :source_topology_id, 1),
            object_id=get(spec, :object_id, 1),
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
    nodes = Dict{Int,PlantGeom.SceneNodeData{Float64}}()
    mtg = MultiScaleTreeGraph.Node(
        MultiScaleTreeGraph.MutableNodeMTG(:/, :Scene, 0, 0),
        Dict{Symbol,Any}(),
    )

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
        area = area1 + area2
        barycenter = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            (p1[3] + p2[3] + p3[3] + p4[3]) / 4,
        )
        source_topology_id = Int(get(spec, :source_topology_id, i))
        nodes[i] = PlantGeom.SceneNodeData(area, barycenter, source_topology_id)
        type_name = String(get(spec, :type, "plate"))
        group = String(get(spec, :group, "plate"))
        mesh = GeometryBasics.Mesh(
            GeometryBasics.Point{3,Float32}[
                GeometryBasics.Point{3,Float32}(Float32(p1[1]), Float32(p1[2]), Float32(p1[3])),
                GeometryBasics.Point{3,Float32}(Float32(p2[1]), Float32(p2[2]), Float32(p2[3])),
                GeometryBasics.Point{3,Float32}(Float32(p3[1]), Float32(p3[2]), Float32(p3[3])),
                GeometryBasics.Point{3,Float32}(Float32(p4[1]), Float32(p4[2]), Float32(p4[3])),
            ],
            PlantGeom.Face3[(1, 2, 3), (1, 3, 4)],
        )
        MultiScaleTreeGraph.Node(
            i,
            mtg,
            MultiScaleTreeGraph.MutableNodeMTG(:+, Symbol(type_name), i, 1),
            Dict{Symbol,Any}(
                :geometry => PlantGeom.Geometry(ref_mesh=PlantGeom.RefMesh("synthetic_$(i)", mesh)),
                :group => group,
                :functional_group => group,
                :type => type_name,
                :object_id => Int(get(spec, :object_id, source_topology_id)),
                :source_topology_id => source_topology_id,
            ),
        )
    end

    PlantGeom.SceneGeometry(
        mtg,
        GeometryBasics.Mesh(points, faces),
        face2node,
        nodes,
        "synthetic_scene_cases",
        (minimum(xs), minimum(ys), maximum(xs), maximum(ys)),
    )
end

function _synthetic_meteo_row(;
    date::Dates.Date=Dates.Date(2020, 6, 21),
    start_time::Dates.Time=Dates.Time(12),
    duration_seconds::Float64=1.0,
    latitude::Float64=0.0,
    ri_par_f::Float64=100.0,
    ri_nir_f::Float64=0.0,
    direct_fraction::Float64=1.0,
    sun_azimut::Float64=180.0,
    sun_elevation::Float64=90.0,
)
    start_dt = Dates.DateTime(date, start_time)
    end_dt = start_dt + Dates.Millisecond(round(Int, duration_seconds * 1000))
    (
        date=date,
        hour_start=Dates.format(Dates.Time(start_dt), Dates.DateFormat("HH:MM:SS")),
        hour_end=Dates.format(Dates.Time(end_dt), Dates.DateFormat("HH:MM:SS")),
        step_duration=duration_seconds,
        latitude=latitude,
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
    )
end

function _max_abs_float_dict_diff(a::Dict{Int,Float64}, b::Dict{Int,Float64})
    maximum(abs(get(a, id, 0.0) - get(b, id, 0.0)) for id in union(keys(a), keys(b)); init=0.0)
end

function _synthetic_exact_check(source_id::String, result)::NamedTuple
    if source_id == "single_plate_direct"
        step = result.step
        pa = get(step.first_order.projected_area_per_node, 1, 0.0)
        q = get(step.budget.incident_energy.initial.par, 1, 0.0)
        f = get(step.budget.incident_flux.initial.par, 1, 0.0)
        ok = isapprox(pa, 1.0; atol=1e-12, rtol=1e-12) &&
            isapprox(q, 100.0; atol=1e-10, rtol=1e-10) &&
            isapprox(f, 100.0; atol=1e-10, rtol=1e-10)
        detail = ok ? "" : "single_plate_direct projected=$(pa) q=$(q) f=$(f)"
        return (ok=ok, detail=detail)
    elseif source_id == "partial_overlap_direct"
        step = result.step
        upper_pa = get(step.first_order.projected_area_per_node, 1, 0.0)
        lower_pa = get(step.first_order.projected_area_per_node, 2, 0.0)
        upper_q = get(step.budget.incident_energy.initial.par, 1, 0.0)
        lower_q = get(step.budget.incident_energy.initial.par, 2, 0.0)
        ok = isapprox(upper_pa, 0.5; atol=1e-10, rtol=1e-10) &&
            isapprox(lower_pa, 0.5; atol=1e-10, rtol=1e-10) &&
            isapprox(upper_q, 50.0; atol=1e-9, rtol=1e-9) &&
            isapprox(lower_q, 50.0; atol=1e-9, rtol=1e-9)
        detail = ok ? "" : "partial_overlap_direct upper_pa=$(upper_pa) lower_pa=$(lower_pa) upper_q=$(upper_q) lower_q=$(lower_q)"
        return (ok=ok, detail=detail)
    elseif source_id == "toricity_wraparound"
        step = result.step
        pa = get(step.first_order.projected_area_per_node, 1, 0.0)
        q = get(step.budget.incident_energy.initial.par, 1, 0.0)
        toric = Bool(get(result.meta, "toricity", false))
        target_pa = toric ? 0.4 : 0.19512196866477877
        target_q = toric ? 40.0 : 19.512196866477876
        ok = isapprox(pa, target_pa; atol=1e-5, rtol=1e-9) &&
            isapprox(q, target_q; atol=1e-5, rtol=1e-9)
        detail = ok ? "" : "toricity_wraparound toric=$(toric) projected=$(pa) q=$(q)"
        return (ok=ok, detail=detail)
    elseif source_id == "stacked_scattering"
        step = result.step
        lower_0 = get(step.budget.incident_energy.initial.par, 2, 0.0)
        lower_q = get(step.budget.incident_energy.total.par, 2, 0.0)
        scat_q = isnothing(step.scattering) ? 0.0 : get(step.scattering.added_power.par, 2, 0.0)
        scattering = Bool(get(result.meta, "scattering", false))
        ok =
            if scattering
                scat_q > 0.0 &&
                    isapprox(lower_0, 0.0; atol=1e-12, rtol=1e-12) &&
                    lower_q > 0.0
            else
                isapprox(lower_0, 0.0; atol=1e-12, rtol=1e-12) &&
                    isapprox(lower_q, 0.0; atol=1e-12, rtol=1e-12)
            end
        detail = ok ? "" : "stacked_scattering scattering=$(scattering) lower_0=$(lower_0) lower_q=$(lower_q) scat_q=$(scat_q)"
        return (ok=ok, detail=detail)
    elseif source_id == "cached_series_parity"
        diffs = result.diffs
        ok = all(v == 0.0 for v in values(diffs))
        detail = ok ? "" : "cached parity diffs=$(diffs)"
        return (ok=ok, detail=detail)
    end
    return (ok=true, detail="")
end
