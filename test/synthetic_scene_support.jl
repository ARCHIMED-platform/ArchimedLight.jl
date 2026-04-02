_incident_par_initial_flux(budget) = budget.incident_flux.initial.par
_incident_par_initial_energy(budget) = budget.incident_energy.initial.par
_incident_par_energy(budget) = budget.incident_energy.total.par

function _synthetic_options(;
    sectors::Int=1,
    all_in_turtle::Bool=false,
    scattering::Bool=false,
    pixel_size::Float64=0.01,
    toricity::Bool=true,
    cache_radiation::Bool=false,
)
    ArchimedLight.LightOptions(
        turtle_sectors=sectors,
        all_in_turtle=all_in_turtle,
        scattering=scattering,
        pixel_size=pixel_size,
        toricity=toricity,
        cache_radiation=cache_radiation,
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

function _virtual_sensor_models()
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
        ArchimedLight.GroupModel(
            "sensor";
            types=OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(model="VirtualSensor"),
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
            group=get(spec, :group, "plate"),
            type=get(spec, :type, "plate"),
            source_topology_id=get(spec, :source_topology_id, 1),
            object_id=get(spec, :object_id, 1),
        )
    end
    _synthetic_quad_scene(quad_specs)
end

function _synthetic_quad_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float32}[]
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    nodes = Dict{Int,ArchimedLight.SceneNodeData{Float64}}()

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
        append!(faces, GeometryBasics.TriangleFace{Int}[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
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
        object_id = Int(get(spec, :object_id, source_topology_id))
        group = String(get(spec, :group, "plate"))
        type_name = String(get(spec, :type, "plate"))
        nodes[i] = ArchimedLight.SceneNodeData(area, barycenter, group, type_name, source_topology_id, object_id)
    end

    ArchimedLight.SceneGeometry(
        nothing,
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

function _run_direct(scene, sky, options; models=_default_synthetic_models())
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    budget = ArchimedLight.integrate_light(scene, models, first, nothing, options; step_duration_seconds=1.0)
    (; turtle, fluxes, first, budget)
end
