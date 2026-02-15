import CSV
import PlantGeom
import PlantMeteo
import GeometryBasics
import StaticArrays
import YAML
import Tables
import MultiScaleTreeGraph
import LinearAlgebra: norm, cross

function _to_string_dict(x)
    out = Dict{String,Any}()
    for (k, v) in x
        ks = string(k)
        out[ks] = v isa AbstractDict ? _to_string_dict(v) : v
    end
    out
end

function _cfg_get(d::Dict{String,Any}, keys::Vector{String}, default)
    for k in keys
        haskey(d, k) && return d[k]
    end
    return default
end

function _as_bool(x, default::Bool)
    x === nothing && return default
    x isa Bool && return x
    x isa Number && return x != 0
    x isa String && return lowercase(strip(x)) in ("1", "true", "yes", "y", "on")
    return default
end

function _as_int(x, default::Int)
    x === nothing && return default
    x isa Integer && return Int(x)
    x isa Number && return round(Int, x)
    x isa String && return parse(Int, x)
    return default
end

function _as_float(x, default::Float64)
    x === nothing && return default
    x isa Number && return Float64(x)
    x isa String && return parse(Float64, x)
    return default
end

function _join_if_relative(base::AbstractString, p::AbstractString)
    isabspath(p) ? p : normpath(joinpath(base, p))
end

function read_light_config(path::AbstractString)
    raw = YAML.load_file(path)
    d = raw isa AbstractDict ? _to_string_dict(raw) : Dict{String,Any}()
    base = dirname(path)

    scene_rel = _cfg_get(d, ["scene", "scene_file"], "")
    meteo_rel = _cfg_get(d, ["meteo", "meteo_file"], "")
    scene = _join_if_relative(base, string(scene_rel))
    meteo = _join_if_relative(base, string(meteo_rel))

    props = get(d, "prop", Dict{String,Any}())
    propsd = props isa AbstractDict ? _to_string_dict(props) : Dict{String,Any}()

    all_in_turtle = _as_bool(_cfg_get(propsd, ["all_in_turtle"], get(d, "all_in_turtle", nothing)), false)
    turtle_sectors = _as_int(_cfg_get(propsd, ["turtle_nb_sectors", "turtle_sectors"], get(d, "turtle_sectors", nothing)), 46)
    pixel_size = _as_float(_cfg_get(propsd, ["pixel_size"], get(d, "pixel_size", nothing)), 0.25)
    area_ratio = _as_bool(_cfg_get(propsd, ["area_ratio"], get(d, "area_ratio", nothing)), true)
    scattering = _as_bool(_cfg_get(propsd, ["scattering"], get(d, "scattering", nothing)), true)
    scattering_max_iter = _as_int(
        _cfg_get(propsd, ["scattering_max_iter", "scat_max_iter"], get(d, "scattering_max_iter", nothing)),
        20
    )
    scattering_stop_ratio = _as_float(
        _cfg_get(propsd, ["scattering_stop_ratio", "scat_stop_ratio"], get(d, "scattering_stop_ratio", nothing)),
        0.01
    )
    scattering_coeff_par = _as_float(_cfg_get(propsd, ["scattering_coeff_par"], get(d, "scattering_coeff_par", nothing)), 0.15)
    scattering_coeff_nir = _as_float(_cfg_get(propsd, ["scattering_coeff_nir"], get(d, "scattering_coeff_nir", nothing)), 0.30)
    cache_radiation = _as_bool(_cfg_get(propsd, ["cache_radiation"], get(d, "cache_radiation", nothing)), false)

    return LightConfig(
        scene,
        meteo,
        all_in_turtle,
        turtle_sectors,
        pixel_size,
        area_ratio,
        scattering,
        scattering_max_iter,
        scattering_stop_ratio,
        scattering_coeff_par,
        scattering_coeff_nir,
        cache_radiation,
        d
    )
end

function _triangle_area3d(p1, p2, p3)
    v1 = StaticArrays.SVector{3,Float64}(p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = StaticArrays.SVector{3,Float64}(p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
    0.5 * norm(cross(v1, v2))
end

function _build_scene_geometry(mtg, source_path::AbstractString)
    merged_mesh, face2node = PlantGeom.build_merged_mesh_with_map(mtg; filter_fun=node -> !isnothing(node.geometry))

    verts = GeometryBasics.decompose(PlantGeom.Point3, merged_mesh)
    faces = GeometryBasics.decompose(PlantGeom.Face3, merged_mesh)
    node_area = Dict{Int,Float64}()
    for (i, f) in enumerate(faces)
        n = face2node[i]
        area = _triangle_area3d(verts[f[1]], verts[f[2]], verts[f[3]])
        node_area[n] = get(node_area, n, 0.0) + area
    end
    SceneGeometry(mtg, merged_mesh, face2node, node_area, String(source_path))
end

function read_scene(path::AbstractString; plantgeom_backend=:auto)
    ext = lowercase(splitext(path)[2])
    mtg =
        if ext == ".ops"
            PlantGeom.read_ops(path)
        elseif ext == ".opf"
            PlantGeom.read_opf(path)
        elseif ext == ".gwa"
            PlantGeom.read_gwa(path)
        else
            error("Unsupported scene extension: $ext")
        end
    _build_scene_geometry(mtg, path)
end

function _rows_to_namedtuples(table)
    Tables.rowtable(table) |> collect
end

function read_meteo(path::AbstractString)
    weather, meta =
        try
            w = PlantMeteo.read_weather(path)
            m = try
                PlantMeteo.metadata(w)
            catch
                (; file=path)
            end
            (w, m)
        catch
            f = CSV.File(path; comment="#", delim=';', ignorerepeated=true)
            (f, (; file=path))
        end
    rows = _rows_to_namedtuples(weather)
    MeteoTable(rows, NamedTuple(meta))
end
