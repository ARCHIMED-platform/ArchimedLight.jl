import GeometryBasics
import StaticArrays
import LinearAlgebra: cross, dot, norm
import PlantGeom

function _triangle_area_and_normal(p1, p2, p3)
    v1 = StaticArrays.SVector{3,Float64}(p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    v2 = StaticArrays.SVector{3,Float64}(p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3])
    n = cross(v1, v2)
    na = norm(n)
    if na == 0.0
        return 0.0, StaticArrays.SVector{3,Float64}(0.0, 0.0, 0.0)
    end
    return 0.5 * na, n / na
end

function _view_basis(direction)
    d = direction / max(norm(direction), eps(Float64))
    a = abs(d[3]) < 0.95 ? StaticArrays.SVector{3,Float64}(0.0, 0.0, 1.0) : StaticArrays.SVector{3,Float64}(0.0, 1.0, 0.0)
    u = cross(a, d)
    u = u / max(norm(u), eps(Float64))
    v = cross(d, u)
    return u, v, d
end

@inline function _edge(x1, y1, x2, y2, x, y)
    (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
end

function _project_vertices(vertices, u, v, d)
    n = length(vertices)
    px = Vector{Float64}(undef, n)
    py = Vector{Float64}(undef, n)
    pz = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        p = vertices[i]
        s = StaticArrays.SVector{3,Float64}(p[1], p[2], p[3])
        px[i] = dot(s, u)
        py[i] = dot(s, v)
        pz[i] = dot(s, d)
    end
    return px, py, pz
end

function _rasterize_direction(vertices, faces, face_area, face_normal, direction, pixel_size::Float64, area_ratio::Bool)
    u, v, d = _view_basis(direction)
    px, py, pz = _project_vertices(vertices, u, v, d)

    xmin = minimum(px)
    xmax = maximum(px)
    ymin = minimum(py)
    ymax = maximum(py)

    xspan = max(xmax - xmin, pixel_size)
    yspan = max(ymax - ymin, pixel_size)
    pix = max(pixel_size, eps(Float64))

    nx = max(1, ceil(Int, xspan / pix))
    ny = max(1, ceil(Int, yspan / pix))

    # Keep memory bounded and deterministic.
    max_pixels = 2_000_000
    if nx * ny > max_pixels
        pix = sqrt((xspan * yspan) / max_pixels)
        nx = max(1, ceil(Int, xspan / pix))
        ny = max(1, ceil(Int, yspan / pix))
    end

    pixel_area = pix * pix
    depth = fill(Inf, nx * ny)
    face_at_pixel = fill(0, nx * ny)
    projected_area = zeros(Float64, length(faces))
    hits = zeros(Int, length(faces))

    @inbounds for fi in eachindex(faces)
        nrm = face_normal[fi]
        cos_inc = dot(nrm, -d)
        if cos_inc <= 0.0 || face_area[fi] == 0.0
            continue
        end

        f = faces[fi]
        i1, i2, i3 = f[1], f[2], f[3]
        x1, y1, z1 = px[i1], py[i1], pz[i1]
        x2, y2, z2 = px[i2], py[i2], pz[i2]
        x3, y3, z3 = px[i3], py[i3], pz[i3]

        a2 = _edge(x1, y1, x2, y2, x3, y3)
        abs(a2) <= eps(Float64) && continue

        xminf = min(x1, min(x2, x3))
        xmaxf = max(x1, max(x2, x3))
        yminf = min(y1, min(y2, y3))
        ymaxf = max(y1, max(y2, y3))

        ix0 = clamp(floor(Int, (xminf - xmin) / pix) + 1, 1, nx)
        ix1 = clamp(ceil(Int, (xmaxf - xmin) / pix) + 1, 1, nx)
        iy0 = clamp(floor(Int, (yminf - ymin) / pix) + 1, 1, ny)
        iy1 = clamp(ceil(Int, (ymaxf - ymin) / pix) + 1, 1, ny)

        for iy in iy0:iy1
            yc = ymin + (iy - 0.5) * pix
            for ix in ix0:ix1
                xc = xmin + (ix - 0.5) * pix
                w1 = _edge(x2, y2, x3, y3, xc, yc)
                w2 = _edge(x3, y3, x1, y1, xc, yc)
                w3 = _edge(x1, y1, x2, y2, xc, yc)

                inside = (w1 >= 0.0 && w2 >= 0.0 && w3 >= 0.0) || (w1 <= 0.0 && w2 <= 0.0 && w3 <= 0.0)
                inside || continue

                λ1 = w1 / a2
                λ2 = w2 / a2
                λ3 = w3 / a2
                z = λ1 * z1 + λ2 * z2 + λ3 * z3

                idx = (iy - 1) * nx + ix
                if z < depth[idx]
                    depth[idx] = z
                    face_at_pixel[idx] = fi
                end
            end
        end
    end

    @inbounds for idx in eachindex(face_at_pixel)
        fi = face_at_pixel[idx]
        fi == 0 && continue
        hits[fi] += 1
    end

    @inbounds for fi in eachindex(hits)
        if hits[fi] == 0
            continue
        end
        if area_ratio
            projected_area[fi] = hits[fi] * pixel_area
        else
            projected_area[fi] = face_area[fi] * max(dot(face_normal[fi], -d), 0.0)
        end
    end

    return projected_area, hits
end

function compute_first_order(
    scene::SceneGeometry,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    cfg::LightConfig;
    backend=:raster_cpu
)
    backend == :raster_cpu || error("Unsupported backend: $backend")

    vertices = GeometryBasics.decompose(PlantGeom.Point3, scene.merged_mesh)
    faces = GeometryBasics.decompose(PlantGeom.Face3, scene.merged_mesh)
    nfaces = length(faces)

    face_area = zeros(Float64, nfaces)
    face_normal = Vector{StaticArrays.SVector{3,Float64}}(undef, nfaces)
    @inbounds for i in 1:nfaces
        f = faces[i]
        area, nrm = _triangle_area_and_normal(vertices[f[1]], vertices[f[2]], vertices[f[3]])
        face_area[i] = area
        face_normal[i] = nrm
    end

    node_ids = collect(keys(scene.total_area_per_node))
    projected_area_per_node = Dict(id => 0.0 for id in node_ids)
    incident_par_power_per_node = Dict(id => 0.0 for id in node_ids)
    incident_nir_power_per_node = Dict(id => 0.0 for id in node_ids)
    hits_per_node = Dict(id => 0 for id in node_ids)

    for (k, sector) in enumerate(turtle.sectors)
        par_flux = fluxes.par[k]
        nir_flux = fluxes.nir[k]
        if par_flux == 0.0 && nir_flux == 0.0
            continue
        end

        projected_area, hits = _rasterize_direction(
            vertices,
            faces,
            face_area,
            face_normal,
            sector.direction,
            cfg.pixel_size,
            cfg.area_ratio
        )

        @inbounds for fi in 1:nfaces
            pa = projected_area[fi]
            h = hits[fi]
            if pa <= 0.0 || h == 0
                continue
            end

            nid = scene.face2node[fi]
            projected_area_per_node[nid] = get(projected_area_per_node, nid, 0.0) + pa
            incident_par_power_per_node[nid] = get(incident_par_power_per_node, nid, 0.0) + par_flux * pa
            incident_nir_power_per_node[nid] = get(incident_nir_power_per_node, nid, 0.0) + nir_flux * pa
            hits_per_node[nid] = get(hits_per_node, nid, 0) + h
        end
    end

    FirstOrderResult(
        projected_area_per_node,
        incident_par_power_per_node,
        incident_nir_power_per_node,
        hits_per_node
    )
end
