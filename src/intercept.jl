using Meshes
import Unitful
using Unitful: oneunit
using LinearAlgebra: dot, cross

strip_quantity(x::Unitful.Quantity) = Unitful.ustrip(x)
strip_quantity(x::Real) = Float64(x)

# convert mesh/point Vec/Point (possibly Unitful) to Float64 triple (meters)
point_to_tuple(p) = (
    strip_quantity(to(p)[1]),
    strip_quantity(to(p)[2]),
    strip_quantity(to(p)[3])
)

# convert a direction Vec (possibly Unitful) to Float64 triple
vec_to_tuple(v) = (
    strip_quantity(v[1]),
    strip_quantity(v[2]),
    strip_quantity(v[3])
)

@inline dot3(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline cross3(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = (
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
)
@inline function normalize3(v::NTuple{3,Float64})
    n2 = dot3(v, v)
    if n2 == 0.0
        return (0.0, 0.0, 1.0)
    end
    invn = inv(sqrt(n2))
    return (v[1]*invn, v[2]*invn, v[3]*invn)
end

@inline function ray_box_intersect(o::NTuple{3,Float64}, d::NTuple{3,Float64}, lo::NTuple{3,Float64}, hi::NTuple{3,Float64})
    tmin = -Inf
    tmax = Inf
    for i in 1:3
        oi = o[i]; di = d[i]
        li = lo[i]; hi_ = hi[i]
        if abs(di) < 1e-12
            if oi < li || oi > hi_
                return false
            else
                continue
            end
        end
        invd = 1.0 / di
        t1 = (li - oi) * invd
        t2 = (hi_ - oi) * invd
        t_near = min(t1, t2)
        t_far = max(t1, t2)
        tmin = max(tmin, t_near)
        tmax = min(tmax, t_far)
        if tmax < tmin
            return false
        end
    end
    return tmax >= max(tmin, 0.0)
end

struct TriangleData
    v0::NTuple{3,Float64}
    edge1::NTuple{3,Float64}
    edge2::NTuple{3,Float64}
    normal::NTuple{3,Float64}
end

@inline function triangle_data(v0::NTuple{3,Float64}, v1::NTuple{3,Float64}, v2::NTuple{3,Float64})
    e1 = (v1[1] - v0[1], v1[2] - v0[2], v1[3] - v0[3])
    e2 = (v2[1] - v0[1], v2[2] - v0[2], v2[3] - v0[3])
    n = normalize3(cross3(e1, e2))
    return TriangleData(v0, e1, e2, n)
end

# Möller–Trumbore triangle intersection returning t (distance) or nothing
@inline function tri_intersect(orig::NTuple{3,Float64}, dir::NTuple{3,Float64}, tri::TriangleData)
    eps = 1e-9
    pvec = cross3(dir, tri.edge2)
    det = dot3(tri.edge1, pvec)
    if abs(det) < eps
        return nothing
    end
    invdet = 1.0 / det
    tvec = (orig[1] - tri.v0[1], orig[2] - tri.v0[2], orig[3] - tri.v0[3])
    u = dot3(tvec, pvec) * invdet
    if u < 0.0 || u > 1.0
        return nothing
    end
    qvec = cross3(tvec, tri.edge1)
    v = dot3(dir, qvec) * invdet
    if v < 0.0 || u + v > 1.0
        return nothing
    end
    t = dot3(tri.edge2, qvec) * invdet
    if t > eps
        return t
    end
    return nothing
end

struct InterceptionResult
    absorbed::Dict{String,Dict{Band,Float64}}
    incident::Dict{String,Dict{Band,Float64}}
    areas::Dict{String,Float64}
    hit_counts::Dict{String,Int}
end

Base.getindex(res::InterceptionResult, key::String) = res.absorbed[key]
Base.haskey(res::InterceptionResult, key::String) = haskey(res.absorbed, key)
Base.length(res::InterceptionResult) = length(res.absorbed)
Base.iterate(res::InterceptionResult) = iterate(res.absorbed)
Base.iterate(res::InterceptionResult, state) = iterate(res.absorbed, state)
Base.keys(res::InterceptionResult) = keys(res.absorbed)
Base.values(res::InterceptionResult) = values(res.absorbed)
Base.pairs(res::InterceptionResult) = pairs(res.absorbed)

"""
    compute_interception(scene, sky::SkyConfig, cfg::InterceptionConfig)

Computes first-order (no multiple scattering) irradiance and absorptance per component and band.
Returns an `InterceptionResult` collecting incident (`Ri`) and absorbed (`Ra`) fluxes in W m⁻².
"""
function compute_interception(scene::Scene, sky::SkyConfig, cfg::InterceptionConfig)
    lo, up = extrema(scene.bbox)
    sx = strip_quantity(to(up)[1] - to(lo)[1])
    sy = strip_quantity(to(up)[2] - to(lo)[2])
    nx = max(1, ceil(Int, sx / cfg.pixel_size))
    ny = max(1, ceil(Int, sy / cfg.pixel_size))

    components = scene.components
    ncomp = length(components)

    comp_names = Vector{String}(undef, ncomp)
    comp_triangles = Vector{Vector{TriangleData}}(undef, ncomp)
    comp_bbox_lo = Vector{NTuple{3,Float64}}(undef, ncomp)
    comp_bbox_hi = Vector{NTuple{3,Float64}}(undef, ncomp)
    comp_transparency = Vector{Float64}(undef, ncomp)
    comp_optics = Vector{OpticalProps}(undef, ncomp)
    comp_areas = Vector{Float64}(undef, ncomp)

    absorbed = Dict{String,Dict{Band,Float64}}()
    incident = Dict{String,Dict{Band,Float64}}()
    areas = Dict{String,Float64}()
    hit_counts = Dict{String,Int}()

    band_order = (PAR, NIR, TIR)
    nbands = length(band_order)
    power_accum = zeros(Float64, ncomp, nbands)

    # precompute component data
    for (idx, comp) in pairs(components)
        comp_names[idx] = comp.name

        tris = TriangleData[]
        verts = vertices(comp.mesh)
        for tri_conn in topology(comp.mesh)
            ids = indices(tri_conn)
            v0 = point_to_tuple(verts[ids[1]])
            v1 = point_to_tuple(verts[ids[2]])
            v2 = point_to_tuple(verts[ids[3]])
            push!(tris, triangle_data(v0, v1, v2))
        end
        comp_triangles[idx] = tris

        bbox = boundingbox(comp.mesh)
        comp_bbox_lo[idx] = point_to_tuple(bbox.min)
        comp_bbox_hi[idx] = point_to_tuple(bbox.max)

        props = get_optics(comp.name)
        comp_optics[idx] = props
        comp_transparency[idx] = clamp(props.transparency, 0.0, 1.0)

        absorbed[comp.name] = Dict{Band,Float64}(PAR => 0.0, NIR => 0.0, TIR => 0.0)
        incident[comp.name] = Dict{Band,Float64}(PAR => 0.0, NIR => 0.0, TIR => 0.0)
        comp_areas[idx] = strip_quantity(area(comp.mesh))
        areas[comp.name] = comp_areas[idx]
        hit_counts[comp.name] = 0
    end

    lo_tuple = (strip_quantity(to(lo)[1]), strip_quantity(to(lo)[2]), strip_quantity(to(lo)[3]))
    up_tuple = (strip_quantity(to(up)[1]), strip_quantity(to(up)[2]), strip_quantity(to(up)[3]))
    span_x = up_tuple[1] - lo_tuple[1]
    span_y = up_tuple[2] - lo_tuple[2]
    cell_w = span_x / nx
    cell_h = span_y / ny
    pixel_area = max(cell_w * cell_h, 1e-9)

    sector_dirs = [vec_to_tuple(sector.dir) for sector in sky.sectors]
    cos_weight_sum = 0.0
    for (sector, dir) in zip(sky.sectors, sector_dirs)
        if dir[3] < 0.0
            cos_weight_sum += sector.weight * abs(dir[3])
        end
    end
    cos_weight_sum = cos_weight_sum <= 0 ? 1.0 : cos_weight_sum
    sector_flux = [Tuple((sector.weight * (dir[3] < 0 ? abs(dir[3]) : 0.0) / cos_weight_sum) * get(sky.irradiance, band, 0.0)
                         for band in band_order)
                   for (sector, dir) in zip(sky.sectors, sector_dirs)]

    hits = Vector{Tuple{Float64,Int}}()
    best_ts = Vector{Float64}(undef, ncomp)
    hit_mask = BitVector(undef, ncomp)
    sector_area = zeros(Float64, ncomp)
    sector_power = zeros(Float64, ncomp, nbands)

    up_z = up_tuple[3]
    lo_z = lo_tuple[3]
    span_z = up_z - lo_z
    launch_z = up_z + (iszero(span_z) ? 1.0 : span_z)

    for (sector_idx, sector) in enumerate(sky.sectors)
        d = sector_dirs[sector_idx]
        if d[3] >= 0.0
            continue
        end
        fill!(sector_area, 0.0)
        fill!(sector_power, 0.0)
        sector_weights = sector_flux[sector_idx]

        for iy in 1:ny, ix in 1:nx
            x = lo_tuple[1] + (ix - 0.5) * cell_w
            y = lo_tuple[2] + (iy - 0.5) * cell_h
            o = (x, y, launch_z)

            fill!(best_ts, Inf)
            fill!(hit_mask, false)

            for comp_idx in 1:ncomp
                if !ray_box_intersect(o, d, comp_bbox_lo[comp_idx], comp_bbox_hi[comp_idx])
                    continue
                end
                best_t = Inf
                for tri in comp_triangles[comp_idx]
                    t = tri_intersect(o, d, tri)
                    if t !== nothing && t < best_t
                        best_t = t
                    end
                end
                if best_t < Inf
                    best_ts[comp_idx] = best_t
                    hit_mask[comp_idx] = true
                end
            end

            if !any(hit_mask)
                continue
            end

            empty!(hits)
            for idx in 1:ncomp
                if hit_mask[idx]
                    push!(hits, (best_ts[idx], idx))
                end
            end
            sort!(hits, by = x -> x[1])

            incident_area = pixel_area
            for (t, comp_idx) in hits
                incident_area <= 1e-12 && break

                transparency = comp_transparency[comp_idx]
                passing = incident_area * transparency
                intercepted_area = incident_area - passing
                if intercepted_area <= 0.0
                    incident_area = passing
                    continue
                end

                sector_area[comp_idx] += intercepted_area
                for (band_idx, flux) in enumerate(sector_weights)
                    sector_power[comp_idx, band_idx] += intercepted_area * flux
                end

                hit_counts[comp_names[comp_idx]] += 1

                incident_area = passing
            end
        end

        for comp_idx in 1:ncomp
            area_sum = sector_area[comp_idx]
            if area_sum <= 0.0
                continue
            end
            ratio = comp_areas[comp_idx] / area_sum
            for band_idx in 1:nbands
                power_accum[comp_idx, band_idx] += sector_power[comp_idx, band_idx] * ratio
            end
        end
    end

    # consolidate incident/absorbed per band
    for (idx, name) in enumerate(comp_names)
        area_val = areas[name]
        props = comp_optics[idx]
        area_val <= 0 && continue
        for (band_idx, band) in enumerate(band_order)
            Ri = power_accum[idx, band_idx] / area_val
            incident[name][band] = Ri
            α = absorptance(props, band)
            absorbed[name][band] = Ri * α
        end
    end

    return InterceptionResult(absorbed, incident, areas, hit_counts)
end
