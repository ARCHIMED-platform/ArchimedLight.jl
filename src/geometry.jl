using Meshes
using LinearAlgebra: norm

const Point3 = Point
const Vec3 = Vec{3}
const Ray3 = Ray
const Triangle3 = Triangle
const SimpleMesh3 = SimpleMesh

"""
    first_hit(ray::Ray, mesh::SimpleMesh)

Return (hit::Bool, tmin) for the nearest triangle intersection, with tmin
having the same length units as the mesh/ray coordinates.
"""
function first_hit(ray::Ray3, mesh::SimpleMesh3)
    v = vertices(mesh)
    tmin = measure(ray) # typemax length in appropriate units
    hit_found = false
    for tri_conn in topology(mesh)
        idx = indices(tri_conn)
        tri = Triangle3(v[idx[1]], v[idx[2]], v[idx[3]])
        I = intersection(ray, tri)
        if type(I) != NotIntersecting
            p = get(I) # intersection point
            t = norm(p - ray(0))
            if t > zero(t) && t < tmin
                tmin = t
                hit_found = true
            end
        end
    end
    return hit_found, tmin
end
