module SkySets

export sky_sectors_for

using ..ArchimedLight: Sector
using LinearAlgebra: norm
using Meshes

"""
    sky_sectors_for(n)

Return canonical Archimed “turtle” sky sectors for diffuse sky:
- 6: exact seed set from Archimed common (downward directions).
- others: fallback to uniform hemisphere (placeholder until tables are added).
"""
function sky_sectors_for(n::Int)
    if n == 6
        elevation_deg = (90.00, 26.57, 26.57, 26.57, 26.57, 26.57)
        azimuth_deg   = (0.0,   180.0, 252.0, 324.0, 36.0,  108.0)
        dirs = Vector{Vec{3}}()
        for k in 1:6
            el = deg2rad(elevation_deg[k])
            az = deg2rad(azimuth_deg[k] + 180.0)
            x = cos(el)
            y = cos(el)
            x *= sin(az)
            y *= cos(az)
            z = sin(el)
            # Java turtle is upward (ground->sky). For ray casting downward, flip sign.
            d = Vec(-x, -y, -z)
            push!(dirs, Vec(d / norm(d)))
        end
        w = 1.0 / length(dirs)
        return [Sector(d, w) for d in dirs]
    else
        m = max(1, round(Int, sqrt(n)))
        sectors = Sector[]
        total = 0.0
        for i in 1:m, j in 1:m
            u = (i-0.5)/m
            v = (j-0.5)/m
            θ = acos(1 - u)
            ϕ = 2π * v
            d = Vec(sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), -abs(cos(θ)))
            if d[3] < zero(d[3])
                push!(sectors, Sector(Vec(d / norm(d)), 1.0))
                total += 1.0
            end
        end
        return [Sector(s.dir, s.weight / total) for s in sectors]
    end
end

deg2rad(x) = x * (π / 180.0)

end # module
