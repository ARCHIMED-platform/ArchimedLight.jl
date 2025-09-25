using Meshes
using LinearAlgebra: norm

struct Sector
    dir::Vec3      # unit direction (toward scene)
    weight::Float64 # relative weight, sum to 1 over sky
end

struct SkyConfig
    sectors::Vector{Sector}
    irradiance::Dict{Band,Float64} # W m^-2 per band (diffuse only)
end

function uniform_hemisphere(n::Int)
    # naive stratified sampling over phi, theta for downward hemisphere (z<0)
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
    # normalize weights (Sector is immutable)
    return [Sector(s.dir, s.weight / total) for s in sectors]
end

function SkyConfig(; count::Int=16, PAR_Wm2::Float64=400.0, NIR_Wm2::Float64=400.0)
    sectors = uniform_hemisphere(count)
    irr = Dict{Band,Float64}(PAR=>PAR_Wm2, NIR=>NIR_Wm2, TIR=>0.0)
    return SkyConfig(sectors, irr)
end

struct InterceptionConfig
    pixel_size::Float64
    scattering::Bool
    all_in_turtle::Bool
    radiation_timestep::Float64
end

InterceptionConfig(; pixel_size::Float64=0.1, scattering::Bool=false, all_in_turtle::Bool=false, radiation_timestep::Float64=0.0) = InterceptionConfig(pixel_size, scattering, all_in_turtle, radiation_timestep)
