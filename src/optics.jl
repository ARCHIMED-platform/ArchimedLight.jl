@enum Band PAR NIR TIR

struct OpticalProps
    # per-band reflectance and transmittance; absorptance = 1 - ρ - τ
    rho::Dict{Band,Float64}
    tau::Dict{Band,Float64}
    transparency::Float64 # 0..1 additional straight-through factor
    scattering::Dict{Band,Float64}   # per-band scattering factor [0..1]
end

function OpticalProps(; transparency=0.0)
    rho = Dict{Band,Float64}(PAR => 0.1, NIR => 0.4, TIR => 0.0)
    tau = Dict{Band,Float64}(PAR => 0.05, NIR => 0.4, TIR => 0.0)
    scat = Dict{Band,Float64}(PAR => 0.0, NIR => 0.0, TIR => 0.0)
    return OpticalProps(rho, tau, transparency, scat)
end

absorptance(op::OpticalProps, b::Band) = max(0.0, 1.0 - get(op.rho, b, 0.0) - get(op.tau, b, 0.0))

# optics registry per component name
const _OPTICS = Dict{String,OpticalProps}()

set_optics!(component_name::String, props::OpticalProps) = (_OPTICS[component_name] = props)
get_optics(name::String) = get(_OPTICS, name, OpticalProps())
