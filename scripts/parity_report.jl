#!/usr/bin/env julia

using ArchimedLight

include(joinpath(dirname(@__DIR__), "test", "java_parity_harness.jl"))

function _reldiff(a, b)
    (a === nothing || b === nothing) && return nothing
    da = Float64(a)
    db = Float64(b)
    den = max(abs(db), eps(Float64))
    abs(da - db) / den
end

function _fmt(x; digits=4)
    x === nothing && return "NA"
    x isa Integer && return string(x)
    x isa AbstractFloat && return string(round(x, digits=digits))
    return string(x)
end

function _print_row(name, metric, julia_val, expected_val)
    rd = _reldiff(julia_val, expected_val)
    println(rpad(name, 28), "  ", rpad(metric, 20), "  julia=", lpad(_fmt(julia_val), 12), "  expected=", lpad(_fmt(expected_val), 12), "  rel_err=", _fmt(rd))
end

function main()
    fixtures = Dict(f.name => f for f in light_parity_fixtures())

    println("ArchimedLight parity report (baseline)")
    println("====================================")

    for (name, all_turtle) in [("test-hitcount2", false), ("test-hitcount2", true), ("test-weighted-sun", nothing), ("test-scattering-one-plate", nothing)]
        fx = fixtures[name]
        rep = fixture_parity_report(fx; all_in_turtle=all_turtle)
        suffix = all_turtle === nothing ? "" : (all_turtle ? " [all_in_turtle=true]" : " [all_in_turtle=false]")
        println()
        println(name * suffix)
        println(repeat("-", length(name * suffix)))
        _print_row(name, "total_hits", rep.snapshot.total_hits, rep.expected_hitcount_total)
        _print_row(name, "sun_azimuth_deg", rep.snapshot.sun_azimuth, rep.expected_sun_azimuth)
        _print_row(name, "sun_elevation_deg", rep.snapshot.sun_elevation, rep.expected_sun_elevation)
        _print_row(name, "scattering_sum", rep.snapshot.total_par - rep.snapshot.total_par0, rep.expected_scattering_total)
    end
end

main()

