#!/usr/bin/env julia

using ArchimedLight
using CSV
using Tables

include(joinpath(dirname(@__DIR__), "test", "java_parity_harness.jl"))

_parse_int(x, d) = x === nothing ? d : something(tryparse(Int, String(x)), d)

function _fixture_by_name(name::String)
    for fx in light_parity_fixtures()
        fx.name == name && return fx
    end
    return nothing
end

function _csv_columns(rows)
    collect(String.(propertynames(first(rows))))
end

function _obs_component_rows(scene, series, cfg, meteo_rows, cols::Vector{String})
    out = Dict{String,Any}[]
    for i in eachindex(series)
        append!(
            out,
            ArchimedLight.component_values_table(
                scene,
                series[i],
                cfg;
                meteo_row=meteo_rows[i],
                step_number=i - 1,
                columns=cols,
            ).rows,
        )
    end
    sort!(out; by=r -> (Int(r["step_number"]), Int(r["item_id"]), Int(r["component_id"])))
    out
end

function _obs_summary_rows(scene, series, cfg, meteo_rows, cols::Vector{String})
    out = ArchimedLight.summary_values_table(
        scene,
        series,
        cfg;
        meteo_rows=meteo_rows,
        start_step_number=0,
        columns=cols,
    ).rows
    sort!(
        out;
        by=r -> (
            Int(get(r, "step_number", 0)),
            string(get(r, "group", "")),
            string(get(r, "type", "")),
            Int(get(r, "item_id", 0)),
        ),
    )
    out
end

function _float_or_string(row, col::String)
    if row isa AbstractDict
        v = row[col]
    else
        v = getproperty(row, Symbol(col))
    end
    return v isa Number ? Float64(v) : string(v)
end

function _compare_rows(expected_rows, observed_rows, key_cols::Vector{String}, value_cols::Vector{String}; atol=1e-3, rtol=3e-3)
    exp = Dict(Tuple(_float_or_string(r, c) for c in key_cols) => r for r in expected_rows)
    obs = Dict(Tuple(_float_or_string(r, c) for c in key_cols) => r for r in observed_rows)
    keys_all = union(keys(exp), keys(obs))
    missing_keys = [k for k in keys_all if !haskey(obs, k)]
    extra_keys = [k for k in keys_all if !haskey(exp, k)]
    mismatches = NamedTuple[]

    for k in keys_all
        haskey(exp, k) && haskey(obs, k) || continue
        er = exp[k]
        orow = obs[k]
        for c in value_cols
            ev = _float_or_string(er, c)
            ov = _float_or_string(orow, c)
            if ev isa Float64 && ov isa Float64
                if !isapprox(ov, ev; atol=atol, rtol=rtol)
                    push!(mismatches, (key=k, column=c, expected=ev, observed=ov, absdiff=abs(ov - ev)))
                end
            elseif ov != ev
                push!(mismatches, (key=k, column=c, expected=ev, observed=ov, absdiff=NaN))
            end
        end
    end

    sort!(missing_keys)
    sort!(extra_keys)
    sort!(mismatches; by=m -> isfinite(m.absdiff) ? m.absdiff : -1.0, rev=true)
    return missing_keys, extra_keys, mismatches
end

function main(args)
    isempty(args) && error("Usage: julia --project=. scripts/compare_fixture_mismatches.jl <fixture> [top_n=20]")
    fixture_name = String(args[1])
    top_n = length(args) >= 2 ? _parse_int(args[2], 20) : 20

    fx = _fixture_by_name(fixture_name)
    fx === nothing && error("Unknown fixture: $(fixture_name)")
    fx.expected_dir === nothing && error("Fixture has no expected dir: $(fixture_name)")

    cfg = ArchimedLight.read_light_config(fx.config_path)
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    selected = ArchimedLight.prepare_meteo(meteo, cfg)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)

    expected_comp = read_java_csv(joinpath(fx.expected_dir, "component_values.csv"))
    comp_cols = _csv_columns(expected_comp)
    observed_comp = _obs_component_rows(scene, series, cfg, selected.rows, comp_cols)
    comp_missing, comp_extra, comp_mismatches = _compare_rows(
        expected_comp,
        observed_comp,
        ["step_number", "item_id", "component_id"],
        [c for c in comp_cols if !(c in ("step_number", "item_id", "component_id"))],
    )

    println("fixture=", fixture_name)
    println("component missing_keys=", length(comp_missing), " extra_keys=", length(comp_extra), " mismatches=", length(comp_mismatches))
    for m in Iterators.take(comp_mismatches, top_n)
        println("  ", m)
    end

    summary_path = joinpath(fx.expected_dir, "summary.csv")
    if isfile(summary_path)
        expected_summary = read_java_csv(summary_path)
        if !isempty(expected_summary)
            summary_cols = _csv_columns(expected_summary)
            observed_summary = _obs_summary_rows(scene, series, cfg, selected.rows, summary_cols)
            sum_missing, sum_extra, sum_mismatches = _compare_rows(
                expected_summary,
                observed_summary,
                [c for c in ("step_number", "group", "type", "item_id") if c in summary_cols],
                [c for c in summary_cols if !(c in ("step_number", "group", "type", "item_id"))],
            )
            println("summary missing_keys=", length(sum_missing), " extra_keys=", length(sum_extra), " mismatches=", length(sum_mismatches))
            for m in Iterators.take(sum_mismatches, top_n)
                println("  ", m)
            end
        end
    end
end

main(ARGS)
