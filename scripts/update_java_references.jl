#!/usr/bin/env julia

"""Run ARCHIMED Java reference simulations for ArchimedLight tests.

This script executes the legacy Java implementation for each integration
fixture under `test/java` and copies the generated outputs under
`test/java_reference`. Tests can then rely on these stored artefacts instead of
spawning the Java simulator on every run.

Usage:
    julia --project scripts/update_java_references.jl

The script expects the ARCHIMED reference jar to be discoverable via
`JavaTestUtils.locate_java_jar()` (e.g. `ARCHIMED_JAR` env var or the default
build tree).
"""

using Printf
using Dates

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const TEST_ROOT = joinpath(REPO_ROOT, "test")
const CASE_ROOT = joinpath(TEST_ROOT, "java")
const REF_ROOT = joinpath(TEST_ROOT, "java_reference")

include(joinpath(TEST_ROOT, "java_utils.jl"))
using .JavaTestUtils

function write_metadata(path::AbstractString, meta::Dict{String,Any})
    order = ["case", "config", "extra_args", "success", "runs", "generated_at"]
    open(path, "w") do io
        for key in order
            value = get(meta, key, nothing)
            value === nothing && continue
            print(io, key, " = ")
            if isa(value, Bool)
                println(io, value ? "true" : "false")
            elseif isa(value, AbstractString)
                println(io, '"', replace(value, "\"" => "\\\""), '"')
            elseif isa(value, AbstractVector)
                items = join([string('"', replace(String(v), "\"" => "\\\""), '"') for v in value], ", ")
                println(io, "[", items, "]")
            else
                println(io, string(value))
            end
        end
    end
end

function copy_reference(case_name::String, config::String; extra_args::Vector{String}=String[])
    jar = locate_java_jar()
    jar === nothing && error("ARCHIMED jar not found. Set ARCHIMED_JAR or build the Java project.")

    case_dir = joinpath(CASE_ROOT, case_name)
    isdir(case_dir) || error("Case directory not found: $case_dir")

    tag = JavaTestUtils.reference_tag(config, extra_args)
    dest_dir = joinpath(REF_ROOT, case_name, tag)
    if isdir(dest_dir)
        rm(dest_dir; recursive=true, force=true)
    end
    mkpath(dest_dir)

    result = JavaTestUtils.with_java_case(case_dir, config; jar=jar, extra_args=extra_args) do success, log_text, run_dirs
        # Copy run directories (if any)
        copied = String[]
        for run_dir in run_dirs
            name = basename(run_dir)
            target = joinpath(dest_dir, name)
            cp(run_dir, target; force=true)
            push!(copied, name)
        end

        # Write log text and metadata
        write(joinpath(dest_dir, "log.txt"), log_text)
        meta = Dict{String,Any}(
            "case" => case_name,
            "config" => config,
            "extra_args" => extra_args,
            "success" => success,
            "runs" => copied,
            "generated_at" => string(Dates.now())
        )
        write_metadata(joinpath(dest_dir, "metadata.toml"), meta)

        return success
    end
    return result
end

function main()
    mkpath(REF_ROOT)

    jobs = Dict(
        "test-absorb" => [
            (config="config.yml", extra_args=String[], expect_success=true)
        ],
        "test-pixelsize" => [
            (config="config.yml", extra_args=String[], expect_success=true),
            (config="config.yml", extra_args=["--prop", "pixel_size=49.9999"], expect_success=true),
            (config="config.yml", extra_args=["--prop", "pixel_size=50.0001"], expect_success=false)
        ],
        "test-ops-output" => [
            (config="config1.yml", extra_args=String[], expect_success=true),
            (config="config2.yml", extra_args=String[], expect_success=true),
            (config="config3.yml", extra_args=String[], expect_success=true),
            (config="config4.yml", extra_args=String[], expect_success=true),
            (config="config5.yml", extra_args=String[], expect_success=false),
            (config="config6.yml", extra_args=String[], expect_success=false)
        ],
        "test-lightsource1" => [
            (config="config1.yml", extra_args=["-d", "--no-upperhit-pixtable"], expect_success=true),
            (config="config2.yml", extra_args=["-d"], expect_success=true)
        ]
    )

    failures = String[]
    for (case, specs) in sort(collect(jobs))
        println("[case] ", case)
        for spec in specs
            config = spec.config
            extra_args = spec.extra_args
            expect_success = spec.expect_success
            tag = JavaTestUtils.reference_tag(config, extra_args)
            print(@sprintf("  %-30s", tag))
            try
                success = copy_reference(case, config; extra_args=extra_args)
                ok = (success == expect_success)
                if ok
                    println(expect_success ? " ✓" : " ✓ (expected failure)")
                else
                    msg = expect_success ? "failed run" : "unexpected success"
                    println(" ✗ (", msg, ")")
                    push!(failures, "$case/$tag")
                end
            catch err
                println(" ✗ (error: $(err))")
                push!(failures, "$case/$tag")
            end
        end
    end

    if !isempty(failures)
        println()
        println("Failures:")
        for f in failures
            println("  - ", f)
        end
        exit(1)
    end
    println("\nAll reference runs completed.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
