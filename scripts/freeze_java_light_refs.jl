#!/usr/bin/env julia

using Dates
using ArchimedLight

const _DEFAULT_JAR = joinpath(dirname(@__DIR__), "example", "archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar")
const _TESTS_ROOT = joinpath(dirname(@__DIR__), "java_implementation", "archimed-lib-2018", "tests")

const _FIXTURES = Dict(
    "test-compare-cafeier1" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-cafeier2" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-cafeier3" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-cafeier4" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-cafeier5" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-simpleplant" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-two_coffee" => (config_rel="archimed2018/config.yml", java_args=String[]),
    "test-compare-timestep" => (config_rel="archimed2018/config.yml", java_args=String[]),
)

const _ALL_FIXTURES = sort!(collect(keys(_FIXTURES)))
const _FROZEN_FILES = ["component_values.csv", "scene_values.csv", "summary.csv", "log-sun-position.csv", "log-iteration-scat-par.csv", "meteo.csv"]

function _usage()
    println(
        """
        Usage:
          julia --project=. scripts/freeze_java_light_refs.jl [--jar <path>] [--fixtures <name1,name2,...>] [--force] [--keep-output]

        Defaults:
          --jar      $(_DEFAULT_JAR)
          --fixtures $(join(_ALL_FIXTURES, ","))

        Notes:
          - This script writes frozen Java references to <fixture>/expected/.
          - Existing expected folders are not overwritten unless --force is provided.
        """,
    )
end

function _parse_args(args::Vector{String})
    jar = _DEFAULT_JAR
    fixtures = copy(_ALL_FIXTURES)
    force = false
    keep_output = false

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--jar"
            i += 1
            i <= length(args) || error("missing value for --jar")
            jar = abspath(args[i])
        elseif a == "--fixtures"
            i += 1
            i <= length(args) || error("missing value for --fixtures")
            fixtures = filter(!isempty, strip.(split(args[i], ',')))
        elseif a == "--force"
            force = true
        elseif a == "--keep-output"
            keep_output = true
        elseif a in ("-h", "--help")
            _usage()
            exit(0)
        else
            error("unknown argument: $(a)")
        end
        i += 1
    end

    return (jar=jar, fixtures=fixtures, force=force, keep_output=keep_output)
end

function _resolve_output_run_dir(output_base::String)
    if isfile(joinpath(output_base, "component_values.csv"))
        return output_base
    end

    dirs = String[]
    for name in readdir(output_base)
        p = joinpath(output_base, name)
        isdir(p) || continue
        isfile(joinpath(p, "component_values.csv")) || continue
        push!(dirs, p)
    end
    isempty(dirs) && error("no simulation output found under $(output_base)")
    sort!(dirs)
    return dirs[end]
end

function _copy_frozen_files!(src_dir::String, expected_dir::String)
    copied = String[]
    for f in _FROZEN_FILES
        src = joinpath(src_dir, f)
        isfile(src) || continue
        dst = joinpath(expected_dir, f)
        cp(src, dst; force=true)
        push!(copied, f)
    end
    isempty(copied) && error("no known CSV outputs found in $(src_dir)")
    copied
end

function _write_manifest(expected_dir::String, fixture::String, cfg_path::String, jar::String, copied::Vector{String})
    path = joinpath(expected_dir, "FREEZE_INFO.txt")
    ts = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS")
    open(path, "w") do io
        println(io, "fixture=$(fixture)")
        println(io, "config=$(cfg_path)")
        println(io, "jar=$(jar)")
        println(io, "timestamp_utc=$(ts)")
        println(io, "files=$(join(copied, ','))")
    end
end

function _freeze_one_fixture(fixture_name::AbstractString, jar::String; force::Bool=false, keep_output::Bool=false)
    fixture = String(fixture_name)
    haskey(_FIXTURES, fixture) || error("unknown fixture $(repr(fixture)); available: $(join(_ALL_FIXTURES, ','))")
    spec = _FIXTURES[fixture]
    fixture_dir = joinpath(_TESTS_ROOT, fixture)
    cfg_path = normpath(joinpath(fixture_dir, spec.config_rel))
    isfile(cfg_path) || error("missing config file: $(cfg_path)")
    cfg_dir = dirname(cfg_path)

    cfg = ArchimedLight.read_light_config(cfg_path)
    out_base = ArchimedLight.output_directory(cfg)
    expected_dir = joinpath(fixture_dir, "expected")

    if isdir(expected_dir)
        force || error("expected directory already exists for $(fixture): $(expected_dir). Use --force to overwrite.")
        rm(expected_dir; recursive=true, force=true)
    end
    mkpath(expected_dir)

    isdir(out_base) && rm(out_base; recursive=true, force=true)

    cmd = Cmd(`java -jar $jar $(basename(cfg_path)) $(spec.java_args...)`; dir=cfg_dir)
    println(">>> running $(fixture)")
    run(cmd)

    out_run_dir = _resolve_output_run_dir(out_base)
    copied = _copy_frozen_files!(out_run_dir, expected_dir)
    _write_manifest(expected_dir, fixture, cfg_path, jar, copied)

    if !keep_output && isdir(out_base)
        rm(out_base; recursive=true, force=true)
    end

    println("    frozen $(length(copied)) files -> $(expected_dir)")
    nothing
end

function main(args::Vector{String})
    opts = _parse_args(args)
    isfile(opts.jar) || error("Java jar not found: $(opts.jar)")

    fixtures = opts.fixtures
    isempty(fixtures) && error("empty fixture list")

    println("Java jar: $(opts.jar)")
    println("Fixtures: $(join(fixtures, ", "))")
    println("Force overwrite: $(opts.force)")
    println("Keep Java output folders: $(opts.keep_output)")
    println()

    t0 = time()
    for fx in fixtures
        _freeze_one_fixture(fx, opts.jar; force=opts.force, keep_output=opts.keep_output)
    end
    elapsed = round(time() - t0; digits=3)
    println()
    println("Done in $(elapsed)s")
end

main(ARGS)
