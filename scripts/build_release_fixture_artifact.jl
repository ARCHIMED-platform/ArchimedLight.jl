#!/usr/bin/env julia

using Pkg.Artifacts
using SHA
using Dates

function usage()
    println(
        """
        Build an artifact tarball for heavy release-only fixture regression.

        Usage:
          julia --project=. scripts/build_release_fixture_artifact.jl [options]

        Options:
          --artifact-name <name>   Artifact name (default: archimedlight-release-fixtures)
          --test-root <path>       Source test root (default: ./test)
          --artifacts-toml <path>  Artifacts.toml to update (default: ./Artifacts.toml)
          --tarball <path>         Output tar.gz path (default: ./release-fixtures.tar.gz)
          --url <url>              Download URL used to bind artifact in Artifacts.toml
          --no-prune-output        Keep fixture output/ directories in bundle
          -h, --help               Show this help
        """,
    )
end

function parse_args(args)
    repo_root = normpath(joinpath(dirname(@__DIR__)))
    opt = Dict{String,Any}(
        "artifact_name" => "archimedlight-release-fixtures",
        "test_root" => joinpath(repo_root, "test"),
        "artifacts_toml" => joinpath(repo_root, "Artifacts.toml"),
        "tarball" => joinpath(repo_root, "release-fixtures.tar.gz"),
        "url" => "",
        "prune_output" => true,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            usage()
            exit(0)
        elseif a == "--artifact-name"
            i += 1
            opt["artifact_name"] = args[i]
        elseif a == "--test-root"
            i += 1
            opt["test_root"] = args[i]
        elseif a == "--artifacts-toml"
            i += 1
            opt["artifacts_toml"] = args[i]
        elseif a == "--tarball"
            i += 1
            opt["tarball"] = args[i]
        elseif a == "--url"
            i += 1
            opt["url"] = args[i]
        elseif a == "--no-prune-output"
            opt["prune_output"] = false
        else
            error("Unknown argument: $(repr(a))")
        end
        i += 1
    end
    return opt
end

function copy_tree_pruned!(src::AbstractString, dst::AbstractString; prune_output::Bool)
    mkpath(dst)
    for name in readdir(src)
        src_path = joinpath(src, name)
        dst_path = joinpath(dst, name)
        if isdir(src_path)
            if prune_output && name == "output"
                continue
            end
            copy_tree_pruned!(src_path, dst_path; prune_output=prune_output)
        else
            cp(src_path, dst_path; force=true)
        end
    end
end

function copy_required_files!(dst::AbstractString, test_root::AbstractString; prune_output::Bool)
    required = [
        "fixtures",
        "references",
        "reference_images",
        "fixtures_manifest.toml",
        "julia_fixture_harness.jl",
        "test_fixtures.jl",
    ]
    for rel in required
        src = joinpath(test_root, rel)
        ispath(src) || error("Missing required path in test root: $(src)")
        dst_path = joinpath(dst, rel)
        if rel == "fixtures"
            copy_tree_pruned!(src, dst_path; prune_output=prune_output)
        else
            cp(src, dst_path; force=true)
        end
    end
end

function write_release_runtests!(dst::AbstractString)
    path = joinpath(dst, "runtests_release.jl")
    script = """
    using Test
    using ArchimedLight

    old = get(ENV, "ARCHIMEDLIGHT_TEST_PROFILE", "")
    ENV["ARCHIMEDLIGHT_TEST_PROFILE"] = "fixtures"
    include(joinpath(@__DIR__, "julia_fixture_harness.jl"))
    include(joinpath(@__DIR__, "test_fixtures.jl"))
    if isempty(old)
        delete!(ENV, "ARCHIMEDLIGHT_TEST_PROFILE")
    else
        ENV["ARCHIMEDLIGHT_TEST_PROFILE"] = old
    end
    """
    open(path, "w") do io
        write(io, script)
    end
    return path
end

function main(args)
    opt = parse_args(args)
    test_root = normpath(String(opt["test_root"]))
    artifacts_toml = normpath(String(opt["artifacts_toml"]))
    tarball = normpath(String(opt["tarball"]))
    artifact_name = String(opt["artifact_name"])
    url = strip(String(opt["url"]))
    prune_output = Bool(opt["prune_output"])

    isdir(test_root) || error("test root does not exist: $(test_root)")
    mkpath(dirname(tarball))

    hash = create_artifact() do artdir
        copy_required_files!(artdir, test_root; prune_output=prune_output)
        write_release_runtests!(artdir)
    end

    archive_artifact(hash, tarball)
    sha = bytes2hex(open(sha256, tarball))

    println("Artifact name: ", artifact_name)
    println("git-tree-sha1: ", hash)
    println("tarball: ", tarball)
    println("sha256: ", sha)
    println("timestamp_utc: ", Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))

    if !isempty(url)
        bind_artifact!(
            artifacts_toml,
            artifact_name,
            hash;
            download_info=[(url, sha)],
            force=true,
            lazy=true,
        )
        println("Updated Artifacts.toml: ", artifacts_toml)
    else
        println("No URL provided; Artifacts.toml was not updated.")
        println("Once uploaded, rerun with --url <download-url> to bind the artifact.")
    end
end

main(ARGS)
