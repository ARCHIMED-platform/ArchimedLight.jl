#!/usr/bin/env julia

using Dates
using Downloads
using Pkg.Artifacts
using SHA

const DEFAULT_SOURCE_URL = "https://zenodo.org/records/21134328/files/benchmark_scenes.zip?download=1"

function usage()
    println(
        """
        Build a Julia artifact tarball for benchmark OPS scenes.

        Julia artifacts unpack tar archives, not zip files. This script downloads or reuses the
        Zenodo benchmark_scenes.zip source, removes macOS metadata files, archives the cleaned
        dataset as a tar.gz artifact, and binds the artifact hash in Artifacts.toml.

        Usage:
          julia --project=. scripts/build_benchmark_scenes_artifact.jl [options]

        Options:
          --artifact-name <name>   Artifact name (default: archimedlight-benchmark-scenes)
          --artifacts-toml <path>  Artifacts.toml to update (default: ./Artifacts.toml)
          --source-url <url>       Source zip URL (default: Zenodo benchmark_scenes.zip)
          --zip <path>             Reuse an existing source zip instead of downloading
          --tarball <path>         Output tar.gz path (default: ./benchmark-scenes.tar.gz)
          --url <url>              Download URL for the generated tar.gz, used in Artifacts.toml
          --keep-metadata          Keep .DS_Store, __MACOSX, and AppleDouble files
          --no-bind                Do not update Artifacts.toml
          -h, --help               Show this help
        """,
    )
end

function parse_args(args)
    repo_root = normpath(joinpath(dirname(@__DIR__)))
    opt = Dict{String,Any}(
        "artifact_name" => "archimedlight-benchmark-scenes",
        "artifacts_toml" => joinpath(repo_root, "Artifacts.toml"),
        "source_url" => DEFAULT_SOURCE_URL,
        "zip" => "",
        "tarball" => joinpath(repo_root, "benchmark-scenes.tar.gz"),
        "url" => "",
        "keep_metadata" => false,
        "bind" => true,
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
        elseif a == "--artifacts-toml"
            i += 1
            opt["artifacts_toml"] = args[i]
        elseif a == "--source-url"
            i += 1
            opt["source_url"] = args[i]
        elseif a == "--zip"
            i += 1
            opt["zip"] = args[i]
        elseif a == "--tarball"
            i += 1
            opt["tarball"] = args[i]
        elseif a == "--url"
            i += 1
            opt["url"] = args[i]
        elseif a == "--keep-metadata"
            opt["keep_metadata"] = true
        elseif a == "--no-bind"
            opt["bind"] = false
        else
            error("Unknown argument: $(repr(a))")
        end
        i += 1
    end
    return opt
end

function _should_skip(name::AbstractString; keep_metadata::Bool)
    keep_metadata && return false
    return name == ".DS_Store" || name == "__MACOSX" || startswith(name, "._")
end

function _copy_tree_filtered!(src::AbstractString, dst::AbstractString; keep_metadata::Bool)
    mkpath(dst)
    for name in readdir(src)
        _should_skip(name; keep_metadata=keep_metadata) && continue
        src_path = joinpath(src, name)
        dst_path = joinpath(dst, name)
        if isdir(src_path)
            _copy_tree_filtered!(src_path, dst_path; keep_metadata=keep_metadata)
        else
            cp(src_path, dst_path; force=true)
        end
    end
end

function _source_zip!(opt::Dict{String,Any}, workdir::AbstractString)
    existing = strip(String(opt["zip"]))
    if !isempty(existing)
        path = normpath(existing)
        isfile(path) || error("zip file does not exist: $(path)")
        return path
    end

    source_url = String(opt["source_url"])
    zip_path = joinpath(workdir, "benchmark_scenes.zip")
    println("Downloading source zip: ", source_url)
    Downloads.download(source_url, zip_path)
    return zip_path
end

function _extract_zip(zip_path::AbstractString, workdir::AbstractString)
    unzip = Sys.which("unzip")
    unzip === nothing && error("Could not find `unzip` on PATH; pass --zip after extracting manually or install unzip.")
    extract_root = joinpath(workdir, "extracted")
    mkpath(extract_root)
    run(Cmd([unzip, "-q", zip_path, "-d", extract_root]))
    source_dir = joinpath(extract_root, "benchmark_scenes")
    isdir(source_dir) || error("zip did not contain expected benchmark_scenes/ directory: $(zip_path)")
    return source_dir
end

function main(args)
    opt = parse_args(args)
    artifact_name = String(opt["artifact_name"])
    artifacts_toml = normpath(String(opt["artifacts_toml"]))
    tarball = normpath(String(opt["tarball"]))
    artifact_url = strip(String(opt["url"]))
    keep_metadata = Bool(opt["keep_metadata"])
    bind = Bool(opt["bind"])

    mkpath(dirname(tarball))

    workdir = mktempdir()
    zip_path = _source_zip!(opt, workdir)
    source_zip_sha = bytes2hex(open(sha256, zip_path))
    source_dir = _extract_zip(zip_path, workdir)

    hash = create_artifact() do artdir
        _copy_tree_filtered!(
            source_dir,
            joinpath(artdir, "benchmark_scenes");
            keep_metadata=keep_metadata,
        )
    end

    println("Archiving artifact to tarball...")
    archive_artifact(hash, tarball)
    tarball_sha = bytes2hex(open(sha256, tarball))

    if bind
        if isempty(artifact_url)
            bind_artifact!(artifacts_toml, artifact_name, hash; force=true, lazy=true)
        else
            bind_artifact!(
                artifacts_toml,
                artifact_name,
                hash;
                download_info=[(artifact_url, tarball_sha)],
                force=true,
                lazy=true,
            )
        end
    end

    println("Artifact name: ", artifact_name)
    println("git-tree-sha1: ", hash)
    println("source zip: ", zip_path)
    println("source zip sha256: ", source_zip_sha)
    println("tarball: ", tarball)
    println("tarball sha256: ", tarball_sha)
    println("timestamp_utc: ", Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    if bind
        println("Updated Artifacts.toml: ", artifacts_toml)
        isempty(artifact_url) && println("No tarball URL provided; Artifacts.toml was bound without download info.")
    end
end

main(ARGS)
