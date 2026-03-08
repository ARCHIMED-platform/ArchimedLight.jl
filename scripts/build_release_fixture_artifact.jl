#!/usr/bin/env julia

using Pkg.Artifacts
using SHA
using Dates
using TOML

function usage()
    println(
        """
        Build an artifact tarball for heavy release-only fixture regression (data only).

        Usage:
          julia --project=. scripts/build_release_fixture_artifact.jl [options]

        Options:
          --artifact-name <name>   Artifact name (default: archimedlight-release-fixtures)
          --test-root <path>       Source release dataset root (default: ./test)
          --artifacts-toml <path>  Artifacts.toml to update (default: ./Artifacts.toml)
          --tarball <path>         Output tar.gz path (default: ./release-fixtures.tar.gz)
          --url <url>              Download URL used to bind artifact in Artifacts.toml
          --no-prune-output        Keep fixture output/ directories in bundle
          --no-prune-unused        Keep disabled/unreferenced fixtures, references, and images
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
        "prune_unused" => true,
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
        elseif a == "--no-prune-unused"
            opt["prune_unused"] = false
        else
            error("Unknown argument: $(repr(a))")
        end
        i += 1
    end
    return opt
end

function _copy_tree_filtered!(
    src::AbstractString,
    dst::AbstractString;
    prune_output::Bool,
    skip_top_dirs::Set{String}=Set{String}(),
)
    mkpath(dst)
    for name in readdir(src)
        name == ".DS_Store" && continue
        src_path = joinpath(src, name)
        dst_path = joinpath(dst, name)
        if isdir(src_path)
            if prune_output && name == "output"
                continue
            end
            if name in skip_top_dirs
                continue
            end
            _copy_tree_filtered!(src_path, dst_path; prune_output=prune_output)
        else
            cp(src_path, dst_path; force=true)
        end
    end
end

function _fixture_dir_from_config(config_rel::AbstractString)
    rel = normpath(config_rel)
    parts = split(rel, '/')
    length(parts) >= 2 || error("Invalid fixture config path in manifest: $(config_rel)")
    parts[1] == "fixtures" || error("Fixture config path must start with fixtures/: $(config_rel)")
    return parts[2]
end

function _first_fixture_subdir_from_config(config_rel::AbstractString, fixture_dir::AbstractString)
    root = normpath(joinpath("fixtures", fixture_dir))
    rel = normpath(relpath(normpath(config_rel), root))
    startswith(rel, "..") && return ""
    parts = split(rel, '/')
    length(parts) <= 1 && return ""
    return parts[1]
end

function _copy_fixture!(
    src_root::AbstractString,
    dst_root::AbstractString,
    fixture_dir::AbstractString,
    config_rel::AbstractString;
    prune_output::Bool,
)
    src = joinpath(src_root, "fixtures", fixture_dir)
    isdir(src) || error("Missing fixture directory: $(src)")
    dst = joinpath(dst_root, "fixtures", fixture_dir)
    mkpath(dst)

    preferred = _first_fixture_subdir_from_config(config_rel, fixture_dir)
    has2016 = isdir(joinpath(src, "archimed2016"))
    has2018 = isdir(joinpath(src, "archimed2018"))

    if has2016 && has2018 && preferred in ("archimed2016", "archimed2018")
        skip = preferred == "archimed2018" ? Set(["archimed2016"]) : Set(["archimed2018"])
        _copy_tree_filtered!(src, dst; prune_output=prune_output, skip_top_dirs=skip)
    else
        _copy_tree_filtered!(src, dst; prune_output=prune_output)
    end
end

function _copy_selected_references!(src_root::AbstractString, dst_root::AbstractString, fixture_ids::Vector{String})
    src_ref_root = joinpath(src_root, "references", "fixtures")
    dst_ref_root = joinpath(dst_root, "references", "fixtures")
    src_img_root = joinpath(src_root, "reference_images", "fixtures")
    dst_img_root = joinpath(dst_root, "reference_images", "fixtures")
    mkpath(dst_ref_root)
    mkpath(dst_img_root)

    for id in fixture_ids
        src_ref = joinpath(src_ref_root, id)
        if isdir(src_ref)
            _copy_tree_filtered!(src_ref, joinpath(dst_ref_root, id); prune_output=false)
        end
        src_img = joinpath(src_img_root, "$(id)_montage.png")
        if isfile(src_img)
            cp(src_img, joinpath(dst_img_root, "$(id)_montage.png"); force=true)
        end
    end
end

function _copy_all_references!(src_root::AbstractString, dst_root::AbstractString)
    for rel in ("references", "reference_images")
        src = joinpath(src_root, rel)
        isdir(src) || continue
        _copy_tree_filtered!(src, joinpath(dst_root, rel); prune_output=false)
    end
end

function _load_manifest(src_root::AbstractString)
    manifest_path = joinpath(src_root, "fixtures_manifest.toml")
    isfile(manifest_path) || error("Missing fixtures manifest: $(manifest_path)")
    return TOML.parsefile(manifest_path)
end

function _enabled_fixture_entries(manifest::Dict{String,Any})
    fixtures = get(manifest, "fixtures", Any[])
    return Any[f for f in fixtures if get(f, "enabled", true)]
end

function _write_manifest!(dst_root::AbstractString, manifest::Dict{String,Any}, fixtures_entries)
    out = Dict{String,Any}()
    if haskey(manifest, "defaults")
        out["defaults"] = manifest["defaults"]
    end
    out["fixtures"] = fixtures_entries
    open(joinpath(dst_root, "fixtures_manifest.toml"), "w") do io
        TOML.print(io, out; sorted=true)
    end
end

function copy_release_dataset!(dst::AbstractString, src_root::AbstractString; prune_output::Bool, prune_unused::Bool)
    manifest = _load_manifest(src_root)
    all_entries = get(manifest, "fixtures", Any[])
    enabled_entries = _enabled_fixture_entries(manifest)
    kept_entries = prune_unused ? enabled_entries : all_entries

    fixture_ids = String[String(item["id"]) for item in kept_entries if haskey(item, "id")]
    unique_ids = unique(fixture_ids)

    mkpath(joinpath(dst, "fixtures"))
    copied_fixture_dirs = Set{String}()
    for item in kept_entries
        haskey(item, "id") || continue
        haskey(item, "config") || continue
        config_rel = String(item["config"])
        fixture_dir = _fixture_dir_from_config(config_rel)
        fixture_dir in copied_fixture_dirs && continue
        _copy_fixture!(src_root, dst, fixture_dir, config_rel; prune_output=prune_output)
        push!(copied_fixture_dirs, fixture_dir)
    end

    if prune_unused
        _copy_selected_references!(src_root, dst, unique_ids)
    else
        _copy_all_references!(src_root, dst)
    end

    _write_manifest!(dst, manifest, kept_entries)

    println("Fixtures in source manifest: ", length(all_entries))
    println("Fixtures kept in artifact: ", length(kept_entries))
    println("Unique fixture ids kept: ", length(unique_ids))
    println("Unique fixture directories kept: ", length(copied_fixture_dirs))

    if prune_unused
        src_fixture_root = joinpath(src_root, "fixtures")
        src_fixture_dirs = isdir(src_fixture_root) ? Set(filter(d -> isdir(joinpath(src_fixture_root, d)), readdir(src_fixture_root))) : Set{String}()
        removed_fixture_dirs = sort!(collect(setdiff(src_fixture_dirs, copied_fixture_dirs)))

        src_ref_root = joinpath(src_root, "references", "fixtures")
        src_ref_dirs = isdir(src_ref_root) ? Set(filter(d -> isdir(joinpath(src_ref_root, d)), readdir(src_ref_root))) : Set{String}()
        removed_ref_dirs = sort!(collect(setdiff(src_ref_dirs, Set(unique_ids))))

        src_img_root = joinpath(src_root, "reference_images", "fixtures")
        src_img_ids =
            if isdir(src_img_root)
                Set(replace(f, "_montage.png" => "") for f in readdir(src_img_root) if endswith(f, "_montage.png"))
            else
                Set{String}()
            end
        removed_image_ids = sort!(collect(setdiff(src_img_ids, Set(unique_ids))))

        println("Removed unused fixture directories: ", length(removed_fixture_dirs))
        println("Removed unused reference directories: ", length(removed_ref_dirs))
        println("Removed unused reference images: ", length(removed_image_ids))
        !isempty(removed_fixture_dirs) && println("  fixture dirs sample: ", join(first(removed_fixture_dirs, min(10, length(removed_fixture_dirs))), ", "))
    end
end

function main(args)
    opt = parse_args(args)
    test_root = normpath(String(opt["test_root"]))
    artifacts_toml = normpath(String(opt["artifacts_toml"]))
    tarball = normpath(String(opt["tarball"]))
    artifact_name = String(opt["artifact_name"])
    url = strip(String(opt["url"]))
    prune_output = Bool(opt["prune_output"])
    prune_unused = Bool(opt["prune_unused"])

    isdir(test_root) || error("test root does not exist: $(test_root)")
    mkpath(dirname(tarball))

    hash = create_artifact() do artdir
        copy_release_dataset!(artdir, test_root; prune_output=prune_output, prune_unused=prune_unused)
    end

    println("Archiving artifact to tarball (this can take several minutes)...")
    archive_artifact(hash, tarball)
    println("Computing tarball sha256...")
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
