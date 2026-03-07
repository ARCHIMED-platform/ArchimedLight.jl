using Artifacts: artifact_hash, artifact_path
using Pkg.Artifacts: ensure_artifact_installed

repo_root = normpath(joinpath(@__DIR__, ".."))
artifacts_toml = joinpath(repo_root, "Artifacts.toml")
artifact_name = get(ENV, "ARCHIMEDLIGHT_RELEASE_ARTIFACT_NAME", "archimedlight-release-fixtures")
dataset_root = strip(get(ENV, "ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR", ""))

if isempty(dataset_root) && isfile(artifacts_toml)
    hash = artifact_hash(artifact_name, artifacts_toml)
    if hash !== nothing
        ensure_artifact_installed(artifact_name, artifacts_toml)
        dataset_root = artifact_path(hash)
    end
end

isempty(dataset_root) && error(
    "Release dataset not found. Set ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR or register artifact $(repr(artifact_name)) in Artifacts.toml.",
)
isdir(dataset_root) || error("Release dataset directory does not exist: $(dataset_root)")

release_runtests = joinpath(dataset_root, "runtests_release.jl")
isfile(release_runtests) || error(
    "Missing release entrypoint: $(release_runtests). Expected release dataset to provide runtests_release.jl.",
)

include(release_runtests)
