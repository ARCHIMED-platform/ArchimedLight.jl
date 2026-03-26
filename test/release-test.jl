@testmodule ReleaseHarness begin
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

    manifest_path = joinpath(dataset_root, "fixtures_manifest.toml")
    isfile(manifest_path) || error("Missing fixtures manifest in release dataset: $(manifest_path)")

    old_data_root = get(ENV, "ARCHIMEDLIGHT_RELEASE_DATA_ROOT", nothing)
    ENV["ARCHIMEDLIGHT_RELEASE_DATA_ROOT"] = dataset_root
    try
        include(joinpath(@__DIR__, "release", "harness.jl"))
        include(joinpath(@__DIR__, "release", "runner.jl"))
    finally
        if old_data_root === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RELEASE_DATA_ROOT")
        else
            ENV["ARCHIMEDLIGHT_RELEASE_DATA_ROOT"] = old_data_root
        end
    end
end

@testitem "Release fixture test-absorb" tags=[:release, :slow, :test_absorb] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-absorb")
end

@testitem "Release fixture test-absorb2" tags=[:release, :slow, :test_absorb2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-absorb2")
end

@testitem "Release fixture test-area_ratio" tags=[:release, :slow, :test_area_ratio] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-area_ratio")
end

@testitem "Release fixture test-area_ratio2" tags=[:release, :slow, :test_area_ratio2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-area_ratio2")
end

@testitem "Release fixture test-area_ratio3" tags=[:release, :slow, :test_area_ratio3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-area_ratio3")
end

@testitem "Release fixture test-area_ratio4" tags=[:release, :slow, :test_area_ratio4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-area_ratio4")
end

@testitem "Release fixture test-area_ratio5" tags=[:release, :slow, :test_area_ratio5] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-area_ratio5")
end

@testitem "Release fixture test-cached-radiation" tags=[:release, :slow, :test_cached_radiation] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cached-radiation")
end

@testitem "Release fixture test-cached-radiation2" tags=[:release, :slow, :test_cached_radiation2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cached-radiation2")
end

@testitem "Release fixture test-cached-radiation3" tags=[:release, :slow, :test_cached_radiation3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cached-radiation3")
end

@testitem "Release fixture test-cached-radiation4" tags=[:release, :slow, :test_cached_radiation4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cached-radiation4")
end

@testitem "Release fixture test-cafeier" tags=[:release, :slow, :test_cafeier] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier")
end

@testitem "Release fixture test-cafeier2" tags=[:release, :slow, :test_cafeier2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier2")
end

@testitem "Release fixture test-cafeier_sensor" tags=[:release, :slow, :test_cafeier_sensor] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier_sensor")
end

@testitem "Release fixture test-cafeier_sensor2" tags=[:release, :slow, :test_cafeier_sensor2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier_sensor2")
end

@testitem "Release fixture test-cafeier_sensor3" tags=[:release, :slow, :test_cafeier_sensor3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier_sensor3")
end

@testitem "Release fixture test-cafeier_sensor4" tags=[:release, :slow, :test_cafeier_sensor4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-cafeier_sensor4")
end

@testitem "Release fixture test-compare-cafeier1" tags=[:release, :slow, :test_compare_cafeier1] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier1")
end

@testitem "Release fixture test-compare-cafeier1-lowsun" tags=[:release, :slow, :test_compare_cafeier1_lowsun] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier1-lowsun")
end

@testitem "Release fixture test-compare-cafeier1-lowsun1800" tags=[:release, :slow, :test_compare_cafeier1_lowsun1800] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier1-lowsun1800")
end

@testitem "Release fixture test-compare-cafeier2" tags=[:release, :slow, :test_compare_cafeier2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier2")
end

@testitem "Release fixture test-compare-cafeier3" tags=[:release, :slow, :test_compare_cafeier3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier3")
end

@testitem "Release fixture test-compare-cafeier4" tags=[:release, :slow, :test_compare_cafeier4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier4")
end

@testitem "Release fixture test-compare-cafeier5" tags=[:release, :slow, :test_compare_cafeier5] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-cafeier5")
end

@testitem "Release fixture test-compare-simpleplant" tags=[:release, :slow, :test_compare_simpleplant] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-simpleplant")
end

@testitem "Release fixture test-compare-sky06" tags=[:release, :slow, :test_compare_sky06] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-sky06")
end

@testitem "Release fixture test-compare-sky16" tags=[:release, :slow, :test_compare_sky16] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-sky16")
end

@testitem "Release fixture test-compare-sky46" tags=[:release, :slow, :test_compare_sky46] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-sky46")
end

@testitem "Release fixture test-compare-timestep" tags=[:release, :slow, :test_compare_timestep] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-timestep")
end

@testitem "Release fixture test-compare-two_coffee" tags=[:release, :slow, :test_compare_two_coffee] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-compare-two_coffee")
end

@testitem "Release fixture test-customband" tags=[:release, :slow, :test_customband] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-customband")
end

@testitem "Release fixture test-hitcount" tags=[:release, :slow, :test_hitcount] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-hitcount")
end

@testitem "Release fixture test-hitcount2" tags=[:release, :slow, :test_hitcount2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-hitcount2")
end

@testitem "Release fixture test-hitcount3" tags=[:release, :slow, :test_hitcount3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-hitcount3")
end

@testitem "Release fixture test-independant-steps" tags=[:release, :slow, :test_independant_steps] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-independant-steps")
end

@testitem "Release fixture test-independant-steps2" tags=[:release, :slow, :test_independant_steps2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-independant-steps2")
end

@testitem "Release fixture test-links" tags=[:release, :slow, :test_links] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links")
end

@testitem "Release fixture test-links-pixeltable" tags=[:release, :slow, :test_links_pixeltable] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links-pixeltable")
end

@testitem "Release fixture test-links-pixeltable2" tags=[:release, :slow, :test_links_pixeltable2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links-pixeltable2")
end

@testitem "Release fixture test-links-sensor-plates" tags=[:release, :slow, :test_links_sensor_plates] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links-sensor-plates")
end

@testitem "Release fixture test-links-sensor-plates2" tags=[:release, :slow, :test_links_sensor_plates2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links-sensor-plates2")
end

@testitem "Release fixture test-links-stats" tags=[:release, :slow, :test_links_stats] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links-stats")
end

@testitem "Release fixture test-links2" tags=[:release, :slow, :test_links2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links2")
end

@testitem "Release fixture test-links3" tags=[:release, :slow, :test_links3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links3")
end

@testitem "Release fixture test-links4" tags=[:release, :slow, :test_links4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-links4")
end

@testitem "Release fixture test-meteo-stepduration" tags=[:release, :slow, :test_meteo_stepduration] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-meteo-stepduration")
end

@testitem "Release fixture test-parallel" tags=[:release, :slow, :test_parallel] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-parallel")
end

@testitem "Release fixture test-parallel2" tags=[:release, :slow, :test_parallel2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-parallel2")
end

@testitem "Release fixture test-save_on_disk1" tags=[:release, :slow, :test_save_on_disk1] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk1")
end

@testitem "Release fixture test-save_on_disk2" tags=[:release, :slow, :test_save_on_disk2] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk2")
end

@testitem "Release fixture test-save_on_disk3" tags=[:release, :slow, :test_save_on_disk3] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk3")
end

@testitem "Release fixture test-save_on_disk4" tags=[:release, :slow, :test_save_on_disk4] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk4")
end

@testitem "Release fixture test-save_on_disk5" tags=[:release, :slow, :test_save_on_disk5] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk5")
end

@testitem "Release fixture test-save_on_disk6" tags=[:release, :slow, :test_save_on_disk6] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-save_on_disk6")
end

@testitem "Release fixture test-scattering-one-plate" tags=[:release, :slow, :test_scattering_one_plate] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-scattering-one-plate")
end

@testitem "Release fixture test-scattering-two-plates" tags=[:release, :slow, :test_scattering_two_plates] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-scattering-two-plates")
end

@testitem "Release fixture test-scattering-two-simpleplants" tags=[:release, :slow, :test_scattering_two_simpleplants] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-scattering-two-simpleplants")
end

@testitem "Release fixture test-skytir" tags=[:release, :slow, :test_skytir] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-skytir")
end

@testitem "Release fixture test-timestep1-manysteps" tags=[:release, :slow, :test_timestep1_manysteps] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep1-manysteps")
end

@testitem "Release fixture test-timestep1-onestep" tags=[:release, :slow, :test_timestep1_onestep] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep1-onestep")
end

@testitem "Release fixture test-timestep2-manysteps" tags=[:release, :slow, :test_timestep2_manysteps] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep2-manysteps")
end

@testitem "Release fixture test-timestep2-onestep" tags=[:release, :slow, :test_timestep2_onestep] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep2-onestep")
end

@testitem "Release fixture test-timestep3-manysteps" tags=[:release, :slow, :test_timestep3_manysteps] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep3-manysteps")
end

@testitem "Release fixture test-timestep3-onestep" tags=[:release, :slow, :test_timestep3_onestep] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-timestep3-onestep")
end

@testitem "Release fixture test-toricity-two-simpleplants-border" tags=[:release, :slow, :test_toricity_two_simpleplants_border] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-toricity-two-simpleplants-border")
end

@testitem "Release fixture test-weighted-sun" tags=[:release, :slow, :test_weighted_sun] setup=[ReleaseHarness] begin
    ReleaseHarness.run_release_fixture!("test-weighted-sun")
end
