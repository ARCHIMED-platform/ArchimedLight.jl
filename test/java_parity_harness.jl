import CSV
import Tables

struct ParityFixture
    name::String
    config_path::String
    scene_path::String
    meteo_path::String
    expected_dir::Union{Nothing,String}
end

function _repo_root()
    normpath(joinpath(@__DIR__, ".."))
end

function _fixture_root()
    joinpath(_repo_root(), "java_implementation", "archimed-lib-2018", "tests")
end

function _fixture_expected_dir(test_dir::String)
    d = joinpath(test_dir, "expected")
    isdir(d) ? d : nothing
end

function _fixture_row(name::String)
    test_dir = joinpath(_fixture_root(), name)
    cfg = joinpath(test_dir, "config.yml")
    scene = ""
    meteo = ""
    if isfile(cfg)
        c = ArchimedLight.read_light_config(cfg)
        scene = c.scene
        meteo = c.meteo
    end
    ParityFixture(name, cfg, scene, meteo, _fixture_expected_dir(test_dir))
end

function light_parity_fixtures()
    names = [
        "test-hitcount",
        "test-hitcount2",
        "test-hitcount3",
        "test-weighted-sun",
        "test-skytir",
        "test-area_ratio5",
        "test-scattering-one-plate",
        "test-scattering-two-plates",
        "test-scattering-divergence",
        "test-links",
        "test-links-pixeltable",
        "test-cached-radiation2",
        "test-cached-radiation3",
        "test-customband",
    ]
    [_fixture_row(n) for n in names]
end

function read_java_csv(path::AbstractString)
    Tables.rowtable(CSV.File(path; delim=';', ignorerepeated=true, comment="#")) |> collect
end

