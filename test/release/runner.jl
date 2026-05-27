using Test
using ReferenceTests
using Dates

const _CACHE_PAIR_FIXTURES = Set(["test-cached-radiation", "test-cached-radiation3"])

function _cached_pair_fixture(fx::JuliaFixture)
    cfg2 = joinpath(dirname(fx.config_path), "config2.yml")
    isfile(cfg2) || error("$(fx.id): missing paired cached-radiation config2.yml")
    return JuliaFixture(
        fx.id * "-config2",
        cfg2,
        fx.visual_metric,
        true,
        fx.scene_override,
        fx.meteo_override,
        fx.force_scattering,
    )
end

function run_cached_radiation_pair_fixture!(fx::JuliaFixture; index::Int=1, total::Int=1)
    @testset "$(fx.id)" begin
        data1 = fixture_runtime_data(fx)
        data2 = fixture_runtime_data(_cached_pair_fixture(fx))
        out1 = mktempdir()
        out2 = mktempdir()
        observed1 = write_fixture_observed_outputs!(fx, out1; data=data1)
        observed2 = write_fixture_observed_outputs!(fx, out2; data=data2)

        for name in ("component_values.csv",)
            @test haskey(observed1, name)
            @test haskey(observed2, name)
            cmp = compare_csv_reference(observed1[name], observed2[name]; label="$(fx.id):$(name)")
            if !cmp.ok
                @info "Cached radiation pair mismatch" fixture = fx.id file = name missing = cmp.missing extra = cmp.extra mismatch = cmp.mismatch detail = cmp.detail
            end
            @test cmp.ok
        end
    end
    return nothing
end

function run_release_fixture!(fx::JuliaFixture; index::Int=1, total::Int=1)
    if fx.id in _CACHE_PAIR_FIXTURES
        return run_cached_radiation_pair_fixture!(fx; index=index, total=total)
    end

    started = now(UTC)
    @testset "$(fx.id)" begin
        data = fixture_runtime_data(fx)
        outdir = mktempdir()
        observed = write_fixture_observed_outputs!(fx, outdir; data=data)

        refs = fixture_numeric_reference_paths(fx; existing_only=true)
        @test !isempty(refs)
        for required in ("component_values.csv", "scene_values.csv", "meteo.csv")
            @test haskey(refs, required)
        end
        for (name, ref_path) in refs
            @test isfile(ref_path)
            obs_path = get(observed, name, "")
            cmp =
                if isempty(obs_path) || !isfile(obs_path)
                    (ok=false, missing=0, extra=0, mismatch=0, detail="observed file missing: $(name)")
                else
                    compare_csv_reference(ref_path, obs_path; label="$(fx.id):$(name)")
                end
            if !cmp.ok
                @info "Fixture numeric mismatch" fixture = fx.id file = name missing = cmp.missing extra = cmp.extra mismatch = cmp.mismatch detail = cmp.detail
            end
            @test cmp.ok
        end

        img_ref = fixture_reference_image_path(fx)
        @test isfile(img_ref)
        fig = render_fixture_montage(fx; data=data)
        @test_reference img_ref fig by = ReferenceTests.psnr_equality(35)
    end
    elapsed = round(time() - datetime2unix(started); digits=3)
    return nothing
end

function run_release_fixture!(id::AbstractString)
    fx = fixture_by_id(id)
    fx === nothing && error("Unknown release fixture: $(repr(id))")
    return run_release_fixture!(fx)
end

function run_release_fixture_regression!(fixtures::Vector{JuliaFixture}=julia_fixtures())
    @testset "Julia fixture regression (release dataset)" begin
        @test !isempty(fixtures)

        total = length(fixtures)
        for (i, fx) in enumerate(fixtures)
            run_release_fixture!(fx; index=i, total=total)
        end
    end
    return nothing
end
