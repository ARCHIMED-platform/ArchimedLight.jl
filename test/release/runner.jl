using Test
using ReferenceTests
using Dates

function run_release_fixture!(fx::JuliaFixture; index::Int=1, total::Int=1)
    started = now(UTC)
    @info "Release fixture start" index=index total=total fixture=fx.id config=fx.config_path started_utc=string(started)
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
                @info "Fixture numeric mismatch" fixture=fx.id file=name missing=cmp.missing extra=cmp.extra mismatch=cmp.mismatch detail=cmp.detail
            end
            @test cmp.ok
        end

        img_ref = fixture_reference_image_path(fx)
        @test isfile(img_ref)
        fig = render_fixture_montage(fx; data=data)
        @test_reference img_ref fig by=ReferenceTests.psnr_equality(35)
    end
    elapsed = round(time() - datetime2unix(started); digits=3)
    @info "Release fixture done" index=index total=total fixture=fx.id elapsed_seconds=elapsed
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
