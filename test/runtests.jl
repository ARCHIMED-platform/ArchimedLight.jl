using TestItemRunner

const _DEFAULT_TEST_FILTER = ti -> !(:release in ti.tags)

TestItemRunner.run_tests(joinpath(@__DIR__, ".."); filter=_DEFAULT_TEST_FILTER, verbose=true)

# Infos:
# Run this for generating figures from the release fixtures VS the reference ones:
# julia --project=test scripts/render_release_fixture_diffs.jl --all
