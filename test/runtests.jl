using TestItemRunner

const _DEFAULT_TEST_FILTER = ti -> !(:release in ti.tags)

TestItemRunner.run_tests(joinpath(@__DIR__, ".."); filter=_DEFAULT_TEST_FILTER, verbose=true)
