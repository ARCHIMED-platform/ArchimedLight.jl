# Release Process

This checklist is for maintainers preparing a registered Julia release of
ArchimedLight.jl. Run it from the repository root unless noted otherwise.

Use a clean branch and keep a short release log with the version, date, commit
SHA, Julia version, platform, fixture artifact SHA256, and the result of each
gate below.

## v0.1.3 Registry Baseline

`0.1.3` is the compatibility baseline for ArchimedLight's first release in
Julia's General registry. Earlier GitHub `0.1.x` tags were development
snapshots and are outside this stability promise.

The v0.1.3 reference refresh covers the corrected sunset hour angle and its
downstream radiation fields, plus two already-reviewed light-transfer
corrections included in the release commit: transparency in the fast upper-hit
path and authoritative explicit sun/partition inputs in directional substeps.
The reviewed release-fixture manifest is:

| Fixture | Affected step(s) | Accepted reason |
| --- | --- | --- |
| `test-cached-radiation` | 0, 12:00-12:30 | the fast upper-hit path now applies the model's `transparency: 0.1`; paired cache paths must still agree |
| `test-cached-radiation3` | 1, 05:30-06:00; sunlit steps 1-4 | corrected dawn overlap plus the fast upper-hit transparency correction; paired cache paths must still agree |
| `test-cached-radiation4` | 1, 05:30-06:00 | corrected dawn overlap |
| `test-compare-cafeier1` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier1-lowsun` | 1, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier1-lowsun1800` | 0, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier2` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier3` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier4` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-compare-cafeier5` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-links-pixeltable2` | 20, 18:00-18:30 | corrected sunset overlap |
| `test-timestep2-manysteps` | 3, 05:30-05:35 | corrected interval is before dawn and therefore has zero radiation |
| `test-timestep2-onestep` | 0, 05:15-05:45 | corrected dawn overlap |
| `test-timestep3-manysteps` | 2, 18:25-18:30 | corrected interval is after sunset and therefore has zero radiation |
| `test-timestep3-onestep` | 0, 05:15-05:45 | corrected dawn overlap |
| `test-hitcount` | 0, 11:50-11:51 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-hitcount2` | 0, 12:00-12:30 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-links` | 0, 08:00-08:30 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-links2` | 0, 08:00-08:30 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-links3` | 0, 08:00-08:30 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-links4` | 0, 08:00-08:30 | the fast upper-hit path now applies the model's `transparency: 0.1` |
| `test-toricity-two-simpleplants-border` | 0, 09:00-09:30 | explicit sun coordinates and direct fraction are now retained by directional substeps |
| `test-weighted-sun` | 2 and 28 | corrected dawn/sunset overlap and duration-weighted representative sun |

For the solar and explicit-partition fixtures, only solar-position, radiation,
derived flux/energy, and radiative image coloring may change. Meteo inputs,
geometry, topology, camera, silhouettes, component area, barycentres, and sky
fractions must remain stable. For the eight transparency fixtures, incident
flux/energy, `sky_fraction`, and radiative coloring may change by the model's
10% transparent share; sky forcing, meteo inputs, geometry, topology, camera,
silhouettes, component area, and barycentres must remain stable.

## 1. Choose The Version

1. Pick the SemVer version to register, for example `0.1.3`.
2. Update `Project.toml`.
3. Update citation metadata:
   - `CITATION.cff`
   - `CITATION.bib`
4. Check release-facing URLs and repository metadata:
   - `README.md`
   - `docs/make.jl`
   - `.lychee.toml`
   - `Artifacts.toml`, if a new release fixture artifact will be published

The git tag produced by Julia registration/TagBot must match the version, for
example `v0.1.3`.

## 2. Run The Fast Gate

```bash
julia --project=test test/runtests.jl
```

This runs the default `@testitem` suite and excludes only tests tagged
`:release`. The fast suite includes the package quality checks in
`test/aqua-test.jl` and the public API JET smoke checks in `test/jet-test.jl`.

Useful focused runs:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :fast in ti.tags, verbose=true)'
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :aqua in ti.tags, verbose=true)'
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :jet in ti.tags, verbose=true)'
```

## 3. Run The Documentation Gate

Generate the README and documentation videos on a local machine with enough
memory:

```bash
julia docs/make_video.jl
```

This creates `archimedlight_day_cycle_1.mp4`,
`archimedlight_day_cycle_2.mp4`, and their poster PNGs in
`docs/src/assets/`. Review the four generated files and commit them with the
release changes. Video rendering is intentionally a manual release step rather
than part of the `Docs` workflow because the agrivoltaic scene requires more
memory than the CI runner provides.

Build the manual:

```bash
julia --project=docs docs/make.jl
```

Run doctests:

```bash
julia --project=docs -e 'using Documenter: DocMeta, doctest; using ArchimedLight; DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true); doctest(ArchimedLight)'
```

Run the link checker when `lychee` is available:

```bash
lychee --config .lychee.toml README.md docs/src
```

The GitHub `Docs` workflow also instantiates the docs environment, runs
doctests, then deploys with `julia-actions/julia-docdeploy`. It runs on pushes
to `main`, tags, relevant pull requests, and manual dispatch. `DOCUMENTER_KEY`
must be configured for deployment and for TagBot.

## 4. Prepare Release Fixtures

Heavy full-fixture data stays outside the package repository and is exposed
through the lazy artifact in `Artifacts.toml`.

By default, keep using the current release fixture artifact. Only rebuild and
publish a new dataset when the release tests or regression matrix reveal changes
that are expected, reviewed, and accepted as the new reference behavior.

Run release tests against the artifact recorded in `Artifacts.toml`:

```bash
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :release in ti.tags, verbose=true)'
```

If the artifact URL is not available yet, or if you want to test a local fixture
checkout directly, point the harness at an extracted dataset:

```bash
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test -e 'using TestItemRunner; TestItemRunner.run_tests("test"; filter=ti -> :release in ti.tags, verbose=true)'
```

To inspect visual drift:

```bash
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test scripts/render_release_fixture_diffs.jl --all
```

If the observed changes are not intended, fix the code and keep the existing
fixture artifact. If the changes are intended and the heavy fixture dataset
already contains the accepted references, build a new data-only artifact tarball
from that checkout:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures-v0.1.3.tar.gz \
  --url https://github.com/VEZY/ArchimedLight.jl/releases/download/v0.1.3/archimedlight-release-fixtures-v0.1.3.tar.gz
```

If the current Julia outputs are the intended new references, let the script
refresh them before packaging and binding:

```bash
julia --project=. scripts/build_release_fixture_artifact.jl \
  --test-root /path/to/heavy-test-root \
  --tarball /tmp/archimedlight-release-fixtures-v0.1.3.tar.gz \
  --refresh-references \
  --url https://github.com/VEZY/ArchimedLight.jl/releases/download/v0.1.3/archimedlight-release-fixtures-v0.1.3.tar.gz
```

The script prints the artifact `git-tree-sha1`, tarball `sha256`, and a UTC
timestamp. Record these in the release log. The URL does not need to exist yet,
but the tarball uploaded later must match the SHA256 recorded in
`Artifacts.toml`.

Upload the exact verified tarball asset named in the `Artifacts.toml` URL to the
GitHub release before General registration. Do not rebuild or refresh references
at that point: the uploaded bytes must be the same bytes tested locally.

## 5. Run The Regression Matrix

Run the normal matrix:

```bash
julia --project=test test/regression_matrix/runtests.jl
```

Run the release profile against the recorded artifact:

```bash
ARCHIMEDLIGHT_REGRESSION_PROFILE=release \
julia --project=test test/regression_matrix/runtests.jl
```

Or run it against a local extracted fixture dataset:

```bash
ARCHIMEDLIGHT_REGRESSION_PROFILE=release \
ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR=/path/to/release-fixtures \
julia --project=test test/regression_matrix/runtests.jl
```

The matrix writes `regression_report.csv` under
`test/regression_matrix/reports/latest` by default. To keep release evidence
outside the working tree:

```bash
ARCHIMEDLIGHT_REGRESSION_REPORT_DIR=/tmp/archimedlight-regression-v0.1.3 \
julia --project=test test/regression_matrix/runtests.jl
```

Only strict cases fail the test run on drift; non-strict cases are report-only.
Use the report to decide whether observed drift is expected before registering.

Refresh in-repo baselines only when the drift is intended:

```bash
ARCHIMEDLIGHT_REGRESSION_UPDATE=true julia --project=test test/regression_matrix/runtests.jl
```

## 6. Run Benchmarks

Dispatch the GitHub `Benchmarks` workflow, or run the benchmark project locally
if needed:

```bash
julia --project=benchmark benchmark/benchmarks.jl
```

Record the first baseline for the release, especially one-step, cached series,
and scattering-heavy cases. The GitHub workflow uses
`MilesCranmer/AirspeedVelocity.jl@action-v1` and can be run manually.

## 7. Review Automation And Secrets

Before registering, verify:

- `Register Package` workflow exists and can be manually dispatched with the
  intended version.
- `TagBot` workflow has access to `GITHUB_TOKEN` and `DOCUMENTER_KEY`.
- `Docs` workflow has access to `DOCUMENTER_KEY`.
- `Test` workflow is green on the release commit, including tags.
- `CODECOV_TOKEN` is configured if coverage upload is expected.

Repository automation relevant to release:

- `.github/workflows/Register.yml`: manual Julia General registration action.
- `.github/workflows/TagBot.yml`: creates the git tag and GitHub release after
  the General registry PR is merged.
- `.github/workflows/Docs.yml`: builds, doctests, and deploys documentation.
- `.github/workflows/Test.yml`: runs the test matrix on `main`, tags, and
  manual dispatch.
- `.github/workflows/Benchmarks.yml`: runs AirspeedVelocity benchmarks.

## 8. Register The Release

1. Commit the final release changes and record the exact commit SHA.
2. Push the release commit to `main` and wait for `Test` and `Docs` to pass.
3. Tag that exact verified commit as `v0.1.3`, push the tag, and create the
   GitHub release.
4. Upload the already-tested release fixture tarball to that release without
   rebuilding it.
5. From an empty Julia depot, install the tagged package, download the artifact
   through `Artifacts.toml`, and run the smoke test.
6. Dispatch the `Register Package` workflow with version `0.1.3` only after the
   tagged install and artifact download succeed.
7. Wait for the Julia General registry PR to pass and merge. If TagBot runs,
   confirm it recognizes the existing exact tag/release rather than replacing
   them.
8. Confirm `Test` and `Docs` pass on the tag.
9. Confirm the stable docs URL resolves:
    `https://VEZY.github.io/ArchimedLight.jl/stable`.

## 9. Post-Release Checks

- Install the registered version in a fresh temporary environment.
- Run a small smoke example from the README or Getting Started docs.
- Confirm the GitHub release notes, citation metadata, and docs all mention the
  same version.
- If release fixture assets were uploaded after TagBot, verify the GitHub asset
  checksum still matches `Artifacts.toml`.
