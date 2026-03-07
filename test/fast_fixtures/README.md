# Fast Fixtures

This folder contains small, explicit test fixtures meant to be readable as usage examples.

Each case has:

- `input/`: regular simulation inputs (`config.yml`, `meteo.csv`, scene/model files).
- `expected/`: committed references used by tests.

Cases:

- `sky_06_direct`: sky-only stages, 6 sectors, `all_in_turtle=false`.
- `sky_16_turtle`: sky-only stages, 16 sectors, `all_in_turtle=true`.
- `sky_46_direct`: sky-only stages, 46 sectors, `all_in_turtle=false`.
- `simpleplant_16_notoric`: one simple plant, 16 sectors, toricity disabled.
- `simpleplant_16_toric`: one simple plant, 16 sectors, toricity enabled.
