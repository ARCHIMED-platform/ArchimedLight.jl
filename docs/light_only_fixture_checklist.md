# Light-Only Fixture Checklist

This checklist tracks Java test fixtures not yet fully covered by frozen Java snapshot parity in Julia.

## Wired as hard assertions in Julia

- `test-parallel`: cache on/off invariance (`Ri_PAR_0_f`, `Ri_NIR_0_f`, `Ri_PAR_f`, `Ri_NIR_f`) is asserted in `/Users/rvezy/Documents/dev/ArchimedLight/test/test_phase0_phase1.jl`.
- `test-parallel2`: same invariance as above, for both `all_in_turtle=true` and `all_in_turtle=false`.
- `test-rhvpd`: RH/VPD meteo `use` consistency now enforced and asserted (success/failure cases).
- `test-compare-cafeier1`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-cafeier2`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-cafeier3`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-cafeier4`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-cafeier5`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-simpleplant`: frozen Java snapshots added under `expected/`; component/scene/summary parity asserted.
- `test-compare-two_coffee`: frozen Java snapshots added under `expected/`; component/scene parity asserted.
- `test-compare-timestep`: frozen Java snapshots added under `expected/`; component/scene parity asserted.

## Light-only fixtures still pending frozen Java snapshots

- `test-cached-radiation5`: timing/cache behavior fixture uses legacy `app.properties` flow; no in-repo frozen `expected/` snapshots yet.

## Non-light scope (kept out of ArchimedLight parity)

- `test-autoparams` (mixed auto-activation of scattering/photosynthesis/energy balance)
- `test-emissivity` (thermal/emissivity)
- `test-skytir` (TIR meteo channel)
- `test-meteo-atmco2` (photosynthesis driver)
- `test-enbalance-flag1`
- `test-enbalance-flag2`
- `test-enbalance-flag3`
- `test-enbalance-params`
- `test-farquharc3-params`
- `test-nfacestomata1`
- `test-nfacestomata2`
- `test-nfacestomata3`
- `test-compare-photo-farqhuar-coffee`
- `test-compare-photo-farqhuar-coffee2`
- `test-compare-photo-farqhuar-coffee3`
- `test-compare-photo-farqhuar-enbalance-coffee`
- `test-compare-photo-farqhuar-enbalance-coffee2`
- `test-compare-photo-farqhuar-enbalance-coffee3`
- `test-compare-photo-farqhuar-enbalance-palm`
- `test-compare-photo-nrh`
- `test-compare-photo-nrh2`
- `test-compare-photo-nrh3`

## Infrastructure-only fixture

- `test-dir`: command-line relative config path behavior (not light physics parity).
