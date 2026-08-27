function _component_values_fixture(
    ncomponents::Int;
    components_per_owner::Int=4,
)
    ncomponents > 0 || throw(ArgumentError("ncomponents must be positive"))
    components_per_owner > 0 ||
        throw(ArgumentError("components_per_owner must be positive"))

    node_ids = collect(1:ncomponents)
    node_values(scale) = Dict(
        id => scale * (1.0 + (id % 97) / 97)
        for id in node_ids
    )
    spectral(par_scale, nir_scale) = ArchimedLight.SpectralNodeValues(
        node_values(par_scale),
        node_values(nir_scale),
    )
    initial_total(par_initial, nir_initial, par_total, nir_total) =
        ArchimedLight.InitialTotalSpectralNodeValues(
            spectral(par_initial, nir_initial),
            spectral(par_total, nir_total),
        )

    source_owner = [
        PlantGeom.SourceOwnerKey(1, fld(component - 1, components_per_owner) + 1)
        for component in node_ids
    ]
    radiative_area = [0.5 + (component % 13) / 13 for component in node_ids]
    metadata = ArchimedLight.LightComponentMetadata(
        node_ids,
        source_owner,
        radiative_area,
    )
    budget = ArchimedLight.LightBudget(
        initial_total(10.0, 20.0, 11.0, 21.0),
        initial_total(2.0, 3.0, 4.0, 5.0),
        initial_total(6.0, 7.0, 8.0, 9.0),
        initial_total(1.0, 2.0, 3.0, 4.0),
        Dict{String,Dict{Int,Float64}}(),
        Dict{String,Dict{Int,Float64}}(),
    )
    empty_spectral = ArchimedLight.SpectralNodeValues(
        Dict{Int,Float64}(),
        Dict{Int,Float64}(),
    )
    first_order = ArchimedLight.FirstOrderResult(
        Dict{Int,Float64}(),
        empty_spectral,
        Dict{Int,Int}(),
    )
    step = ArchimedLight.LightStepResult(
        ArchimedLight.SkyState(180.0, 45.0, 100.0, 50.0, 1.0, 0.0),
        ArchimedLight.TurtleGrid(ArchimedLight.TurtleSector[]),
        ArchimedLight.DirectionalFluxes(Int[], Float64[], Float64[]),
        first_order,
        nothing,
        budget,
        Dict{String,Float64}(),
        nothing,
        nothing,
        nothing,
        metadata,
    )
    aggregation = ArchimedLight.CompiledComponentAggregation(step)
    table = ArchimedLight.component_values(step, aggregation)
    return (step=step, aggregation=aggregation, table=table)
end

const COMPONENT_VALUES_FIXTURES = (
    "1k components / 250 owners" => _component_values_fixture(1_000),
    "10k components / 2.5k owners" => _component_values_fixture(10_000),
)

SUITE["Component publication"] = BenchmarkGroup()
for (label, fixture) in COMPONENT_VALUES_FIXTURES
    group = SUITE["Component publication"][label] = BenchmarkGroup()
    step = fixture.step
    aggregation = fixture.aggregation
    table = fixture.table

    group["compile owner aggregation"] =
        @benchmarkable ArchimedLight.CompiledComponentAggregation($step) evals = 1
    group["component table"] =
        @benchmarkable ArchimedLight.component_values($step; level=:component) evals = 1
    group["owner table allocated"] =
        @benchmarkable ArchimedLight.component_values($step, $aggregation) evals = 1
    group["owner table refill"] =
        @benchmarkable ArchimedLight.component_values!($table, $step, $aggregation)
end
