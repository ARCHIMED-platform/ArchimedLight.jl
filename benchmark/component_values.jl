const _COMPONENT_VALUES_PER_OWNER = 4

function _component_values_node_values(ids::Vector{Int}, scale::Float64)
    values = Dict{Int,Float64}()
    sizehint!(values, length(ids))
    @inbounds for id in ids
        values[id] = scale * (1.0 + (id % 97) / 97)
    end
    return values
end

function _component_values_spectral(ids::Vector{Int}, par_scale::Float64, nir_scale::Float64)
    return ArchimedLight.SpectralNodeValues(
        _component_values_node_values(ids, par_scale),
        _component_values_node_values(ids, nir_scale),
    )
end

function _component_values_initial_total(
    ids::Vector{Int},
    par_initial::Float64,
    nir_initial::Float64,
    par_total::Float64,
    nir_total::Float64,
)
    return ArchimedLight.InitialTotalSpectralNodeValues(
        _component_values_spectral(ids, par_initial, nir_initial),
        _component_values_spectral(ids, par_total, nir_total),
    )
end

function _component_values_step(
    ncomponents::Int;
    components_per_owner::Int=_COMPONENT_VALUES_PER_OWNER,
)
    ncomponents > 0 || throw(ArgumentError("ncomponents must be positive"))
    components_per_owner > 0 || throw(ArgumentError("components_per_owner must be positive"))

    node_ids = collect(1:ncomponents)
    source_owner = Vector{PlantGeom.SourceOwnerKey}(undef, ncomponents)
    radiative_area = Vector{Float64}(undef, ncomponents)
    @inbounds for component in eachindex(node_ids)
        owner = fld(component - 1, components_per_owner) + 1
        source_owner[component] = PlantGeom.SourceOwnerKey(1, owner)
        radiative_area[component] = 0.5 + (component % 13) / 13
    end
    metadata = ArchimedLight.LightComponentMetadata(
        node_ids,
        source_owner,
        radiative_area,
    )

    budget = ArchimedLight.LightBudget(
        _component_values_initial_total(node_ids, 10.0, 20.0, 11.0, 21.0),
        _component_values_initial_total(node_ids, 2.0, 3.0, 4.0, 5.0),
        _component_values_initial_total(node_ids, 6.0, 7.0, 8.0, 9.0),
        _component_values_initial_total(node_ids, 1.0, 2.0, 3.0, 4.0),
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
    return ArchimedLight.LightStepResult(
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
end

function _component_values_fixture(ncomponents::Int)
    step = _component_values_step(ncomponents)
    aggregation = ArchimedLight.CompiledComponentAggregation(step)
    table = ArchimedLight.component_values(step, aggregation)
    return (step=step, aggregation=aggregation, table=table)
end

# Keep this call behind a concrete function barrier. Measuring the same call at
# global scope includes REPL/compiler bookkeeping and does not represent the
# steady-state path used by a scene-level model.
@noinline function _refill_component_values!(table, step, aggregation)
    return ArchimedLight.component_values!(table, step, aggregation)
end

function _component_values_allocation_gate(fixture)
    table = fixture.table
    step = fixture.step
    aggregation = fixture.aggregation
    _refill_component_values!(table, step, aggregation)
    _refill_component_values!(table, step, aggregation)
    bytes = @allocated _refill_component_values!(table, step, aggregation)
    bytes == 0 || error(
        "component_values! allocated $bytes bytes after warmup for " *
        "$(length(step.component_metadata.node_id)) components",
    )
    return bytes
end

function _raw_component_metric_checksum(step)
    budget = step.budget
    metadata = step.component_metadata
    metadata === nothing && throw(ArgumentError("component metadata is required"))
    sources = (
        budget.incident_flux.initial.par,
        budget.incident_flux.initial.nir,
        budget.incident_flux.total.par,
        budget.incident_flux.total.nir,
        budget.incident_energy.initial.par,
        budget.incident_energy.initial.nir,
        budget.incident_energy.total.par,
        budget.incident_energy.total.nir,
        budget.absorbed_flux.initial.par,
        budget.absorbed_flux.initial.nir,
        budget.absorbed_flux.total.par,
        budget.absorbed_flux.total.nir,
        budget.absorbed_energy.initial.par,
        budget.absorbed_energy.initial.nir,
        budget.absorbed_energy.total.par,
        budget.absorbed_energy.total.nir,
    )
    checksum = 0.0
    @inbounds for source in sources
        for id in metadata.node_id
            checksum += get(source, id, 0.0)
        end
    end
    return checksum
end

const COMPONENT_VALUES_1K = _component_values_fixture(1_000)
const COMPONENT_VALUES_10K = _component_values_fixture(10_000)

# This is an executable benchmark-suite gate, not a timing threshold. Runtime
# ratios are reported by BenchmarkTools; allocations must remain exactly zero.
const COMPONENT_VALUES_ALLOCATION_GATE = (
    components_1k=_component_values_allocation_gate(COMPONENT_VALUES_1K),
    components_10k=_component_values_allocation_gate(COMPONENT_VALUES_10K),
)

SUITE["Component publication"] = BenchmarkGroup()
for (label, fixture) in (
    "1k components / 250 owners" => COMPONENT_VALUES_1K,
    "10k components / 2.5k owners" => COMPONENT_VALUES_10K,
)
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
        @benchmarkable _refill_component_values!($table, $step, $aggregation)
    group["raw 16-map read lower bound"] =
        @benchmarkable _raw_component_metric_checksum($step)
end
