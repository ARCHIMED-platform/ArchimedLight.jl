@testmodule ComponentValuesHelper begin
    using ArchimedLight
    import PlantGeom

    const COMPONENT_IDS = [11, 12, 13]
    const OWNER_A = PlantGeom.SourceOwnerKey(4, 101)
    const OWNER_B = PlantGeom.SourceOwnerKey(7, 202)
    const OWNERS = [OWNER_A, OWNER_A, OWNER_B]
    const RADIATIVE_AREAS = [1.0, 3.0, 2.0]

    _node_values(values) = Dict(COMPONENT_IDS .=> Float64.(values))

    _spectral(par, nir) = ArchimedLight.SpectralNodeValues(
        _node_values(par),
        _node_values(nir),
    )

    _initial_total(par_initial, nir_initial, par_total, nir_total) =
        ArchimedLight.InitialTotalSpectralNodeValues(
            _spectral(par_initial, nir_initial),
            _spectral(par_total, nir_total),
        )

    function controlled_budget(; factor=1.0)
        scale(values) = factor .* values
        return ArchimedLight.LightBudget(
            _initial_total(
                scale([10, 30, 50]),
                scale([20, 40, 60]),
                scale([11, 31, 51]),
                scale([21, 41, 61]),
            ),
            _initial_total(
                scale([2, 7, 11]),
                scale([3, 8, 12]),
                scale([4, 9, 13]),
                scale([5, 10, 14]),
            ),
            _initial_total(
                scale([6, 16, 26]),
                scale([7, 17, 27]),
                scale([8, 18, 28]),
                scale([9, 19, 29]),
            ),
            _initial_total(
                scale([1, 5, 9]),
                scale([2, 6, 10]),
                scale([3, 7, 11]),
                scale([4, 8, 12]),
            ),
            Dict{String,Dict{Int,Float64}}(),
            Dict{String,Dict{Int,Float64}}(),
        )
    end

    function controlled_step(
        metadata=ArchimedLight.LightComponentMetadata(
            copy(COMPONENT_IDS),
            copy(OWNERS),
            copy(RADIATIVE_AREAS),
        );
        factor=1.0,
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
            controlled_budget(; factor=factor),
            Dict{String,Float64}(),
            nothing,
            nothing,
            nothing,
            metadata,
        )
    end

    const FLUX_COLUMNS = (
        :Ri_PAR_0_f,
        :Ri_NIR_0_f,
        :Ri_PAR_f,
        :Ri_NIR_f,
        :Ra_PAR_0_f,
        :Ra_NIR_0_f,
        :Ra_PAR_f,
        :Ra_NIR_f,
    )

    const ENERGY_COLUMNS = (
        :Ri_PAR_0_q,
        :Ri_NIR_0_q,
        :Ri_PAR_q,
        :Ri_NIR_q,
        :Ra_PAR_0_q,
        :Ra_NIR_0_q,
        :Ra_PAR_q,
        :Ra_NIR_q,
    )

    function expected_component_columns(budget)
        return (
            Ri_PAR_0_f=[budget.incident_flux.initial.par[id] for id in COMPONENT_IDS],
            Ri_NIR_0_f=[budget.incident_flux.initial.nir[id] for id in COMPONENT_IDS],
            Ri_PAR_f=[budget.incident_flux.total.par[id] for id in COMPONENT_IDS],
            Ri_NIR_f=[budget.incident_flux.total.nir[id] for id in COMPONENT_IDS],
            Ri_PAR_0_q=[budget.incident_energy.initial.par[id] for id in COMPONENT_IDS],
            Ri_NIR_0_q=[budget.incident_energy.initial.nir[id] for id in COMPONENT_IDS],
            Ri_PAR_q=[budget.incident_energy.total.par[id] for id in COMPONENT_IDS],
            Ri_NIR_q=[budget.incident_energy.total.nir[id] for id in COMPONENT_IDS],
            Ra_PAR_0_f=[budget.absorbed_flux.initial.par[id] for id in COMPONENT_IDS],
            Ra_NIR_0_f=[budget.absorbed_flux.initial.nir[id] for id in COMPONENT_IDS],
            Ra_PAR_f=[budget.absorbed_flux.total.par[id] for id in COMPONENT_IDS],
            Ra_NIR_f=[budget.absorbed_flux.total.nir[id] for id in COMPONENT_IDS],
            Ra_PAR_0_q=[budget.absorbed_energy.initial.par[id] for id in COMPONENT_IDS],
            Ra_NIR_0_q=[budget.absorbed_energy.initial.nir[id] for id in COMPONENT_IDS],
            Ra_PAR_q=[budget.absorbed_energy.total.par[id] for id in COMPONENT_IDS],
            Ra_NIR_q=[budget.absorbed_energy.total.nir[id] for id in COMPONENT_IDS],
        )
    end
end

@testitem "Zero-area scene facets stay outside component coupling metadata" tags = [:component_values, :fast] begin
    using ArchimedLight

    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    sim, meteo = read_simulation(config; plot_paving_override=0)
    update_options!(
        sim,
        LightOptions(
            sim.options;
            turtle_sectors=1,
            pixel_size=0.1,
            scattering=false,
            cache_pixel_table=false,
        ),
    )

    raw_area = ArchimedLight._interception_area_per_node_local(
        sim.scene,
        sim.models,
        sim.options,
    )
    zero_area_ids = sort!(Int[id for (id, area) in raw_area if iszero(area)])
    positive_area_ids = sort!(Int[id for (id, area) in raw_area if area > 0.0])
    @test !isempty(zero_area_ids)
    @test !isempty(positive_area_ids)

    step = run_light(sim, first(meteo))
    metadata = step.component_metadata
    @test metadata !== nothing
    @test collect(metadata.node_id) == positive_area_ids
    @test collect(metadata.radiative_area) == [raw_area[id] for id in positive_area_ids]
    @test isempty(intersect(zero_area_ids, collect(metadata.node_id)))
    @test all(>(0.0), metadata.radiative_area)
end

@testitem "Typed component light values follow LightBudget" tags = [:component_values, :fast] setup = [ComponentValuesHelper] begin
    using ArchimedLight
    import PlantGeom
    import Tables

    step = ComponentValuesHelper.controlled_step()
    table = component_values(step; par_energy_to_photon=5.0)
    expected = ComponentValuesHelper.expected_component_columns(step.budget)

    @test Tables.istable(typeof(table))
    @test Tables.columnaccess(typeof(table))
    @test table.node_id == ComponentValuesHelper.COMPONENT_IDS
    @test table.source_owner == ComponentValuesHelper.OWNERS
    @test table.radiative_mesh_area == ComponentValuesHelper.RADIATIVE_AREAS
    @test :area ∉ propertynames(table)
    @test eltype(table.node_id) === Int
    @test eltype(table.source_owner) === PlantGeom.SourceOwnerKey
    @test eltype(table.radiative_mesh_area) === Float64
    @test_throws Base.CanonicalIndexError setindex!(table.node_id, 99, 1)
    @test_throws Base.CanonicalIndexError setindex!(table.source_owner, ComponentValuesHelper.OWNER_B, 1)
    @test_throws Base.CanonicalIndexError setindex!(table.radiative_mesh_area, 99.0, 1)

    identity_columns = (
        (table.node_id, 99),
        (table.source_owner, ComponentValuesHelper.OWNER_B),
        (table.radiative_mesh_area, 99.0),
    )
    for (column, replacement) in identity_columns
        before = copy(column)
        @test parent(column) === column
        @test_throws Base.CanonicalIndexError setindex!(parent(column), replacement, 1)
        @test column == before
    end

    for name in keys(expected)
        @test getproperty(table, name) == getproperty(expected, name)
        @test eltype(getproperty(table, name)) === Float64
    end
    @test table.Ra_SW_f == table.Ra_PAR_f .+ table.Ra_NIR_f
    @test table.aPPFD == 5.0 .* table.Ra_PAR_f
    @test_throws ArgumentError component_values(step; level=:plant)
    @test_throws ArgumentError component_values(step; par_energy_to_photon=0.0)
    @test_throws ArgumentError component_values(step; par_energy_to_photon=NaN)

    raw = light_metric_values(step, :Ra_PAR_f)
    @test table.Ra_PAR_f == [raw[id] for id in table.node_id]

    no_metadata = ArchimedLight.LightStepResult(
        step.sky,
        step.turtle,
        step.fluxes,
        step.first_order,
        step.scattering,
        step.budget,
        step.extra_band_irradiance,
    )
    @test_throws ArgumentError component_values(no_metadata)
    @test_throws DimensionMismatch ArchimedLight.LightComponentMetadata(
        [1, 2],
        PlantGeom.SourceOwnerKey[ComponentValuesHelper.OWNER_A],
        [1.0, 2.0],
    )
    @test_throws ArgumentError ArchimedLight.LightComponentMetadata(
        [1, 1],
        PlantGeom.SourceOwnerKey[ComponentValuesHelper.OWNER_A, ComponentValuesHelper.OWNER_B],
        [1.0, 2.0],
    )
    for invalid_area in (0.0, -1.0, NaN, Inf)
        @test_throws ArgumentError ArchimedLight.LightComponentMetadata(
            [1],
            PlantGeom.SourceOwnerKey[ComponentValuesHelper.OWNER_A],
            [invalid_area],
        )
    end

    ids = [1]
    owners = PlantGeom.SourceOwnerKey[ComponentValuesHelper.OWNER_A]
    areas = [1.0]
    retained = ArchimedLight.LightComponentMetadata(ids, owners, areas)
    ids[1] = 2
    owners[1] = ComponentValuesHelper.OWNER_B
    areas[1] = 2.0
    @test retained.node_id == [1]
    @test retained.source_owner == [ComponentValuesHelper.OWNER_A]
    @test retained.radiative_area == [1.0]
end

@testitem "Component aggregation is owner keyed and conservative" tags = [:component_values, :fast] setup = [ComponentValuesHelper] begin
    using ArchimedLight

    step = ComponentValuesHelper.controlled_step()
    component_table = component_values(step)
    aggregation = CompiledComponentAggregation(step)
    owner_table = component_values(step, aggregation; par_energy_to_photon=5.0)
    implicit_owner_table = component_values(step; level=:source_owner, par_energy_to_photon=5.0)

    @test owner_table == implicit_owner_table
    @test owner_table.source_owner == [
        ComponentValuesHelper.OWNER_A,
        ComponentValuesHelper.OWNER_B,
    ]
    @test owner_table.radiative_mesh_area == [4.0, 2.0]

    for name in ComponentValuesHelper.ENERGY_COLUMNS
        component_values_column = getproperty(component_table, name)
        owner_values_column = getproperty(owner_table, name)
        @test owner_values_column == [
            component_values_column[1] + component_values_column[2],
            component_values_column[3],
        ]
        @test sum(owner_values_column) == sum(component_values_column)
    end

    for name in ComponentValuesHelper.FLUX_COLUMNS
        component_values_column = getproperty(component_table, name)
        owner_values_column = getproperty(owner_table, name)
        @test owner_values_column == [
            (component_values_column[1] + 3.0 * component_values_column[2]) / 4.0,
            component_values_column[3],
        ]
        @test sum(owner_values_column .* owner_table.radiative_mesh_area) ==
              sum(component_values_column .* component_table.radiative_mesh_area)
    end

    @test owner_table.Ra_SW_f == owner_table.Ra_PAR_f .+ owner_table.Ra_NIR_f
    @test owner_table.aPPFD == 5.0 .* owner_table.Ra_PAR_f
end

@testitem "Compiled aggregation validates generation and reuses buffers" tags = [:component_values, :fast] setup = [ComponentValuesHelper] begin
    using ArchimedLight
    import PlantGeom

    @noinline function refill_component_values!(table, step, aggregation)
        return component_values!(table, step, aggregation)
    end

    step = ComponentValuesHelper.controlled_step()
    aggregation = CompiledComponentAggregation(step)
    @test_throws Base.CanonicalIndexError setindex!(aggregation.owner_keys, ComponentValuesHelper.OWNER_B, 1)
    @test_throws Base.CanonicalIndexError setindex!(aggregation.component_to_owner, 0, 1)
    @test_throws Base.CanonicalIndexError setindex!(aggregation.owner_radiative_area, 0.0, 1)
    @test_throws ArgumentError CompiledComponentAggregation(
        step.component_metadata,
        collect(aggregation.owner_keys),
        [0, 1, 2],
        collect(aggregation.owner_radiative_area),
    )
    table = component_values(step, aggregation)
    baseline = component_values(step, aggregation)
    column_names = (
        ComponentValuesHelper.FLUX_COLUMNS...,
        ComponentValuesHelper.ENERGY_COLUMNS...,
        :Ra_SW_f,
        :aPPFD,
    )
    result_columns = map(name -> getproperty(table, name), column_names)

    scaled_step = ComponentValuesHelper.controlled_step(step.component_metadata; factor=2.0)
    returned = component_values!(table, scaled_step, aggregation)
    @test returned === table
    @test all(
        getproperty(table, name) === result_columns[i]
        for (i, name) in pairs(column_names)
    )
    @test all(
        getproperty(table, name) == 2.0 .* getproperty(baseline, name)
        for name in column_names
    )

    refill_component_values!(table, scaled_step, aggregation)
    refill_component_values!(table, scaled_step, aggregation)
    @test @allocated(refill_component_values!(table, scaled_step, aggregation)) == 0

    same_values_new_generation = ArchimedLight.LightComponentMetadata(
        copy(step.component_metadata.node_id),
        copy(step.component_metadata.source_owner),
        copy(step.component_metadata.radiative_area),
    )
    other_step = ComponentValuesHelper.controlled_step(same_values_new_generation)
    @test_throws ArgumentError component_values(other_step, aggregation)
    @test_throws ArgumentError component_values!(table, other_step, aggregation)

    foreign_identity_buffer = merge(table, (source_owner=copy(table.source_owner),))
    @test_throws ArgumentError component_values!(foreign_identity_buffer, scaled_step, aggregation)
    foreign_area_buffer = merge(table, (radiative_mesh_area=copy(table.radiative_mesh_area),))
    @test_throws ArgumentError component_values!(foreign_area_buffer, scaled_step, aggregation)

    baseline_after_valid_fill = copy(table.Ra_PAR_f)
    short_column_buffer = merge(table, (aPPFD=zeros(Float64, 1),))
    @test_throws DimensionMismatch component_values!(short_column_buffer, scaled_step, aggregation)
    integer_column_buffer = merge(table, (aPPFD=zeros(Int, length(table.aPPFD)),))
    @test_throws ArgumentError component_values!(integer_column_buffer, scaled_step, aggregation)
    read_only_column_buffer = merge(
        table,
        (aPPFD=ArchimedLight._ReadOnlyVector(copy(table.aPPFD)),),
    )
    @test_throws ArgumentError component_values!(read_only_column_buffer, scaled_step, aggregation)
    @test_throws ArgumentError component_values!(
        table,
        scaled_step,
        aggregation;
        par_energy_to_photon=NaN,
    )
    @test table.Ra_PAR_f == baseline_after_valid_fill
end

@testitem "Component identity snapshots survive scene replacement" tags = [:component_values, :fast] begin
    using ArchimedLight
    using CSV
    using GeometryBasics
    using PlantGeom
    import Tables

    triangle(z) = GeometryBasics.Mesh(
        GeometryBasics.Point3f[(0, 0, z), (1, 0, z), (0, 1, z)],
        GeometryBasics.TriangleFace{Int}[GeometryBasics.TriangleFace{Int}(1, 2, 3)],
    )
    scene(id, z; units=SceneUnits()) = make_scene(
        domain=(0.0, 0.0, 1.0, 1.0),
        units=units,
    ) do builder
        add_object!(builder, triangle(z); group="probe", type="Leaf", id=id)
    end

    models = models_for(
        "probe" => ("Leaf" => translucent(par=0.15, nir=0.30),),
    )
    options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=false,
        scattering=false,
        toricity=false,
        pixel_size=0.1,
    )
    centimetre_scene = scene(
        9,
        0.0;
        units=SceneUnits(length=PlantGeom.Unitful.u"cm"),
    )
    @test_throws ArgumentError ArchimedLight.prepare_light_cache(
        centimetre_scene,
        models,
        options,
    )

    sim = LightSimulation(scene(1, 0.0), models; options=options)
    sky = SkyState(180.0, 90.0, 100.0, 50.0, 1.0, 0.0)
    old_step = run_light(sim, sky; step_duration_seconds=60.0)
    old_metadata = old_step.component_metadata
    old_table = component_values(old_step)
    old_aggregation = CompiledComponentAggregation(old_step)

    @test old_metadata !== nothing
    @test old_table.radiative_mesh_area == [0.5]
    @test old_table.Ra_PAR_f == [old_step.budget.absorbed_flux.total.par[only(old_table.node_id)]]

    legacy_path = joinpath(mktempdir(), "component_values.csv")
    write_component_values(legacy_path, sim, old_step)
    legacy_row = only(Tables.rowtable(CSV.File(legacy_path; delim=';', normalizenames=false)))
    @test Int(legacy_row.node_id) == only(old_table.node_id)
    @test Float64(legacy_row.Ri_PAR_0_f) == only(old_table.Ri_PAR_0_f)
    @test Float64(legacy_row.Ra_PAR_0_q) == only(old_table.Ra_PAR_0_q)

    update_scene!(sim, scene(2, 1.0))
    @test old_step.component_metadata === old_metadata
    @test component_values(old_step) == old_table
    @test component_values(old_step, old_aggregation).source_owner == old_table.source_owner

    new_step = run_light(sim, sky; step_duration_seconds=60.0)
    @test new_step.component_metadata !== old_metadata
    # Independently assembled scenes may reuse the same scene-local owner key;
    # generation safety is carried by the retained metadata snapshot identity.
    @test new_step.component_metadata.source_owner == old_metadata.source_owner
    @test_throws ArgumentError component_values(new_step, old_aggregation)

    new_table = component_values(new_step)
    raw = light_metric_values(new_step, :Ri_PAR_0_f)
    @test new_table.Ri_PAR_0_f == [raw[id] for id in new_table.node_id]
end
