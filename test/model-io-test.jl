@testitem "Model IO YAML and typed models agree" tags=[:model_io, :fast] begin
    using OrderedCollections: OrderedDict

    tmp = mktempdir()
    model_path = joinpath(tmp, "plant.yml")
    write(
        model_path,
        join(
            [
                "Group: plant",
                "Type:",
                "  Leaf:",
                "    Interception:",
                "      model: Translucent",
                "      transparency: 0.1",
                "      optical_properties:",
                "        PAR: 0.15",
                "        NIR: 0.30",
            ],
            "\n",
        ) * "\n",
    )

    read_back = ArchimedLight.read_models(model_path)
    prepared = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "plant";
            types=OrderedDict(
                "Leaf" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.1,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])

    @test collect(keys(read_back.groups)) == collect(keys(prepared.groups))
    @test read_back["plant"].types["Leaf"].interception.model == "Translucent"
    @test isapprox(read_back["plant"].types["Leaf"].interception.optical_properties.par, 0.15; atol=1e-12)
    @test isapprox(read_back["plant"].types["Leaf"].interception.transparency, 0.1; atol=1e-12)
end

@testitem "Model IO missing NIR falls back to defaults" tags=[:model_io, :fast] begin
    tmp = mktempdir()
    model_path = joinpath(tmp, "plant.yml")
    write(
        model_path,
        join(
            [
                "Group: plant",
                "Type:",
                "  Leaf:",
                "    Interception:",
                "      model: Translucent",
                "      optical_properties:",
                "        PAR: 0.15",
            ],
            "\n",
        ) * "\n",
    )

    read_back = ArchimedLight.read_models(model_path)
    props = read_back["plant"].types["Leaf"].interception.optical_properties
    @test isapprox(props.par, 0.15; atol=1e-12)
    @test isapprox(props.nir, 0.0; atol=1e-12)
    coeffs = ArchimedLight._group_optical_coeffs(read_back)
    @test coeffs[("plant", "Leaf")]["PAR"] == 0.15
    @test !haskey(coeffs[("plant", "Leaf")], "NIR")
end

@testitem "Custom optical coefficients drive scattering and sensors stay passive" tags=[:model_io, :fast] begin
    tmp = mktempdir()
    model_path = joinpath(tmp, "plant.yml")
    write(
        model_path,
        join(
            [
                "Group: plant",
                "Type:",
                "  Leaf:",
                "    Interception:",
                "      model: Translucent",
                "      optical_properties:",
                "        PAR: 0.20",
                "        NIR: 0.30",
                "        custom: 0.80",
                "        TIR: 0.95",
                "  Sensor:",
                "    Interception:",
                "      model: VirtualSensor",
            ],
            "\n",
        ) * "\n",
    )

    models = ArchimedLight.read_models(model_path)
    props = models["plant"].types["Leaf"].interception.optical_properties
    @test props.extras["custom"] == 0.80

    group_coeffs = ArchimedLight._group_optical_coeffs(models)
    @test group_coeffs[("plant", "Leaf")]["CUSTOM"] == 0.80
    @test !haskey(group_coeffs[("plant", "Leaf")], "TIR")

    node_ids = [1, 2, 3]
    node_group = Dict(nid => "plant" for nid in node_ids)
    node_type = Dict(1 => "Leaf", 2 => "Leaf", 3 => "Sensor")
    options = ArchimedLight.LightOptions(scattering_max_iter=4, scattering_stop_ratio=1e-12)
    coeff_par, coeff_nir = ArchimedLight._coeff_maps_by_node(
        node_ids,
        node_group,
        node_type,
        group_coeffs,
        options,
    )
    graph = ArchimedLight.ScatteringTransferGraph(
        Dict((2, 1) => 1, (2, 3) => 1),
        Dict(nid => 1 for nid in node_ids),
        node_ids,
        node_group,
        node_type,
        group_coeffs,
        coeff_par,
        coeff_nir,
        options.scattering_coeff_par,
        options.scattering_coeff_nir,
    )
    first = ArchimedLight.FirstOrderResult(
        Dict(nid => 1.0 for nid in node_ids),
        ArchimedLight.SpectralNodeValues(
            Dict(nid => 0.0 for nid in node_ids),
            Dict(nid => 0.0 for nid in node_ids),
        ),
        Dict(nid => 1 for nid in node_ids),
    )
    custom = ArchimedLight.compute_scattering_band(
        graph,
        first,
        options;
        band="custom",
        initial_power_per_node=Dict(1 => 100.0, 2 => 0.0, 3 => 100.0),
    )

    # The leaf contributes 100 * 0.80 / 2 = 40. The sensor contributes zero
    # for this custom band rather than falling back to the PAR default.
    @test isapprox(custom.added_power_per_node[2], 40.0; atol=1e-12)
end

@testitem "Thermal irradiance is not an extra light band" tags=[:model_io, :fast] begin
    meteo_row = (RI_PAR_f=100.0, RI_UV_f=25.0, RI_TIR_f=450.0)
    extras = ArchimedLight._extra_band_irradiance(meteo_row)

    @test extras == Dict("UV" => 25.0)
    @test !haskey(extras, "TIR")
end

@testitem "Scattering convergence retains Java virtual-sensor observations" tags=[:model_io, :fast] begin
    node_ids = [1, 2, 3]
    node_group = Dict(1 => "plant", 2 => "sensor", 3 => "plant")
    node_type = Dict(1 => "Leaf", 2 => "Sensor", 3 => "Leaf")
    group_coeffs = Dict(
        ("plant", "Leaf") => Dict("PAR" => 0.15),
        ("sensor", "Sensor") => Dict("PAR" => 0.0, "__ALL_BANDS__" => 0.0),
    )
    coeff_par = Dict(1 => 0.15, 2 => 0.0, 3 => 0.15)
    coeff_nir = Dict(1 => 0.30, 2 => 0.0, 3 => 0.30)
    graph = ArchimedLight.ScatteringTransferGraph(
        Dict((1, 3) => 1, (2, 3) => 1, (3, 1) => 1, (2, 1) => 1),
        Dict(nid => 1 for nid in node_ids),
        node_ids,
        node_group,
        node_type,
        group_coeffs,
        coeff_par,
        coeff_nir,
        0.15,
        0.30,
    )
    first = ArchimedLight.FirstOrderResult(
        Dict(nid => 1.0 for nid in node_ids),
        ArchimedLight.SpectralNodeValues(
            Dict(1 => 0.0, 2 => 0.0, 3 => 100.0),
            Dict(nid => 0.0 for nid in node_ids),
        ),
        Dict(nid => 1 for nid in node_ids),
    )
    options = ArchimedLight.LightOptions(
        scattering_max_iter=10,
        scattering_stop_ratio=0.01,
    )

    result = ArchimedLight.compute_scattering_band(graph, first, options; band="PAR")

    # Java's scene-wide stopping total includes sensor observations. They do
    # not re-emit (coefficient zero), but their observed share keeps the third
    # iteration above the previous iteration's stopping decision.
    @test result.iterations == 3
    @test result.converged
    @test isapprox(result.added_power_per_node[2], 8.1046875; atol=1e-12, rtol=1e-12)
    @test isapprox(result.added_power_per_node[1], 7.5421875; atol=1e-12, rtol=1e-12)
    @test isapprox(result.added_power_per_node[3], 0.5625; atol=1e-12, rtol=1e-12)
end

@testitem "Model IO scene write round-trip keeps attached attrs" tags=[:model_io, :fast] begin
    import MultiScaleTreeGraph

    include(joinpath(@__DIR__, "support.jl"))

    fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
    row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options))
    step = ArchimedLight.run_light_step(fixture.scene, fixture.models, row, fixture.options)
    ArchimedLight.attach_light_step!(fixture.scene, step; fields=[:incident_par_initial_energy])

    outdir = mktempdir()
    out_opf = joinpath(outdir, "scene.opf")
    ArchimedLight.write_scene(out_opf, fixture.scene)

    opf = ArchimedLight.read_scene(out_opf).mtg
    found_attr = Ref(false)
    MultiScaleTreeGraph.traverse!(opf) do node
        if haskey(node, :Ri_PAR_0_q)
            found_attr[] = true
            return false
        end
        return true
    end
    @test found_attr[]
end
