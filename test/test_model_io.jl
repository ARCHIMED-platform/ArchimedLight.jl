import MultiScaleTreeGraph

@testset "Model IO" begin
    @testset "YAML and typed models agree" begin
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

    @testset "Missing NIR optical property falls back to runtime defaults" begin
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

    @testset "Scene write round-trip keeps attached attrs" begin
        fixture = load_fixture_inputs(joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input"))
        row = first(ArchimedLight.prepare_meteo(fixture.meteo, fixture.options).rows)
        step = ArchimedLight.run_light_step(fixture.scene, fixture.models, row, fixture.options)
        ArchimedLight.attach_light_step!(fixture.scene, step; fields=[:incident_par_initial_energy])

        outdir = mktempdir()
        out_opf = joinpath(outdir, "scene.opf")
        ArchimedLight.write_scene(out_opf, fixture.scene)

        opf = ArchimedLight.read_scene(out_opf).mtg
        found_attr = false
        MultiScaleTreeGraph.traverse!(opf) do node
            if haskey(node, :Ri_PAR_0_q)
                found_attr = true
                return false
            end
            return true
        end
        @test found_attr
    end
end
