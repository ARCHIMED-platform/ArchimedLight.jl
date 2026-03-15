using OrderedCollections: OrderedDict
import MultiScaleTreeGraph

@testset "Model IO" begin
    @testset "Ordered YAML round-trip" begin
        tmp = mktempdir()
        mkpath(joinpath(tmp, "models"))

        cfg_path = joinpath(tmp, "config.yml")
        model_path = joinpath(tmp, "models", "plant.yml")
        meteo_path = joinpath(tmp, "meteo.csv")

        write(
            cfg_path,
            join(
                [
                    "scene: scene/exported.ops",
                    "models:",
                    "  - models/plant.yml",
                    "meteo: meteo.csv",
                    "pixel_size: 1.0",
                    "scattering: false",
                ],
                "\n",
            ) * "\n",
        )
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
                ],
                "\n",
            ) * "\n",
        )
        write(meteo_path, "date;hour_start;hour_end;Rg;Tac;latitude;longitude;altitude\n2016/07/01;08:00:00;09:00:00;0;25;0;0;0\n")

        cfg = ArchimedLight.read_light_config(cfg_path)
        @test collect(keys(cfg.raw))[1:5] == ["scene", "models", "meteo", "pixel_size", "scattering"]

        ArchimedLight.set_parameter!(cfg, 2.5, "config", "pixel_size")
        ArchimedLight.set_parameter!(cfg, 0.25, "models", 1, "Type", "Leaf", "Interception", "transparency")

        @test isapprox(cfg.pixel_size, 0.025; atol=1e-12)

        outdir = joinpath(tmp, "written")
        out_cfg = ArchimedLight.write_light_inputs(outdir, cfg; scene_rel="scene/exported.ops")
        out_cfg_txt = read(out_cfg, String)
        out_model_txt = read(joinpath(outdir, "models", "plant.yml"), String)

        @test findfirst("scene:", out_cfg_txt) < findfirst("models:", out_cfg_txt) < findfirst("meteo:", out_cfg_txt)
        @test occursin("pixel_size: 2.5", out_cfg_txt)
        @test findfirst("Group:", out_model_txt) < findfirst("Type:", out_model_txt)
        @test occursin("transparency: 0.25", out_model_txt)
    end

    @testset "Scene export writes OPF attrs" begin
        case_root = joinpath(@__DIR__, "fast_fixtures", "simpleplant_16_notoric", "input")
        cfg = ArchimedLight.read_light_config(joinpath(case_root, "config.yml"))
        ArchimedLight.set_parameter!(cfg, true, "config", "export_ops")
        ArchimedLight.set_parameter!(cfg, OrderedDict("Ri_PAR_0_q" => true), "config", "opf_variables")

        scene = ArchimedLight.read_scene(cfg.scene)
        meteo = ArchimedLight.read_meteo(cfg.meteo)
        selected = ArchimedLight.prepare_meteo(meteo, cfg)
        step = only(ArchimedLight.run_light_series(scene, meteo, cfg))

        outdir = mktempdir()
        files = ArchimedLight.write_light_outputs(
            scene,
            step,
            cfg;
            meteo_row=only(selected.rows),
            outdir=outdir,
            write_component=false,
            write_scene=false,
            write_summary=false,
            write_sun_position_log=false,
            write_scattering_log=false,
        )

        @test haskey(files, "ops_step_00001")
        @test haskey(files, "config")

        exported_opf = joinpath(outdir, "opf", "simple_OPF_shapes.opf")
        @test isfile(exported_opf)

        opf = PlantGeom.read_opf(exported_opf, attr_type=Dict)
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
