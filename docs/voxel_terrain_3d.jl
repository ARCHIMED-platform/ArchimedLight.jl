const _VOXEL_DOCS_SOIL_DRY = SoilOpticalProperties(
    par_reflectance=0.12,
    nir_reflectance=0.32,
)
const _VOXEL_DOCS_SOIL_MOIST = SoilOpticalProperties(
    par_reflectance=0.24,
    nir_reflectance=0.42,
)
const _VOXEL_DOCS_SOIL_LITTER = SoilOpticalProperties(
    par_reflectance=0.38,
    nir_reflectance=0.52,
)

function _voxel_docs_elevation(x, y)
    broad_slope = 0.055x + 0.025y
    central_mound = 0.52exp(-((x - 4.8)^2 + 0.7(y - 4.2)^2) / 5.5)
    shallow_hollow = 0.28exp(-((x - 2.0)^2 + (y - 6.0)^2) / 2.0)
    ripples = 0.13sinpi(x / 4) * cospi(y / 4)
    # Keep the full terrain shape, but below the first trunk voxel so the
    # terrain reads as the ground the tree emerges from rather than a surface
    # hiding its base.
    return 0.60 * (
        0.42 + broad_slope + central_mound - shallow_hollow + ripples
    )
end

function _voxel_docs_tree_pad(nx, ny, nz)
    pad = zeros(nx, ny, nz)
    tree_x, tree_y = 4.35, 4.65
    for k in 1:nz, j in 1:ny, i in 1:nx
        x, y, z = i - 0.5, j - 0.5, k - 0.5
        ground = _voxel_docs_elevation(x, y)

        # A one-voxel-wide lower column makes the tree silhouette explicit.
        trunk = (x - tree_x)^2 + (y - tree_y)^2 <= 0.42^2 &&
                z > ground + 0.10 && z <= 4.60

        # An ellipsoid with a deterministic, slightly ruffled boundary gives
        # the crown an organic outline without introducing randomness.
        crown_metric = ((x - tree_x) / 2.55)^2 +
                       ((y - tree_y) / 2.25)^2 +
                       ((z - 6.30) / 2.40)^2
        crown_irregularity = 0.10sin(1.7x + 0.7y + 0.4z) +
                             0.07cos(1.2y - 1.1z)
        boundary_hole = crown_metric > 0.62 && mod(5i + 3j + 7k, 17) == 0
        crown = crown_metric + crown_irregularity <= 1.0 && !boundary_hole

        crown_density = crown ?
            1.02 + 0.68clamp(1.0 - crown_metric, 0.0, 1.0) : 0.0
        trunk_density = trunk ? 0.34 : 0.0
        # Give the trunk precedence where it enters the crown, otherwise its
        # uppermost voxel inherits the green crown colour and shortens the
        # apparent stem.
        pad[i, j, k] = trunk ? trunk_density : crown_density
    end
    return pad
end

function _voxel_docs_scene()
    nx, ny, nz = 8, 8, 9
    grid = VoxelGrid((0, 0, 0), (8, 8, 9), _voxel_docs_tree_pad(nx, ny, nz))

    x = collect(range(0.0, 8.0; length=17))
    y = collect(range(0.0, 8.0; length=17))
    elevation = [_voxel_docs_elevation(xi, yi) for xi in x, yi in y]
    material_ids = [
        xi + yi < 7 ? 1 : (xi > 10 && yi > 8 ? 3 : 2) for
        xi in 1:(length(x) - 1), yi in 1:(length(y) - 1)
    ]
    terrain = HeightFieldTerrain(
        x,
        y,
        elevation,
        _VOXEL_DOCS_SOIL_DRY;
        material_ids=material_ids,
        materials=Dict(
            1 => _VOXEL_DOCS_SOIL_DRY,
            2 => _VOXEL_DOCS_SOIL_MOIST,
            3 => _VOXEL_DOCS_SOIL_LITTER,
        ),
    )

    upward = [
        (0.0, 0.0, 1.0),
        (0.8, 0.0, 0.6),
        (-0.8, 0.0, 0.6),
        (0.0, 0.8, 0.6),
        (0.0, -0.8, 0.6),
    ]
    directions = vcat(upward, [(-d[1], -d[2], -d[3]) for d in upward])
    quadrature = VoxelScatteringQuadrature(directions, fill(0.1, length(directions)))
    rays = [
        (origin=(2.0, 3.1, 9.0), direction=(0.12, 0.08, -1.0)),
        (origin=(4.5, 4.4, 9.0), direction=(-0.08, 0.12, -1.0)),
        (origin=(6.7, 6.3, 9.0), direction=(-0.14, -0.10, -1.0)),
    ]
    return grid, terrain, quadrature, rays
end

function _voxel_docs_layer_visibility!(session, checkbox, layers...)
    for layer in layers
        Bonito.onjs(session, checkbox.value, js"""async visible => {
            const meshes = await $(layer);
            for (const mesh of meshes) {
                mesh.visible = visible;
                if (mesh.plot_object) {
                    mesh.plot_object.plot_data.visible = visible;
                }
            }
        }""")
    end
    return nothing
end

function _voxel_docs_control(label, checkbox)
    return Bonito.DOM.label(
        checkbox,
        Bonito.DOM.span(label),
        style=Bonito.Styles(
            "display" => "flex",
            "align-items" => "center",
            "gap" => "0.55rem",
            "white-space" => "nowrap",
            "font-size" => "0.9rem",
            "cursor" => "pointer",
        ),
    )
end

"""Build the self-contained interactive voxel/terrain figure used by Documenter."""
function voxel_terrain_app()
    return Bonito.App() do session
        grid, terrain, quadrature, rays = _voxel_docs_scene()
        figure = WGLMakie.Figure(
            size=(960, 620),
            backgroundcolor=WGLMakie.RGBf(0.965, 0.975, 0.985),
        )
        axis = WGLMakie.LScene(
            figure[1, 1];
            show_axis=true,
            scenekw=(
                backgroundcolor=WGLMakie.RGBf(0.965, 0.975, 0.985),
                clear=true,
            ),
        )
        plot = voxelplot!(
            axis,
            grid,
            terrain;
            rays=rays,
            quadrature=quadrature,
            cutaway=:none,
            transparency=false,
            pad_threshold=0.2,
            # Preserve crown translucency in the colormap, but keep the brown
            # trunk opaque so the terrain cannot visually bleed through it.
            voxel_alpha=1.0,
            voxel_gap=0.10,
            voxel_colormap=[
                WGLMakie.RGBAf(0.43, 0.25, 0.10, 1.00),
                WGLMakie.RGBAf(0.56, 0.36, 0.12, 1.00),
                WGLMakie.RGBAf(0.57, 0.74, 0.27, 0.68),
                WGLMakie.RGBAf(0.27, 0.58, 0.19, 0.64),
                WGLMakie.RGBAf(0.07, 0.35, 0.12, 0.62),
            ],
            voxel_colorrange=(0.30, 1.70),
            lattice_color=(WGLMakie.RGBf(0.40, 0.49, 0.56), 0.18),
            lattice_linewidth=0.7,
            terrain_alpha=1.0,
            terrain_colorrange=(0.08, 0.42),
            terrain_wireframe_color=(:black, 0.14),
            show_normals=false,
        )
        WGLMakie.cam3d!(axis)
        WGLMakie.update_cam!(
            axis.scene,
            (12.0, -14.0, 12.0),
            (4.0, 4.0, 4.4),
            (0.0, 0.0, 1.0),
        )

        full_voxels,
        front_voxels,
        lattice,
        terrain_surface,
        terrain_wireframe,
        incoming_rays,
        terrain_hits,
        reflected_rays,
        surface_normals = plot.plots

        show_voxels = Bonito.Checkbox(true)
        show_lattice = Bonito.Checkbox(true)
        show_terrain = Bonito.Checkbox(true)
        show_rays = Bonito.Checkbox(true)
        show_normals = Bonito.Checkbox(false)
        front_cutaway = Bonito.Checkbox(false)

        _voxel_docs_layer_visibility!(session, show_lattice, lattice)
        _voxel_docs_layer_visibility!(
            session,
            show_terrain,
            terrain_surface,
            terrain_wireframe,
        )
        _voxel_docs_layer_visibility!(
            session,
            show_rays,
            incoming_rays,
            terrain_hits,
            reflected_rays,
        )
        _voxel_docs_layer_visibility!(session, show_normals, surface_normals)

        Bonito.evaljs(session, js"""
        (() => {
            const show = $(show_voxels.value);
            const cutaway = $(front_cutaway.value);
            const fullPromise = $(full_voxels);
            const frontPromise = $(front_voxels);

            function setVisible(meshes, visible) {
                for (const mesh of meshes) {
                    mesh.visible = visible;
                    if (mesh.plot_object) {
                        mesh.plot_object.plot_data.visible = visible;
                    }
                }
            }

            async function synchronizeVoxels() {
                const [full, front] = await Promise.all([fullPromise, frontPromise]);
                setVisible(full, show.value && !cutaway.value);
                setVisible(front, show.value && cutaway.value);
            }

            Bonito.onany([show, cutaway], synchronizeVoxels);
            synchronizeVoxels();
        })();
        """)

        controls = Bonito.DOM.div(
            _voxel_docs_control("Tree voxels", show_voxels),
            _voxel_docs_control("Voxel lattice", show_lattice),
            _voxel_docs_control("Terrain surface", show_terrain),
            _voxel_docs_control("Incoming and reflected rays", show_rays),
            _voxel_docs_control("Surface normals", show_normals),
            _voxel_docs_control("Front cutaway", front_cutaway);
            style=Bonito.Styles(
                "display" => "flex",
                "flex-wrap" => "wrap",
                "gap" => "0.65rem 1.15rem",
                "padding" => "0.8rem 1rem",
                "border-bottom" => "1px solid rgba(90, 110, 130, 0.25)",
                "background" => "rgba(248, 250, 252, 0.96)",
            ),
        )
        legend = Bonito.DOM.div(
            "Brown: stem scaffold · Green: foliage PAD · Earth tones: terrain · Gold: incoming path · Red: soil hit · Purple/orange: Lambertian weights · Blue: local normal";
            style=Bonito.Styles(
                "padding" => "0.55rem 1rem 0.7rem",
                "font-size" => "0.82rem",
                "color" => "#475569",
                "background" => "rgba(248, 250, 252, 0.96)",
            ),
        )
        canvas = Bonito.DOM.div(
            WGLMakie.WithConfig(figure; resize_to=:parent);
            style=Bonito.Styles(
                "width" => "100%",
                "height" => "clamp(440px, 67vw, 640px)",
                "min-height" => "440px",
            ),
        )
        return Bonito.DOM.div(
            controls,
            canvas,
            legend;
            id="archimed-voxelplot-docs",
            style=Bonito.Styles(
                "width" => "100%",
                "overflow" => "hidden",
                "border" => "1px solid rgba(90, 110, 130, 0.28)",
                "border-radius" => "0.55rem",
                "background" => "#f8fafc",
            ),
        )
    end
end
