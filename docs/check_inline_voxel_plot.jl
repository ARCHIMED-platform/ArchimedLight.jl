"""Fail the documentation build if the interactive voxel figure was not exported safely."""
function check_inline_voxel_plot(build_directory)
    candidates = (
        joinpath(build_directory, "voxel_terrain.html"),
        joinpath(build_directory, "voxel_terrain", "index.html"),
    )
    page_path = findfirst(isfile, candidates)
    page_path === nothing && error(
        "could not find the generated voxel terrain documentation page in $build_directory",
    )
    html = read(candidates[page_path], String)

    required_fragments = (
        "archimed-voxelplot-docs",
        "Tree voxels",
        "Terrain surface",
        "Incoming and reflected rays",
        "Surface normals",
        "Front cutaway",
    )
    for fragment in required_fragments
        occursin(fragment, html) || error(
            "interactive voxel documentation is missing the fragment: $fragment",
        )
    end

    forbidden_fragments = (
        "voxel_terrain_3d.html",
        "<iframe",
        "ws://localhost",
        "wss://localhost",
        "ws://127.0.0.1",
        "wss://127.0.0.1",
    )
    for fragment in forbidden_fragments
        occursin(fragment, html) && error(
            "interactive voxel documentation contains a forbidden fragment: $fragment",
        )
    end
    return candidates[page_path]
end
