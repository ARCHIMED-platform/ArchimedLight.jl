function _voxel_header_triplet(lines, key::AbstractString, parser)
    expected = lowercase(key)
    for line in lines
        stripped = strip(line)
        startswith(stripped, "#") || continue
        words = split(replace(stripped, ':' => ' '))
        isempty(words) && continue
        lowercase(words[1]) == expected || continue
        length(words) >= 4 || throw(ArgumentError("invalid $key voxel header"))
        return ntuple(i -> parser(words[i + 1]), 3)
    end
    throw(ArgumentError("missing $key voxel header"))
end

"""
    read_voxel_grid(path)::VoxelGrid

Read a historical ARCHIMED `VOXEL SPACE` file. The grid must contain exactly
one row for every 0-based `(i, j, k)` coordinate. `PAD` and `PadBVTotal` are
accepted as plant-area-density column names; additional columns are ignored.
"""
function read_voxel_grid(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty voxel file: $path"))
    strip(lines[1]) == "VOXEL SPACE" ||
        throw(ArgumentError("invalid voxel file: missing VOXEL SPACE identifier"))

    minimum = _voxel_header_triplet(lines, "#min_corner", x -> parse(Float64, x))
    maximum = _voxel_header_triplet(lines, "#max_corner", x -> parse(Float64, x))
    dimensions = _voxel_header_triplet(lines, "#split", x -> parse(Int, x))
    all(>(0), dimensions) || throw(ArgumentError("voxel split values must be positive"))

    header_index = findfirst(eachindex(lines)) do index
        index == 1 && return false
        stripped = strip(lines[index])
        !isempty(stripped) && !startswith(stripped, "#")
    end
    header_index === nothing && throw(ArgumentError("missing voxel column header"))
    header = split(strip(lines[header_index]))
    length(header) >= 4 || throw(ArgumentError("invalid voxel column header"))
    lowercase.(header[1:3]) == ["i", "j", "k"] ||
        throw(ArgumentError("voxel columns must start with i j k"))

    pad_column = findfirst(name -> lowercase(name) in ("pad", "padbvtotal"), header)
    pad_column === nothing && throw(ArgumentError("voxel file has no PAD or PadBVTotal column"))

    pad = zeros(Float64, dimensions...)
    seen = falses(dimensions...)
    for line_index in (header_index + 1):length(lines)
        stripped = strip(lines[line_index])
        (isempty(stripped) || startswith(stripped, "#")) && continue
        words = split(stripped)
        length(words) >= length(header) ||
            throw(ArgumentError("voxel row $line_index has too few columns"))
        i = parse(Int, words[1]) + 1
        j = parse(Int, words[2]) + 1
        k = parse(Int, words[3]) + 1
        checkbounds(Bool, pad, i, j, k) ||
            throw(ArgumentError("voxel row $line_index has out-of-bounds index"))
        seen[i, j, k] && throw(ArgumentError("duplicate voxel index at row $line_index"))
        pad[i, j, k] = parse(Float64, words[pad_column])
        seen[i, j, k] = true
    end
    all(seen) || throw(ArgumentError("voxel file does not define every grid cell"))
    return VoxelGrid(minimum, maximum, pad)
end

function _check_voxel_result_size(grid::VoxelGrid, result::VoxelFirstOrderResult)
    expected = size(grid)
    for values in (result.par_energy, result.nir_energy, result.par_flux, result.nir_flux)
        size(values) == expected || throw(DimensionMismatch("voxel result does not match grid dimensions"))
    end
    return nothing
end

"""
    write_voxel_grid(path, grid, result=nothing)

Write a regular grid in historical `.vox` layout. When `result` is supplied,
the file also contains first-order PAR/NIR energy (`J voxel⁻¹ step⁻¹`) and
plant-area-normalized flux (`W m⁻² plant`).
"""
function write_voxel_grid(
    path::AbstractString,
    grid::VoxelGrid,
    result::Union{Nothing,VoxelFirstOrderResult}=nothing,
)
    result === nothing || _check_voxel_result_size(grid, result)
    open(path, "w") do io
        println(io, "VOXEL SPACE")
        println(io, "#min_corner: ", join(grid.minimum, ' '))
        println(io, "#max_corner: ", join(grid.maximum, ' '))
        println(io, "#split: ", join(size(grid), ' '))
        println(io, "#LAD_TYPE: Spherical")
        if result === nothing
            println(io, "i j k PAD")
        else
            println(io, "i j k PAD Ri_PAR_0_q Ri_NIR_0_q Ri_PAR_0_f Ri_NIR_0_f")
        end
        nx, ny, nz = size(grid)
        for i in 1:nx, j in 1:ny, k in 1:nz
            values = Any[i - 1, j - 1, k - 1, grid.pad[i, j, k]]
            if result !== nothing
                append!(values, (
                    result.par_energy[i, j, k],
                    result.nir_energy[i, j, k],
                    result.par_flux[i, j, k],
                    result.nir_flux[i, j, k],
                ))
            end
            println(io, join(values, ' '))
        end
    end
    return path
end
