#!/usr/bin/env julia

using CairoMakie

const REPO_ROOT = dirname(@__DIR__)
const DEFAULT_INPUT =
    joinpath(REPO_ROOT, "benchmark", "results", "artifact_scene_matrix_cpu.tsv")
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "docs", "src", "assets")
const DEFAULT_SUMMARY =
    joinpath(REPO_ROOT, "benchmark", "results", "artifact_scene_matrix_summary.md")
const OVERVIEW_FILENAME = "benchmark_scene_matrix_overview.png"
const EFFECTS_FILENAME = "benchmark_scene_matrix_effects.png"

const EXPECTED_SCENES = [
    "simple_plant",
    "coffee",
    "wheat",
    "agrivoltaics_wheat",
    "elaeis",
    "elaeis_two_plants",
]
const EXPECTED_DIRECTIONS = [16, 46, 136]
const EXPECTED_PIXELS_CM = [0.1, 0.5, 1.0, 10.0]
const EXPECTED_BOOL = [false, true]
const EXPECTED_STEPS = [1, 12, 24, 48]
const PLANNED_ROWS_PER_SCENE =
    length(EXPECTED_DIRECTIONS) *
    length(EXPECTED_PIXELS_CM) *
    length(EXPECTED_BOOL)^2 *
    length(EXPECTED_STEPS)

const SCENE_LABELS = Dict(
    "simple_plant" => "Simple plant",
    "coffee" => "Coffee",
    "wheat" => "Wheat",
    "agrivoltaics_wheat" => "Agrivoltaics + wheat",
    "elaeis" => "Oil palm",
    "elaeis_two_plants" => "Two oil palms",
)
const SCENE_COLORS = Dict(
    "simple_plant" => :dodgerblue3,
    "coffee" => :darkorange2,
    "wheat" => :seagreen3,
    "agrivoltaics_wheat" => :mediumpurple3,
    "elaeis" => :firebrick3,
    "elaeis_two_plants" => :slateblue3,
)

struct BenchmarkRow
    scene::String
    backend::String
    directions::Int
    pixel_size_cm::Float64
    scattering::Bool
    cache_radiation::Bool
    steps::Int
    samples::Int
    nodes::Int
    faces::Int
    objects::Int
    domain_x::Float64
    domain_y::Float64
    pixel_cells_estimate::Int
    status::String
    median_ms::Float64
    min_ms::Float64
    max_ms::Float64
    median_alloc_mb::Float64
    resolved_backend::String
    geometry_mode::String
end

function usage()
    println(
        """
        Visualize the CPU artifact-scene benchmark matrix.

        Usage:
          julia --project=benchmark benchmark/visualize_artifact_scene_matrix.jl [options]

        Options:
          --input=<path>       Input TSV (default: benchmark/results/artifact_scene_matrix_cpu.tsv)
          --output-dir=<path>  Figure directory (default: docs/src/assets)
          --summary=<path>     Generated Markdown summary path
          --no-summary         Do not write the Markdown summary file
          -h, --help           Show this help
        """,
    )
end

function parse_args(args)
    input = DEFAULT_INPUT
    output_dir = DEFAULT_OUTPUT_DIR
    summary = DEFAULT_SUMMARY
    write_summary = true
    for arg in args
        if arg in ("-h", "--help")
            usage()
            return nothing
        elseif startswith(arg, "--input=")
            input = normpath(split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = normpath(split(arg, '='; limit=2)[2])
        elseif startswith(arg, "--summary=")
            summary = normpath(split(arg, '='; limit=2)[2])
        elseif arg == "--no-summary"
            write_summary = false
        else
            error("Unknown argument: $(repr(arg)). Run with --help for usage.")
        end
    end
    return (
        input=input,
        output_dir=output_dir,
        summary=summary,
        write_summary=write_summary,
    )
end

function _parse_bool(value::AbstractString, column::AbstractString, line_number::Int)
    value == "true" && return true
    value == "false" && return false
    error("Expected true or false in column $column at line $line_number, got $(repr(value)).")
end

function _required_columns()
    return [
        "scene",
        "backend",
        "directions",
        "pixel_size_cm",
        "scattering",
        "cache_radiation",
        "steps",
        "samples",
        "nodes",
        "faces",
        "objects",
        "domain_x",
        "domain_y",
        "pixel_cells_estimate",
        "status",
        "median_ms",
        "min_ms",
        "max_ms",
        "median_alloc_mb",
        "resolved_backend",
        "geometry_mode",
    ]
end

function read_results(path::AbstractString)
    isfile(path) || error("Benchmark TSV does not exist: $path")
    rows = BenchmarkRow[]
    open(path, "r") do io
        eof(io) && error("Benchmark TSV is empty: $path")
        columns = split(readline(io), '\t'; keepempty=true)
        column_index = Dict(name => i for (i, name) in pairs(columns))
        missing_columns = setdiff(_required_columns(), collect(keys(column_index)))
        isempty(missing_columns) ||
            error("Benchmark TSV is missing columns: $(join(sort!(missing_columns), ", ")).")

        for (offset, line) in enumerate(eachline(io))
            line_number = offset + 1
            isempty(line) && continue
            fields = split(line, '\t'; keepempty=true)
            length(fields) == length(columns) || error(
                "Line $line_number has $(length(fields)) fields; expected $(length(columns)).",
            )
            value = name -> fields[column_index[name]]
            push!(
                rows,
                BenchmarkRow(
                    String(value("scene")),
                    String(value("backend")),
                    parse(Int, value("directions")),
                    parse(Float64, value("pixel_size_cm")),
                    _parse_bool(value("scattering"), "scattering", line_number),
                    _parse_bool(value("cache_radiation"), "cache_radiation", line_number),
                    parse(Int, value("steps")),
                    parse(Int, value("samples")),
                    parse(Int, value("nodes")),
                    parse(Int, value("faces")),
                    parse(Int, value("objects")),
                    parse(Float64, value("domain_x")),
                    parse(Float64, value("domain_y")),
                    parse(Int, value("pixel_cells_estimate")),
                    String(value("status")),
                    parse(Float64, value("median_ms")),
                    parse(Float64, value("min_ms")),
                    parse(Float64, value("max_ms")),
                    parse(Float64, value("median_alloc_mb")),
                    String(value("resolved_backend")),
                    String(value("geometry_mode")),
                ),
            )
        end
    end
    return rows
end

function validate_results(rows)
    isempty(rows) && error("Benchmark TSV contains no data rows.")
    backends = unique(String[row.backend for row in rows])
    backends == ["normal_cpu"] ||
        error("Expected only the normal_cpu backend, found: $(join(sort!(backends), ", ")).")

    successful = filter(row -> row.status == "ok", rows)
    isempty(successful) && error("Benchmark TSV contains no successful rows.")
    all(row -> row.samples > 0, successful) ||
        error("Successful rows must report a positive sample count.")
    all(row -> isfinite(row.median_ms) && row.median_ms > 0.0, successful) ||
        error("Successful rows must have positive finite timings.")
    all(row -> isfinite(row.median_alloc_mb) && row.median_alloc_mb > 0.0, successful) ||
        error("Successful rows must have positive finite allocation measurements.")
    all(row -> row.resolved_backend == "RasterCPUBackend", successful) ||
        error("Successful rows must resolve to RasterCPUBackend.")
    return rows
end

_is_success(row::BenchmarkRow) =
    row.status == "ok" &&
    isfinite(row.median_ms) &&
    row.median_ms > 0.0 &&
    isfinite(row.median_alloc_mb) &&
    row.median_alloc_mb > 0.0

function _observed_scenes(rows)
    return [scene for scene in EXPECTED_SCENES if any(row -> row.scene == scene, rows)]
end

function _quantile(values, probability::Float64)
    isempty(values) && error("Cannot summarize an empty collection.")
    0.0 <= probability <= 1.0 || error("Probability must be in [0, 1].")
    sorted_values = sort!(Float64[Float64(value) for value in values])
    index = 1.0 + (length(sorted_values) - 1) * probability
    lower = floor(Int, index)
    upper = ceil(Int, index)
    lower == upper && return sorted_values[lower]
    weight = index - lower
    return (1.0 - weight) * sorted_values[lower] + weight * sorted_values[upper]
end

function ratio_summary(effects)
    ratios = Float64[effect.ratio for effect in effects if isfinite(effect.ratio)]
    isempty(ratios) && error("No matched ratios are available for this comparison.")
    return (
        median=_quantile(ratios, 0.5),
        q25=_quantile(ratios, 0.25),
        q75=_quantile(ratios, 0.75),
        n=length(ratios),
    )
end

function _pair_key(row::BenchmarkRow, varied::Symbol)
    varied == :cache_radiation &&
        return (
            row.scene,
            row.directions,
            row.pixel_size_cm,
            row.scattering,
            row.steps,
        )
    varied == :scattering &&
        return (
            row.scene,
            row.directions,
            row.pixel_size_cm,
            row.cache_radiation,
            row.steps,
        )
    varied == :directions &&
        return (
            row.scene,
            row.pixel_size_cm,
            row.scattering,
            row.cache_radiation,
            row.steps,
        )
    varied == :pixel_size_cm &&
        return (
            row.scene,
            row.directions,
            row.scattering,
            row.cache_radiation,
            row.steps,
        )
    varied == :steps &&
        return (
            row.scene,
            row.directions,
            row.pixel_size_cm,
            row.scattering,
            row.cache_radiation,
        )
    error("Unsupported paired comparison field: $varied")
end

function _metric(row::BenchmarkRow, metric::Symbol)
    metric == :median_ms && return row.median_ms
    metric == :median_alloc_mb && return row.median_alloc_mb
    error("Unsupported metric: $metric")
end

function paired_ratios(
    rows,
    varied::Symbol,
    numerator_value,
    denominator_value;
    metric::Symbol=:median_ms,
)
    groups = Dict{Any,Dict{Any,BenchmarkRow}}()
    for row in rows
        _is_success(row) || continue
        grouped_values = get!(groups, _pair_key(row, varied), Dict{Any,BenchmarkRow}())
        field_value = getfield(row, varied)
        haskey(grouped_values, field_value) &&
            error("Duplicate successful benchmark configuration for $varied=$(field_value).")
        grouped_values[field_value] = row
    end

    effects = NamedTuple[]
    for grouped_values in values(groups)
        haskey(grouped_values, numerator_value) || continue
        haskey(grouped_values, denominator_value) || continue
        numerator = grouped_values[numerator_value]
        denominator = grouped_values[denominator_value]
        ratio = _metric(numerator, metric) / _metric(denominator, metric)
        push!(
            effects,
            (
                scene=numerator.scene,
                directions=numerator.directions,
                pixel_size_cm=numerator.pixel_size_cm,
                scattering=numerator.scattering,
                cache_radiation=numerator.cache_radiation,
                steps=numerator.steps,
                ratio=ratio,
            ),
        )
    end
    return effects
end

function coverage_table(rows)
    return [
        begin
            scene_rows = filter(row -> row.scene == scene, rows)
            ok = count(row -> row.status == "ok", scene_rows)
            skipped = count(row -> row.status == "skipped_huge_case", scene_rows)
            (
                scene=scene,
                rows=length(scene_rows),
                ok=ok,
                skipped=skipped,
                missing=max(PLANNED_ROWS_PER_SCENE - length(scene_rows), 0),
            )
        end for scene in EXPECTED_SCENES
    ]
end

function _find_row(rows; kwargs...)
    index = findfirst(rows) do row
        all(getfield(row, key) == value for (key, value) in kwargs)
    end
    index === nothing && return nothing
    return rows[index]
end

function scene_table(rows)
    summaries = NamedTuple[]
    for scene in _observed_scenes(rows)
        scene_rows = filter(row -> row.scene == scene && _is_success(row), rows)
        isempty(scene_rows) && continue
        representative = _find_row(
            scene_rows;
            directions=16,
            pixel_size_cm=1.0,
            scattering=false,
            cache_radiation=false,
            steps=1,
        )
        fine_grid = _find_row(scene_rows; pixel_size_cm=0.1)
        one_cm_grid = _find_row(scene_rows; pixel_size_cm=1.0)
        representative === nothing && continue
        fine_grid === nothing && continue
        one_cm_grid === nothing && continue
        push!(
            summaries,
            (
                scene=scene,
                nodes=representative.nodes,
                faces=representative.faces,
                domain_x=representative.domain_x,
                domain_y=representative.domain_y,
                cells_01cm=fine_grid.pixel_cells_estimate,
                cells_1cm=one_cm_grid.pixel_cells_estimate,
                representative_ms=representative.median_ms,
            ),
        )
    end
    return summaries
end

function detect_step_anomalies(rows; factor::Float64=20.0)
    groups = Dict{Any,Vector{BenchmarkRow}}()
    for row in rows
        _is_success(row) || continue
        key = (
            row.scene,
            row.directions,
            row.pixel_size_cm,
            row.scattering,
            row.cache_radiation,
        )
        push!(get!(groups, key, BenchmarkRow[]), row)
    end

    anomalies = Dict{Any,NamedTuple}()
    for group in values(groups)
        sort!(group; by=row -> row.steps)
        for i in 1:(length(group) - 1)
            earlier = group[i]
            later = group[i + 1]
            step_ratio = later.steps / earlier.steps
            runtime_ratio = later.median_ms / earlier.median_ms
            if runtime_ratio > factor * step_ratio
                key = (
                    later.scene,
                    later.directions,
                    later.pixel_size_cm,
                    later.scattering,
                    later.cache_radiation,
                    later.steps,
                )
                anomalies[key] = (
                    row=later,
                    reason="runtime increased $(round(runtime_ratio; sigdigits=4))× for $(round(step_ratio; sigdigits=4))× more steps",
                )
            elseif earlier.median_ms / later.median_ms > factor
                drop_ratio = earlier.median_ms / later.median_ms
                key = (
                    earlier.scene,
                    earlier.directions,
                    earlier.pixel_size_cm,
                    earlier.scattering,
                    earlier.cache_radiation,
                    earlier.steps,
                )
                anomalies[key] = (
                    row=earlier,
                    reason="runtime fell $(round(drop_ratio; sigdigits=4))× when the step count increased",
                )
            end
        end
    end
    return sort!(
        collect(values(anomalies));
        by=item -> (
            item.row.scene,
            item.row.directions,
            item.row.pixel_size_cm,
            item.row.steps,
        ),
    )
end

_fmt(value::Integer) = string(value)
_fmt(value::AbstractFloat; sigdigits::Int=5) =
    isfinite(value) ? string(round(value; sigdigits=sigdigits)) : string(value)

function render_summary(rows, input_path::AbstractString)
    io = IOBuffer()
    println(io, "# Artifact scene benchmark summary")
    println(io)
    println(
        io,
        "Generated by benchmark/visualize_artifact_scene_matrix.jl from ",
        basename(input_path),
        ".",
    )
    println(io)

    println(io, "## Coverage")
    println(io)
    println(io, "| Scene | Rows | OK | Skipped | Missing from 192 |")
    println(io, "| --- | ---: | ---: | ---: | ---: |")
    for item in coverage_table(rows)
        println(
            io,
            "| ",
            item.scene,
            " | ",
            item.rows,
            " | ",
            item.ok,
            " | ",
            item.skipped,
            " | ",
            item.missing,
            " |",
        )
    end
    println(io)

    println(io, "## Recorded scene complexity")
    println(io)
    println(
        io,
        "| Scene | Nodes | Faces | Domain (m) | Estimated cells at 0.1 cm | Estimated cells at 1 cm | Representative runtime (ms) |",
    )
    println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for item in scene_table(rows)
        println(
            io,
            "| ",
            get(SCENE_LABELS, item.scene, item.scene),
            " | ",
            item.nodes,
            " | ",
            item.faces,
            " | ",
            _fmt(item.domain_x),
            " × ",
            _fmt(item.domain_y),
            " | ",
            item.cells_01cm,
            " | ",
            item.cells_1cm,
            " | ",
            _fmt(item.representative_ms),
            " |",
        )
    end
    println(io)
    println(
        io,
        "Representative runtime uses 16 directions, 1 cm pixels, one step, scattering off, and cache off.",
    )
    println(io)

    cache_effects = paired_ratios(rows, :cache_radiation, false, true)
    cache_alloc_effects =
        paired_ratios(rows, :cache_radiation, false, true; metric=:median_alloc_mb)
    scattering_effects = paired_ratios(rows, :scattering, true, false)
    direction_46 = paired_ratios(rows, :directions, 46, 16)
    direction_136 = paired_ratios(rows, :directions, 136, 16)

    println(io, "## Matched effects")
    println(io)
    println(io, "| Effect | Group | Median ratio | IQR | Matched pairs |")
    println(io, "| --- | --- | ---: | ---: | ---: |")
    for step in EXPECTED_STEPS
        summary = ratio_summary(filter(effect -> effect.steps == step, cache_effects))
        println(
            io,
            "| Cache speedup (uncached / cached) | ",
            step,
            " step",
            step == 1 ? "" : "s",
            " | ",
            _fmt(summary.median),
            "× | ",
            _fmt(summary.q25),
            "–",
            _fmt(summary.q75),
            "× | ",
            summary.n,
            " |",
        )
    end
    for cache_state in (false, true)
        summary = ratio_summary(
            filter(effect -> effect.cache_radiation == cache_state, scattering_effects),
        )
        println(
            io,
            "| Scattering overhead (on / off) | cache ",
            cache_state ? "on" : "off",
            " | ",
            _fmt(summary.median),
            "× | ",
            _fmt(summary.q25),
            "–",
            _fmt(summary.q75),
            "× | ",
            summary.n,
            " |",
        )
    end
    for (label, effects) in (("46 / 16", direction_46), ("136 / 16", direction_136))
        summary = ratio_summary(effects)
        println(
            io,
            "| Direction cost | ",
            label,
            " | ",
            _fmt(summary.median),
            "× | ",
            _fmt(summary.q25),
            "–",
            _fmt(summary.q75),
            "× | ",
            summary.n,
            " |",
        )
    end
    println(io)

    println(io, "## Cache allocation reduction")
    println(io)
    println(io, "| Steps | Median uncached / cached allocation | IQR | Matched pairs |")
    println(io, "| ---: | ---: | ---: | ---: |")
    for step in filter(!=(1), EXPECTED_STEPS)
        summary = ratio_summary(filter(effect -> effect.steps == step, cache_alloc_effects))
        println(
            io,
            "| ",
            step,
            " | ",
            _fmt(summary.median),
            "× | ",
            _fmt(summary.q25),
            "–",
            _fmt(summary.q75),
            "× | ",
            summary.n,
            " |",
        )
    end
    println(io)

    anomalies = detect_step_anomalies(rows)
    println(io, "## Timing anomalies")
    println(io)
    if isempty(anomalies)
        println(io, "No obvious cross-step timing anomalies were detected.")
    else
        println(io, "| Scene | Directions | Pixel (cm) | Scattering | Cache | Steps | Runtime (ms) | Reason |")
        println(io, "| --- | ---: | ---: | --- | --- | ---: | ---: | --- |")
        for item in anomalies
            row = item.row
            println(
                io,
                "| ",
                row.scene,
                " | ",
                row.directions,
                " | ",
                _fmt(row.pixel_size_cm),
                " | ",
                row.scattering,
                " | ",
                row.cache_radiation,
                " | ",
                row.steps,
                " | ",
                _fmt(row.median_ms),
                " | ",
                item.reason,
                " |",
            )
        end
    end
    println(io)
    println(
        io,
        "Ratios compare configurations that differ in exactly one parameter. IQR describes variation across configurations, not repeated-run uncertainty.",
    )
    return String(take!(io))
end

function _controlled_overview_rows(rows)
    selected = filter(rows) do row
        _is_success(row) &&
            row.directions == 16 &&
            row.steps == 12 &&
            !row.scattering &&
            !row.cache_radiation
    end
    isempty(selected) && error("No rows match the controlled overview slice.")
    return selected
end

function _positive_limits(values; padding::Float64=1.25)
    minimum_value = minimum(values)
    maximum_value = maximum(values)
    minimum_value > 0.0 || error("Log-scale values must be positive.")
    return (minimum_value / padding, maximum_value * padding)
end

function plot_overview(rows, output_path::AbstractString)
    selected = _controlled_overview_rows(rows)
    scenes = _observed_scenes(selected)
    CairoMakie.activate!(type="png")
    figure = Figure(size=(1560, 900), backgroundcolor=:white)
    Label(
        figure[0, 1:length(scenes)],
        "Controlled CPU cost: 16 directions, 12 steps, scattering off, cache off";
        fontsize=24,
        font=:bold,
    )

    all_x = Float64[row.pixel_cells_estimate for row in selected]
    all_runtime = Float64[row.median_ms for row in selected]
    all_allocation = Float64[row.median_alloc_mb for row in selected]
    x_limits = _positive_limits(all_x)
    runtime_limits = _positive_limits(all_runtime; padding=1.5)
    allocation_limits = _positive_limits(all_allocation; padding=1.5)

    for (column, scene) in pairs(scenes)
        scene_rows = sort!(
            filter(row -> row.scene == scene, selected);
            by=row -> row.pixel_cells_estimate,
        )
        length(scene_rows) == length(EXPECTED_PIXELS_CM) ||
            @warn "Controlled overview slice is incomplete" scene rows=length(scene_rows)
        color = get(SCENE_COLORS, scene, :steelblue)
        x = Float64[row.pixel_cells_estimate for row in scene_rows]
        runtime = Float64[row.median_ms for row in scene_rows]
        allocation = Float64[row.median_alloc_mb for row in scene_rows]

        runtime_axis = Axis(
            figure[1, column];
            title=get(SCENE_LABELS, scene, scene),
            xscale=log10,
            yscale=log10,
            ylabel=column == 1 ? "Runtime (ms)" : "",
        )
        lines!(runtime_axis, x, runtime; color, linewidth=3)
        scatter!(runtime_axis, x, runtime; color, markersize=12)
        xlims!(runtime_axis, x_limits...)
        ylims!(runtime_axis, runtime_limits...)

        allocation_axis = Axis(
            figure[2, column];
            xscale=log10,
            yscale=log10,
            xlabel="Estimated horizontal cells",
            ylabel=column == 1 ? "Allocated (MiB)" : "",
        )
        lines!(allocation_axis, x, allocation; color, linewidth=3)
        scatter!(allocation_axis, x, allocation; color, markersize=12)
        xlims!(allocation_axis, x_limits...)
        ylims!(allocation_axis, allocation_limits...)
    end

    mkpath(dirname(output_path))
    save(output_path, figure; px_per_unit=1.5)
    return output_path
end

function _interval_panel!(
    position,
    summaries,
    labels;
    title::AbstractString,
    ylabel::AbstractString="Runtime ratio",
    colors=fill(:steelblue, length(summaries)),
)
    length(summaries) == length(labels) ||
        error("Each interval summary needs one category label.")
    axis = Axis(
        position;
        title,
        ylabel,
        yscale=log10,
        xticks=(collect(eachindex(labels)), labels),
        xticklabelrotation=pi / 10,
    )
    hlines!(axis, [1.0]; color=:gray55, linestyle=:dash, linewidth=1.5)
    for (index, summary) in pairs(summaries)
        color = colors[index]
        lines!(
            axis,
            [index, index],
            [summary.q25, summary.q75];
            color,
            linewidth=6,
        )
        scatter!(
            axis,
            [index],
            [summary.median];
            color,
            markersize=15,
            strokecolor=:white,
            strokewidth=1.5,
        )
    end
    lower = min(0.8, minimum(summary.q25 for summary in summaries) / 1.35)
    upper = max(1.25, maximum(summary.q75 for summary in summaries) * 1.35)
    ylims!(axis, lower, upper)
    return axis
end

function plot_effects(rows, output_path::AbstractString)
    cache_effects = paired_ratios(rows, :cache_radiation, false, true)
    scattering_effects = paired_ratios(rows, :scattering, true, false)
    direction_46 = paired_ratios(rows, :directions, 46, 16)
    direction_136 = paired_ratios(rows, :directions, 136, 16)
    fine_grid_effects = paired_ratios(rows, :pixel_size_cm, 0.1, 1.0)

    cache_summaries = [
        ratio_summary(filter(effect -> effect.steps == step, cache_effects)) for
        step in EXPECTED_STEPS
    ]
    scattering_summaries = NamedTuple[]
    for cache_state in (false, true)
        matched = filter(
            effect -> effect.cache_radiation == cache_state,
            scattering_effects,
        )
        summary = ratio_summary(matched)
        summary.median < 100.0 || error(
            "Implausible scattering ratio for cache=$(cache_state): $(summary.median)",
        )
        push!(scattering_summaries, summary)
    end
    direction_summaries = [ratio_summary(direction_46), ratio_summary(direction_136)]
    scenes = _observed_scenes(rows)
    fine_grid_summaries = [
        ratio_summary(filter(effect -> effect.scene == scene, fine_grid_effects)) for
        scene in scenes
    ]

    CairoMakie.activate!(type="png")
    figure = Figure(size=(1450, 1000), backgroundcolor=:white)
    Label(
        figure[0, 1:2],
        "Matched CPU parameter effects";
        fontsize=24,
        font=:bold,
    )
    _interval_panel!(
        figure[1, 1],
        cache_summaries,
        string.(EXPECTED_STEPS);
        title="Cache speedup by series length\n(uncached / cached)",
        colors=[:gray45, :dodgerblue3, :dodgerblue3, :dodgerblue3],
    )
    _interval_panel!(
        figure[1, 2],
        scattering_summaries,
        ["cache off", "cache on"];
        title="Scattering overhead\n(scattering on / off)",
        colors=[:darkorange2, :mediumpurple3],
    )
    _interval_panel!(
        figure[2, 1],
        direction_summaries,
        ["46 / 16", "136 / 16"];
        title="Directional-resolution cost",
        colors=[:seagreen3, :seagreen4],
    )
    _interval_panel!(
        figure[2, 2],
        fine_grid_summaries,
        [get(SCENE_LABELS, scene, scene) for scene in scenes];
        title="Fine-grid cost by scene\n(0.1 cm / 1 cm)",
        colors=[get(SCENE_COLORS, scene, :steelblue) for scene in scenes],
    )

    Label(
        figure[3, 1:2],
        "Point: median matched ratio    Thick line: interquartile range    Dashed line: no effect (1×)";
        fontsize=15,
        color=:gray35,
    )
    mkpath(dirname(output_path))
    save(output_path, figure; px_per_unit=1.5)
    return output_path
end

function main(args=ARGS)
    options = parse_args(args)
    options === nothing && return nothing
    rows = validate_results(read_results(options.input))
    summary = render_summary(rows, options.input)

    overview_path = joinpath(options.output_dir, OVERVIEW_FILENAME)
    effects_path = joinpath(options.output_dir, EFFECTS_FILENAME)
    plot_overview(rows, overview_path)
    plot_effects(rows, effects_path)

    if options.write_summary
        mkpath(dirname(options.summary))
        write(options.summary, summary)
    end

    print(summary)
    println("Wrote ", overview_path)
    println("Wrote ", effects_path)
    options.write_summary && println("Wrote ", options.summary)
    return (
        overview=overview_path,
        effects=effects_path,
        summary=options.write_summary ? options.summary : nothing,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
