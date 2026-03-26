#!/usr/bin/env julia

using Artifacts: artifact_hash, artifact_path
using Pkg.Artifacts: ensure_artifact_installed
using ArchimedLight
using CairoMakie
using CSV
using Tables

function _release_dataset_root()
    dataset_root = strip(get(ENV, "ARCHIMEDLIGHT_RELEASE_FIXTURES_DIR", ""))
    if isempty(dataset_root)
        repo_root = normpath(joinpath(@__DIR__, ".."))
        artifacts_toml = joinpath(repo_root, "Artifacts.toml")
        artifact_name = get(ENV, "ARCHIMEDLIGHT_RELEASE_ARTIFACT_NAME", "archimedlight-release-fixtures")
        if isfile(artifacts_toml)
            hash = artifact_hash(artifact_name, artifacts_toml)
            if hash !== nothing
                ensure_artifact_installed(artifact_name, artifacts_toml)
                dataset_root = artifact_path(hash)
            end
        end
    end
    isempty(dataset_root) && error("Release dataset not found.")
    return normpath(dataset_root)
end

const _RELEASE_ROOT = _release_dataset_root()
ENV["ARCHIMEDLIGHT_RELEASE_DATA_ROOT"] = _RELEASE_ROOT

include(joinpath(dirname(@__DIR__), "test", "release", "harness.jl"))

function _reference_metric_map(fx::JuliaFixture, data, metric::String; step_number::Int=0)
    ref_path = joinpath(fixture_reference_dir(fx), "component_values.csv")
    rows = collect(Tables.rowtable(CSV.File(ref_path; delim=';', normalizenames=false)))
    rows_by_key = Dict{Tuple{Int,Int,Int},Any}()
    for row in rows
        key = (Int(row.step_number), Int(row.item_id), Int(row.component_id))
        rows_by_key[key] = row
    end

    out = Dict{Int,Float64}()
    for (nid, (item_id, component_id)) in _simulation_output_keys(data.scene, data.models, data.options)
        row = get(rows_by_key, (step_number, Int(item_id), Int(component_id)), nothing)
        row === nothing && continue
        metric_sym = Symbol(metric)
        metric_sym in propertynames(row) || error(
            "Metric $(repr(metric)) is not available in $(ref_path) for fixture $(repr(fx.id)). " *
            "Available columns: $(join(string.(propertynames(row)), ", ")).",
        )
        value = getproperty(row, metric_sym)
        value === missing && continue
        out[nid] = Float64(value)
    end
    return out
end

function _plot_metric!(slot, geometry, values; title::String, colorrange, colormap=:viridis)
    vv = _vertex_values_for_step(geometry.vertices, geometry.faces, geometry.face2node, values)
    ax = Axis3(
        slot;
        title=title,
        aspect=:data,
        azimuth=1.45,
        elevation=0.35,
        perspectiveness=0.65,
    )
    plot = mesh!(
        ax,
        geometry.vertices,
        geometry.faces;
        color=vv,
        colormap=colormap,
        colorrange=colorrange,
        nan_color=:lightgray,
        shading=true,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    return plot
end

function _render_image_diff(ref_path::AbstractString, obs_path::AbstractString, out_path::AbstractString; title::String)
    ref = CairoMakie.load(ref_path)
    obs = CairoMakie.load(obs_path)
    size(ref) == size(obs) || error("Image sizes differ: $(size(ref)) vs $(size(obs))")

    diff = map(ref, obs) do a, b
        RGBAf(abs(Float32(a.r) - Float32(b.r)), abs(Float32(a.g) - Float32(b.g)), abs(Float32(a.b) - Float32(b.b)), 1.0)
    end

    fig = Figure(size=(1800, 760))
    ax1 = Axis(fig[1, 1], title="Reference")
    ax2 = Axis(fig[1, 2], title="Observed")
    ax3 = Axis(fig[1, 3], title="|Observed - Reference|")
    for ax in (ax1, ax2, ax3)
        hidedecorations!(ax)
        hidespines!(ax)
    end
    image!(ax1, rotr90(ref))
    image!(ax2, rotr90(obs))
    image!(ax3, rotr90(diff))
    Label(fig[0, 1:3], title, fontsize=24)

    mkpath(dirname(out_path))
    save(out_path, fig)
    println(out_path)
    return out_path
end

function render_fixture_diff(fixture_id::String; metric::Union{Nothing,String}=nothing, step_number::Int=0)
    fx = fixture_by_id(fixture_id)
    fx === nothing && error("Unknown fixture $(repr(fixture_id))")

    if metric === nothing
        data = fixture_runtime_data(fx)
        outdir = joinpath(dirname(@__DIR__), "test", "release", "visual_diffs")
        observed_path = joinpath(outdir, "$(fixture_id)_$(fx.visual_metric)_observed.png")
        write_fixture_reference_image!(fx; out_path=observed_path, data=data)
        diff_path = joinpath(outdir, "$(fixture_id)_$(fx.visual_metric)_image_diff.png")
        return _render_image_diff(
            fixture_reference_image_path(fx),
            observed_path,
            diff_path;
            title="$(fixture_id) | metric=$(fx.visual_metric)",
        )
    end

    data = fixture_runtime_data(fx)
    geometry = ArchimedLight._scene_geometry_for_interception(data.scene, data.models, data.options)
    observed = _budget_metric_map(data.series[step_number+1], metric)
    reference = _reference_metric_map(fx, data, metric; step_number=step_number)
    diff = Dict{Int,Float64}()
    for nid in union(keys(observed), keys(reference))
        diff[nid] = abs(get(observed, nid, 0.0) - get(reference, nid, 0.0))
    end

    finite_vals = Float64[]
    append!(finite_vals, values(observed))
    append!(finite_vals, values(reference))
    main_range =
        if isempty(finite_vals)
            (0.0, 1.0)
        else
            lo = minimum(finite_vals)
            hi = maximum(finite_vals)
            lo == hi ? (lo - 1e-12, hi + 1e-12) : (lo, hi)
        end
    diff_hi = isempty(diff) ? 1.0 : max(maximum(values(diff)), 1e-12)

    fig = Figure(size=(1700, 760))
    p1 = _plot_metric!(fig[1, 1], geometry, reference; title="Reference $(metric)", colorrange=main_range)
    _plot_metric!(fig[1, 2], geometry, observed; title="Observed $(metric)", colorrange=main_range)
    p3 = _plot_metric!(fig[1, 3], geometry, diff; title="|Observed - Reference|", colorrange=(0.0, diff_hi), colormap=:magma)
    Label(fig[0, 1:3], "$(fixture_id) | step=$(step_number)", fontsize=24)
    Colorbar(fig[1, 4], p1, label=metric)
    Colorbar(fig[1, 5], p3, label="abs diff")

    outdir = joinpath(dirname(@__DIR__), "test", "release", "visual_diffs")
    mkpath(outdir)
    outpath = joinpath(outdir, "$(fixture_id)_$(metric)_step$(step_number)_diff.png")
    save(outpath, fig)
    println(outpath)
    return outpath
end

function main(args)
    metric = nothing
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
        println("Usage:")
        println("  julia --project=test scripts/render_release_fixture_diffs.jl [fixture-id ...]")
        println("  julia --project=test scripts/render_release_fixture_diffs.jl --all")
        println("  julia --project=test scripts/render_release_fixture_diffs.jl --metric Ri_PAR_f [fixture-id ...]")
        println("")
        println("With no fixture ids, renders the default quick set:")
        println("  test-scattering-one-plate")
        println("  test-scattering-two-plates")
        println("")
        println("By default, each fixture uses its `visual_metric` from the release manifest")
        println("and compares the rendered reference image against the observed rendered image.")
        println("Use --metric to force the old component-metric topology diff mode.")
        return nothing
    elseif arg == "--metric"
            i == length(args) && error("--metric requires a value")
            metric = args[i + 1]
            i += 2
            continue
        else
            push!(positional, arg)
        end
        i += 1
    end

    fixtures =
        if any(arg -> arg == "--all", positional)
            [fx.id for fx in julia_fixtures()]
        elseif isempty(positional)
            ["test-scattering-one-plate", "test-scattering-two-plates"]
        else
            filter(!=("--all"), positional)
        end

    for fx in fixtures
        render_fixture_diff(fx; metric=metric, step_number=0)
    end
    return nothing
end

main(ARGS)
