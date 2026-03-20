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
        value = getproperty(row, Symbol(metric))
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

function render_fixture_diff(fixture_id::String; metric::String="Ri_NIR_q", step_number::Int=0)
    fx = fixture_by_id(fixture_id)
    fx === nothing && error("Unknown fixture $(repr(fixture_id))")
    data = fixture_runtime_data(fx)
    geometry = ArchimedLight._scene_geometry_for_interception(data.scene, data.models, data.options)
    observed = _budget_metric_map(data.series[step_number + 1], metric)
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
    fixtures = isempty(args) ? ["test-scattering-one-plate", "test-scattering-two-plates"] : collect(args)
    for fx in fixtures
        render_fixture_diff(fx; metric="Ri_NIR_q", step_number=0)
    end
end

main(ARGS)
