using CairoMakie

const _REF_BG = RGBf(0.985, 0.984, 0.975)
const _REF_INK = RGBf(0.14, 0.16, 0.18)
const _REF_MUTED = RGBf(0.45, 0.48, 0.50)
const _REF_BLUE = RGBf(0.16, 0.45, 0.72)
const _REF_TEAL = RGBf(0.12, 0.62, 0.58)
const _REF_GOLD = RGBf(0.88, 0.68, 0.18)
const _REF_RED = RGBf(0.76, 0.28, 0.24)
const _REF_GRID = RGBf(0.83, 0.84, 0.82)
const _REF_LEAF = RGBf(0.37, 0.63, 0.28)
const _REF_SENSOR = RGBf(0.77, 0.48, 0.18)

function _hide_axis!(ax)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

function _box!(ax, x, y, w, h; color, strokecolor=_REF_INK, linewidth=2, radius=0.05)
    poly!(
        ax,
        Makie.Polygon(
            Point2f[
                (x, y),
                (x + w, y),
                (x + w, y + h),
                (x, y + h),
            ],
        ),
        color=color,
        strokecolor=strokecolor,
        strokewidth=linewidth,
    )
    if radius > 0
        lines!(ax, [x, x + w, x + w, x, x], [y, y, y + h, y + h, y], color=strokecolor, linewidth=linewidth)
    end
    return nothing
end

function _arrow!(ax, p1, p2; color=_REF_INK, linewidth=3, arrowsize=18)
    arrows2d!(ax, [Point2f(p1)], [Vec2f(p2[1] - p1[1], p2[2] - p1[2])], color=color, lengthscale=1,
        shaftwidth=linewidth * 0.012, tipwidth=arrowsize * 0.012, tiplength=arrowsize * 0.02)
    return nothing
end

function _label!(ax, x, y, s; kwargs...)
    text!(ax, x, y, text=s; color=_REF_INK, align=(:center, :center), kwargs...)
end

function pipeline_overview_figure(path::AbstractString)
    fig = Figure(size=(1200, 360), backgroundcolor=_REF_BG)
    ax = Axis(fig[1, 1], limits=((0, 14), (0, 4)), backgroundcolor=_REF_BG)
    _hide_axis!(ax)

    boxes = [
        (0.4, 1.5, 2.1, 1.0, _REF_BLUE, "Scene +\nmodels + meteo"),
        (3.0, 1.5, 2.1, 1.0, _REF_TEAL, "Sky state\nand turtle fluxes"),
        (5.6, 1.5, 2.2, 1.0, _REF_GOLD, "Projected pixel\nstacks by direction"),
        (8.4, 1.5, 2.1, 1.0, _REF_LEAF, "First-order\ninterception"),
        (11.0, 1.5, 2.0, 1.0, _REF_RED, "Iterative\nscattering"),
    ]
    for (x, y, w, h, c, lbl) in boxes
        _box!(ax, x, y, w, h, color=c)
        _label!(ax, x + w / 2, y + h / 2, lbl, fontsize=22, font=:bold)
    end

    for x in (2.5, 5.1, 7.9, 10.5)
        _arrow!(ax, (x, 2.0), (x + 0.4, 2.0), color=_REF_INK)
    end

    _box!(ax, 8.4, 0.25, 4.6, 0.7, color=RGBAf(_REF_MUTED.r, _REF_MUTED.g, _REF_MUTED.b, 0.10), linewidth=1.5)
    _label!(ax, 10.7, 0.6, "Export Ri / Ra per component, scene, and summary tables", fontsize=20)
    _arrow!(ax, (10.6, 1.5), (10.6, 0.98), color=_REF_MUTED, linewidth=2.5, arrowsize=15)
    _arrow!(ax, (12.0, 1.5), (12.0, 0.98), color=_REF_MUTED, linewidth=2.5, arrowsize=15)

    _label!(ax, 6.9, 3.35, "ARCHIMED light pipeline", fontsize=28, font=:bold)
    _label!(ax, 6.9, 2.95, "Java `MirProcess` and Julia `run_light_step` follow the same sequence", fontsize=18, color=_REF_MUTED)

    save(path, fig)
    return path
end

function projection_area_ratio_figure(path::AbstractString)
    fig = Figure(size=(1200, 520), backgroundcolor=_REF_BG)
    ax1 = Axis(fig[1, 1], limits=((0, 10), (0, 8)), backgroundcolor=_REF_BG)
    ax2 = Axis(fig[1, 2], limits=((0, 10), (0, 8)), backgroundcolor=_REF_BG)
    _hide_axis!(ax1)
    _hide_axis!(ax2)

    _label!(ax1, 5, 7.5, "Projection onto the plot", fontsize=24, font=:bold)
    tri = Point2f[(2.0, 5.5), (7.0, 6.3), (5.4, 3.2)]
    poly!(ax1, Makie.Polygon(tri), color=RGBAf(_REF_LEAF.r, _REF_LEAF.g, _REF_LEAF.b, 0.75), strokecolor=_REF_INK, strokewidth=2)
    _label!(ax1, 5.1, 5.4, "triangle in scene", fontsize=18)
    _arrow!(ax1, (8.6, 6.6), (6.9, 5.3), color=_REF_RED)
    _label!(ax1, 8.7, 6.9, "direction", fontsize=18, color=_REF_RED)

    lines!(ax1, [1.0, 8.8], [1.4, 1.4], color=_REF_MUTED, linewidth=2)
    _label!(ax1, 8.2, 1.0, "plot plane", fontsize=16, color=_REF_MUTED)
    proj = Point2f[(2.2, 2.6), (6.6, 3.4), (5.0, 1.9)]
    poly!(ax1, Makie.Polygon(proj), color=RGBAf(_REF_BLUE.r, _REF_BLUE.g, _REF_BLUE.b, 0.35), strokecolor=_REF_BLUE, strokewidth=2)

    for p in tri
        _arrow!(ax1, (p[1], p[2]), (p[1] + 0.3, p[2] - 2.5), color=RGBAf(_REF_MUTED.r, _REF_MUTED.g, _REF_MUTED.b, 0.7), linewidth=2, arrowsize=13)
    end

    _label!(ax2, 5, 7.5, "Rasterization and area-ratio correction", fontsize=24, font=:bold)
    x0, y0, dx = 1.2, 1.0, 1.2
    for i in 0:5
        lines!(ax2, [x0 + i * dx, x0 + i * dx], [y0, y0 + 6dx/5], color=_REF_GRID, linewidth=1.5)
        lines!(ax2, [x0, x0 + 6dx], [y0 + i * dx / 1.0, y0 + i * dx / 1.0], color=_REF_GRID, linewidth=1.5)
    end
    raster = Point2f[(2.0, 5.2), (6.3, 5.8), (5.4, 2.0), (2.7, 1.9)]
    poly!(ax2, Makie.Polygon(raster), color=RGBAf(_REF_BLUE.r, _REF_BLUE.g, _REF_BLUE.b, 0.14), strokecolor=_REF_BLUE, strokewidth=2)

    hit_cells = [(2, 3), (3, 3), (4, 3), (2, 2), (3, 2), (4, 2), (4, 1)]
    for (i, j) in hit_cells
        rect = Point2f[
            (x0 + i * dx, y0 + j * dx),
            (x0 + (i + 1) * dx, y0 + j * dx),
            (x0 + (i + 1) * dx, y0 + (j + 1) * dx),
            (x0 + i * dx, y0 + (j + 1) * dx),
        ]
        poly!(ax2, Makie.Polygon(rect), color=RGBAf(_REF_GOLD.r, _REF_GOLD.g, _REF_GOLD.b, 0.55), strokecolor=_REF_INK, strokewidth=0.8)
    end

    _label!(ax2, 7.7, 5.9, "covered pixels", fontsize=16, color=_REF_GOLD)
    _label!(ax2, 7.7, 4.9, "projected triangle", fontsize=16, color=_REF_BLUE)
    _label!(ax2, 7.2, 3.45, "area ratio =", fontsize=20, align=(:left, :center))
    text!(ax2, 7.2, 2.6, text="projected mesh area\nprojected hit-pixel area", align=(:left, :center), fontsize=18, color=_REF_INK)
    lines!(ax2, [7.15, 9.55], [2.95, 2.95], color=_REF_INK, linewidth=2)
    _label!(ax2, 7.9, 1.45, "This rescales pixel counts back\ntoward geometric projected area.", fontsize=16, color=_REF_MUTED, align=(:left, :center))

    save(path, fig)
    return path
end

function scattering_transfer_figure(path::AbstractString)
    fig = Figure(size=(1200, 500), backgroundcolor=_REF_BG)
    ax = Axis(fig[1, 1], limits=((0, 14), (0, 9)), backgroundcolor=_REF_BG)
    _hide_axis!(ax)

    _label!(ax, 7.0, 8.35, "Scattering transfer on a pixel stack", fontsize=26, font=:bold)
    _label!(ax, 7.0, 7.85, "The same ordered hits used for interception are reused for iterative energy exchange.", fontsize=18, color=_REF_MUTED)

    _box!(ax, 5.1, 1.0, 3.0, 5.8, color=RGBAf(_REF_MUTED.r, _REF_MUTED.g, _REF_MUTED.b, 0.08), linewidth=1.5)
    _label!(ax, 6.6, 0.55, "one projected pixel", fontsize=17, color=_REF_MUTED)

    components = [
        (5.6, 5.5, 2.0, 0.8, _REF_LEAF, "leaf A"),
        (5.6, 4.2, 2.0, 0.55, _REF_SENSOR, "virtual\nsensor"),
        (5.6, 2.7, 2.0, 0.8, _REF_BLUE, "leaf B"),
        (5.6, 1.4, 2.0, 0.7, _REF_GOLD, "ground"),
    ]
    for (x, y, w, h, c, lbl) in components
        _box!(ax, x, y, w, h, color=c)
        _label!(ax, x + w / 2, y + h / 2, lbl, fontsize=18, font=:bold)
    end

    _arrow!(ax, (4.0, 6.9), (5.35, 6.0), color=_REF_RED)
    _label!(ax, 3.0, 7.2, "incoming\nfirst-order light", fontsize=18, color=_REF_RED)

    _arrow!(ax, (8.1, 5.95), (9.8, 4.2), color=_REF_TEAL)
    _arrow!(ax, (8.1, 3.1), (9.8, 4.8), color=_REF_TEAL)
    _label!(ax, 11.2, 4.5, "adjacent visible\nhits exchange scattered\nenergy upward and\ndownward", fontsize=18, color=_REF_TEAL)

    _arrow!(ax, (6.6, 5.5), (6.6, 4.8), color=_REF_MUTED, linewidth=2.2, arrowsize=14)
    _arrow!(ax, (6.6, 4.75), (6.6, 3.55), color=_REF_MUTED, linewidth=2.2, arrowsize=14)
    _label!(ax, 3.2, 3.9, "virtual sensors stay in the stack\nbut are transparent for transfer logic", fontsize=17, color=_REF_MUTED)

    text!(ax, 9.2, 1.55, text="node scattering pool\n× waveband coefficient\n÷ directional hit count\n÷ 2", align=(:left, :center), fontsize=18, color=_REF_INK)
    lines!(ax, [9.0, 12.7], [2.45, 2.45], color=_REF_GRID, linewidth=2)

    save(path, fig)
    return path
end

function generate_reference_figures(outdir::AbstractString=joinpath(@__DIR__, "src", "assets"))
    mkpath(outdir)
    CairoMakie.activate!(type="svg")
    set_theme!(Theme(font="TeX Gyre Heros", backgroundcolor=_REF_BG))

    pipeline_overview_figure(joinpath(outdir, "archimed_pipeline_overview.svg"))
    projection_area_ratio_figure(joinpath(outdir, "archimed_projection_area_ratio.svg"))
    scattering_transfer_figure(joinpath(outdir, "archimed_scattering_transfer.svg"))
    return outdir
end
