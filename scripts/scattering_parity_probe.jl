using ArchimedLight
using CSV
using Tables
using Printf

function _repo_root()
    normpath(joinpath(@__DIR__, ".."))
end

function _fixture_root()
    joinpath(_repo_root(), "java_implementation", "archimed-lib-2018", "tests")
end

function _fixture_cfg_path(name::String)
    base = joinpath(_fixture_root(), name)
    p = joinpath(base, "archimed2018", "config.yml")
    isfile(p) && return p
    p = joinpath(base, "config.yml")
    isfile(p) && return p
    error("fixture config not found for $(name)")
end

function _expected_component_path(name::String)
    d = joinpath(_fixture_root(), name, "expected")
    for f in ("component_values.csv", "component_values_allinturtle.csv", "component_values_notallinturtle.csv")
        p = joinpath(d, f)
        isfile(p) && return p
    end
    return nothing
end

function _csv_columns(rows)
    isempty(rows) && return String[]
    String[string(n) for n in propertynames(first(rows))]
end

function _row_get_by_col(row, col::String)
    if row isa AbstractDict
        return get(row, col, missing)
    end
    s = Symbol(col)
    return s in propertynames(row) ? getproperty(row, s) : missing
end

function _row_key(row, key_cols::Vector{String})
    Tuple(_row_get_by_col(row, c) for c in key_cols)
end

function _row_float(row, col::String)
    v = _row_get_by_col(row, col)
    v === missing && return nothing
    v isa Number && return Float64(v)
    p = tryparse(Float64, strip(string(v)))
    p === nothing ? nothing : p
end

function _compare_rows(expected_rows, observed_rows; key_cols::Vector{String}, value_cols::Vector{String}, atol::Float64=1e-3, rtol::Float64=3e-3, top_n::Int=20)
    exp_map = Dict{Tuple,Any}(_row_key(r, key_cols) => r for r in expected_rows)
    obs_map = Dict{Tuple,Any}(_row_key(r, key_cols) => r for r in observed_rows)
    keys_all = Set(keys(exp_map))
    union!(keys_all, keys(obs_map))
    mism = NamedTuple[]
    for k in keys_all
        haskey(exp_map, k) || continue
        haskey(obs_map, k) || continue
        e = exp_map[k]
        o = obs_map[k]
        for c in value_cols
            ev = _row_float(e, c)
            ov = _row_float(o, c)
            (ev === nothing || ov === nothing) && continue
            if isnan(ev) && isnan(ov)
                continue
            end
            abs_err = abs(ov - ev)
            rel_err = abs_err / max(abs(ev), eps(Float64))
            if !(abs_err <= atol || rel_err <= rtol)
                push!(mism, (key=k, col=c, expected=ev, observed=ov, abs_err=abs_err, rel_err=rel_err))
            end
        end
    end
    sort!(mism, by=x -> x.rel_err, rev=true)
    length(mism) > top_n && resize!(mism, top_n)
    return mism
end

function _stack_component_rows(scene, series, cfg, meteo_rows, cols::Vector{String})
    out = Dict{String,Any}[]
    for i in eachindex(series)
        rows = ArchimedLight.component_values_table(scene, series[i], cfg; meteo_row=meteo_rows[i], step_number=i - 1, columns=cols).rows
        append!(out, rows)
    end
    out
end

function _java_pairs_from_debug(debug_dir::String)
    files = [joinpath(debug_dir, @sprintf("log-pixeltable-dir%02d-step00.csv", d)) for d in 0:15]
    pairs = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Int}()
    for f in files
        isfile(f) || continue
        rows = Tables.rowtable(CSV.File(f; delim=';'))
        bypix = Dict{Tuple{Int,Int},Vector{Tuple{Int,Tuple{Int,Int}}}}()
        for r in rows
            key = (Int(r.plantId), Int(r.nodeId))
            push!(get!(bypix, (Int(r.x), Int(r.y)), Vector{Tuple{Int,Tuple{Int,Int}}}()), (Int(r.n), key))
        end
        for st in values(bypix)
            sort!(st, by=t -> t[1])
            for i in 1:(length(st) - 1)
                a = st[i][2]
                b = st[i + 1][2]
                pairs[(a, b)] = get(pairs, (a, b), 0) + 1
                pairs[(b, a)] = get(pairs, (b, a), 0) + 1
            end
        end
    end
    pairs
end

function _julia_pairs_for_fixture(scene, cfg, row)
    sky = ArchimedLight.compute_sky(row, cfg)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    pair_counts, _, _, _ = ArchimedLight._pair_counts_for_scattering(scene, turtle, cfg)
    keymap = ArchimedLight._interception_java_keys(scene, cfg)
    pairs = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Int}()
    for ((to, from), c) in pair_counts
        kt = get(keymap, to, (1, to + 1))
        kf = get(keymap, from, (1, from + 1))
        pairs[(kt, kf)] = get(pairs, (kt, kf), 0) + c
    end
    return pairs
end

function _pair_rel_sad(jpairs, zpairs)
    keys_all = Set(keys(jpairs))
    union!(keys_all, keys(zpairs))
    sad = 0
    for k in keys_all
        sad += abs(get(jpairs, k, 0) - get(zpairs, k, 0))
    end
    return sad / max(sum(values(zpairs)), 1)
end

function _java_iter0_total_w(debug_dir::String, band::String)
    p = joinpath(debug_dir, "log-iteration-scat-$(lowercase(band)).csv")
    isfile(p) || return nothing
    rows = Tables.rowtable(CSV.File(p; delim=';'))
    s = 0.0
    for r in rows
        Int(r.iter) == 0 || continue
        s += Float64(r.scat_W)
    end
    return s
end

function _julia_iter0_total_w(scene, cfg, row, band::String)
    sky = ArchimedLight.compute_sky(row, cfg)
    turtle = ArchimedLight.build_turtle(cfg, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    first_order = ArchimedLight.compute_first_order(scene, turtle, fluxes, cfg)
    graph = ArchimedLight.build_scattering_transfer_graph(scene, turtle, first_order, cfg)
    coeff_default = uppercase(band) == "NIR" ? cfg.scattering_coeff_nir : cfg.scattering_coeff_par
    coeff_by_node = ArchimedLight._coeff_by_node(graph, band, coeff_default)
    current =
        if uppercase(band) == "NIR"
            Dict{Int,Float64}(nid => get(first_order.incident_nir_power_per_node, nid, 0.0) for nid in graph.node_ids)
        else
            Dict{Int,Float64}(nid => get(first_order.incident_par_power_per_node, nid, 0.0) for nid in graph.node_ids)
        end

    hit_energy = Dict{Int,Float64}()
    for nid in graph.node_ids
        nh = get(graph.all_hits, nid, 0)
        nh <= 0 && continue
        hit_energy[nid] = get(current, nid, 0.0) * get(coeff_by_node, nid, coeff_default) / nh / 2.0
    end

    total_next = 0.0
    for ((_, from), cnt) in graph.pair_counts
        total_next += cnt * get(hit_energy, from, 0.0)
    end
    return total_next
end

function _print_top_mism(mism)
    for m in mism
        println("  key=", m.key, " col=", m.col, " expected=", m.expected, " observed=", m.observed, " rel_err=", m.rel_err)
    end
end

function main()
    if isempty(ARGS)
        println("Usage: julia --project=. scripts/scattering_parity_probe.jl <fixture-name> [java-debug-output-dir]")
        return
    end
    fixture = ARGS[1]
    cfg_path = _fixture_cfg_path(fixture)
    cfg = ArchimedLight.read_light_config(cfg_path)
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    selected = ArchimedLight.prepare_meteo(meteo, cfg)
    series = ArchimedLight.run_light_series(scene, meteo, cfg)

    exp_comp_path = _expected_component_path(fixture)
    if exp_comp_path !== nothing
        expected = Tables.rowtable(CSV.File(exp_comp_path; delim=';'))
        cols = _csv_columns(expected)
        observed = _stack_component_rows(scene, series, cfg, selected.rows, cols)
        mism = _compare_rows(
            expected,
            observed;
            key_cols=["step_number", "item_id", "component_id"],
            value_cols=[c for c in cols if !(c in ("step_number", "item_id", "component_id"))],
            atol=1e-3,
            rtol=3e-3,
            top_n=15,
        )
        println("component mismatches(top 15) = ", length(mism))
        _print_top_mism(mism)
    else
        println("no frozen expected component_values.csv for fixture ", fixture)
    end

    debug_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(cfg_path), "output", "000002")
    if !isdir(debug_dir)
        println("debug output directory not found: ", debug_dir)
        return
    end

    row = first(selected.rows)
    jpar = _julia_iter0_total_w(scene, cfg, row, "PAR")
    jnir = _julia_iter0_total_w(scene, cfg, row, "NIR")
    zpar = _java_iter0_total_w(debug_dir, "PAR")
    znir = _java_iter0_total_w(debug_dir, "NIR")
    println("iter0 scattered total W: PAR julia=", jpar, " java=", zpar)
    println("iter0 scattered total W: NIR julia=", jnir, " java=", znir)

    jpairs = _julia_pairs_for_fixture(scene, cfg, row)
    zpairs = _java_pairs_from_debug(debug_dir)
    println("pair rel_sad = ", _pair_rel_sad(jpairs, zpairs))
end

main()
