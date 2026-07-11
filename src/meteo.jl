const _METEO_SOLAR_TO_PAR = 0.48
const _METEO_SOLAR_TO_NIR = 0.52
const _METEO_ALIAS_RTOL = 1.0e-10
const _METEO_ALIAS_ATOL = 1.0e-10
const _METEO_ENERGY_RTOL = 0.01
const _METEO_ENERGY_ATOL = 1.0
const _METEO_HIGH_IRRADIANCE_WARNING = 2000.0

const _METEO_ALIASES = (
    ri_sw_f=(:Ri_SW_f, :RI_SW_f, :Rg, :rg, :sw_global, :global),
    ri_par_f=(:Ri_PAR_f, :RI_PAR_f, :PAR, :par),
    ri_nir_f=(:Ri_NIR_f, :RI_NIR_f, :NIR, :nir),
    clearness=(:clearness, :Kt),
    latitude=(:latitude, :lat),
    sun_azimuth=(:sun_azimuth, :sun_azimut),
    sun_elevation=(:sun_elevation,),
    direct_fraction=(:direct_fraction, :fDIR_SW, :Fd),
    ri_sw_direct=(:Ri_SW_f_direct, :RI_SW_f_direct),
    ri_sw_diffuse=(:Ri_SW_f_diffuse, :RI_SW_f_diffuse),
)

const _METEO_CANONICAL_NAMES = (
    ri_sw_f=:Ri_SW_f,
    ri_par_f=:Ri_PAR_f,
    ri_nir_f=:Ri_NIR_f,
    clearness=:clearness,
    latitude=:latitude,
    sun_azimuth=:sun_azimuth,
    sun_elevation=:sun_elevation,
    direct_fraction=:direct_fraction,
    ri_sw_direct=:Ri_SW_f_direct,
    ri_sw_diffuse=:Ri_SW_f_diffuse,
)

"""
    ResolvedMeteoStep

Internal, value-resolved meteorological forcing for one simulation step.

All main-waveband irradiances are finite. Canonical date/start/end values are
optional so explicit solar geometry can be used without astronomical inputs.
Other optional values are `nothing` when they must still be derived by the sky
model. `sources` records the selected raw column/metadata name or the derivation
path used for each logical variable.
"""
struct ResolvedMeteoStep
    row_index::Int
    duration_seconds::Float64
    date::Union{Nothing,Dates.Date}
    start_hour::Union{Nothing,Float64}
    end_hour::Union{Nothing,Float64}
    ri_sw_f::Float64
    ri_par_f::Float64
    ri_nir_f::Float64
    clearness::Union{Nothing,Float64}
    latitude_deg::Union{Nothing,Float64}
    sun_azimuth_deg::Union{Nothing,Float64}
    sun_elevation_deg::Union{Nothing,Float64}
    direct_fraction::Union{Nothing,Float64}
    diffuse_fraction::Union{Nothing,Float64}
    ri_sw_direct::Union{Nothing,Float64}
    ri_sw_diffuse::Union{Nothing,Float64}
    extra_irradiance::Dict{String,Float64}
    radiation_source::Symbol
    solar_geometry_source::Symbol
    partition_source::Symbol
    sources::Dict{Symbol,Symbol}
end

struct MeteoResolutionResult
    step::Union{Nothing,ResolvedMeteoStep}
    errors::Vector{String}
    warnings::Vector{String}
    infos::Vector{String}
end

struct _MeteoLogicalValue
    value::Union{Nothing,Float64}
    source::Union{Nothing,Symbol}
    unavailable::Vector{Pair{Symbol,String}}
end

_meteo_row_label(row_index::Int) = "Meteo row $(row_index)"

function _meteo_row_names(row)
    row isa PlantMeteo.TimeStepRow && return Tuple(Symbol.(Tables.columnnames(row)))
    return Tuple(Symbol.(propertynames(row)))
end

function _meteo_metadata(row)
    row isa PlantMeteo.TimeStepRow || return NamedTuple()
    source = getfield(row, :source)
    hasfield(typeof(source), :metadata) || return NamedTuple()
    return getfield(source, :metadata)
end

function _meteo_parse_number(value)
    value === nothing && return nothing, "nothing"
    value === missing && return nothing, "missing"
    if value isa Number
        parsed = Float64(value)
        return isfinite(parsed) ? (parsed, "") : (nothing, string(parsed))
    end

    text = strip(string(value))
    isempty(text) && return nothing, "blank"
    lowercase(text) in ("missing", "na", "n/a", "null", "nothing") &&
        return nothing, text
    parsed = tryparse(Float64, text)
    parsed === nothing && return nothing, "not numeric ($(repr(value)))"
    return isfinite(parsed) ? (parsed, "") : (nothing, string(parsed))
end

function _meteo_collect_level_values(container, names, aliases)
    finite = Pair{Symbol,Float64}[]
    unavailable = Pair{Symbol,String}[]
    for name in aliases
        name in names || continue
        value, reason = _meteo_parse_number(getproperty(container, name))
        if value === nothing
            push!(unavailable, name => reason)
        else
            push!(finite, name => value)
        end
    end
    return finite, unavailable
end

function _meteo_choose_finite_value!(
    errors::Vector{String},
    row_index::Int,
    logical_name::Symbol,
    canonical_name::Symbol,
    finite::Vector{Pair{Symbol,Float64}},
)
    isempty(finite) && return nothing, nothing
    reference = first(finite).second
    is_irradiance = logical_name in (
        :ri_sw_f,
        :ri_par_f,
        :ri_nir_f,
        :ri_sw_direct,
        :ri_sw_diffuse,
    )
    if any(
        item -> !(
            is_irradiance ?
            _meteo_energy_agrees(item.second, reference) :
            isapprox(
                item.second,
                reference;
                rtol=_METEO_ALIAS_RTOL,
                atol=_METEO_ALIAS_ATOL,
            )
        ),
        @view(finite[2:end]),
    )
        details = join(("$(item.first)=$(repr(item.second))" for item in finite), ", ")
        push!(
            errors,
            "$(_meteo_row_label(row_index)) has conflicting aliases for `$(canonical_name)`: $(details). Keep one spelling or make the values agree.",
        )
        return nothing, nothing
    end

    canonical_index = findfirst(item -> item.first == canonical_name, finite)
    chosen = canonical_index === nothing ? first(finite) : finite[canonical_index]
    return chosen.second, chosen.first
end

function _meteo_resolve_logical_value!(
    errors::Vector{String},
    row,
    row_index::Int,
    logical_name::Symbol,
    row_names,
    metadata,
    metadata_names,
)
    aliases = getproperty(_METEO_ALIASES, logical_name)
    canonical_name = getproperty(_METEO_CANONICAL_NAMES, logical_name)
    row_finite, row_unavailable = _meteo_collect_level_values(row, row_names, aliases)
    if !isempty(row_finite)
        value, source = _meteo_choose_finite_value!(
            errors,
            row_index,
            logical_name,
            canonical_name,
            row_finite,
        )
        return _MeteoLogicalValue(value, source, row_unavailable)
    end

    metadata_finite, metadata_unavailable =
        _meteo_collect_level_values(metadata, metadata_names, aliases)
    value, source = _meteo_choose_finite_value!(
        errors,
        row_index,
        logical_name,
        canonical_name,
        metadata_finite,
    )
    return _MeteoLogicalValue(
        value,
        source,
        vcat(row_unavailable, metadata_unavailable),
    )
end

function _meteo_parse_time_value(value)
    value === nothing && return nothing, "nothing"
    value === missing && return nothing, "missing"
    value isa Dates.Time && return value, ""
    value isa Dates.DateTime && return Dates.Time(value), ""

    text = strip(string(value))
    isempty(text) && return nothing, "blank"
    lowercase(text) in ("missing", "na", "n/a", "null", "nothing") &&
        return nothing, text
    for format in (Dates.DateFormat("HH:MM:SS"), Dates.DateFormat("HH:MM"))
        parsed = try
            Dates.Time(text, format)
        catch
            nothing
        end
        parsed === nothing || return parsed, ""
    end
    return nothing, "not a valid time ($(repr(value)))"
end

_meteo_decimal_hour(time::Dates.Time) =
    Dates.hour(time) + Dates.minute(time) / 60 + Dates.second(time) / 3600 +
    Dates.millisecond(time) / 3.6e6

function _meteo_parse_date_value(value)
    value === nothing && return nothing, nothing, "nothing"
    value === missing && return nothing, nothing, "missing"
    value isa Dates.DateTime &&
        return Dates.Date(value), _meteo_decimal_hour(Dates.Time(value)), ""
    value isa Dates.Date && return value, nothing, ""

    text = strip(string(value))
    isempty(text) && return nothing, nothing, "blank"
    lowercase(text) in ("missing", "na", "n/a", "null", "nothing") &&
        return nothing, nothing, text

    datetime = try
        Dates.DateTime(text)
    catch
        nothing
    end
    datetime === nothing || return (
        Dates.Date(datetime),
        _meteo_decimal_hour(Dates.Time(datetime)),
        "",
    )

    for format in (Dates.DateFormat("yyyy/mm/dd"), Dates.DateFormat("yyyy-mm-dd"))
        parsed = try
            Dates.Date(text, format)
        catch
            nothing
        end
        parsed === nothing || return parsed, nothing, ""
    end
    return nothing, nothing, "not a valid date ($(repr(value)))"
end

function _meteo_parse_integer_value(value, name::Symbol)
    parsed, reason = _meteo_parse_number(value)
    parsed === nothing && return nothing, reason
    rounded = round(Int, parsed)
    isapprox(parsed, rounded; rtol=0.0, atol=0.0) ||
        return nothing, "$(name) must be an integer ($(repr(value)))"
    return rounded, ""
end

function _meteo_resolve_date(row, names)
    reasons = String[]
    if :date in names
        date, datetime_hour, reason = _meteo_parse_date_value(getproperty(row, :date))
        date === nothing || return date, datetime_hour, :date, reasons
        push!(reasons, "`date` is $(reason)")
    end

    if :dayofyear in names
        doy, doy_reason =
            _meteo_parse_integer_value(getproperty(row, :dayofyear), :dayofyear)
        year = 2000
        year_reason = ""
        if :year in names
            parsed_year, year_reason =
                _meteo_parse_integer_value(getproperty(row, :year), :year)
            parsed_year === nothing || (year = parsed_year)
        end
        if doy !== nothing && isempty(year_reason)
            try
                last_doy = Dates.dayofyear(Dates.Date(year, 12, 31))
                if 1 <= doy <= last_doy
                    return Dates.Date(year, 1, 1) + Dates.Day(doy - 1), nothing, :dayofyear, reasons
                end
                doy_reason = "dayofyear=$(repr(doy)) is outside 1:$(last_doy)"
            catch err
                year_reason = "invalid year $(repr(year)): $(sprint(showerror, err))"
            end
        end
        isempty(doy_reason) || push!(reasons, "`dayofyear` is $(doy_reason)")
        isempty(year_reason) || push!(reasons, "`year` is $(year_reason)")
    end

    isempty(reasons) && push!(reasons, "neither `date` nor `dayofyear` is provided")
    return nothing, nothing, nothing, reasons
end

function _meteo_resolve_start_time(row, datetime_hour, names)
    reasons = String[]
    for name in (:hour_start, :hour)
        name in names || continue
        value, reason = _meteo_parse_time_value(getproperty(row, name))
        value === nothing || return _meteo_decimal_hour(value), name, reasons
        push!(reasons, "`$(name)` is $(reason)")
    end
    datetime_hour === nothing || return datetime_hour, :date, reasons
    isempty(reasons) &&
        push!(reasons, "no `hour_start`, `hour`, or DateTime-valued `date` is provided")
    return nothing, nothing, reasons
end

function _meteo_resolve_end_time(row, names)
    :hour_end in names ||
        return nothing, nothing, ["no `hour_end` is provided"]
    value, reason = _meteo_parse_time_value(getproperty(row, :hour_end))
    value === nothing && return nothing, nothing, ["`hour_end` is $(reason)"]
    return _meteo_decimal_hour(value), :hour_end, String[]
end

function _meteo_resolve_explicit_duration(row, names)
    reasons = String[]
    for name in (:step_duration, :duration)
        name in names || continue
        duration = try
            Float64(_duration_seconds_strict(getproperty(row, name); field_name=string(name)))
        catch err
            push!(reasons, "`$(name)` is invalid: $(sprint(showerror, err))")
            nothing
        end
        duration === nothing || return duration, name, reasons
    end
    isempty(reasons) && push!(reasons, "no `step_duration` or `duration` is provided")
    return nothing, nothing, reasons
end

function _meteo_resolve_temporal!(
    errors::Vector{String},
    row,
    row_index::Int,
    sources,
    row_names,
)
    date, datetime_hour, date_source, date_reasons = _meteo_resolve_date(row, row_names)
    start_hour, start_source, start_reasons =
        _meteo_resolve_start_time(row, datetime_hour, row_names)
    explicit_end_hour, explicit_end_source, end_reasons =
        _meteo_resolve_end_time(row, row_names)
    duration_seconds, duration_source, duration_reasons =
        _meteo_resolve_explicit_duration(row, row_names)

    end_before_start =
        start_hour !== nothing && explicit_end_hour !== nothing &&
        explicit_end_hour < start_hour
    end_before_start && push!(
        errors,
        "end is before start at meteo row $(row_index); cross-midnight rollover is not enabled.",
    )

    if duration_seconds !== nothing && start_hour !== nothing &&
       explicit_end_hour !== nothing && !end_before_start
        hours_duration = (explicit_end_hour - start_hour) * 3600.0
        tolerance = max(1.0e-6, 1.0e-9 * max(abs(hours_duration), abs(duration_seconds), 1.0))
        isapprox(duration_seconds, hours_duration; atol=tolerance, rtol=0.0) || push!(
            errors,
            "$(_meteo_row_label(row_index)) has conflicting temporal inputs: `hour_start` + `hour_end` imply $(repr(hours_duration)) seconds, but `$(duration_source)` provides $(repr(duration_seconds)) seconds.",
        )
    end

    if duration_seconds === nothing && start_hour !== nothing &&
       explicit_end_hour !== nothing && !end_before_start
        derived_duration = (explicit_end_hour - start_hour) * 3600.0
        if isfinite(derived_duration) && derived_duration > 0.0
            duration_seconds = derived_duration
            duration_source = :derived_from_hours
        end
    end
    if duration_seconds === nothing
        details = join(unique(vcat(duration_reasons, start_reasons, end_reasons)), "; ")
        push!(
            errors,
            "$(_meteo_row_label(row_index)) needs a positive `duration`/`step_duration` or a valid `hour_start` + `hour_end` interval ($(details)).",
        )
    end

    end_hour = nothing
    end_source = nothing
    if start_hour !== nothing && duration_seconds !== nothing
        end_hour = start_hour + duration_seconds / 3600.0
        end_source = duration_source
    elseif start_hour !== nothing && explicit_end_hour !== nothing && !end_before_start
        end_hour = explicit_end_hour
        end_source = explicit_end_source
    end

    date_source === nothing || (sources[:date] = date_source)
    start_source === nothing || (sources[:start_hour] = start_source)
    end_source === nothing || (sources[:end_hour] = end_source)
    duration_source === nothing || (sources[:duration_seconds] = duration_source)

    interval_reasons = String[]
    date === nothing && append!(interval_reasons, date_reasons)
    start_hour === nothing && append!(interval_reasons, start_reasons)
    end_hour === nothing && append!(interval_reasons, end_reasons)
    return (
        duration_seconds=duration_seconds,
        date=date,
        start_hour=start_hour,
        end_hour=end_hour,
        interval_reasons=unique(interval_reasons),
    )
end

_meteo_resolve_temporal!(errors::Vector{String}, row, row_index::Int, sources) =
    _meteo_resolve_temporal!(errors, row, row_index, sources, _meteo_row_names(row))

function _meteo_require_solar_interval!(errors, temporal, row_index::Int)
    for reason in temporal.interval_reasons
        push!(
            errors,
            "$(_meteo_row_label(row_index)) cannot reconstruct the solar interval because $(reason).",
        )
    end
    return nothing
end

function _meteo_has_solar_interval(temporal)
    temporal.date !== nothing && temporal.start_hour !== nothing &&
        temporal.end_hour !== nothing
end

function _meteo_has_resolved_solar_interval(step::ResolvedMeteoStep)
    step.date !== nothing && step.start_hour !== nothing && step.end_hour !== nothing
end

function _meteo_energy_tolerance(a::Float64, b::Float64)
    max(_METEO_ENERGY_ATOL, _METEO_ENERGY_RTOL * max(abs(a), abs(b)))
end

function _meteo_energy_agrees(a::Float64, b::Float64)
    abs(a - b) <= _meteo_energy_tolerance(a, b)
end

function _meteo_check_range!(
    errors::Vector{String},
    row_index::Int,
    name::Symbol,
    value::Union{Nothing,Float64},
    lower::Float64,
    upper::Float64;
    upper_inclusive::Bool=true,
)
    value === nothing && return
    in_upper = upper_inclusive ? value <= upper : value < upper
    value >= lower && in_upper && return
    bracket = upper_inclusive ? "]" : ")"
    push!(
        errors,
        "$(_meteo_row_label(row_index)) has `$(name)`=$(repr(value)), outside the expected range [$(lower), $(upper)$(bracket). Pass `check_boundaries=false` only if this value is deliberate.",
    )
end

function _meteo_check_nonnegative!(
    errors::Vector{String},
    row_index::Int,
    name::Symbol,
    value::Union{Nothing,Float64},
)
    value === nothing && return
    value >= 0.0 && return
    push!(
        errors,
        "$(_meteo_row_label(row_index)) has `$(name)`=$(repr(value)); irradiance must be non-negative. Pass `check_boundaries=false` only if this value is deliberate.",
    )
end

function _meteo_warn_high_irradiance!(
    warnings::Vector{String},
    row_index::Int,
    name::Symbol,
    value::Union{Nothing,Float64},
)
    value === nothing && return
    value <= _METEO_HIGH_IRRADIANCE_WARNING && return
    push!(
        warnings,
        "$(_meteo_row_label(row_index)) has unusually high `$(name)`=$(repr(value)) W m^-2; check the value and units.",
    )
end

function _meteo_resolve_extra_bands!(
    errors::Vector{String},
    warnings::Vector{String},
    row,
    row_index::Int,
    row_names,
    metadata,
    metadata_names;
    check_boundaries::Bool,
)
    function collect_finite(container, names)
        by_band = Dict{String,Vector{Pair{Symbol,Float64}}}()
        for name in names
            text = String(name)
            upper = uppercase(text)
            startswith(upper, "RI_") || continue
            endswith(upper, "_F") || continue
            upper in ("RI_SW_F", "RI_PAR_F", "RI_NIR_F", "RI_TIR_F") &&
                continue
            band = upper[4:(end - 2)]
            isempty(band) && continue
            value, _ = _meteo_parse_number(getproperty(container, name))
            value === nothing && continue
            push!(get!(by_band, band, Pair{Symbol,Float64}[]), name => value)
        end
        return by_band
    end

    row_values = collect_finite(row, row_names)
    metadata_values = collect_finite(metadata, metadata_names)
    extras = Dict{String,Float64}()
    for band in sort!(collect(union(keys(row_values), keys(metadata_values))))
        finite = get(row_values, band, Pair{Symbol,Float64}[])
        isempty(finite) &&
            (finite = get(metadata_values, band, Pair{Symbol,Float64}[]))
        reference = first(finite).second
        if any(item -> !_meteo_energy_agrees(item.second, reference), @view(finite[2:end]))
            details = join(("$(item.first)=$(repr(item.second))" for item in finite), ", ")
            push!(
                errors,
                "$(_meteo_row_label(row_index)) has conflicting aliases for extra band `$(band)`: $(details).",
            )
            continue
        end
        canonical_index = findfirst(
            item -> startswith(String(item.first), "Ri_") &&
                    endswith(String(item.first), "_f"),
            finite,
        )
        selected = canonical_index === nothing ? first(finite) : finite[canonical_index]
        extras[band] = selected.second
        if check_boundaries
            _meteo_check_nonnegative!(
                errors,
                row_index,
                Symbol("Ri_$(band)_f"),
                selected.second,
            )
            _meteo_warn_high_irradiance!(
                warnings,
                row_index,
                Symbol("Ri_$(band)_f"),
                selected.second,
            )
        end
    end
    return extras
end

function _meteo_resolve_clearness_radiation!(
    errors::Vector{String},
    row_index::Int,
    clearness::Float64,
    latitude::Union{Nothing,Float64},
    temporal,
)
    has_interval = _meteo_has_solar_interval(temporal)
    has_interval || _meteo_require_solar_interval!(errors, temporal, row_index)
    if latitude === nothing || !has_interval
        push!(
            errors,
            "$(_meteo_row_label(row_index)) can only derive radiation from `clearness` when finite latitude, date, start time, and end/duration are available.",
        )
        return nothing
    end
    try
        return _global_wm2_from_clearness(
            clearness,
            deg2rad(latitude),
            Dates.dayofyear(temporal.date),
            temporal.start_hour,
            temporal.end_hour,
        )
    catch err
        push!(
            errors,
            "$(_meteo_row_label(row_index)) could not derive radiation from `clearness`: $(sprint(showerror, err))",
        )
        return nothing
    end
end

function _resolve_meteo_step(
    row,
    options::LightOptions=LightOptions();
    row_index::Int=1,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    errors = String[]
    warnings = String[]
    infos = String[]
    sources = Dict{Symbol,Symbol}()
    row_names = _meteo_row_names(row)
    metadata = _meteo_metadata(row)
    metadata_names = propertynames(metadata)

    resolved = Dict{Symbol,_MeteoLogicalValue}()
    for logical_name in propertynames(_METEO_ALIASES)
        value = _meteo_resolve_logical_value!(
            errors,
            row,
            row_index,
            logical_name,
            row_names,
            metadata,
            metadata_names,
        )
        resolved[logical_name] = value
        value.source === nothing || (sources[logical_name] = value.source)
    end

    sw = resolved[:ri_sw_f].value
    par = resolved[:ri_par_f].value
    nir = resolved[:ri_nir_f].value
    clearness = resolved[:clearness].value
    latitude = resolved[:latitude].value
    sun_azimuth = resolved[:sun_azimuth].value
    sun_elevation = resolved[:sun_elevation].value
    explicit_direct_fraction = resolved[:direct_fraction].value
    direct_sw = resolved[:ri_sw_direct].value
    diffuse_sw = resolved[:ri_sw_diffuse].value
    nir_is_required_source =
        options.nir_interception || (sw === nothing && par === nothing && nir !== nothing)

    temporal = _meteo_resolve_temporal!(errors, row, row_index, sources, row_names)
    duration_seconds = temporal.duration_seconds

    if check_boundaries
        _meteo_check_nonnegative!(errors, row_index, :Ri_SW_f, sw)
        _meteo_check_nonnegative!(errors, row_index, :Ri_PAR_f, par)
        nir_is_required_source &&
            _meteo_check_nonnegative!(errors, row_index, :Ri_NIR_f, nir)
        _meteo_check_nonnegative!(errors, row_index, :Ri_SW_f_direct, direct_sw)
        _meteo_check_nonnegative!(errors, row_index, :Ri_SW_f_diffuse, diffuse_sw)
        _meteo_check_range!(errors, row_index, :clearness, clearness, 0.0, 1.0)
        _meteo_check_range!(errors, row_index, :latitude, latitude, -90.0, 90.0)
        _meteo_check_range!(
            errors,
            row_index,
            :sun_azimuth,
            sun_azimuth,
            0.0,
            360.0;
            upper_inclusive=false,
        )
        _meteo_check_range!(
            errors,
            row_index,
            :sun_elevation,
            sun_elevation,
            -90.0,
            90.0,
        )
        _meteo_check_range!(
            errors,
            row_index,
            :direct_fraction,
            explicit_direct_fraction,
            0.0,
            1.0,
        )
        for (name, value) in (
            (:Ri_SW_f, sw),
            (:Ri_PAR_f, par),
            (:Ri_NIR_f, nir_is_required_source ? nir : nothing),
            (:Ri_SW_f_direct, direct_sw),
            (:Ri_SW_f_diffuse, diffuse_sw),
        )
            _meteo_warn_high_irradiance!(warnings, row_index, name, value)
        end
    end

    total_candidates = Pair{Symbol,Float64}[]
    sw === nothing || push!(total_candidates, :ri_sw_f => sw)
    if par !== nothing && nir !== nothing &&
       (options.nir_interception || sw === nothing)
        push!(total_candidates, :spectral_sum => par + nir)
    end
    if direct_sw !== nothing && diffuse_sw !== nothing
        push!(total_candidates, :direct_diffuse_sum => direct_sw + diffuse_sw)
    end
    if direct_sw !== nothing && explicit_direct_fraction !== nothing &&
       0.0 < explicit_direct_fraction <= 1.0
        push!(
            total_candidates,
            :direct_and_fraction => direct_sw / explicit_direct_fraction,
        )
    end
    if diffuse_sw !== nothing && explicit_direct_fraction !== nothing &&
       0.0 <= explicit_direct_fraction < 1.0
        push!(
            total_candidates,
            :diffuse_and_fraction => diffuse_sw / (1.0 - explicit_direct_fraction),
        )
    end

    total_conflict = false
    if length(total_candidates) > 1
        reference = first(total_candidates)
        for candidate in @view(total_candidates[2:end])
            _meteo_energy_agrees(reference.second, candidate.second) && continue
            tolerance = _meteo_energy_tolerance(reference.second, candidate.second)
            push!(
                errors,
                "$(_meteo_row_label(row_index)) has conflicting total shortwave estimates: `$(reference.first)`=$(repr(reference.second)) and `$(candidate.first)`=$(repr(candidate.second)) differ by more than $(repr(tolerance)) W m^-2.",
            )
            total_conflict = true
        end
    end

    radiation_source = :unresolved
    if !isempty(total_candidates)
        selected = first(total_candidates)
        sw = selected.second
        radiation_source = selected.first
        if selected.first != :ri_sw_f || !haskey(sources, :ri_sw_f)
            sources[:ri_sw_f] = selected.first
        end
    elseif par !== nothing
        sw = par / _METEO_SOLAR_TO_PAR
        radiation_source = :par_only
        sources[:ri_sw_f] = :derived_from_par
    elseif nir !== nothing
        sw = nir / _METEO_SOLAR_TO_NIR
        radiation_source = :nir_only
        sources[:ri_sw_f] = :derived_from_nir
    elseif clearness !== nothing
        sw = _meteo_resolve_clearness_radiation!(
            errors,
            row_index,
            clearness,
            latitude,
            temporal,
        )
        radiation_source = :clearness
        sw === nothing || (sources[:ri_sw_f] = :derived_from_clearness)
    else
        unavailable = Pair{Symbol,String}[]
        for name in (:ri_sw_f, :ri_par_f, :ri_nir_f, :clearness)
            append!(unavailable, resolved[name].unavailable)
        end
        detail = isempty(unavailable) ? "all sources are absent" : join(
            ("$(item.first)=$(item.second)" for item in unavailable),
            ", ",
        )
        push!(
            errors,
            "$(_meteo_row_label(row_index)) cannot resolve incident radiation ($(detail)). Provide a finite `Ri_SW_f`, `Ri_PAR_f`, `Ri_NIR_f`, or `clearness` value.",
        )
    end

    if sw !== nothing && !total_conflict
        if options.nir_interception
            if par !== nothing && nir !== nothing
                # Both measured bands are preserved after the consistency check above.
            elseif par !== nothing
                residual = sw - par
                if residual < -_meteo_energy_tolerance(sw, par)
                    push!(
                        errors,
                        "$(_meteo_row_label(row_index)) has `Ri_PAR_f`=$(repr(par)) greater than resolved `Ri_SW_f`=$(repr(sw)).",
                    )
                else
                    nir = max(residual, 0.0)
                    sources[:ri_nir_f] = :derived_as_sw_minus_par
                end
            elseif nir !== nothing
                residual = sw - nir
                if residual < -_meteo_energy_tolerance(sw, nir)
                    push!(
                        errors,
                        "$(_meteo_row_label(row_index)) has `Ri_NIR_f`=$(repr(nir)) greater than resolved `Ri_SW_f`=$(repr(sw)).",
                    )
                else
                    par = max(residual, 0.0)
                    sources[:ri_par_f] = :derived_as_sw_minus_nir
                end
            else
                par = _METEO_SOLAR_TO_PAR * sw
                nir = _METEO_SOLAR_TO_NIR * sw
                sources[:ri_par_f] = :derived_from_sw
                sources[:ri_nir_f] = :derived_from_sw
            end
        else
            if par === nothing
                par = _METEO_SOLAR_TO_PAR * sw
                sources[:ri_par_f] = :derived_from_sw
            elseif par > sw + _meteo_energy_tolerance(sw, par)
                push!(
                    errors,
                    "$(_meteo_row_label(row_index)) has `Ri_PAR_f`=$(repr(par)) greater than resolved `Ri_SW_f`=$(repr(sw)).",
                )
            end
            nir = 0.0
            sources[:ri_nir_f] = :disabled
        end
    end

    if check_boundaries && sw !== nothing && par !== nothing && nir !== nothing
        _meteo_check_nonnegative!(errors, row_index, :Ri_SW_f, sw)
        _meteo_check_nonnegative!(errors, row_index, :Ri_PAR_f, par)
        options.nir_interception &&
            _meteo_check_nonnegative!(errors, row_index, :Ri_NIR_f, nir)
        _meteo_warn_high_irradiance!(warnings, row_index, :Ri_SW_f, sw)
        _meteo_warn_high_irradiance!(warnings, row_index, :Ri_PAR_f, par)
        options.nir_interception &&
            _meteo_warn_high_irradiance!(warnings, row_index, :Ri_NIR_f, nir)
    end

    if sw !== nothing
        if direct_sw !== nothing && diffuse_sw === nothing
            residual = sw - direct_sw
            if residual < -_meteo_energy_tolerance(sw, direct_sw)
                push!(
                    errors,
                    "$(_meteo_row_label(row_index)) has direct shortwave $(repr(direct_sw)) greater than total shortwave $(repr(sw)).",
                )
            else
                diffuse_sw = max(residual, 0.0)
                sources[:ri_sw_diffuse] = :derived_as_sw_minus_direct
            end
        elseif diffuse_sw !== nothing && direct_sw === nothing
            residual = sw - diffuse_sw
            if residual < -_meteo_energy_tolerance(sw, diffuse_sw)
                push!(
                    errors,
                    "$(_meteo_row_label(row_index)) has diffuse shortwave $(repr(diffuse_sw)) greater than total shortwave $(repr(sw)).",
                )
            else
                direct_sw = max(residual, 0.0)
                sources[:ri_sw_direct] = :derived_as_sw_minus_diffuse
            end
        end
    end

    component_fraction = nothing
    if direct_sw !== nothing && diffuse_sw !== nothing
        component_total = direct_sw + diffuse_sw
        component_fraction = component_total > 0.0 ? direct_sw / component_total : 0.0
    end
    if explicit_direct_fraction !== nothing && component_fraction !== nothing &&
       sw !== nothing && sw > 0.0 &&
       !isapprox(explicit_direct_fraction, component_fraction; rtol=0.0, atol=0.01)
        push!(
            errors,
            "$(_meteo_row_label(row_index)) has `direct_fraction`=$(repr(explicit_direct_fraction)), but direct/diffuse components imply $(repr(component_fraction)).",
        )
    end

    direct_fraction =
        explicit_direct_fraction === nothing ? component_fraction : explicit_direct_fraction
    diffuse_fraction = direct_fraction === nothing ? nothing : 1.0 - direct_fraction
    partition_source =
        explicit_direct_fraction !== nothing ? :explicit_fraction :
        component_fraction !== nothing ? :direct_diffuse_components : :automatic

    can_auto_sun = latitude !== nothing && _meteo_has_solar_interval(temporal)
    solar_geometry_source = :automatic
    if sun_azimuth !== nothing && sun_elevation !== nothing
        solar_geometry_source = :explicit
    elseif sun_azimuth !== nothing || sun_elevation !== nothing
        if can_auto_sun
            provided = sun_azimuth === nothing ? :sun_elevation : :sun_azimuth
            push!(
                warnings,
                "$(_meteo_row_label(row_index)) provides only `$(provided)`; the partial explicit position is ignored and both sun coordinates will be reconstructed.",
            )
            sun_azimuth = nothing
            sun_elevation = nothing
        else
            _meteo_require_solar_interval!(errors, temporal, row_index)
            missing_name = sun_azimuth === nothing ? :sun_azimuth : :sun_elevation
            push!(
                errors,
                "$(_meteo_row_label(row_index)) is missing `$(missing_name)` and lacks finite latitude/date/interval data needed to reconstruct the sun position.",
            )
        end
    elseif !can_auto_sun && !(sw !== nothing && iszero(sw))
        _meteo_require_solar_interval!(errors, temporal, row_index)
        columns = join(string.(_meteo_row_names(row)), ", ")
        push!(
            errors,
            "$(_meteo_row_label(row_index)) needs finite `latitude` metadata/column unless both `sun_azimuth` and `sun_elevation` are provided. Available columns: $(columns).",
        )
    elseif !can_auto_sun && sw !== nothing && iszero(sw)
        sun_azimuth = 0.0
        sun_elevation = -90.0
        solar_geometry_source = :zero_forcing_default
        push!(
            infos,
            "$(_meteo_row_label(row_index)) uses a below-horizon default sun position because resolved radiation is zero.",
        )
    end

    if solar_geometry_source == :explicit && !can_auto_sun && sw !== nothing && sw > 0.0 &&
       direct_fraction === nothing && clearness === nothing
        push!(
            errors,
            "$(_meteo_row_label(row_index)) needs finite `direct_fraction`, direct/diffuse components, or `clearness` when explicit sun coordinates are used without date/time/latitude data.",
        )
    end

    extras = _meteo_resolve_extra_bands!(
        errors,
        warnings,
        row,
        row_index,
        row_names,
        metadata,
        metadata_names;
        check_boundaries=check_boundaries,
    )

    step = nothing
    if sw !== nothing && par !== nothing && nir !== nothing && duration_seconds !== nothing &&
       !total_conflict
        step = ResolvedMeteoStep(
            row_index,
            duration_seconds,
            temporal.date,
            temporal.start_hour,
            temporal.end_hour,
            sw,
            par,
            nir,
            clearness,
            latitude,
            sun_azimuth,
            sun_elevation,
            direct_fraction,
            diffuse_fraction,
            direct_sw,
            diffuse_sw,
            extras,
            radiation_source,
            solar_geometry_source,
            partition_source,
            sources,
        )
    end

    return MeteoResolutionResult(step, unique(errors), unique(warnings), unique(infos))
end

function _resolved_meteo_step_or_error(
    row,
    options::LightOptions=LightOptions();
    row_index::Int=1,
    check_boundaries::Bool=options.check_meteo_boundaries,
)
    result = _resolve_meteo_step(
        row,
        options;
        row_index=row_index,
        check_boundaries=check_boundaries,
    )
    if !isempty(result.errors) || result.step === nothing
        details = isempty(result.errors) ? "Meteo resolution did not produce a step." :
                  join(("- " * error for error in result.errors), "\n")
        error("Invalid meteo input:\n$(details)")
    end
    return result.step
end
