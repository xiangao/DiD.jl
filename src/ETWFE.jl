module ETWFE

using DataFrames
using LinearAlgebra
using Statistics

export att_gt, emfx, dataset

# ---------------------------------------------------------------------------
# Dataset loader
# ---------------------------------------------------------------------------

"""
    dataset("mpdta")

Return a balanced panel dataset mimicking Callaway & Sant'Anna's mpdta:
500 counties × 4 years (2004–2007), with three treated cohorts (2004, 2006, 2007)
and a never-treated group.  Includes `first_treat_str` and `year_str` string columns
required for ETWFE cohort×time interactions via StatsModels.
"""
function dataset(name::String)
    if name == "mpdta"
        years     = [2004, 2005, 2006, 2007]
        n_per_cohort = 125          # 125 counties each → 500 total
        cohorts   = vcat(fill(2004, n_per_cohort),
                         fill(2006, n_per_cohort),
                         fill(2007, n_per_cohort),
                         fill(0,    n_per_cohort))   # 0 = never treated
        n_counties = length(cohorts)

        countyreal  = repeat(1:n_counties, inner = length(years))
        year_vec    = repeat(years,        outer = n_counties)
        first_treat = repeat(cohorts,      inner = length(years))

        # Treatment indicator: 1 iff county is in a treated cohort AND year ≥ cohort year
        treat = [g > 0 && t >= g ? 1 : 0 for (g, t) in zip(first_treat, year_vec)]

        # Synthetic outcome: unit FE + time trend + treatment effect ≈ 0.3 + noise
        # Deterministic noise for reproducibility in the textbook
        unit_fe   = repeat([sin(i * 0.37) for i in 1:n_counties], inner = length(years))
        time_fe   = repeat([0.0, 0.08, 0.16, 0.24], outer = n_counties)
        te_true   = 0.30 .* treat
        eps       = [cos(i * 1.37) * 0.4 for i in 1:(n_counties * length(years))]
        lemp      = unit_fe .+ time_fe .+ te_true .+ eps
        lpop      = repeat([10.0 + sin(i * 0.53) for i in 1:n_counties],
                           inner = length(years))

        # String columns needed for StatsModels categorical dummies in ETWFE
        first_treat_str = string.(first_treat)
        year_str        = string.(year_vec)

        return DataFrame(
            countyreal      = countyreal,
            year            = year_vec,
            year_str        = year_str,
            first_treat     = first_treat,
            first_treat_str = first_treat_str,
            treat           = treat,
            lemp            = lemp,
            lpop            = lpop,
        )
    end
    error("Unknown dataset: '$name'")
end

# ---------------------------------------------------------------------------
# att_gt stub (Callaway & Sant'Anna API — not implemented)
# ---------------------------------------------------------------------------

function att_gt(formula, data; gname, tname, idname = nothing, xformla = nothing,
                panel = true, control_group = "nevertreated")
    @warn "att_gt is a stub. For ETWFE estimation use feols(...) from Panelest + emfx()."
    return (estimates = zeros(5),
            group     = [2004, 2004, 2006, 2006, 2006],
            time      = [2004, 2005, 2004, 2005, 2006])
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Normal quantile for common confidence levels
function _z_value(level::Real)
    level ≈ 90.0 && return 1.64485
    level ≈ 95.0 && return 1.95996
    level ≈ 99.0 && return 2.57583
    return 1.95996  # fallback: 95 %
end

# Weighted average + SE from VCV sub-matrix (equal weights)
function _aggregate(β::Vector, V::Matrix, idxs::Vector{Int})
    n  = length(idxs)
    w  = fill(1.0 / n, n)
    τ  = dot(w, β[idxs])
    se = sqrt(max(0.0, w' * V[idxs, idxs] * w))
    return τ, se
end

# ---------------------------------------------------------------------------
# emfx — aggregate ETWFE coefficients into ATTs
# ---------------------------------------------------------------------------

"""
    emfx(model; type="simple", post_only=true, gvar="first_treat_str",
          tvar="year_str", level=95.0)

Aggregate cohort×time ETWFE coefficients from a Panelest model into
average treatment effects.

# Arguments
- `model`    : a `PanelestModel` from `Panelest.feols`, estimated with an
               ETWFE formula such as
               `@formula(y ~ x + first_treat_str & year_str + fe(id) + fe(t))`.
- `type`     : aggregation type — `"simple"` (overall ATT), `"event"` (event study),
               or `"calendar"` (by calendar year).
- `post_only`: if `true` (default), restrict to post-treatment cells (T ≥ G).
- `gvar`     : name of the cohort column used in the interaction (default `"first_treat_str"`).
- `tvar`     : name of the time column used in the interaction (default `"year_str"`).
- `level`    : confidence level in percent (default 95.0).

# Returns
A `DataFrame` with columns `estimate`, `std_error`, `conf_low`, `conf_high`
and a group column (`type`, `event`, or `calendar`).
"""
function emfx(model;
              type      = "simple",
              post_only = true,
              gvar      = "first_treat_str",
              tvar      = "year_str",
              level     = 95.0)

    names = model.coefnames
    β     = model.beta
    V     = model.vcov
    z     = _z_value(level)

    # StatsModels interaction coefnames:  "gvar: G & tvar: T"  (or reversed)
    pat1 = Regex("^" * gvar * ": (\\S+) & " * tvar * ": (\\S+)\$")
    pat2 = Regex("^" * tvar * ": (\\S+) & " * gvar * ": (\\S+)\$")

    # Collect (coefficient index, cohort G, period T) for each interaction coef
    cells = NamedTuple{(:idx, :G, :T), Tuple{Int,Int,Int}}[]
    for (i, n) in enumerate(names)
        m = match(pat1, n)
        if !isnothing(m)
            G = parse(Int, m.captures[1])
            T = parse(Int, m.captures[2])
            push!(cells, (idx = i, G = G, T = T))
            continue
        end
        m = match(pat2, n)
        if !isnothing(m)
            T = parse(Int, m.captures[1])
            G = parse(Int, m.captures[2])
            push!(cells, (idx = i, G = G, T = T))
        end
    end

    # Exclude never-treated cohort (G = 0)
    filter!(c -> c.G > 0, cells)
    # Optionally restrict to post-treatment cells (T ≥ G)
    post_only && filter!(c -> c.T >= c.G, cells)

    if isempty(cells)
        error(
            "emfx: no cohort×time interaction coefficients found.\n" *
            "Expected coefnames like '$gvar: 2006 & $tvar: 2007'.\n" *
            "Check that gvar='$gvar' and tvar='$tvar' match the column names in your formula."
        )
    end

    if type == "simple"
        idxs     = [c.idx for c in cells]
        τ, se    = _aggregate(β, V, idxs)
        return DataFrame(
            type      = ["ATT"],
            estimate  = [τ],
            std_error = [se],
            conf_low  = [τ - z * se],
            conf_high = [τ + z * se],
        )

    elseif type == "event"
        events = sort(unique(c.T - c.G for c in cells))
        df = DataFrame(event     = Int[],
                       estimate  = Float64[],
                       std_error = Float64[],
                       conf_low  = Float64[],
                       conf_high = Float64[])
        for e in events
            idxs  = [c.idx for c in cells if c.T - c.G == e]
            τ, se = _aggregate(β, V, idxs)
            push!(df, (event = e, estimate = τ, std_error = se,
                       conf_low = τ - z * se, conf_high = τ + z * se))
        end
        return df

    elseif type == "calendar"
        cal_years = sort(unique(c.T for c in cells))
        df = DataFrame(calendar  = Int[],
                       estimate  = Float64[],
                       std_error = Float64[],
                       conf_low  = Float64[],
                       conf_high = Float64[])
        for t in cal_years
            idxs  = [c.idx for c in cells if c.T == t]
            τ, se = _aggregate(β, V, idxs)
            push!(df, (calendar = t, estimate = τ, std_error = se,
                       conf_low = τ - z * se, conf_high = τ + z * se))
        end
        return df

    else
        error("emfx: unknown type '$type'. Choose \"simple\", \"event\", or \"calendar\".")
    end
end

end  # module DiD
