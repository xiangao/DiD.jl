# ETWFE.jl — DEPRECATED

**This package is deprecated. Use [`Panelest.jl`](https://github.com/xiangao/Panelest.jl)'s
`etwfe()` + `emfx()` instead — it implements the same Wooldridge (2021)
extended two-way fixed effects (ETWFE) estimator, with a correct
identification check and test coverage that this package never had.**

```julia
using Panelest, DataFrames

df = Panelest.dataset("mpdta")
m  = Panelest.etwfe(df, @formula(lemp ~ lpop); gvar = :first_treat, tvar = :year)

Panelest.emfx(m, type = "simple")   # overall ATT
Panelest.emfx(m, type = "event")    # event-study
Panelest.emfx(m, type = "calendar") # by calendar year
```

## Why this package is deprecated

`att_gt` in this package was always a stub — it never estimated anything and
returned hardcoded placeholder numbers (see `src/ETWFE.jl`). The only real
functionality here, `emfx()` parsing raw `Panelest.feols` coefficient names,
had no test coverage exercising that code path and is superseded by
`Panelest.jl`'s own `etwfe()` + `emfx(::ETWFEResult)`, which wraps the
regression construction, validates that every cohort has a pre-treatment
period (dropping and warning about cohorts that don't, since their effect
isn't identified), and is covered by regression tests.

`att_gt` now raises an error directing you to `Panelest.jl` instead of
silently returning fake numbers.

## References

- Callaway, B. & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
- Wooldridge, J. M. (2021). Two-way fixed effects, the two-way Mundlak regression, and difference-in-differences estimators.
