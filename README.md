# ETWFE.jl

This is a small Julia implementation of the extended two-way fixed effects
(ETWFE) estimator for staggered difference-in-differences. The setup follows
Wooldridge (2021) and is close in spirit to the `did` workflow of
Callaway and Sant'Anna (2021).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/xiangao/DiD.jl")
```

## Usage

```julia
using ETWFE, DataFrames

# Load the built-in mpdta panel dataset.
df = dataset("mpdta")

# Estimate cohort-by-time ATT(g,t) cells.
res = att_gt(df, :lemp, :year, :countyreal, :first_treat)

# Aggregate the cells.
es  = emfx(res, type="event")     # event-study relative to treatment
cal = emfx(res, type="calendar")  # by calendar year
ov  = emfx(res, type="overall")   # single overall ATT
```

## Functions

- `att_gt(data, y, t, id, g; ...)` — estimate ATT(g,t) for each cohort × time cell
- `emfx(model; type, level)` — aggregate ATT(g,t) cells to event-study, calendar, or overall
- `dataset(name)` — load built-in datasets (`"mpdta"`)

## References

- Callaway, B. & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
- Wooldridge, J. M. (2021). Two-way fixed effects, the two-way Mundlak regression, and difference-in-differences estimators.
