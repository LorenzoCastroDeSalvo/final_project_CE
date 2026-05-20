# Optimal Monetary Policy According to HANK

CCA CompEcon term project by Lorenzo Castro De Salvo and Paolo Diotallevi.

This repository is a Julia replication-light project for selected computational
objects from Sushant Acharya, Edouard Challe, and Keshav Dogra, "Optimal
Monetary Policy According to HANK." The project focuses on the tractable HANK
model's closed-form household policies, the inequality-related monetary-policy
weight `Omega`, and policy/value-function exhibits generated from the analytical
solution.

- Paper: [Optimal Monetary Policy According to HANK](https://www.aeaweb.org/articles?id=10.1257/aer.20200239)
- Replication package: [openICPSR project 184261](https://doi.org/10.3886/E184261V1)
- Course: CCA CompEcon

## Installation

From the repository root:

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Run

The package entry point is `run()`:

```julia
using HANKPolicies
run()
```

This writes reproducible plots and data to `images/`, including the comparative
statics for `Omega`, household policy functions, flow utility, and an evaluated
value function under the analytical policy rules.

You can also run the same entry point from the shell:

```powershell
julia --project=. scripts/run_all.jl
```

## Tests

```julia
import Pkg
Pkg.test()
```

The tests check that the package loads, core analytical relationships hold, and
`run()` executes on a small test grid.

## Report and Online Documentation

The project report is in `report.qmd`. Render it locally with:

```powershell
quarto render report.qmd
```

The GitHub Actions workflow in `.github/workflows/publish.yml` renders the
Quarto report and publishes it to GitHub Pages on pushes to `main`. In the
repository settings, set Pages source to "GitHub Actions."

Generated figures are stored in `images/`. The original notebook-export script
used during development is kept in `scripts/copia-notebook-complete-savefigs-oswald.jl`,
and its saved paper-style figures are preserved in `figures/`.
