# CLAUDE.md — DRDeed.jl Core Package

**Parent:** `../CLAUDE.md` (Project Orchestrator)

**Status:** ✅ **v0.5.0** — All features implemented

**Repository:** `https://github.com/mmogib/DRDeed.jl`

---

## Role

Core Julia optimization package implementing DEED, DR-DEED, and Game-Theoretic DEED solvers.

---

## Package Structure

```
src/
├── DRDeed.jl          # Module entry
├── types.jl           # Data/solution types
├── data.jl            # Data generation (PJM + Saudi)
├── utils.jl           # File I/O
└── models/deedmodel/
    ├── deed.jl        # Basic DEED
    ├── drdeed.jl      # DR-DEED
    └── gtdrdeed.jl    # GT-DEED (main)
```

---

## Models

| Model | Variables | Description |
|-------|-----------|-------------|
| **DEED** | q[j,t] | Basic economic emission dispatch |
| **DR-DEED** | q, χ, ω | + demand response curtailment |
| **GT-DEED** | q, p, s, χ, ω, x, f, h, y | Full game-theoretic (9 classes) |

---

## Public API

```julia
# Model constructors:
deed(generators, periods)
drdeed(customers, generators, periods)
gtdrdeed(customers, generators, periods)

# Data:
getGTDRDeedData(customers, generators, periods)  # PJM
getSaudiGTDeedData()                              # Saudi (v0.5.0)

# Solving:
model = gtdrdeed(5, 6, 24)
result = model[:scalarized](; w=[1/3, 1/3, 1/3])

# Results:
result.solution.Cost, .Emission, .Utility, .Losst
result.solution.q, .χ, .ω, .p, .s, .x, .f, .h
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| JuMP | Optimization modeling |
| Ipopt | Interior-point NLP solver |
| DataFrames, XLSX | Data I/O |

---

## Related Sub-Contexts

- `../papercode/CLAUDE.md` — experiment code using this package
- `../papertext/CLAUDE.md` — paper describing these models
