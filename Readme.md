# DRDeed.jl

[![Version](https://img.shields.io/badge/version-0.5.0-blue.svg)](https://github.com/mmogib/DRDeed.jl)
[![Julia](https://img.shields.io/badge/Julia-1.10+-purple.svg)](https://julialang.org/)

A Julia package for solving Dynamic Economic Emission Dispatch (DEED) problems with game-theoretic demand response in smart grids.

## Overview

This package implements the optimization models described in:

> **A Stackelberg Game for Multi-Objective Demand Response in Dynamic Economic Emission Dispatch**
>
> Norah Almuraysil, Mohammed Alshahrani, Slim Belhaiza

The package provides three model variants:

| Model | Description |
|-------|-------------|
| **DEED** | Basic dynamic economic emission dispatch |
| **DR-DEED** | DEED with demand response (curtailment only) |
| **SGSD-DEED** | Full Stackelberg game with 9 decision variable classes |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mmogib/DRDeed.jl.git")
```

## Quick Start

```julia
using DRDeed

# Create model: 5 customers, 6 generators, 24 periods
model = gtdrdeed(5, 6, 24)

# Solve with balanced weights [cost, emission, utility]
result = model[:scalarized](; w=[1/3, 1/3, 1/3])

# Access results
println("Cost: ", result.solution.Cost)
println("Emission: ", result.solution.Emission)
println("Utility: ", result.solution.Utility)

# Convert to DataFrame for analysis
using DataFrames
solrep = solution2DataFrame(result.solution)
```

## Available Models

```julia
# Basic DEED (no demand response)
deed(generators, periods)

# DR-DEED with curtailment
drdeed(customers, generators, periods)

# Full game-theoretic model (SGSD-DEED)
gtdrdeed(customers, generators, periods)
```

## Data Generation

```julia
# PJM-based data (default)
data = getGTDRDeedData(customers, generators, periods)

# Saudi Eastern Province data
data = getSaudiGTDeedData()
```

## Dependencies

- **JuMP** - Optimization modeling
- **Ipopt** - Interior-point NLP solver (recommended: v3.14 with MA27)
- **DataFrames** - Data manipulation
- **XLSX** - Excel I/O

## Paper Reproduction

To reproduce all paper results, see the companion repository:
[https://github.com/mmogib/NMSPaperMAY25](https://github.com/mmogib/NMSPaperMAY25)

## License

MIT License
