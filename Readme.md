# DRDeed package

## Install

```julia
add https://github.com/mmogib/DRDeed.jl.git
```

## Usage

```julia
using DRDeed
sol = gtdrdeed(5, 6, 24)(0.3 * ones(3))
solrep = solution2DataFrame(sol.solution)
vscodedisplay.(map(l -> solrep[l], collect(keys(solrep))))
```
