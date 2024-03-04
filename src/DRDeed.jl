module DRDeed
include("dependents.jl")
include("types.jl")
include("data.jl")
include("models/models.jl")
include("presentations.jl")


# solveit,
export getDRDeedData,
  DeedData,
  DeedSolution,
  deed,
  drdeed,
  getLambda,
  getDemand,
  getCostCoefficients,
  getProvidersPowerLimits,
  getMatrixB,
  getCustomersFunctionData,
  getBudgetLimit,
  slimmodel,
  gtdrdeed,
  data2DataFrame,
  solution2DataFrame,
  getGTDRDeedData
end # module DRDeed
