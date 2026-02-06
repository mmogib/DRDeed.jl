module DRDeed
include("dependents.jl")
include("types.jl")
include("data.jl")
include("models/models.jl")
include("presentations.jl")
include("utils.jl")


# solveit,
export getDRDeedData,
  DeedData,
  DeedSolution,
  getLambda,
  getDemand,
  getCostCoefficients,
  getProvidersPowerLimits,
  getMatrixB,
  getCustomersFunctionData,
  getBudgetLimit,
  data2DataFrame,
  solution2DataFrame,
  getGTDRDeedData,
  slim_data,
  saveModel,
  loadModel,
  solveModel,
  SuccessResult,
  FailResult,
  BranchData,
  NetworkData,
  getIEEE30Data,
  getSaudiDeedData,
  getSaudiGTDeedData
end # module DRDeed
