function data2DataFrame(data::DeedData)
  mxgc = max(data.customers, data.generators)
  df0 = DataFrame(
    customers = repeat([data.customers], mxgc),
    generators = repeat([data.generators], mxgc),
    periods = repeat([data.periods], mxgc),
    UB = repeat([data.UB], mxgc),
  )
  df1 = DataFrame(
    a = data.a,
    b = data.b,
    c = data.c,
    e = data.e,
    f = data.f,
    g = data.g,
    pjmin = data.pjmin,
    pjmax = data.pjmax,
    DR = data.DR,
    UR = data.UR,
  )
  dfB = DataFrame(Dict(zip(["B_generator_$i" for i = 1:data.generators], eachcol(data.B))))
  df2 = DataFrame(CM = data.CM, K1 = data.K1, K2 = data.K2, theta = data.θ)
  dfDt = hcat(
    DataFrame(demand = data.Dt),
    DataFrame(Dict(zip(["lambda_$i" for i = 1:data.customers], eachrow(data.λ)))),
  )
  Dict(
    :Dtotal => dfDt,
    :Bmatrix => dfB,
    :generators => hcat(df0[1:data.generators, :], df1),
    :customers => hcat(df0[1:data.customers, :], df2),
  )

end

function data2DataFrame(data::GTDeedData)
  deedDf = data2DataFrame(data.deedData)
  df1 = DataFrame(
    Dict(zip(["customer_cost_a", "customer_cost_b", "customer_cost_c"], eachcol(data.costCoef))),
  )
  df2 = DataFrame(
    Dict(zip(["customer_$(i)_demand" for i = 1:data.deedData.customers], eachcol(data.CDemandt))),
  )
  df3 = DataFrame(pimin = data.pimin, pimax = data.pimax)
  Dict(deedDf..., :customers_coef => df1, :customers_demand => df2, :customers_powerlimits => df3)
end
"""
  deedSolution::DRDeedSolution
  p::Matrix{Float64}
  x::Matrix{Float64}
  y::Array{Float64,3}
  CustomerCost::Vector{Float64}
"""
function solution2DataFrame(solution::GTDRDeedSolution)
  drdeed = solution2DataFrame(solution.deedSolution)
  p = solution.p
  s = solution.s
  h = solution.h
  f = solution.f
  x = solution.x
  y = solution.y
  customers, generators, periods = size(y)
  yallj = Matrix{Float64}(undef, customers, periods)

  foreach(1:customers) do i
    yallj[i, :] = sum(y[i, j, :] for j = 1:generators)
  end
  cost = solution.customerCost
  pdf = DataFrame(Dict(zip(["customer_$(i)_p" for i = 1:customers], eachrow(p))))
  sdf = DataFrame(Dict(zip(["customer_$(i)_s" for i = 1:customers], eachrow(s[:, 2:end]))))
  hdf = DataFrame(Dict(zip(["customer_$(i)_h" for i = 1:customers], eachrow(h[:, 1:end-1]))))
  fdf = DataFrame(Dict(zip(["customer_$(i)_f" for i = 1:customers], eachrow(f))))
  cdf = DataFrame(Dict(zip(["customer_$(i)_utility_value" for i = 1:customers], eachrow(cost))))
  xdf = DataFrame(Dict(zip(["customer_$(i)_x" for i = 1:customers], eachrow(x))))
  ydf = DataFrame(Dict(zip(["customer_$(i)_y_allj" for i = 1:customers], eachrow(yallj))))
  df = hcat(cdf, pdf, sdf, hdf, fdf, xdf, ydf)
  Dict(drdeed..., :customers => df)
end
function solution2DataFrame(solution::DRDeedSolution)
  χ = solution.χ
  ω = solution.ω
  q = solution.q
  Cost = solution.Cost
  Emission = solution.Emission
  Utility = solution.Utility
  Loss = solution.Losst
  PowerGenerated = sum(q)
  customers = size(χ, 1)
  generators = size(q, 1)
  df1 = DataFrame(
    :cost => [Cost],
    :emission => [Emission],
    :utility => [Utility],
    :loss => sum(Loss),
    :power => [PowerGenerated],
  )
  dfloss = DataFrame(:losst => Loss)
  df2 = DataFrame(Dict(zip(["customer_$(i)_chi" for i = 1:customers], eachrow(χ))))
  df3 = DataFrame(Dict(zip(["customer_$(i)_w" for i = 1:customers], eachrow(ω))))

  df4 = DataFrame(Dict(zip(["generator_$(j)_q" for j = 1:generators], eachrow(q))))

  Dict(:summary => df1, :customers => hcat(dfloss, df2, df3), :generators => hcat(dfloss, df4))
end