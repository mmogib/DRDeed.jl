function data2DataFrame(data::SlimData)
  # customers::Int
  # generators::Int
  # periods::Int
  # Caith::Matrix{Float64}
  # Caitp::Matrix{Float64} # rows: customers, cols: periods
  # Caits::Matrix{Float64}
  # Caitf::Matrix{Float64}
  # Cajtp::Matrix{Float64}
  # Cajts::Matrix{Float64}
  # Caijtx::Array{Float64,3}
  # Caijty::Array{Float64,3}
  # DemandDt::Matrix{Float64}
  # CProductionLimitst::Matrix{Float64} # Lit
  # CPowerUpperLimitst::Matrix{Float64} # Kit
  # CStorageUpperLimitst::Matrix{Float64} # Cit
  # PProductionLimitst::Matrix{Float64} # Ljt
  # PPowerUpperLimitst::Matrix{Float64} # Kjt
  # PStorageUpperLimitst::Matrix{Float64} # Cjt
  mxgc = max(data.customers, data.generators)
  df0 = DataFrame(
    num_customers = repeat([data.customers], mxgc),
    num_generators = repeat([data.generators], mxgc),
    periods = repeat([data.periods], mxgc),
  )

  dfch = DataFrame(Dict(zip(["customer_Cai$(i)h" for i = 1:data.periods], eachcol(data.Caith))))
  dfcp = DataFrame(Dict(zip(["customer_Cai$(i)p" for i = 1:data.periods], eachcol(data.Caitp))))
  dfcs = DataFrame(Dict(zip(["customer_Cai$(i)s" for i = 1:data.periods], eachcol(data.Caits))))
  dfcf = DataFrame(Dict(zip(["customer_Cai$(i)f" for i = 1:data.periods], eachcol(data.Caitf))))
  dfcDt = DataFrame(Dict(zip(["customer_Di$(i)" for i = 1:data.periods], eachcol(data.DemandDt))))
  dfcLt = DataFrame(
    Dict(zip(["customer_Li$(i)" for i = 1:data.periods], eachcol(data.CProductionLimitst))),
  )
  dfcKt = DataFrame(
    Dict(zip(["customer_Ki$(i)" for i = 1:data.periods], eachcol(data.CPowerUpperLimitst))),
  )

  dfgp = DataFrame(Dict(zip(["generators_Caj$(i)p" for i = 1:data.periods], eachcol(data.Cajtp))))
  dfgs = DataFrame(Dict(zip(["generators_Caj$(i)s" for i = 1:data.periods], eachcol(data.Cajts))))
  dfgLt = DataFrame(
    Dict(zip(["generators_Lj$(i)" for i = 1:data.periods], eachcol(data.PProductionLimitst))),
  )
  dfgKt = DataFrame(
    Dict(zip(["generators_Kj$(i)" for i = 1:data.periods], eachcol(data.PPowerUpperLimitst))),
  )
  dfgCt = DataFrame(
    Dict(zip(["generators_Cj$(i)" for i = 1:data.periods], eachcol(data.PStorageUpperLimitst))),
  )
  dfcC = map(1:data.customers) do i
    df1 = map(1:data.generators) do j
      d1 = DataFrame(
        Dict(zip(["generators_Ca$(i)$(j)$(t)x" for t = 1:data.periods], data.Caijtx[i, j, :])),
      )
      d2 = DataFrame(
        Dict(zip(["generators_Ca$(i)$(j)$(t)y" for t = 1:data.periods], data.Caijty[i, j, :])),
      )
      hcat(d1, d2)
    end
    hcat(df1...)
  end
  # cxdf = 
  Dict(
    :slim_customers => hcat(df0[1:data.customers, :], dfch, dfcp, dfcs, dfcf, dfcDt, dfcLt, dfcKt),
    :slim_generators => hcat(df0[1:data.generators, :], dfgp, dfgs, dfgLt, dfgKt, dfgCt),
    Dict(zip([Symbol("Cxy_$i") for i = 1:data.customers], dfcC))...,
  )

end

function data2DataFrame(data::DeedData)
  mxgc = max(data.customers, data.generators)
  df0 = DataFrame(
    num_customers = repeat([data.customers], mxgc),
    num_generators = repeat([data.generators], mxgc),
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
  customers = data.deedData.customers
  # generators = data.deedData.generators
  # periods = data.deedData.periods
  deedDf = data2DataFrame(data.deedData)
  fld_names = fieldnames(GTDeedData) |> x -> filter(y -> y != :deedData, x)
  customer_matrix_fld_names = filter(fld_names) do fld
    d = getfield(data, fld)
    if isa(d, Matrix)
      r, _ = size(d)
      r == customers
    else
      false
    end
  end
  customer_matrix_df = map(customer_matrix_fld_names) do fld
    d = getfield(data, fld)
    DataFrame(Dict(zip([String(fld) for i = 1:customers], eachcol(d))))
  end
  customer_vector_fld_names = filter(fld_names) do fld
    d = getfield(data, fld)
    if isa(d, Vector)
      length(d) == customers
    else
      false
    end
  end
  customer_vector_df = map(customer_vector_fld_names) do fld
    d = getfield(data, fld)
    DataFrame(String(fld) => d)
  end

  df3 = DataFrame(pimin = data.pimin, pimax = data.pimax)
  Dict(
    deedDf...,
    :customers_data => hcat(customer_vector_df..., customer_matrix_df...),
    :customers_powerlimits => df3,
  )
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
  gtdr_summary =
    hcat(drdeed[:summary], DataFrame(customer_power = [sum(p)], customer_cost = [sum(cost)]))
  Dict(drdeed..., :gtrdcustomers => df, :gtdr_summary => gtdr_summary)
end
function solution2DataFrame(solution::DRDeedSolution)
  χ = solution.χ
  ω = solution.ω
  q = solution.q
  w = solution.w
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
    :weights => join(w, ","),
  )
  dfloss = DataFrame(:losst => Loss)
  df2 = DataFrame(Dict(zip(["customer_$(i)_chi" for i = 1:customers], eachrow(χ))))
  df3 = DataFrame(Dict(zip(["customer_$(i)_w" for i = 1:customers], eachrow(ω))))

  df4 = DataFrame(Dict(zip(["generator_$(j)_q" for j = 1:generators], eachrow(q))))

  Dict(:summary => df1, :drcustomers => hcat(dfloss, df2, df3), :drgenerators => hcat(dfloss, df4))
end

function solution2DataFrame(solution::DeedSolution)
  P = solution.P
  w = solution.w
  Cost = solution.Cost
  Emission = solution.Emission
  Loss = solution.Loss
  PowerGenerated = sum(P)
  generators = size(P, 1)
  df1 = DataFrame(
    :cost => [Cost],
    :emission => [Emission],
    :loss => sum(Loss),
    :power => [PowerGenerated],
    :weights => join(w, ","),
  )
  dfloss = DataFrame(:losst => Loss)
  df4 = DataFrame(Dict(zip(["generator_$(j)_P" for j = 1:generators], eachrow(P))))

  Dict(:summary => df1, :generators => hcat(dfloss, df4))
end