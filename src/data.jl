
# slim data
function matrix2DF(M::Array{<:Number,2}, name::String)
  DataFrame(Dict(zip(["$(name)_$i" for i = 1:size(M, 2)], eachcol(M))))
end

function slim_random(mn, mx)
  if mn > mx
    println("mn $mn mx $mx")
    println("Error Random Generator")
    return mn
  end
  return mn + rand() * (mx - mn)
end

function slim_data(n::Int, m::Int, T::Int)
  dfit = DataFrame()
  dfjt = DataFrame()
  dfijx = Vector{DataFrame}(undef, T)
  dfijy = Vector{DataFrame}(undef, T)

  aith = zeros(n, T)
  for i = 1:n
    for t = 1:T
      aith[i, t] = slim_random(30, 80) / 100
    end
  end
  dfit = hcat(dfit, matrix2DF(aith, "aith"))
  aitp = zeros(Float64, n, T)
  for i = 1:n
    for t = 1:T
      aitp[i, t] = slim_random(round(aith[i, t] * 100), 100) / 100
    end

  end
  dfit = hcat(dfit, matrix2DF(aitp, "aitp"))

  aits = copy(aitp)
  dfit = hcat(dfit, matrix2DF(aits, "aits"))


  aijtx = zeros(Float64, (n, m, T))
  max_aijtx = 0.0

  for i = 1:n
    for j = 1:m
      for t = 1:T
        aijtx[i, j, t] = slim_random((aits[i, t] * 100), 100) / 100
        if max_aijtx < aijtx[i, j, t]
          max_aijtx = aijtx[i, j, t]
        end
      end
    end
  end
  dfijx = map(t -> matrix2DF(aijtx[:, :, t], "aijtx"), 1:T)

  aitf = zeros(Float64, n, T)
  for i = 1:n
    for t = 1:T
      if t >= 2
        aitf[i, t] =
          slim_random(max(aitf[i, t-1] - aits[i, t-1] - aith[i, t], max_aijtx * 100), 100) / 100
      else
        aitf[i, t] = slim_random(max_aijtx * 100, 100) / 100
      end
    end

  end
  dfit = hcat(dfit, matrix2DF(aitf, "aitf"))
  ajts = zeros(Float64, m, T)
  for j = 1:m
    for t = 1:T
      ajts[j, t] = slim_random(30, 80) / 100
    end
  end
  dfjt = hcat(dfjt, matrix2DF(ajts, "ajts"))
  ajtp = copy(ajts)

  dfjt = hcat(dfjt, matrix2DF(ajtp, "ajtp"))

  aijty = zeros(Float64, (n, m, T))
  for i = 1:n
    for j = 1:m
      for t = 1:T
        aijty[i, j, t] = 1 / aijtx[i, j, t]
      end
    end
  end
  dfijy = map(t -> matrix2DF(aijty[:, :, t], "aijty"), 1:T)

  D = zeros(n, T)
  for i = 1:n
    for t = 1:T
      D[i, t] = slim_random(3, 7)
    end
  end

  Ki = zeros(n, T)
  for i = 1:n
    for t = 1:T
      Ki[i, t] = slim_random(1, 3)
    end
  end

  Ci = zeros(n, T)
  for i = 1:n
    for t = 1:T
      Ci[i, t] = slim_random(1, 4)
    end
  end

  Li = zeros(n, T)
  for i = 1:n
    for t = 1:T
      Li[i, t] = slim_random(2, 5)
    end
  end
  dfit = hcat(
    dfit,
    matrix2DF(D, "Dit"),
    matrix2DF(Ki, "Kit"),
    matrix2DF(Ci, "Cit"),
    matrix2DF(Li, "Lit"),
  )

  Kj = zeros(m, T)
  for j = 1:m
    for t = 1:T
      Kj[j, t] = trunc(Int, slim_random(0, 13))
    end
  end

  Cj = zeros(m, T)
  for j = 1:m
    for t = 1:T
      Cj[j, t] = slim_random(7, 13)
    end
  end

  Lj = zeros(m, T)
  for j = 1:m
    for t = 1:T
      Lj[j, t] = slim_random(8, 14)
    end
  end
  dfjt = hcat(dfjt, matrix2DF(Kj, "Kjt"), matrix2DF(Cj, "Cjt"), matrix2DF(Lj, "Ljt"))
  return Dict(:it => dfit, :jt => dfjt, :ijx => dfijx, :ijy => dfijy)
end

function getSlimData(customers::Int = 3, generators::Int = 1, periods::Int = 2)
  dfs = slim_data(customers, generators, periods)
  aith = select(dfs[:it], ["aith_$i" for i = 1:periods]) |> Matrix
  aitp = select(dfs[:it], ["aitp_$i" for i = 1:periods]) |> Matrix
  aits = select(dfs[:it], ["aits_$i" for i = 1:periods]) |> Matrix
  aitf = select(dfs[:it], ["aitf_$i" for i = 1:periods]) |> Matrix
  ajtp = select(dfs[:jt], ["ajtp_$i" for i = 1:periods]) |> Matrix
  ajts = select(dfs[:jt], ["ajts_$i" for i = 1:periods]) |> Matrix
  aijtx = Array{Float64,3}(undef, customers, generators, periods)
  dfx = dfs[:ijx]
  for t = 1:periods
    aijtx[:, :, t] = select(dfx[t], ["aijtx_$i" for i = 1:generators]) |> Matrix
  end
  aijty = Array{Float64,3}(undef, customers, generators, periods)
  dfy = dfs[:ijy]
  for t = 1:periods
    aijty[:, :, t] = select(dfy[t], ["aijty_$i" for i = 1:generators]) |> Matrix
  end
  Dit = select(dfs[:it], ["Dit_$i" for i = 1:periods]) |> Matrix
  Lit = select(dfs[:it], ["Lit_$i" for i = 1:periods]) |> Matrix
  Kit = select(dfs[:it], ["Kit_$i" for i = 1:periods]) |> Matrix
  Cit = select(dfs[:it], ["Cit_$i" for i = 1:periods]) |> Matrix
  Ljt = select(dfs[:jt], ["Ljt_$i" for i = 1:periods]) |> Matrix
  Kjt = select(dfs[:jt], ["Kjt_$i" for i = 1:periods]) |> Matrix
  Cjt = select(dfs[:jt], ["Cjt_$i" for i = 1:periods]) |> Matrix
  # if customers == 3 && generators == 1 && periods == 2
  #   aith = [0.3 0.58; 0.39 0.7; 0.59 0.53]
  #   aitp = [0.54 0.95; 0.89 0.92; 0.65 0.93]
  #   aits = [0.54 0.95; 0.89 0.92; 0.65 0.93]
  #   aitf = [0.97 0.97; 0.99 0.98; 0.97 0.97]

  #   ajts = [0.3 0.48]
  #   ajtp = [0.3 0.48]

  #   aijtx = Array{Float64,3}(undef, 3, 1, 2)
  #   aijtx[:, 1, :] = [0.86 0.97; 0.91 0.92; 0.67 0.95]
  #   aijty = Array{Float64,3}(undef, 3, 1, 2)
  #   aijty[:, 1, :] = [1.163 1.031; 1.099 1.087; 1.493 1.053]

  #   Dit = [5 5; 5 5; 3 5]
  #   Kit = [1 1; 1 2; 2 2]
  #   Cit = [2 1; 3 3; 3 3]
  #   Lit = [3 2; 3 2; 4 2]
  #   Kjt = [10 10]
  #   Cjt = [12 12]
  #   Ljt = [11 10]

  #   return SlimData(
  #     customers,
  #     generators,
  #     periods,
  #     aith,
  #     aitp,
  #     aits,
  #     aitf,
  #     ajtp,
  #     ajts,
  #     aijtx,
  #     aijty,
  #     Dit,
  #     Lit,
  #     Kit,
  #     Cit,
  #     Ljt,
  #     Kjt,
  #     Cjt,
  #   )
  # else
  #   # 0< aith <= aitp =aits <= aijtx <= aitf
  #   # aijtx <= aijty <= aitf
  #   d = Normal()
  #   d37 = truncated(d, 0.299, 0.709)
  #   d35 = truncated(d, 0.299, 0.509)
  #   d1 = truncated(d, 0.15, 0.5)
  #   # d5 = truncated(d, 0.5, 0.95)
  #   # d9 = truncated(d, 0.9, 1)
  #   aith = rand(d37, customers, periods)
  #   aitp = aith .+ (d.σ / 2)
  #   aits = aitp
  #   aitf = aits .+ (2 * d.σ)

  #   ajts = rand(d35, generators, periods)
  #   ajtp = ajts

  #   aijty = Array{Float64,3}(undef, customers, generators, periods)
  #   foreach(1:periods) do t
  #     tmp = aitp[:, t] .+ rand(d35, customers, generators)
  #     aijty[:, :, t] = min.(tmp, aitf[:, t])
  #   end

  #   aijtx = Array{Float64,3}(undef, customers, generators, periods)
  #   foreach(1:periods) do t
  #     tmp = aitp[:, t] .+ rand(d1, customers, generators)
  #     aijtx[:, :, t] = min.(tmp, aijty[:, :, t])
  #   end


  #   Dit = rand(3:10, customers, periods)
  #   Kit = rand(1:5, customers, periods)
  #   Cit = rand(1:5, customers, periods)
  #   Lit = rand(1:5, customers, periods)
  #   Kjt = rand(10:30, generators, periods)
  #   Cjt = rand(10:30, generators, periods)
  #   Ljt = rand(10:30, generators, periods)

  #   return SlimData(
  #     customers,
  #     generators,
  #     periods,
  #     aith,
  #     aitp,
  #     aits,
  #     aitf,
  #     ajtp,
  #     ajts,
  #     aijtx,
  #     aijty,
  #     Dit,
  #     Lit,
  #     Kit,
  #     Cit,
  #     Ljt,
  #     Kjt,
  #     Cjt,
  #   )
  # end
  return SlimData(
    customers,
    generators,
    periods,
    aith,
    aitp,
    aits,
    aitf,
    ajtp,
    ajts,
    aijtx,
    aijty,
    Dit,
    Lit,
    Kit,
    Cit,
    Ljt,
    Kjt,
    Cjt,
  )
end

# deed data


function getBudgetLimit(generators::Int)
  U = if generators == 6
    50_000
  elseif generators == 10
    100_000
  else
    10_000 * generators
  end
  U
end
function getCustomersFunctionData(customers::Int)
  # K1j, K2j θj CMj
  K5 = [
    1.847 11.64 0 200
    1.378 11.63 0.1734 280
    1.079 11.32 0.4828 410
    0.9124 11.5 0.7208 500
    0.8794 11.21 1 700
  ]
  K7 = [
    1.847 11.64 0 180
    1.378 11.63 0.14 230
    1.079 11.32 0.26 310
    0.9124 11.5 0.37 390
    0.8794 11.21 0.55 440
    1.378 11.63 0.84 530
    1.5231 11.5 1 600
  ]
  K = if customers == 5
    K5
  elseif customers == 7
    K7
  else
    nd = map(enumerate(zip(mean.(eachcol(K7)), std.(eachcol(K7))))) do (i, (m, s))
      i == 3 ? truncated(Normal(), 0, 1) : Normal(m, s)
    end
    tmp = Matrix{Float64}(undef, customers, 4)
    foreach(1:4) do i
      r = rand(nd[i], customers)
      tmp[:, i] = i == 3 ? vcat(0, sort(r)[2:end-1], 1) : r
    end
    tmp
  end
  K

end
function getMatrixB(generators::Int)
  B1 =
    1e-4 * [
      0.420 0.051 0.045 0.057 0.078 0.066
      0.051 0.180 0.039 0.048 0.045 0.060
      0.045 0.039 0.195 0.051 0.072 0.057
      0.057 0.048 0.051 0.213 0.090 0.075
      0.078 0.045 0.072 0.090 0.207 0.096
      0.066 0.060 0.057 0.075 0.096 0.255
    ]
  B2 =
    1e-5 * [
      4.9 1.4 1.5 1.5 1.6 1.7 1.7 1.8 1.9 2.0
      1.4 4.5 1.6 1.6 1.7 1.5 1.5 1.6 1.8 1.8
      1.5 1.6 3.9 1.0 1.2 1.2 1.4 1.4 1.6 1.6
      1.5 1.6 1.0 4.0 1.4 1.0 1.1 1.2 1.4 1.5
      1.6 1.7 1.2 1.4 3.5 1.1 1.3 1.3 1.5 1.6
      1.7 1.5 1.2 1.0 1.1 3.6 1.3 1.2 1.4 1.5
      1.7 1.5 1.4 1.1 1.3 1.2 3.8 1.6 1.6 1.8
      1.8 1.6 1.4 1.2 1.3 1.2 1.6 4.0 1.5 1.6
      1.9 1.8 1.6 1.4 1.5 1.4 1.6 1.5 4.2 1.9
      2.0 1.8 1.6 1.5 1.6 1.5 1.8 1.6 1.9 4.4
    ]

  B = if generators == 6
    B1
  elseif generators == 10
    B2
  else
    d = truncated(Normal(0, 1), minimum(B2), maximum(B1))
    M = rand(d, generators, generators)
    M + M' + diagm(rand(d, generators))
  end
  B
end

function getLambda(customers::Int, periods::Int = 24)
  # $/MW 
  # source Kim, H.-J., & Kim, M.-K. (2019). Multi-Objective Based Optimal Energy Management of Grid-Connected Microgrid Considering Advanced Demand Response. Energies, 12(21), Article 21. https://doi.org/10.3390/en12214142

  l = [
    27.61 29.41 28.24 26.69 29.01 33.96 83.97 81.10 110.60 74.12 78.95 66.85 47.98 66.82 48.50 49.21 66.65 61.49 56.19 57.92 49.16 54.00 34.37 30.30
    28.30 30.07 28.87 28.76 32.24 36.67 89.46 82.88 112.93 75.43 80.19 67.55 48.58 67.74 49.35 50.28 69.36 66.57 57.67 59.38 49.86 54.38 34.67 30.71
    28.79 30.53 29.28 29.28 32.64 37.15 90.65 83.79 114.11 76.09 80.65 67.76 48.63 68.07 49.69 50.87 70.29 67.19 58.25 59.98 50.36 54.84 34.96 31.00
    26.93 28.79 27.66 27.66 31.20 35.38 85.71 79.06 107.72 72.40 77.29 65.75 47.10 65.55 47.41 49.94 66.05 59.69 54.48 55.58 48.31 53.46 33.98 29.89
    27.60 29.44 28.33 28.32 31.66 35.99 87.70 81.06 110.44 73.95 78.93 66.67 47.93 66.74 48.47 49.19 67.71 66.24 56.53 57.98 48.96 53.63 34.21 30.20
  ]
  d = if customers == 5
    l
  else
    d = truncated(Normal(0, 60), 26, 115)
    rand(d, customers, periods)
  end
  d
end
function getProvidersPowerLimits(generators::Int)
  pl_6u = [
    100 500 120 80
    50 200 90 50
    80 300 100 65
    50 150 90 50
    50 200 90 50
    50 120 90 50
  ]
  pl_10u = [
    150 470 80 80
    135 460 80 80
    73 340 80 80
    60 300 50 50
    73 243 50 50
    57 160 50 50
    20 130 30 30
    47 120 30 30
    20 80 30 30
    55 55 30 30
  ]
  pl = if generators == 6
    pl_6u
  elseif generators == 10
    pl_10u
  else
    tmp = Matrix{Float64}(undef, generators, 4)
    foreach(1:4) do i
      tmp[:, i] = rand(vcat(pl_6u[:, i], pl_10u[:, i]), generators)
    end
    foreach(1:generators) do i
      tmp[:, 1:2] = sort(tmp[:, 1:2], dims = 1)
    end
    tmp
  end
  pl
end
function getCostCoefficients(generators::Int)
  # a +bP +cP^2 $/h $/MWh $/MW^2h
  # e +fP +fP^2 lb/h lb/MWh lb/MW^2h

  cost_6units = [
    240 7 0.007 13.8593 0.32767 0.00419
    200 10 0.0095 13.8593 0.32767 0.00419
    220 8.5 0.009 40.2669 -0.54551 0.00683
    200 11 0.009 40.2669 -0.54551 0.00683
    220 10.5 0.008 42.8955 -0.51116 0.00461
    190 12 0.0075 42.8955 -0.51116 0.00461
  ]
  cost_10units = [
    958.2 21.6 0.00043 360.0012 -3.9864 0.04702
    1313.6 21.05 0.00063 350.0056 -3.9524 0.04652
    604.97 20.81 0.00039 330.0056 -3.9023 0.04652
    471.6 23.9 0.0007 330.0056 -3.9023 0.04652
    480.29 21.62 0.00079 13.8593 0.3277 0.0042
    601.75 17.87 0.00056 13.8593 0.3277 0.0042
    502.7 16.51 0.00211 40.2669 -0.5455 0.0068
    639.4 23.23 0.0048 40.2669 -0.5455 0.0068
    455.6 19.58 0.10908 42.8955 -0.5112 0.0046
    692.4 22.54 0.00951 42.8955 -0.5112 0.0046
  ]
  cs = if (generators == 6)
    cost_6units[:, 1:6]
  elseif generators == 10
    cost_10units[:, 1:6]
  else
    c6umin = minimum(cost_6units[:, 1:6], dims = 1) ./ 6
    c6umax = maximum(cost_6units[:, 1:6], dims = 1) ./ 6
    c10umin = minimum(cost_10units[:, 1:6], dims = 1) ./ 10
    c10umax = maximum(cost_10units[:, 1:6], dims = 1) ./ 10
    cmin = 0.5 * (c6umin + c10umin) * generators
    cmax = 0.5 * (c6umax + c10umax) * generators
    nd = Normal.(cmax, 0.5 * (cmax - cmin))
    C = Matrix{Float64}(undef, generators, 6)
    foreach(1:6) do i
      C[:, i] = rand(nd[i], 1, generators)
    end
    C
  end

  cs
end
function getCustomersDemand(initialDt::Vector{Float64}, customers::Int, periods::Int)
  if periods > 24
    # @error "The data will be provided for one day only"
    throw("The data will be provided for one day only")
  end
  d = fit(Normal, initialDt)
  rand(d, customers, periods)
end
function getDemand(generators::Int, periods::Int = 24)
  # MW
  if periods > 24
    # @error "The data will be provided for one day only"
    throw("The data will be provided for one day only")
  end
  Dt = if generators == 6
    # 6 unit - source Elaiw, A. M., Xia, X., & Shehata, A. M. (2012). Application of model predictive control to optimal dynamic dispatch of generation with emission limitations. Electric Power Systems Research, 84(1), 31–44. https://doi.org/10.1016/j.epsr.2011.09.024
    # 10 unit source Xia, X., Zhang, J., & Elaiw, A. (2011). An application of model predictive control to the dynamic economic dispatch of power generation. Control Engineering Practice, 19(6), 638–648. https://doi.org/10.1016/j.conengprac.2011.03.001

    [
      955,
      942,
      935,
      930,
      935,
      963,
      989,
      1023,
      1126,
      1150,
      1201,
      1235,
      1190,
      1251,
      1263,
      1250,
      1221,
      1202,
      1159,
      1092,
      1023,
      984,
      975,
      960,
    ]
  elseif generators == 10
    [
      1036,
      1110,
      1258,
      1406,
      1480,
      1628,
      1702,
      1776,
      1924,
      2072,
      2146,
      2220,
      2072,
      1924,
      1776,
      1554,
      1480,
      1628,
      1776,
      2072,
      1924,
      1628,
      1332,
      1184,
    ]
  else
    mn, mx = 985, 1750
    d = Normal(mx, 0.5 * (mx - mn))
    midday = rand(13:16)
    round.(Int, vcat(sort(rand(d, midday)), sort(rand(d, 24 - midday), rev = true)))
  end
  D = if periods == 24
    Dt
  else
    interval_length = floor(24 / periods) |> Int
    vcat([sum(Dt[i:min(24, i * interval_length)]) for i = 1:periods]...) |> shuffle

  end
  D
end
function getGTDRDeedData(customers::Int, generators::Int, periods::Int)
  # ahdot_it::Matrix{Float64}
  # adot_it::Matrix{Float64}
  # bdot_it::Matrix{Float64}
  # cdot_it::Matrix{Float64}
  # ddot_it::Matrix{Float64}
  # edot_it::Matrix{Float64}
  # nudot_it::Matrix{Float64}
  # afdot_it::Matrix{Float64}
  deed_data = getDRDeedData(customers, generators, periods)
  # 0< aith <= aitp =aits <= aijtx <= aitf
  means_of_generators_coef =
    repeat(mean.([deed_data.a deed_data.b deed_data.c]) ./ customers, customers)
  dists = fit.(Normal, eachcol(means_of_generators_coef))
  adot_it, bdot_it, cdot_it = rand(dists[1], customers, periods),
  rand(dists[2], customers, periods),
  rand(dists[3], customers, periods)

  ddot_it, edot_it = bdot_it, cdot_it
  ahdot_it = adot_it .- 0.25 * dists[1].σ
  nudot_it = adot_it .+ 0.25 * dists[1].σ
  afdot_it = nudot_it .+ 0.125 * dists[1].σ
  # 1-3, 4-5, x: 6, f:7, h:8
  pertCustomerDt = deed_data.Dt ./ (2 * customers)
  custpmersDt = getCustomersDemand(vec(collect(pertCustomerDt')), customers, periods)
  pmin = rand(2:5, customers)
  pmax = pmin + rand(1:3, customers)
  return GTDeedData(
    deed_data,
    ahdot_it,
    adot_it,
    bdot_it,
    cdot_it,
    ddot_it,
    edot_it,
    nudot_it,
    afdot_it,
    custpmersDt,
    pmin,
    pmax,
  )
end
# function getGTDRData(customers::Int, generators::Int, periods::Int)
#   household_consumption = [
#     10 8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10
#     8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8
#     7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7
#     6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7 6 9 12 15 14 10 8 7 6
#   ]

#   # Define the coefficients matrix
#   # Columns: $\dot{a_i}$, $\dot{b_i}$, $\dot{c_i}$, $\dot{d_i}$, $\dot{e_i}$, $\dot{\nu_1}$, $\dot{\nu_2}$
#   # Rows: Households 1-4 (Morning, Evening)
#   coefficients_matrix = [
#     100 2 0.1 5 1 3 2
#     120 2.5 0.15 6 1.2 3.5 2.2
#     90 1.8 0.08 4 0.8 2.5 1.8
#     110 2.2 0.12 5.5 1.1 3.2 2.1
#     100 2 0.1 5 1 3 2
#     120 2.5 0.15 6 1.2 3.5 2.2
#     90 1.8 0.08 4 0.8 2.5 1.8
#     110 2.2 0.12 5.5 1.1 3.2 2.1
#   ]
#   ad = coefficients_matrix[1:n, 1]
#   bd = coefficients_matrix[1:n, 2]
#   cd = coefficients_matrix[1:n, 3]
#   dd = coefficients_matrix[1:n, 4]
#   ed = coefficients_matrix[1:n, 5]
#   νd = coefficients_matrix[1, 6:7]

#   # Define the coefficients matrix
#   # Columns: $\dot{\alpha_i}$, $\dot{\beta_i}$, $\dot{\gamma_i}$, $\dot{\delta_i}$
#   # Rows: Households 1-4
#   emission_coefficients_matrix = [
#     50 1 0.05 2
#     60 1.2 0.06 2.5
#     45 0.9 0.04 1.8
#     55 1.1 0.055 2.2
#   ]
#   αd = emission_coefficients_matrix[:, 1]
#   βd = emission_coefficients_matrix[:, 2]
#   γd = emission_coefficients_matrix[:, 3]
#   δd = emission_coefficients_matrix[:, 4]
#   # ddd, edd, νdd
#   power_plant_coefficients_matrix = [
#     5 1 3 2 4 1.5
#     5 1.1 3.5 2.5 4.2 1.8
#   ]
#   ddd = power_plant_coefficients_matrix[:, 1]
#   edd = power_plant_coefficients_matrix[:, 2]
#   νdd = power_plant_coefficients_matrix[1, 3:end]
#   δdd = [
#     0.15
#     0.12
#   ]
#   demand_data = [
#     955,
#     942,
#     935,
#     930,
#     935,
#     963,
#     989,
#     1023,
#     1126,
#     1150,
#     1201,
#     1235,
#     1190,
#     1251,
#     1263,
#     1250,
#     1221,
#     1202,
#     1159,
#     1092,
#     1023,
#     984,
#     975,
#     960,
#   ]
#   lambda1 = [
#     1 1.70 3.70 2.70
#     2 1.40 2.70 1.90
#     3 2.20 3.20 1.80
#     4 3.70 2.60 1.90
#     5 4.50 3.80 2.30
#     6 4.70 1.70 0.70
#     7 5.10 2.30 1.40
#     8 5.30 1.50 0.50
#     9 6.70 4.30 2.90
#     10 6.60 4.60 1.60
#     11 6.80 3.50 4.30
#     12 6.20 4.20 4.80
#     13 7.30 4.30 5.10
#     14 7.80 6.30 5.40
#     15 0.50 3.50 5.50
#     16 5.20 5.30 6.10
#     17 6.80 5.30 5.60
#     18 5.70 6.10 6.30
#     19 4.80 2.60 4.50
#     20 3.90 3.60 4.20
#     21 3.80 4.20 3.90
#     22 3.10 3.80 3.20
#     23 2.50 2.30 2.80
#     24 1.90 3.80 4.20
#   ]
#   lambda2 = [
#     1 2.10 3.20 2.50
#     2 1.80 2.50 2.10
#     3 2.50 3.00 2.30
#     4 3.20 2.80 1.60
#     5 4.10 3.40 2.70
#     6 4.40 1.80 1.10
#     7 5.20 2.60 1.80
#     8 5.50 1.90 1.20
#     9 6.90 4.10 3.10
#     10 6.20 4.80 2.30
#     11 6.50 3.80 4.50
#     12 5.90 4.50 5.00
#     13 7.10 4.80 5.30
#     14 7.50 6.60 5.70
#     15 0.80 3.80 5.80
#     16 5.80 5.50 6.20
#     17 6.50 5.50 5.70
#     18 5.40 6.30 6.40
#     19 4.50 2.80 4.70
#     20 3.60 3.70 4.40
#     21 3.50 4.40 4.10
#     22 3.00 4.00 3.40
#     23 2.40 2.80 3.30
#     24 1.80 4.00 4.40
#   ]
#   λ24 = hcat(lambda1[:, 2:end], lambda2[:, 2:end])
#   λ5 = [
#     27 30 27 26 30 37 88 81 110 75 81 70 50 70 50 51 69 62 55 57 51 52 32 27
#     28 31 28 27 31 38 89 82 111 76 82 71 51 71 51 52 70 63 56 58 52 53 33 28
#     29 32 29 28 32 39 90 83 112 77 83 72 52 72 52 53 71 64 57 59 53 54 34 29
#     30 33 30 29 33 40 91 84 113 78 84 73 53 73 53 54 72 65 58 60 54 55 35 30
#     31 34 31 30 34 41 92 85 114 79 85 74 54 74 54 55 73 66 59 61 55 56 36 31
#   ]
#   λT = Matrix([mean(λ24[1:12, :], dims = 1); mean(λ24[13:end, :], dims = 1)]')
#   yc = 21.540
#   Dit = 12 * 0.0025
#   Sdi = [5; 4; 5; 4]
#   Sddj = [24; 24]
#   Ldi = [3 2; 3 2; 3 2; 4 2]
#   Lddj = [11 10; 11 10]


# end
function getDRDeedData(customers::Int, generators::Int, periods::Int)
  generators_cost_coof = getCostCoefficients(generators)
  generators_power_limits = getProvidersPowerLimits(generators)
  customers_function_data = getCustomersFunctionData(customers)
  return DeedData(
    customers,
    generators,
    periods,
    getDemand(generators, periods),
    getMatrixB(generators),
    generators_power_limits[1:generators, 1],
    generators_power_limits[1:generators, 2],
    generators_power_limits[1:generators, 3],
    generators_power_limits[1:generators, 4],
    customers_function_data[:, 1],
    customers_function_data[:, 2],
    customers_function_data[:, 3],
    customers_function_data[:, 4],
    getBudgetLimit(generators),
    generators_cost_coof[1:generators, 1],
    generators_cost_coof[1:generators, 2],
    generators_cost_coof[1:generators, 3],
    generators_cost_coof[1:generators, 4],
    generators_cost_coof[1:generators, 5],
    generators_cost_coof[1:generators, 6],
    getLambda(customers, periods),
  )
end

