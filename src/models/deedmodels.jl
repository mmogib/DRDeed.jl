include("deedmodel/deed.jl")
include("deedmodel/drdeed.jl")
include("deedmodel/gtdrdeed.jl")
include("deedmodel/gamethmodels.jl")

# function solveit(factor = 1)
#   n, m, periods = 4, 2, 2
#   # Define the 4x24 matrix for household electricity consumption
#   # kW/h

#   model = Model(Ipopt.Optimizer)
#   @variable(model, p[1:customers, 1:periods] >= 0)
#   @variable(model, s[1:customers, 0:periods] >= 0)
#   @variable(model, y[1:customers, 1:generators, 1:periods] >= 0)
#   @NLexpression(
#     model,
#     Ci[i in 1:customers, t in 1:periods],
#     ad[i] +
#     bd[i] * p[i, t] +
#     cd[i] * p[i, t]^2 +
#     dd[i] * s[i, t] +
#     ed[i] * (s[i, t])^2 +
#     sum(νd[j] * y[i, j, t] for j = 1:generators)
#   )
#   @NLexpression(
#     model,
#     Ei[i in 1:customers, t in 1:periods],
#     αd[i] + βd[i] * p[i, t] + γd[i] * p[i, t]^2 - δd[i] * s[i, t]
#   )

#   @variable(model, q[1:generators, 1:periods] >= 0)
#   @variable(model, r[1:generators, 0:periods] >= 0)
#   @variable(model, x[1:customers, 1:generators, 1:periods] >= 0)
#   @NLexpression(
#     model,
#     Cj[j in 1:generators, t in 1:periods],
#     add[j] + bdd[j] * q[j, t] + cd[j] * q[j, t]^2 + dd[j] * r[j, t] + ed[j] * (r[j, t])^2 -
#     sum(νdd[j] * x[i, j, t] for i = 1:generators)
#   )
#   @NLexpression(
#     model,
#     Ej[j in 1:generators, t in 1:periods],
#     αdd[j] + βdd[j] * q[j, t] + γdd[j] * q[j, t]^2 - δdd[j] * r[j, t]
#   )

#   @NLexpression(
#     model,
#     C,
#     sum(Ci[i, t] for i = 1:customers for t = 1:periods) + sum(Cj[j, t] for j = 1:generators for t = 1:periods)
#   )
#   @NLexpression(
#     model,
#     E,
#     sum(Ei[i, t] for i = 1:customers for t = 1:periods) + sum(Ej[j, t] for j = 1:generators for t = 1:periods)
#   )

#   @variable(model, χ[1:customers, 1:periods] >= 0)
#   @variable(model, ω[1:customers, 1:periods] >= 0)
#   @NLexpression(model, utility, sum(λperiods[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods))

#   @NLexpression(
#     model,
#     losst[t in 1:periods],
#     sum(p[i, t] * BIJ[i, k] * p[k, t] for i = 1:customers for k = 1:customers) +
#     sum(q[j, t] * BIJ[n+j-1, n+k-1] * q[k, t] for j = 1:generators for k = 1:generators) +
#     2 * sum(p[i, t] * BIJ[i, j+n-1] * q[j, t] for i = 1:customers for j = 1:generators)
#   )
#   # @NLconstraint(
#   #   model,
#   #   powerbalance[t in 1:periods],
#   #   sum(p[i, t] for i = 1:customers) + sum(q[j, t] for j = 1:generators) ==
#   #   Dit * n - sum(χ[i, t] for i = 1:customers) + losst[t]
#   # )

#   @NLconstraint(
#     model,
#     powerbalance[t in 1:periods],
#     sum(p[j, t] for j = 1:customers) == Dit * n + sum(χ[i, t] for i = 1:customers) - losst[t],
#   )
#   @constraint(model, [t in 1:periods], qjmin .<= q[:, t] .<= qjmax)
#   @constraint(model, [t in 1:(periods-1)], -DRj .<= (q[:, t+1] - q[:, t]) .<= URj)
#   @NLconstraint(
#     model,
#     benefit[i in 1:customers],
#     sum(ω[i, t] - χ[i, t] * (K1[i] * χ[i, t] + (1 - θ[i]) * K2[i]) for t = 1:periods) >= 0
#   )
#   @NLconstraint(
#     model,
#     benefit2[i in 2:customers],
#     sum(ω[i, t] - χ[i, t] * (K1[i] * χ[i, t] + (1 - θ[i]) * K2[i]) for t = 1:periods) >=
#     sum(ω[i-1, t] - χ[i-1, t] * (K1[i-1] * χ[i-1, t] + (1 - θ[i-1]) * K2[i-1]) for t = 1:periods)
#   )

#   @constraint(model, budget, sum(ω[i, t] for i = 1:customers for t = 1:periods) <= UB)
#   @constraint(model, load[i in 1:customers], sum(χ[i, t] for t = 1:periods) <= CM[i])

#   @constraint(model, storage_limit_i[i in 1:customers], sum(s[i, t] for t = 1:periods) <= Sdi[i])
#   @constraint(model, storage_limit_j[j in 1:generators], sum(r[j, t] for t = 1:periods) <= Sddj[j])

#   @variable(model, h[1:customers, 1:(periods+1)] >= 0)
#   @constraint(model, storage_bvs0_i[i in 1:customers], s[i, 0] == 0)
#   @constraint(model, storage_bvsperiods_i[i in 1:customers], s[i, periods] == 0)
#   @constraint(model, shifted_load_bv_i[i in 1:customers], h[i, 1] == 0)

#   @constraint(model, storage_bvr0_j[j in 1:generators], r[j, 0] == 0)
#   @constraint(model, storage_bvrperiods_j[j in 1:generators], r[j, periods] == 0)

#   @variable(model, f[1:customers, 1:periods] >= 0)
#   @constraint(
#     model,
#     conservation[i in 1:customers, t in 1:periods],
#     f[i, t] + s[i, t] - s[i, t-1] - p[i, t] + sum(x[i, j, t] for j = 1:generators) ==
#     sum(y[i, j, t] for j = 1:generators)
#   )

#   @constraint(model, balance[i in 1:customers, t in 1:periods], h[i, t+1] + f[i, t] == h[i, t] + Dit)

#   @constraint(model, power_sent[i in 1:customers, t in 1:periods], sum(x[i, j, t] for j = 1:generators) <= Ldi[i, t])
#   @constraint(model, power_received[j in 1:generators, t in 1:periods], sum(y[i, j, t] for i = 1:customers) <= Lddj[j, t])

#   @constraint(
#     model,
#     conservation_j[j in 1:generators, t in 1:periods],
#     r[j, t] - r[j, t-1] - q[j, t] + sum(y[i, j, t] for i = 1:customers) == sum(x[i, j, t] for i = 1:customers)
#   )
#   w1 = w2 = w3 = 1
#   @NLobjective(model, Min, w1 * C + w2 * E - w3 * utility)
#   optimize!(model)
#   smry = solution_summary(model)
#   return smry,
#   Dict(zip(["Cost", "Emission", "Utility", "q", "p"], map(x -> value.(x), [C, E, utility, q, p])))
# end