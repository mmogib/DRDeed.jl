@variable(model, q[1:generators, 1:periods] >= 0)
@variable(model, p[1:customers, 1:periods] >= 0)
@variable(model, s[1:customers, 0:periods] >= 0)
@variable(model, χ[1:customers, 1:periods] >= 0)
@variable(model, ω[1:customers, 1:periods] >= 0)
@variable(model, x[1:customers, 1:periods] >= 0)
@variable(model, f[1:customers, 1:periods] >= 0)
@variable(model, h[1:customers, 1:periods+1] >= 0)
@variable(model, y[1:customers, 1:generators, 1:periods] >= 0)

@NLexpression(
  model,
  Ci[i in 1:customers, t in 1:periods],
  gtdata.adot_it[i, t] +
  gtdata.ahdot_it[i, t] * h[i, t] +
  gtdata.bdot_it[i, t] * p[i, t] +
  gtdata.cdot_it[i, t] * p[i, t]^2 +
  gtdata.edot_it[i, t] * s[i, t] +
  gtdata.ddot_it[i, t] * s[i, t]^2 +
  gtdata.nudot_it[i, t] * x[i, t] +
  gtdata.afdot_it[i, t] * f[i, t]
)
@NLexpression(
  model,
  Cj[j in 1:generators, t in 1:periods],
  data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2
)
# +
# sum(gtdata.nudot_it[i, t] * x[i, t] for i = 1:customers)
@NLexpression(
  model,
  Ej[j in 1:generators, t in 1:periods],
  data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2
)

@NLexpression(model, Ccust, sum(Ci[i, t] for i = 1:customers for t = 1:periods))
@NLexpression(model, C, sum(Cj[j, t] for j = 1:generators for t = 1:periods))
@NLexpression(model, E, sum(Ej[j, t] for j = 1:generators for t = 1:periods))

@NLexpression(
  model,
  utility,
  sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods) + Ccust
)

@NLexpression(
  model,
  losst[t in 1:periods],
  sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
)
# Balance <-- 
@constraint(
  model,
  powerbalance1[i in 1:customers, t in 1:periods],
  f[i, t] + s[i, t] - s[i, t-1] - p[i, t] + x[i, t] == sum(y[i, j, t] for j = 1:generators)
)
@constraint(
  model,
  powerbalance2[i in 1:customers, t in 1:periods],
  gtdata.CDemandt[i, t] + h[i, t] == h[i, t+1] + f[i, t]
)
@constraint(
  model,
  powerbalance3[j in 1:generators, t in 1:periods],
  sum(y[i, j, t] for i = 1:customers) == q[j, t]
)
@NLconstraint(
  model,
  powerbalance4[t in 1:periods],
  sum(x[i, t] for i = 1:customers) + sum(q[j, t] for j = 1:generators) ==
  sum(gtdata.CDemandt[i, t] for i = 1:customers) + losst[t] - sum(χ[i, t] for i = 1:customers)
)
## <-- Balance

# limits -->
@constraint(model, [t in 1:periods-1], gtdata.pimin .<= s[:, t] .<= gtdata.pimax)
@constraint(model, [t in 1:periods], gtdata.pimin .<= p[:, t] .<= gtdata.pimax)
@constraint(model, [t in 1:periods], sum(x[:, t]) .<= gtdata.pimax)
@constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
@constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
# <-- limits

@NLconstraint(
  model,
  benefit[i in 1:customers, t in 1:periods],
  ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >= 0
)
@NLconstraint(
  model,
  benefit2[i in 2:customers, t in 1:periods],
  ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >=
  ω[i-1, t] - χ[i-1, t] * (data.K1[i-1] * χ[i-1, t] + (1 - data.θ[i-1]) * data.K2[i-1])
)

@constraint(model, budget, sum(ω[i, t] for i = 1:customers for t = 1:periods) <= data.UB)
@constraint(model, load[i in 1:customers], sum(χ[i, t] for t = 1:periods) <= data.CM[i])
@constraint(model, storage1[i in 1:customers], s[i, 0] == 0)
@constraint(model, storage2[i in 1:customers], s[i, periods] == 0.0)
@constraint(model, shifted_load[i in 1:customers], h[i, 1] == 0)


@constraint(model, balance_chi[i in 1:customers, t in 1:periods], χ[i, t] <= h[i, t] + s[i, t])
@NLobjective(model, Min, w[1] * C + w[2] * E - w[3] * utility)
optimize!(model)