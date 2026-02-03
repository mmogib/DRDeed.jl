function gtdrdeed(
  customers::Int = 5,
  generators::Int = 6,
  periods::Int = 24;
  solver::Union{Nothing,MOI.AbstractOptimizer} = nothing,
  data::Union{Nothing,GTDeedData} = nothing,
)
  gtdata = isnothing(data) ? getGTDRDeedData(customers, generators, periods) : data
  solver =
    isnothing(solver) ?
    optimizer_with_attributes(
      Ipopt.Optimizer,
      MOI.Silent() => true,
      "sb" => "yes",
      "max_iter" => 10_000,
    ) : solver



  data = gtdata.deedData
  function fun(x::Vector{<:Number})
    """
    @variable(model, q[1:generators, 1:periods] >= 0)
    @variable(model, p[1:customers, 1:periods] >= 0)
    @variable(model, s[1:customers, 0:periods] >= 0)
    @variable(model, χ[1:customers, 1:periods] >= 0)
    @variable(model, ω[1:customers, 1:periods] >= 0)
    @variable(model, x[1:customers, 1:periods] >= 0)
    @variable(model, f[1:customers, 1:periods] >= 0)
    @variable(model, h[1:customers, 1:periods+1] >= 0)
    @variable(model, y[1:customers, 1:generators, 1:periods] >= 0)
    """
    GtimeP = generators * periods
    CtimesP = customers * periods
    CtimesP1 = customers * (periods + 1)

    # CtimesGtomesP = customers * GtimeP
    q = reshape(x[1:GtimeP], generators, periods)
    p = reshape(x[GtimeP+1:GtimeP+CtimesP], customers, periods)
    s = reshape(x[GtimeP+CtimesP+1:GtimeP+CtimesP+CtimesP1], customers, periods + 1)
    χ = reshape(x[GtimeP+CtimesP+CtimesP1+1:GtimeP+2CtimesP+CtimesP1], customers, periods)
    ω = reshape(x[GtimeP+2CtimesP+CtimesP1+1:GtimeP+3CtimesP+CtimesP1], customers, periods)
    xx = reshape(x[GtimeP+3CtimesP+CtimesP1+1:GtimeP+4CtimesP+CtimesP1], customers, periods)
    f = reshape(x[GtimeP+4CtimesP+CtimesP1+1:GtimeP+5CtimesP+CtimesP1], customers, periods)
    h = reshape(x[GtimeP+5CtimesP+CtimesP1+1:GtimeP+5CtimesP+2CtimesP1], customers, periods + 1)
    y = reshape(x[GtimeP+5CtimesP+2CtimesP1+1:end], customers, generators, periods)
    Ci = sum(
      gtdata.adot_it[i, t] +
      gtdata.ahdot_it[i, t] * h[i, t] +
      gtdata.bdot_it[i, t] * p[i, t] +
      gtdata.cdot_it[i, t] * p[i, t]^2 +
      gtdata.edot_it[i, t] * s[i, t] +
      gtdata.ddot_it[i, t] * s[i, t]^2 +
      gtdata.nudot_it[i, t] * xx[i, t] +
      gtdata.afdot_it[i, t] * f[i, t] for i = 1:customers for t = 1:periods
    )
    Cj = sum(
      data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2 for j = 1:generators for
      t = 1:periods
    )
    Ej = sum(
      data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2 for j = 1:generators for
      t = 1:periods
    )
    utility = sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods) + Ci
    return vcat(Cj, Ej, -utility)
  end
  function loss_fun(x::Vector{<:Number})
    gtiumeperiod = generators * periods
    q = reshape(x[1:gtiumeperiod], generators, periods)
    sum(
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators) for
      t = 1:periods
    )
  end
  function ws(w::Vector{Float64} = ones(3))
    model = Model(solver)
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

    if (has_values(model))
      deedSol = DRDeedSolution(
        value.(χ),
        value.(ω),
        value.(q),
        value(C),
        value(E),
        value(utility),
        value.(losst),
        w,
      )
      solution = GTDRDeedSolution(
        deedSol,
        value.(p),
        value.(s),
        value.(h),
        value.(f),
        value.(x),
        value.(y),
        value.(Ci),
      )
      return SuccessResult(model, gtdata, solution)
    else
      return FailResult(model, gtdata, termination_status(model))
    end
  end
  function ps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, z)
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
    @NLconstraint(model, (v[1] + z * d[1] - C) >= 0)
    @NLconstraint(model, (v[2] + z * d[2] - E) >= 0)
    @NLconstraint(model, (v[3] + z * d[3] + utility) >= 0)
    @objective(model, Min, z)
    optimize!(model)

    if (has_values(model))
      # x = value.(vcat(q[:], p[:], s[:]))
      sv = value.(s).data
      x = value.(vcat(q[:], p[:], sv[:], χ[:], ω[:], x[:], f[:], h[:], y[:]))
      return x, value.(z)
    else
      return nothing
    end

  end
  function dps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
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
    @variable(model, t)
    @variable(model, w[1:3] >= 0)
    @constraint(model, dot(w, d) == 1)
    @NLconstraint(model, C * w[1] + E * w[2] - w[3] * utility >= t)
    @objective(model, Max, t - dot(w, v))
    optimize!(model)
    if (has_values(model))
      return value.(w)
    else
      return nothing
    end
  end

  return Dict(:ws => ws, :ps => ps, :dps => dps, :f => fun, :loss => loss_fun)
end

export gtdrdeed
"""
function solvefeasible(; w::Vector{Float64} = ones(3))
    @NLobjective(model, Max, 0)
    optimize!(model)

    if (has_values(model))
      deedSol = DRDeedSolution(
        value.(χ),
        value.(ω),
        value.(q),
        value(C),
        value(E),
        value(utility),
        value.(losst),
        w,
      )
      solution = GTDRDeedSolution(
        deedSol,
        value.(p),
        value.(s),
        value.(h),
        value.(f),
        value.(x),
        value.(y),
        value.(Ci),
      )
      return SuccessResult(model, gtdata, solution)
    else
      return FailResult(model, gtdata, termination_status(model))
    end
  end
  function solvehybrid(C0, E0, Util0; w::Vector{Float64} = ones(3))
    @NLconstraint(model, C <= 0.5 * C0)
    @NLconstraint(model, E <= 0.5 * E0)
    @NLconstraint(model, utility >= -Util0)
    @NLobjective(model, Min, w[1] * C + w[2] * E - w[3] * utility)
    optimize!(model)

    if (has_values(model))
      deedSol = DRDeedSolution(
        value.(χ),
        value.(ω),
        value.(q),
        value(C),
        value(E),
        value(utility),
        value.(losst),
        w,
      )
      solution = GTDRDeedSolution(
        deedSol,
        value.(p),
        value.(s),
        value.(h),
        value.(f),
        value.(x),
        value.(y),
        value.(Ci),
      )
      return SuccessResult(model, gtdata, solution)
    else
      return FailResult(model, gtdata, termination_status(model))
    end
  end
  """