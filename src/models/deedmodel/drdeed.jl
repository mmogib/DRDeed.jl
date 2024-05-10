
function drdeed(
  customers::Int = 5,
  generators::Int = 6,
  periods::Int = 24;
  solver::Union{Nothing,MOI.AbstractOptimizer} = nothing,
  data::Union{Nothing,DeedData} = nothing,
)
  data = isnothing(data) ? getDRDeedData(customers, generators, periods) : data
  solver =
    isnothing(solver) ?
    optimizer_with_attributes(
      Ipopt.Optimizer,
      MOI.Silent() => true,
      "sb" => "yes",
      "max_iter" => 10_000,
    ) : solver
  function drdeed_fun(x::Vector{<:Number})
    gtiumeperiod = generators * periods
    ctimesp = customers * periods
    q = reshape(x[1:gtiumeperiod], generators, periods)
    χ = reshape(x[gtiumeperiod+1:gtiumeperiod+ctimesp], customers, periods)
    ω = reshape(x[gtiumeperiod+ctimesp+1:gtiumeperiod+2*ctimesp], customers, periods)
    C = sum(
      data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2 for j = 1:generators for
      t = 1:periods
    )
    E = sum(
      data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2 for j = 1:generators for
      t = 1:periods
    )
    utility = sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods)
    return vcat(C, E, -utility)
  end
  function loss_fun(x::Vector{<:Number})
    gtiumeperiod = generators * periods
    q = reshape(x[1:gtiumeperiod], generators, periods)
    sum(
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators) for
      t = 1:periods
    )
  end
  function ws(w::Vector{<:Number} = ones(3))

    model = Model(solver)
    @variable(model, q[1:generators, 1:periods] >= 0)
    @variable(model, χ[1:customers, 1:periods] >= 0)
    @variable(model, ω[1:customers, 1:periods] >= 0)

    @NLexpression(
      model,
      Cj[j in 1:generators, t in 1:periods],
      data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2
    )
    @NLexpression(
      model,
      Ej[j in 1:generators, t in 1:periods],
      data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2
    )

    @NLexpression(model, C, sum(Cj[j, t] for j = 1:generators for t = 1:periods))
    @NLexpression(model, E, sum(Ej[j, t] for j = 1:generators for t = 1:periods))

    @NLexpression(
      model,
      utility,
      sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods)
    )

    @NLexpression(
      model,
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t] - sum(χ[i, t] for i = 1:customers)
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    # @NLexpression(
    #   model,
    #   bnft[i in 1:customers],
    #   sum(
    #     ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) for t = 1:periods
    #   )
    # )
    @NLconstraint(
      model,
      benefit[t in 1:periods, i = 1:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >= 0
    )
    @NLconstraint(
      model,
      benefit2[t in 1:periods, i = 2:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >=
      ω[i-1, t] - χ[i-1, t] * (data.K1[i-1] * χ[i-1, t] + (1 - data.θ[i-1]) * data.K2[i-1])
    )

    @constraint(model, budget, sum(ω[i, t] for i = 1:customers for t = 1:periods) <= data.UB)
    @constraint(model, load[i in 1:customers], sum(χ[i, t] for t = 1:periods) <= data.CM[i])
    @NLobjective(model, Min, w[1] * C + w[2] * E - w[3] * utility)
    optimize!(model)

    if (has_values(model))
      solution = DRDeedSolution(
        value.(χ),
        value.(ω),
        value.(q),
        value(C),
        value(E),
        value(utility),
        value.(losst),
        w,
      )
      return SuccessResult(model, data, solution)
    else
      return FailResult(model, data, termination_status(model))
    end
  end
  function ps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, z)
    @variable(model, q[1:generators, 1:periods] >= 0)
    @variable(model, χ[1:customers, 1:periods] >= 0)
    @variable(model, ω[1:customers, 1:periods] >= 0)

    @NLexpression(
      model,
      Cj[j in 1:generators, t in 1:periods],
      data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2
    )
    @NLexpression(
      model,
      Ej[j in 1:generators, t in 1:periods],
      data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2
    )

    @NLexpression(model, C, sum(Cj[j, t] for j = 1:generators for t = 1:periods))
    @NLexpression(model, E, sum(Ej[j, t] for j = 1:generators for t = 1:periods))

    @NLexpression(
      model,
      utility,
      sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods)
    )

    @NLexpression(
      model,
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t] - sum(χ[i, t] for i = 1:customers)
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    # @NLexpression(
    #   model,
    #   bnft[i in 1:customers],
    #   sum(
    #     ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) for t = 1:periods
    #   )
    # )
    @NLconstraint(
      model,
      benefit[t in 1:periods, i = 1:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >= 0
    )
    @NLconstraint(
      model,
      benefit2[t in 1:periods, i = 2:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >=
      ω[i-1, t] - χ[i-1, t] * (data.K1[i-1] * χ[i-1, t] + (1 - data.θ[i-1]) * data.K2[i-1])
    )

    @constraint(model, budget, sum(ω[i, t] for i = 1:customers for t = 1:periods) <= data.UB)
    @constraint(model, load[i in 1:customers], sum(χ[i, t] for t = 1:periods) <= data.CM[i])
    @NLconstraint(model, (v[1] + z * d[1] - C) >= 0)
    @NLconstraint(model, (v[2] + z * d[2] - E) >= 0)
    @NLconstraint(model, (v[3] + z * d[3] + utility) >= 0)
    @objective(model, Min, z)
    optimize!(model)
    x = value.(vcat(q[:], χ[:], ω[:]))
    if (has_values(model))
      return x, value.(z)
    else
      return nothing
    end

  end
  function dps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, q[1:generators, 1:periods] >= 0)
    @variable(model, χ[1:customers, 1:periods] >= 0)
    @variable(model, ω[1:customers, 1:periods] >= 0)

    @NLexpression(
      model,
      Cj[j in 1:generators, t in 1:periods],
      data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2
    )
    @NLexpression(
      model,
      Ej[j in 1:generators, t in 1:periods],
      data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2
    )

    @NLexpression(model, C, sum(Cj[j, t] for j = 1:generators for t = 1:periods))
    @NLexpression(model, E, sum(Ej[j, t] for j = 1:generators for t = 1:periods))

    @NLexpression(
      model,
      utility,
      sum(data.λ[i, t] * χ[i, t] - ω[i, t] for i = 1:customers for t = 1:periods)
    )

    @NLexpression(
      model,
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t] - sum(χ[i, t] for i = 1:customers)
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    # @NLexpression(
    #   model,
    #   bnft[i in 1:customers],
    #   sum(
    #     ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) for t = 1:periods
    #   )
    # )
    @NLconstraint(
      model,
      benefit[t in 1:periods, i = 1:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >= 0
    )
    @NLconstraint(
      model,
      benefit2[t in 1:periods, i = 2:customers],
      ω[i, t] - χ[i, t] * (data.K1[i] * χ[i, t] + (1 - data.θ[i]) * data.K2[i]) >=
      ω[i-1, t] - χ[i-1, t] * (data.K1[i-1] * χ[i-1, t] + (1 - data.θ[i-1]) * data.K2[i-1])
    )

    @constraint(model, budget, sum(ω[i, t] for i = 1:customers for t = 1:periods) <= data.UB)
    @constraint(model, load[i in 1:customers], sum(χ[i, t] for t = 1:periods) <= data.CM[i])
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
  return Dict(:ws => ws, :ps => ps, :dps => dps, :f => drdeed_fun, :loss => loss_fun)
end
export drdeed