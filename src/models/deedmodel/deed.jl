function deed(
  generators::Int = 6,
  periods::Int = 24;
  solver::Union{Nothing,MOI.AbstractOptimizer} = nothing,
  data::Union{Nothing,DeedData} = nothing,
)
  data = isnothing(data) ? getDRDeedData(5, generators, periods) : data
  solver =
    isnothing(solver) ?
    optimizer_with_attributes(
      Ipopt.Optimizer,
      MOI.Silent() => true,
      "sb" => "yes",
      "max_iter" => 10_000,
    ) : solver

  function deed_fun(x::Vector{<:Number})
    q = reshape(x, generators, periods)
    C = E = 0
    for j = 1:generators
      for t = 1:periods
        C += data.a[j] + data.b[j] * q[j, t] + data.c[j] * q[j, t]^2
        E += data.e[j] + data.f[j] * q[j, t] + data.g[j] * q[j, t]^2
      end
    end
    vcat(C, E)
  end
  function loss_fun(x::Vector{<:Number})
    q = reshape(x, generators, periods)
    sum(
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators) for
      t = 1:periods
    )
  end
  function ws(w::Vector{Float64})
    model = Model(solver)
    @variable(model, q[1:generators, 1:periods] >= 0)

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
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t]
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    @NLobjective(model, Min, w[1] * C + w[2] * E)
    optimize!(model)

    if (has_values(model))
      solution = DeedSolution(value.(q), value(C), value(E), value.(losst), w)
      return SuccessResult(model, data, solution)
    else
      return FailResult(model, data, termination_status(model))
    end
  end


  function ps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, q[1:generators, 1:periods] >= 0)

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
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t]
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    @variable(model, z)
    @NLconstraint(model, v[1] + z * d[1] - C >= 0)
    @NLconstraint(model, v[2] + z * d[2] - E >= 0)
    @objective(model, Min, z)
    optimize!(model)
    if (has_values(model))
      return value.(q), value.(z)
    else
      return nothing, nothing
    end

  end
  function dps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, q[1:generators, 1:periods] >= 0)

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
      losst[t in 1:periods],
      sum(q[j, t] * data.B[j, k] * q[k, t] for j = 1:generators for k = 1:generators)
    )
    @NLconstraint(
      model,
      powerbalance[t in 1:periods],
      sum(q[j, t] for j = 1:generators) == data.Dt[t] + losst[t]
    )

    @constraint(model, [t in 1:periods], data.pjmin .<= q[:, t] .<= data.pjmax)
    @constraint(model, [t in 1:(periods-1)], -data.DR .<= (q[:, t+1] - q[:, t]) .<= data.UR)
    @variable(model, t)
    @variable(model, w[1:2] >= 0)
    @constraint(model, dot(w, d) == 1)
    @NLconstraint(model, C * w[1] + E * w[2] >= t)
    @objective(model, Max, t - dot(w, v))
    optimize!(model)
    if (has_values(model))
      return value.(w)
    else
      return nothing
    end
  end
  return Dict(:ws => ws, :ps => ps, :dps => dps, :f => deed_fun, :loss => loss_fun)
end

export deed