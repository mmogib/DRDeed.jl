function slimmodel(
  customers::Int = 3,
  generators::Int = 1,
  periods::Int = 2,
  solver::Union{Nothing,MOI.AbstractOptimizer} = nothing,
  data::Union{Nothing,SlimData} = nothing,
)
  # ipopt = optimizer_with_attributes(
  #   Ipopt.Optimizer,
  #   MOI.Silent() => true,
  #   "sb" => "yes",
  #   "max_iter" => 1000_000,
  # )
  # highs = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
  # pavito =
  #   optimizer_with_attributes(Pavito.Optimizer, "mip_solver" => highs, "cont_solver" => ipopt)
  #   M = 1000
  data = isnothing(data) ? getSlimData(customers, generators, periods) : data
  solver =
    isnothing(solver) ?
    optimizer_with_attributes(
      Ipopt.Optimizer,
      MOI.Silent() => true,
      "sb" => "yes",
      "max_iter" => 10_000,
    ) : solver
  function mpmodel(; w::Vector{Float64} = ones(2), obj::Vector{Symbol} = [:Z1])
    model = Model(solver)
    @variable(model, x[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, y[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, f[1:customers, 1:periods] >= 0)     #correct
    @variable(model, h[1:customers, 1:(periods+1)] >= 0) #correct
    @variable(model, ρi[1:customers, 1:periods] >= 0)  #correct
    @variable(model, si[1:customers, 0:periods] >= 0) #correct

    @variable(model, ρj[1:generators, 1:periods] >= 0)
    @variable(model, sj[1:generators, 0:periods] >= 0)

    #variable constraint (initial conditions)
    @constraint(model, [i in 1:customers], si[i, periods] == 0) #correct (6)
    @constraint(model, [i in 1:customers], si[i, 0] == 0)#correct (6)
    @constraint(model, [j in 1:generators], sj[j, 0] == 0) #correct(12)
    @constraint(model, [j in 1:generators], sj[j, periods] == 0) #correct (12)
    @constraint(model, [i in 1:customers], h[i, 1] == 0)#correct (6)
    @constraint(model, [i in 1:customers], h[i, periods+1] == 0) #correct (6)

    # i constrain
    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      f[i, t] + si[i, t] - si[i, t-1] - ρi[i, t] + sum(x[i, :, t]) == sum(y[i, :, t])
    ) #correct (1)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      h[i, t+1] + f[i, t] == h[i, t] + data.DemandDt[i, t]
    ) #correct (2)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      sum(x[i, :, t]) <= data.CProductionLimitst[i, t]
    ) #correct (3)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      ρi[i, t] <= data.CPowerUpperLimitst[i, t]
    )   #correct (4)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      si[i, t] <= data.CStorageUpperLimitst[i, t]
    )#correct (5)


    # j constrain
    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] - sj[j, t-1] - ρj[j, t] + sum(y[:, j, t]) == sum(x[:, j, t])
    )   #correct (8)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sum(y[:, j, t]) <= data.PProductionLimitst[j, t]
    )  #correct (9)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      ρj[j, t] <= data.PPowerUpperLimitst[j, t]
    ) #correct (10)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] <= data.PStorageUpperLimitst[j, t]
    )  #correct 


    @expression(
      model,
      U[i in 1:customers],
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods
      ) + sum(data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods)
    )
    @expression(
      model,
      V[j in 1:generators],
      sum(data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods) +
      sum(data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods)
    )

    solutions = Vector{ModelResult}(undef, length(obj))
    for (i, ob) in enumerate(obj)
      if ob == :Z1
        @objective(model, Max, w[1] * sum(U) + w[2] * sum(V))
      elseif ob == :Z2
        @variable(model, Z2)
        @constraint(model, Z2 .<= U)
        @objective(model, Max, Z2)
      elseif ob == :Z3
        @variable(model, Z3)
        @constraint(model, Z3 .<= V)
        @objective(model, Max, Z3)
      end
      optimize!(model)
      if has_values(model)
        solution = SlimSolution(
          value.(U),
          value.(V),
          value.(x),
          value.(y),
          value.(ρi),
          value.(ρj),
          value.(f),
          value.(h),
          value.(si),
          value.(sj),
        )
        solutions[i] = SuccessResult(model, data, solution)
      else
        solutions[i] = FailResult(model, data, String(Symbol(termination_status(model))))
      end
    end
    solutions
  end
  function ws(w::Vector{<:Number})
    model = Model(solver)
    @variable(model, x[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, y[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, f[1:customers, 1:periods] >= 0)     #correct
    @variable(model, h[1:customers, 1:(periods+1)] >= 0) #correct
    @variable(model, ρi[1:customers, 1:periods] >= 0)  #correct
    @variable(model, si[1:customers, 0:periods] >= 0) #correct
    @variable(model, ρj[1:generators, 1:periods] >= 0)
    @variable(model, sj[1:generators, 0:periods] >= 0)

    #variable constraint (initial conditions)
    @constraint(model, [i in 1:customers], si[i, periods] == 0) #correct (6)
    @constraint(model, [i in 1:customers], si[i, 0] == 0)#correct (6)
    @constraint(model, [j in 1:generators], sj[j, 0] == 0) #correct(12)
    @constraint(model, [j in 1:generators], sj[j, periods] == 0) #correct (12)
    @constraint(model, [i in 1:customers], h[i, 1] == 0)#correct (6)
    @constraint(model, [i in 1:customers], h[i, periods+1] == 0) #correct (6)

    # i constrain
    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      f[i, t] + si[i, t] - si[i, t-1] - ρi[i, t] + sum(x[i, :, t]) == sum(y[i, :, t])
    ) #correct (1)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      h[i, t+1] + f[i, t] == h[i, t] + data.DemandDt[i, t]
    ) #correct (2)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      sum(x[i, :, t]) <= data.CProductionLimitst[i, t]
    ) #correct (3)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      ρi[i, t] <= data.CPowerUpperLimitst[i, t]
    )   #correct (4)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      si[i, t] <= data.CStorageUpperLimitst[i, t]
    )#correct (5)


    # j constrain
    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] - sj[j, t-1] - ρj[j, t] + sum(y[:, j, t]) == sum(x[:, j, t])
    )   #correct (8)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sum(y[:, j, t]) <= data.PProductionLimitst[j, t]
    )  #correct (9)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      ρj[j, t] <= data.PPowerUpperLimitst[j, t]
    ) #correct (10)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] <= data.PStorageUpperLimitst[j, t]
    )  #correct 


    @expression(
      model,
      U[i in 1:customers],
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods
      ) + sum(data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods)
    )
    @expression(
      model,
      V[j in 1:generators],
      sum(data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods) +
      sum(data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods)
    )
    @objective(model, Min, -w[1] * sum(U) - w[2] * sum(V))
    optimize!(model)
    if has_values(model)
      solution = SlimSolution(
        value.(U),
        value.(V),
        value.(x),
        value.(y),
        value.(ρi),
        value.(ρj),
        value.(f),
        value.(h),
        value.(si),
        value.(sj),
      )
      return SuccessResult(model, data, solution)
    else
      return FailResult(model, data, String(Symbol(termination_status(model))))
    end
  end
  function vec2solution(z::Vector{<:Number})
    x, y, f, h, ρi, si, ρj, sj = vec2array(z)
    UV = slim_fun_raw(z)
    U = UV[:U]
    V = UV[:V]
    solution = SlimSolution(U, V, x, y, ρi, ρj, f, h, si, sj)
    SuccessResult(Model(), data, solution)
  end
  function vec2array(z::Vector{<:Number})
    c_g_p = customers * generators * periods
    c_p = customers * periods
    g_p = generators * periods
    x = reshape(z[1:c_g_p], customers, generators, periods)
    y = reshape(z[c_g_p+1:2*c_g_p], customers, generators, periods)
    f = reshape(z[2*c_g_p+1:2*c_g_p+c_p], customers, periods)
    h = reshape(z[2*c_g_p+c_p+1:2*c_g_p+2*c_p+customers], customers, periods + 1)
    ρi = reshape(z[2*c_g_p+2*c_p+customers+1:2*c_g_p+3*c_p+customers], customers, periods)
    si = reshape(z[2*c_g_p+3*c_p+customers+1:2*c_g_p+4*c_p+2*customers], customers, periods + 1)
    ρj = reshape(z[2*c_g_p+4*c_p+2*customers+1:2*c_g_p+4*c_p+2*customers+g_p], generators, periods)
    sj = reshape(
      z[2*c_g_p+4*c_p+2*customers+g_p+1:2*c_g_p+4*c_p+2*customers+2*g_p+generators],
      generators,
      periods + 1,
    )
    return (x, y, f, h, ρi, si, ρj, sj)
  end
  function slim_fun(z::Vector{<:Number})
    x, y, f, h, ρi, si, ρj, sj = vec2array(z)
    # x = value.(vcat(x[:], y[:], f[:], h[:], ρi[:], si[:], ρj[:], sj[:]))
    U =
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods for i = 1:customers
      ) + sum(
        data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods for i = 1:customers
      )

    V =
      sum(
        data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods for
        j = 1:generators
      ) + sum(
        data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods for j = 1:generators
      )
    vcat(U, V)
  end
  function slim_fun_raw(z::Vector{<:Number})
    x, y, f, h, ρi, si, ρj, sj = vec2array(z)
    # 
    U = map(1:customers) do i
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods
      ) + sum(data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods)
    end

    V = map(1:generators) do j
      sum(data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods) +
      sum(data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods)
    end
    Dict(:U => U, :V => V)
  end

  function ps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, x[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, y[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, f[1:customers, 1:periods] >= 0)     #correct
    @variable(model, h[1:customers, 1:(periods+1)] >= 0) #correct
    @variable(model, ρi[1:customers, 1:periods] >= 0)  #correct
    @variable(model, si[1:customers, 0:periods] >= 0) #correct
    @variable(model, ρj[1:generators, 1:periods] >= 0)
    @variable(model, sj[1:generators, 0:periods] >= 0)

    #variable constraint (initial conditions)
    @constraint(model, [i in 1:customers], si[i, periods] == 0) #correct (6)
    @constraint(model, [i in 1:customers], si[i, 0] == 0)#correct (6)
    @constraint(model, [j in 1:generators], sj[j, 0] == 0) #correct(12)
    @constraint(model, [j in 1:generators], sj[j, periods] == 0) #correct (12)
    @constraint(model, [i in 1:customers], h[i, 1] == 0)#correct (6)
    @constraint(model, [i in 1:customers], h[i, periods+1] == 0) #correct (6)

    # i constrain
    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      f[i, t] + si[i, t] - si[i, t-1] - ρi[i, t] + sum(x[i, :, t]) == sum(y[i, :, t])
    ) #correct (1)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      h[i, t+1] + f[i, t] == h[i, t] + data.DemandDt[i, t]
    ) #correct (2)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      sum(x[i, :, t]) <= data.CProductionLimitst[i, t]
    ) #correct (3)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      ρi[i, t] <= data.CPowerUpperLimitst[i, t]
    )   #correct (4)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      si[i, t] <= data.CStorageUpperLimitst[i, t]
    )#correct (5)


    # j constrain
    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] - sj[j, t-1] - ρj[j, t] + sum(y[:, j, t]) == sum(x[:, j, t])
    )   #correct (8)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sum(y[:, j, t]) <= data.PProductionLimitst[j, t]
    )  #correct (9)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      ρj[j, t] <= data.PPowerUpperLimitst[j, t]
    ) #correct (10)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] <= data.PStorageUpperLimitst[j, t]
    )  #correct 


    @expression(
      model,
      U[i in 1:customers],
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods
      ) + sum(data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods)
    )
    @expression(
      model,
      V[j in 1:generators],
      sum(data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods) +
      sum(data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods)
    )
    @variable(model, z)
    @constraint(model, (v[1] + z * d[1] + sum(U)) >= 0)
    @constraint(model, (v[2] + z * d[2] + sum(V)) >= 0)
    @objective(model, Min, z)
    optimize!(model)
    xv = value.(x)
    yv = value.(y)
    fv = value.(f)
    hv = value.(h)
    pv = value.(ρi)
    siv = value.(si)
    qv = value.(ρj)
    sjv = value.(sj)
    xx = vcat(xv[:], yv[:], fv[:], hv[:], pv[:], (siv.data)[:], qv[:], (sjv.data)[:])
    if (has_values(model))
      return xx, value.(z)
    else
      return nothing
    end
  end
  function dps(v::Vector{<:Number}, d::Vector{<:Number})
    model = Model(solver)
    @variable(model, x[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, y[1:customers, 1:generators, 1:periods] >= 0) #correct
    @variable(model, f[1:customers, 1:periods] >= 0)     #correct
    @variable(model, h[1:customers, 1:(periods+1)] >= 0) #correct
    @variable(model, ρi[1:customers, 1:periods] >= 0)  #correct
    @variable(model, si[1:customers, 0:periods] >= 0) #correct
    @variable(model, ρj[1:generators, 1:periods] >= 0)
    @variable(model, sj[1:generators, 0:periods] >= 0)

    #variable constraint (initial conditions)
    @constraint(model, [i in 1:customers], si[i, periods] == 0) #correct (6)
    @constraint(model, [i in 1:customers], si[i, 0] == 0)#correct (6)
    @constraint(model, [j in 1:generators], sj[j, 0] == 0) #correct(12)
    @constraint(model, [j in 1:generators], sj[j, periods] == 0) #correct (12)
    @constraint(model, [i in 1:customers], h[i, 1] == 0)#correct (6)
    @constraint(model, [i in 1:customers], h[i, periods+1] == 0) #correct (6)

    # i constrain
    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      f[i, t] + si[i, t] - si[i, t-1] - ρi[i, t] + sum(x[i, :, t]) == sum(y[i, :, t])
    ) #correct (1)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      h[i, t+1] + f[i, t] == h[i, t] + data.DemandDt[i, t]
    ) #correct (2)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      sum(x[i, :, t]) <= data.CProductionLimitst[i, t]
    ) #correct (3)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      ρi[i, t] <= data.CPowerUpperLimitst[i, t]
    )   #correct (4)

    @constraint(
      model,
      [i in 1:customers, t in 1:periods],
      si[i, t] <= data.CStorageUpperLimitst[i, t]
    )#correct (5)


    # j constrain
    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] - sj[j, t-1] - ρj[j, t] + sum(y[:, j, t]) == sum(x[:, j, t])
    )   #correct (8)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sum(y[:, j, t]) <= data.PProductionLimitst[j, t]
    )  #correct (9)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      ρj[j, t] <= data.PPowerUpperLimitst[j, t]
    ) #correct (10)

    @constraint(
      model,
      [j in 1:generators, t in 1:periods],
      sj[j, t] <= data.PStorageUpperLimitst[j, t]
    )  #correct 


    @expression(
      model,
      U[i in 1:customers],
      sum(
        data.Caitf[i, t] * f[i, t] +
        data.Caits[i, t] * si[i, t] +
        data.Caitp[i, t] * ρi[i, t] +
        data.Caith[i, t] * h[i, t] for t = 1:periods
      ) + sum(data.Caijtx[i, j, t] * x[i, j, t] for j = 1:generators for t = 1:periods)
    )
    @expression(
      model,
      V[j in 1:generators],
      sum(data.Cajtp[j, t] * ρj[j, t] + data.Cajts[j, t] * sj[j, t] for t = 1:periods) +
      sum(data.Caijty[i, j, t] * y[i, j, t] for i = 1:customers for t = 1:periods)
    )
    @variable(model, t)
    @variable(model, w[1:2] >= 0)
    @constraint(model, dot(w, d) == 1)
    @constraint(model, w[1] * sum(U) + sum(V) * w[2] <= t)
    @objective(model, Max, t - dot(w, v))
    optimize!(model)
    if (has_values(model))
      return value.(w)
    else
      return nothing
    end
    optimize!(model)
  end
  return Dict(
    :ws => ws,
    :ps => ps,
    :dps => dps,
    :f => slim_fun,
    :fraw => slim_fun_raw,
    :vec2array => vec2array,
    :vec2solution => vec2solution,
  )
end

function solveModel(::Type{SlimModel}, model::JuMP.Model)
  highs = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
  set_optimizer(model, highs)
  optimize!(model)
  # U = expression_by_name(model, "U")
  value.(model[:x])
end
export slimmodel