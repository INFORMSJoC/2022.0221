mip_piecewise_knapALL = Dict(
                           :jump   => (pwl, object_U, knapsack_cap, object_weight, are_vars_integer, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_knap_jump(pwl,
                                               object_U,
                                               knapsack_cap,
                                               object_weight,
                                               are_vars_integer,
                                               env,
                                               model,
                                               var1D1,
                                               var1D2,
                                               var2D1,
                                               var2D2;
                                               with_solver,
                                               with_solver_parameters,
                                               currlb,
                                               cutoff,
                                               tilim,
                                               print,
                                               warmstart),
                           :grb   => (pwl, object_U, knapsack_cap, object_weight, are_vars_integer, env, model, var1D1, var1D2, w, var2D2, params, currlb, cutoff)-> mip_piecewise_knap_grb(pwl,
                                               object_U,
                                               knapsack_cap,
                                               object_weight,
                                               are_vars_integer,
                                               env,
                                               model,
                                               var1D1,
                                               var1D2,
                                               w,
                                               var2D2,
                                               params,
                                               currlb,
                                               cutoff)
                    )

function mip_piecewise_knap_jump(pwl,
                        object_U,
                        knapsack_cap,
                        object_weight,
                        are_vars_integer,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        var2D1,
                        var2D2,
                        params,
                        currlb,
                        cutoff)


    grbparams = params.with_solver_parameters(currlb, cutoff, params.tilim)
    m = JuMP.Model(optimizer_with_attributes(params.with_solver, grbparams...))
    set_time_limit_sec(m, params.tilim) # it seems that without this line the tilim is not respected ?

    nknap = length(knapsack_cap)
    nobj = length(object_U)
    nperiods = length(pwl)
    if nperiods != nknap*nobj
        println("ERROR incorrect nperiods, should be equal to nknap*nobj !!!!")
    end
    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    @variable(m, w[1 : nknap, 1 : nobj] >= 0, integer = are_vars_integer)
    @variable(m, v[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="v_")
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin, base_name="z_binary_")
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="s_contforbinary_")


    @constraint(m, [i in 1 : nknap], sum(object_weight[j]*w[i, j] for j in 1 : nobj) <= knapsack_cap[i])
    @constraint(m, [j in 1 : nobj], sum(w[i, j] for i in 1 : nknap) <= object_U[j])

    @constraint(m, [t in 1 : nperiods], sum(s[t, j] for j in 1 : length(pwl[t])) == 1)
    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * s[t, j] - v[t, j] <= 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * s[t, j] - v[t, j] >= 0)
    # @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)
    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(v[t, j] for j in 1 : length(pwl[t])) == 0)

    obj = @expression(m, sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    @objective(m, Min, obj)
    #println("SUN before optimize", "number of variables=")
    optimize!(m)
    println("status = $(termination_status(m))")
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        obj0 = 0.0
        return wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end

function mip_piecewise_knap_grb(pwl,
                        object_U,
                        knapsack_cap,
                        object_weight,
                        are_vars_integer,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        w,
                        var2D2,
                        params,
                        currlb,
                        cutoff)

    nknap = length(knapsack_cap)
    nobj = length(object_U)
    nvars = length(pwl)
    if nvars != nknap*nobj
        println("ERROR incorrect nvars, should be equal to nknap*nobj !!!!")
    end

    if !params.warmstartMILP
        model = Gurobi.Optimizer()
    end

    grbparams = gurobi_parameters(currlb, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end

    if params.printon
        MOI.set(model, MOI.Name(), "pbmcnd_"*string(nknap)*"_"*string(nobj)*"_"*string(eps)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable xx_ij
        w=Array{MathOptInterface.VariableIndex, 2}(undef, nknap, nobj)
        for i in 1:nknap
            for j in 1:nobj
                if are_vars_integer
                    w[i,j],_ = MOI.add_constrained_variable(model, MOI.Integer())
                    MOI.set(model, Gurobi.VariableAttribute("LB"), w[i,j], 0.0)
                else
                    w[i,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                end
                MOI.set(model, MOI.VariableName(), w[i,j], "w("*string(i)*","*string(j)*")")
            end
        end

        for i in 1:nknap
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(object_weight[j],w[i, j]) for j in 1:nobj], 0.0),
                       MOI.LessThan(knapsack_cap[i]),
                   );#,    base_name="noflowifnoarc_")
        end

        # for i in 1:nknap - 1
        #     mult = knapsack_cap[i] < knapsack_cap[i + 1] ? 1.0 : -1.0
        #     MOI.add_constraint(
        #                model,
        #                MOI.ScalarAffineFunction(vcat([MOI.ScalarAffineTerm(mult, w[i, j]) for j in 1:nobj], [MOI.ScalarAffineTerm(-mult, w[i + 1, j]) for j in 1:nobj]), 0.0),
        #                MOI.LessThan(0.0),
        #            );#,    base_name="noflowifnoarc_")
        # end
        
        for j in 1:nobj
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0,w[i, j]) for i in 1:nknap], 0.0),
                       MOI.LessThan(object_U[j]),
                   );#,    base_name="noflowifnoarc_")
        end

        #PWL objective coefs
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE) #min because the objective is negative
    end

    for t in 1:nvars
        i = CartesianIndices(w)[t][1]
        j = CartesianIndices(w)[t][2]
        grbXvar = (i-1)*nobj+(j-1)
        breakx, breaky = CCGetBreakpoints(pwl[t], 1e-7)
        nbreaks = length(breakx)
        old_nbreaks = [copy(Cint(nbreaks))]
        ptr_old_nbreaks = Ptr{Cint}(UInt(pointer(old_nbreaks)))
        error = GRBgetpwlobj(model, grbXvar, ptr_old_nbreaks, Ptr{UInt8}(C_NULL), Ptr{UInt8}(C_NULL))
        prevx = zeros(Cdouble, old_nbreaks[1])
        prevy = zeros(Cdouble, old_nbreaks[1])
        if params.warmstartMILP
            error = GRBgetpwlobj(model, grbXvar, ptr_old_nbreaks, prevx, prevy)
            @assert(size(prevx, 1) == old_nbreaks[1], "sizes mismatch! $(size(prevx, 1)) != $(old_nbreaks[1])")
        end
        if nbreaks == old_nbreaks[1] && length(breakx) == length(prevx) && length(breaky) == length(prevy) && maximum(abs.(breakx - prevx)) < 1e-5 && maximum(abs.(breaky - prevy)) < 1e-5
            if print
                println("objective for this variable did not change!")
            end
        else
            error = Gurobi.GRBsetpwlobj(model, grbXvar, nbreaks, breakx, breaky)
            if (error > 0)
                println(GRBgeterrormsg(env))
            end
        end
    end

    if params.printon
        GRBwrite(model, "pbKnap_"*string(nknap)*"_"*string(nobj)*string(eps)*"_grb.lp")
    end

    status = GRBoptimize(model)
    if status != 0
        return [], -1, -1, env, model, var1D1, var1D2, w, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, var1D1, var1D2, w, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), vec(w))
        for i in 1 : nvars
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        # for t in 1 : length(pwl)
        #     pwx = pwl[t](wval[t])
        #     println("lower bound of w[$t] = $(wval[t]) is $pwx")
        # end
        # obj = MOI.get(model, MOI.ObjectiveValue())
        obj = sum(pwl[t](wval[t]).val for t in 1 : length(pwl))
        obj0 = 0.0
        return wval, obj0, obj, env, model, var1D1, var1D2, w, var2D2
    end
end
