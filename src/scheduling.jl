mip_piecewise_schedulingALL = Dict(
                           :jump   => (pwl, durations, demands, time_windows, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_scheduling_jump(pwl,
                                                   durations,
                                                   demands,
                                                   time_windows,
                                                   env,
                                                   model,
                                                   var1D1,
                                                   var1D2,
                                                   var2D1,
                                                   var2D2;
                                                   currlb,
                                                   cutoff,
                                                   tilim,
                                                   print,
                                                   warmstart),
                           :grb   => (pwl, durations, demands, time_windows, env, model, w, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=true, warmstart=false)-> mip_piecewise_scheduling_grb(pwl,
                                                   durations,
                                                   demands,
                                                   time_windows,
                                                   env,
                                                   model,
                                                   w,
                                                   var1D2,
                                                   var2D1,
                                                   var2D2;
                                                   currlb,
                                                   cutoff,
                                                   tilim,
                                                   print,
                                                   warmstart)
                       )

function mip_piecewise_scheduling_jump(pwl,
                        durations,
                        demands,
                        time_windows,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        var2D1,
                        var2D2;
                        with_solver = GLPK.Optimizer,
                        with_solver_parameters::Function = empty_parameters,
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        tilim = 3600,
                        print=false,
                        warmstart=false)

    params = with_solver_parameters(currlb, cutoff, tilim)
    m = JuMP.Model(optimizer_with_attributes(with_solver, params...))

    ntasks = length(durations)
    nperiods = length(pwl)

    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    #@variable(m, x[1 : ntasks, 1 : nperiods] >= 0)
    @variable(m, x[1 : ntasks, 1 : nperiods], Bin)
    @variable(m, w[1 : nperiods] >= 0)
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin)
    @variable(m, y[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)

    @constraint(m, [i in 1 : ntasks], sum(x[i, t] for t in 1 : nperiods if time_windows[i][1] <= t && t < time_windows[i][2]) >= durations[i])
    @constraint(m, [t in 1 : nperiods], w[t] - sum(demands[i] * x[i, t] for i in 1 : ntasks) == 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * y[t, j] - s[t, j] <= 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * y[t, j] - s[t, j] >= 0)
    @constraint(m, [t in 1 : nperiods], sum(y[t, j] for j in 1 : length(pwl[t])) == 1)

    for t in 1 : nperiods - 1
        feast = [i for i in 1 : ntasks if time_windows[i][1] <= t && t < time_windows[i][2]]
        feasnext = [i for i in 1 : ntasks if time_windows[i][1] <= t + 1 && t + 1 < time_windows[i][2]]
        if feast == feasnext
            println("equivalent periods $t and $(t + 1)")
            @constraint(m, sum(x[i, t] - x[i, t + 1] for i in feast) >= 0)
        end
    end

    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(y[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(y[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end
    # @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t]), l in 1 : nlog[t] if binary[t][j][l]], y[t, j] - z[t, l] <= 0)
    # @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t]), l in 1 : nlog[t] if !binary[t][j][l]], y[t, j] + z[t, l] <= 1)

    @constraint(m, [t in 1 : nperiods], w[t] - sum(s[t, j] for j in 1 : length(pwl[t])) == 0)
    @objective(m, Min, sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * y[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    optimize!(m)
    # if termination_status(m) != MOI.OPTIMAL
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
        obj0 = 0
        return wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end


function mip_piecewise_scheduling_grb(pwl,
                        durations,
                        demands,
                        time_windows,
                        env,
                        model,
                        w,
                        var1D2,
                        var2D1,
                        var2D2;
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        tilim = 3600,
                        print=true,
                        warmstart=false)

    ntasks = length(durations)
    nperiods = length(pwl)

    if !warmstart
        model = Gurobi.Optimizer()
    end

    params = gurobi_parameters(0.0, cutoff, tilim, warmstart, Integer(print))
    for p in params
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if print
        MOI.set(model, MOI.Name(), "pbSched_"*string(ntasks)*"_"*string(nperiods)*"_grb.lp")
    end

    if !warmstart
        w=Array{MathOptInterface.VariableIndex, 1}(undef, nperiods)
        for t in 1:nperiods
                w[t],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.set(model, MOI.VariableName(), w[t], "w("*string(t)*")")
                #MOI.add_constraint(model, x[i,j], MOI.GreaterThan(0.0))
        end

        x=Array{MathOptInterface.VariableIndex, 2}(undef, ntasks, nperiods)
        for i in 1:ntasks, t in 1:nperiods
                x[i,t],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
                #x[i,j] = MOI.add_variable(model)
                MOI.set(model, MOI.VariableName(), x[i,t], "x("*string(i)*","*string(t)*")")
                #MOI.add_constraint(model, x[i,j], MOI.GreaterThan(0.0))
        end

        #constraint: sum_{t in PERIODS} x_it >= DURATION[i],  forall i in TASKS
        for i in 1:ntasks
            MOI.add_constraint(
                       model,
                       #MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:nperiods], [x[i,j] for j in 1:nperiods]), 0.0),
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i,t]) for t in 1:nperiods if ((time_windows[i][1] <= t) && (t < time_windows[i][2]))], 0.0),
                       MOI.GreaterThan(1.0*durations[i]),
                   );
        end

        #constraint: w_t - sum_{i in TASKS} DEMANDS[i]*x_it = 0,  forall t in PERIODS
        for t in 1:nperiods
            MOI.add_constraint(
                       model,
                       #MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:nperiods], [x[i,j] for j in 1:nperiods]), 0.0),
                       #MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, w[t]), MOI.ScalarAffineTerm(-1.0*demands[i], x[i,t]) for i in 1:ntasks], 0.0),
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0], [-1.0*demands[i] for i in 1:ntasks]), vcat([w[t]], [x[i,t] for i in 1:ntasks])), 0.0),
                       MOI.EqualTo(0.0),
                   );
        end

        for t in 1 : nperiods - 1
            feast = [i for i in 1 : ntasks if time_windows[i][1] <= t && t < time_windows[i][2]]
            feasnext = [i for i in 1 : ntasks if time_windows[i][1] <= t + 1 && t + 1 < time_windows[i][2]]
            if feast == feasnext
                println("equivalent periods $t and $(t + 1)")
                MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0 for i in feast],[-1.0 for i in feast]), vcat([x[i,t] for i in feast],[x[i,t+1] for i in feast])), 0.0),
                           MOI.GreaterThan(0.0)
                       );
            end
        end

        #PWL objective coefs
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    for t in eachindex(pwl)
        breakx, breaky = CCGetBreakpoints(pwl[t], 1e-7)
        nbreaks = length(breakx)
        old_nbreaks = [copy(Cint(nbreaks))]
        ptr_old_nbreaks = Ptr{Cint}(UInt(pointer(old_nbreaks)))
        error = GRBgetpwlobj(model, grbXvar, ptr_old_nbreaks, Ptr{UInt8}(C_NULL), Ptr{UInt8}(C_NULL))
        prevx = zeros(Cdouble, old_nbreaks[1])
        prevy = zeros(Cdouble, old_nbreaks[1])
        if warmstart
            error = GRBgetpwlobj(model, t - 1, ptr_old_nbreaks, prevx, prevy)
            @assert(size(prevx, 1) == old_nbreaks[1], "sizes mismatch! $(size(prevx, 1)) != $(old_nbreaks[1])")
        end
        if nbreaks == old_nbreaks[1] && length(breakx) == length(prevx) && length(breaky) == length(prevy) && maximum(abs.(breakx - prevx)) < 1e-5 && maximum(abs.(breaky - prevy)) < 1e-5
            if print
                println("objective for this variable did not change!")
            end
        else
            error = Gurobi.GRBsetpwlobj(model, t - 1, nbreaks, breakx, breaky)
            if (error > 0)
                println(GRBgeterrormsg(env))
            end
        end
    end
    #end define objective

    if print
        GRBwrite(model, "pbSched_"*string(ntasks)*"_"*string(nperiods)*"_grb.lp")
    end

    status = GRBoptimize(model)
    if status != 0
        return [], -1, -1, env, model, w, var1D2, var2D1, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, w, var1D2, var2D1, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        obj0 = 0.0
        return wval, obj0, obj, env, model, w, var1D2, var2D1, var2D2
    end
end
