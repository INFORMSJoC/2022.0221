mip_piecewise_transportationALL = Dict(
                           :cpx => (pwl, caps, demands, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)->mip_piecewise_transportation_cpx(pwl,
                                                   caps,
                                                   demands,
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
                                                   print),
                           :jump   => (pwl, caps, demands, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_transportation_jump(pwl,
                                                   caps,
                                                   demands,
                                                   env,
                                                   model,
                                                   var1D1,
                                                   var1D2,
                                                   var2D1,
                                                   var2D2;
                                                   currlb,
                                                   cutoff,
                                                   tilim,
                                                   print),
                           :grb   => (pwl, caps, demands, env, model, var1D1, var1D2, x, var2D2, params, currlb=-1e+20, cutoff=1e+20)-> mip_piecewise_transportation_grb(pwl,
                                                   caps,
                                                   demands,
                                                   env,
                                                   model,
                                                   var1D1,
                                                   var1D2,
                                                   x,
                                                   var2D2,
                                                   params,
                                                   currlb,
                                                   cutoff)
                       )

function mip_piecewise_transportation_jump(pwl,
                        caps,
                        demands,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        var2D1,
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20)

    params = with_solver_parameters(currlb, cutoff, tilim)
    m = JuMP.Model(optimizer_with_attributes(with_solver, params...))

    ndeps = length(caps)
    ncusts = length(demands)
    nperiods = length(pwl)
    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    @variable(m, w[1 : ndeps, 1 : ncusts] >= 0)
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin)
    @variable(m, x[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)

    @constraint(m, [j in 1 : ncusts], sum(w[i, j] for i in 1 : ndeps) >= demands[j])
    @constraint(m, [i in 1 : ndeps], sum(w[i, j] for j in 1 : ncusts) <= caps[i])

    @constraint(m, [t in 1 : nperiods], sum(x[t, j] for j in 1 : length(pwl[t])) == 1)
    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(x[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(x[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * x[t, j] - s[t, j] <= 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * x[t, j] - s[t, j] >= 0)
    # @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)
    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(s[t, j] for j in 1 : length(pwl[t])) == 0)

    obj = @expression(m, sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * x[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    @objective(m, Min, obj)
    write_to_file(m, "modeljump.lp")
    optimize!(m)
    println("status = $(termination_status(m))")
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        println("w=",value.(w))
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

function mip_piecewise_transportation_cpx(pwl,
                        caps,
                        demands,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        var2D1,
                        var2D2;
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        tilim = 3600,
                        print = true,
                        warmstart = false)

    ndeps = length(caps)
    ncusts = length(demands)
    nvars = ndeps*ncusts

    model = CPLEX.Optimizer()
    params = cplex_parameters(0.0, cutoff, tilim)
    for p in params
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if print
        MOI.set(model, CPLEX.PassNames(), true)
        MOI.set(model, MOI.Name(), "pbTransp_"*string(ndeps)*"_"*string(ncusts)*".lp")
    end

    #declare variable x_ij
    x=Array{MathOptInterface.VariableIndex, 2}(undef, ndeps, ncusts)
    for i in 1:ndeps
        for j in 1:ncusts
            x[i,j] = MOI.add_variable(model)
            MOI.set(model, MOI.VariableName(), x[i,j], "x("*string(i)*","*string(j)*")")
            MOI.add_constraint(model, x[i,j], MOI.GreaterThan(0.0))
        end
    end

    #declare variable f_ij, the pwl cost of x_ij
    f=Array{MathOptInterface.VariableIndex, 2}(undef, ndeps, ncusts)
    for i in 1:ndeps
        for j in 1:ncusts
            f[i,j] = MOI.add_variable(model)
            MOI.set(model, MOI.VariableName(), f[i,j], "f("*string(i)*","*string(j)*")")
        end
    end

    #constraint: sum_{j in CUSTOMERS} x_ij <= Capacity[i],  forall i in DEPOTS
    for i in 1:ndeps
        MOI.add_constraint(
                   model,
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:ncusts], [x[i,j] for j in 1:ncusts]), 0.0),
                   MOI.LessThan(1.0*caps[i]),
               );
    end

    #constraint: sum_{i in DEPOTS} x_ij >= Demand[j],  forall j in CUSTOMERS
    for j in 1:ncusts
        MOI.add_constraint(
                   model,
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:ndeps], [x[i,j] for i in 1:ndeps]), 0.0),
                   MOI.GreaterThan(1.0*demands[j]),
               );
    end

    #constraint PWL: link x_ij and f_ij for i in CUSTOMERS, j in DEPOTS
    indexdiff=nvars #difference between the indices of X and F
    for t in 1:length(pwl)
        i = CartesianIndices(x)[t][1]
        j = CartesianIndices(x)[t][2]
        cpxXvar = (i-1)*ncusts+(j-1)
        # breakx, breaky =SUNGetBreakpoints(pwl[t])
        breakx, breaky = CCGetBreakpoints(pwl[t])
        preslope = pwl[t][1].a
        postslope = pwl[t][end].a
        nbreaks = length(breakx)
        pwlname = "pwl("*string(i)*","*string(j)*")"
        CPXaddpwl(model.env, model.lp, cpxXvar+indexdiff, cpxXvar, preslope, postslope, nbreaks, breakx, breaky, pwlname)
    end

    #objective function: min sum_ij f_ij
    MOI.set(
               model,
               MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
               MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(nvars), vec(f)), 0.0),
           );
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    if print
        CPLEX.CPXwriteprob(model.env, model.lp, "pbTransp_"*string(ndeps)*"_"*string(ncusts)*".lp", "LP")
    end

    status = CPLEX.CPXmipopt(model.env, model.lp)
    if status != 0
        return [], -1, -1
    end

    #println("statusmiopt=", status)
    #solstat = CPLEX.CPXgetstat(model.env, model.lp)
    #println("solstat=", status)
    #status = MOI.get(model, MOI.PrimalStatus())

    status = MOI.get(model, MOI.TerminationStatus())
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        xval = MOI.get(model, MOI.VariablePrimal(), vec(x))
        println("xval=",xval)
        for i in 1 : nvars
            if abs(xval[i] - round(Int64, xval[i])) < 1e-9
                xval[i] = round(Int64, xval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        obj0 = 0.0
        return xval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end



function mip_piecewise_transportation_grb(pwl,
                        caps,
                        demands,
                        env,
                        model,
                        var1D1,
                        var1D2,
                        x,
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20)

    ndeps = length(caps)
    ncusts = length(demands)
    nvars = ndeps*ncusts

    if !params.warmstartMILP
        model = Gurobi.Optimizer()
    end

	# println.(pwl)
    grbparams = gurobi_parameters(0.0, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if params.printon
        MOI.set(model, MOI.Name(), "pbTransp_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable x_ij
        x=Array{MathOptInterface.VariableIndex, 2}(undef, ndeps, ncusts)
        for i in 1:ndeps, j in 1:ncusts
                x[i,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                #x[i,j] = MOI.add_variable(model)
                MOI.set(model, MOI.VariableName(), x[i,j], "x("*string(i)*","*string(j)*")")
                #MOI.add_constraint(model, x[i,j], MOI.GreaterThan(0.0))
        end

        #constraint: sum_{j in CUSTOMERS} x_ij <= Capacity[i],  forall i in DEPOTS
        """for i in 1:ndeps
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:ncusts], [x[i,j] for j in 1:ncusts]), 0.0),
                       MOI.LessThan(1.0*caps[i]),
                   );
        end"""
        for i in 1:ndeps
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i,j]) for j in 1:ncusts], 0.0),
                       MOI.LessThan(1.0*caps[i]),
                   );
        end

        #constraint: sum_{i in DEPOTS} x_ij >= Demand[j],  forall j in CUSTOMERS
        """for j in 1:ncusts
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for _ in 1:ndeps], [x[i,j] for i in 1:ndeps]), 0.0),
                       MOI.GreaterThan(1.0*demands[j]),
                   );
        end"""
        for j in 1:ncusts
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i,j]) for i in 1:ndeps], 0.0),
                       MOI.GreaterThan(1.0*demands[j]),
                   );
        end

        #PWL objective coefs
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    for t in eachindex(pwl)
        i = CartesianIndices(x)[t][1]
        j = CartesianIndices(x)[t][2]
        grbXvar = (i-1)*ncusts+(j-1)
        # breakx, breaky = SUNGetBreakpoints(pwl[t])
        breakx, breaky = CCGetBreakpoints(pwl[t])
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
            if params.printon
                println("objective for this variable did not change!")
            end
        else
            error = Gurobi.GRBsetpwlobj(model, grbXvar, nbreaks, breakx, breaky)
            if (error > 0)
                println(GRBgeterrormsg(env))
            end
        end
    end

    # if print
        # GRBwrite(model, "pbTransp_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
#	return  [], -1, -1, env, model, var1D1, var1D2, x, var2D2
    # end

    status = GRBoptimize(model)
    if status != 0
        return [], -1, -1, env, model, var1D1, var1D2, x, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, var1D1, var1D2, x, var2D2
    else
        xval = MOI.get(model, MOI.VariablePrimal(), vec(x))
        # println("xval=",xval)
        for i in 1 : nvars
            if abs(xval[i] - round(Int64, xval[i])) < 1e-9
                xval[i] = round(Int64, xval[i])
            end
        end
	    obj = MOI.get(model, MOI.ObjectiveValue())
        obj0 = 0.0
        return xval, obj0, obj, env, model, var1D1, var1D2, x, var2D2
    end
end
