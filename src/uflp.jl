mip_piecewise_uflp_nlassALL = Dict(
                           :jump   => (pwl, fixed_costs, caps, demands, capcons, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_uflp_nlass_jump(pwl,
                                               fixed_costs,
                                               caps,
                                               demands,
                                               capcons,
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
                           :grb   => (pwl, fixed_costs, caps, demands, capcons, env, model, y, var1D2, w, var2D2, params, currlb, cutoff)-> mip_piecewise_uflp_nlass_grb(pwl,
                                               fixed_costs,
                                               caps,
                                               demands,
                                               capcons,
                                               env,
                                               model,
                                               y,
                                               var1D2,
                                               w,
                                               var2D2,
                                               params,
                                               currlb,
                                               cutoff)
                       )

mip_piecewise_uflp_nldepALL = Dict(
                          :jump   => (pwl, fixed_costs, asscosts, caps, demands, capcons, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_uflp_nldep_jump(pwl,
                                              fixed_costs,
                                              asscosts,
                                              caps,
                                              demands,
                                              capcons,
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
                          :grb   => (pwl, fixed_costs, asscosts, caps, demands, capcons, env, model, w, y, x, var2D2, params, currlb, cutoff)-> mip_piecewise_uflp_nldep_grb(pwl,
                                              fixed_costs,
                                              asscosts,
                                              caps,
                                              demands,
                                              capcons,
                                              env,
                                              model,
                                              w,
                                              y,
                                              x,
                                              var2D2,
                                              params,
                                              currlb,
                                              cutoff)
                      )

function mip_piecewise_uflp_nlass_jump(pwl,
                        fixed_costs,
                        caps,
                        demands,
                        capcons,
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
                        print = false,
                        warmstart = false)

    params = with_solver_parameters(currlb, cutoff, tilim)
    m = JuMP.Model(optimizer_with_attributes(with_solver, params...))

    ndeps = length(caps)
    ncusts = length(demands)
    nperiods = length(pwl)
    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        # println("pad = $(nlog[t]), and log = $(log(2, length(pwl[t]) + 1))")
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    @variable(m, y[1 : ndeps], Bin)
    @variable(m, w[1 : ndeps, 1 : ncusts] >= 0)
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin)
    @variable(m, x[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)

    @constraint(m, [j in 1 : ncusts], sum(w[i, j] for i in 1 : ndeps) >= 1)
    @constraint(m, [i in 1 : ndeps, j in 1 : ncusts], w[i, j] - y[i] <= 0)
    if capcons
        @constraint(m, [i in 1 : ndeps], sum(demands[j] * w[i, j] for j in 1 : ncusts) - caps[i] * y[i] <= 0)
    end

    # @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t]), l in 1 : nlog[t] if binary[t][j][l]], x[t, j] - y[t, l] <= 0)
    # @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t]), l in 1 : nlog[t] if !binary[t][j][l]], x[t, j] + y[t, l] <= 1)

    @constraint(m, [t in 1 : nperiods], sum(x[t, j] for j in 1 : length(pwl[t])) == 1)
    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(x[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(x[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * x[t, j] - s[t, j] <= 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * x[t, j] - s[t, j] >= 0)
    # @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)
    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(s[t, j] for j in 1 : length(pwl[t])) == 0)

    obj = @expression(m, LinearAlgebra.dot(fixed_costs, y) + sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * x[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    @objective(m, Min, obj)
    optimize!(m)
    println("status = $(termination_status(m))")
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        yval = round.(Int64, value.(y))
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        obj0 = LinearAlgebra.dot(fixed_costs, yval)
        return wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end


function mip_piecewise_uflp_nlass_grb(pwl,
                        fixed_costs,
                        caps,
                        demands,
                        capcons,
                        env,
                        model,
                        y,
                        var1D2,
                        w,
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

    grbparams = gurobi_parameters(0.0, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if params.printon
        MOI.set(model, MOI.Name(), "pbUFLPnlass_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable w_ij
        w=Array{MathOptInterface.VariableIndex, 2}(undef, ndeps, ncusts)
        for i in 1:ndeps
            for j in 1:ncusts
                w[i,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.set(model, MOI.VariableName(), w[i,j], "w("*string(i)*","*string(j)*")")
            end
        end
        y=Array{MathOptInterface.VariableIndex, 1}(undef, ndeps)
        for i in 1:ndeps
            y[i],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            MOI.set(model, MOI.VariableName(), y[i], "y("*string(i)*")")
        end

        for j in 1:ncusts
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, w[i,j]) for i in 1:ndeps], 0.0),
                       MOI.GreaterThan(1.0)
                   );
            for i in 1:ndeps
                cons = MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, w[i,j]), MOI.ScalarAffineTerm(-1.0, y[i]) ], 0.0),
                           MOI.LessThan(0.0)
                       )
                MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), cons, 1)
            end
        end

        for i in 1:ndeps
            if capcons
                MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([-1.0*caps[i]], [demands[j] for j in 1:ncusts]), vcat([y[i]],[w[i,j] for j in 1:ncusts])), 0.0),
                           MOI.LessThan(0.0),
                       );
            end
        end

        # obj = @expression(m, LinearAlgebra.dot(fixed_costs, y) + LinearAlgebra.dot(asscosts, x) + sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
        #start define objective
        MOI.set(
                   model,
                   MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(fixed_costs, y), 0.0),
               );
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    #PWL objective coefs
    for t in eachindex(pwl)
        i = CartesianIndices(w)[t][1]
        j = CartesianIndices(w)[t][2]
        grbXvar = (i-1)*ncusts+(j-1)
        
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
            # if print 
            #     println("objective for this variable did not change!")
            # end
        else
            error = Gurobi.GRBsetpwlobj(model, grbXvar, nbreaks, breakx, breaky)
            if (error > 0)
                println(GRBgeterrormsg(env))
            end
        end
    end
    #end define objective

    # if print
    #     Gurobi.GRBwrite(model, "pbUFLPnlass_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
    # end

    status = Gurobi.GRBoptimize(model)
    if status != 0
        return [], -1, -1, env, model, y, var1D2, w, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, y, var1D2, w, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), vec(w))
        yval = round.(Int64, MOI.get(model, MOI.VariablePrimal(), y))
        for i in 1 : nvars
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        obj0 = LinearAlgebra.dot(fixed_costs, yval)
        return wval, obj0, obj, env, model, y, var1D2, w, var2D2
    end
end

function mip_piecewise_uflp_nldep_jump(pwl,
                        fixed_costs,
                        asscosts,
                        caps,
                        demands,
                        capcons,
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
                        print = false,
                        warmstart = false)

    params = with_solver_parameters(currlb, cutoff, tilim)
    m = JuMP.Model(optimizer_with_attributes(with_solver, params...))

    ndeps = length(caps)
    ncusts = length(demands)
    nperiods = length(pwl)

    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        # println("pad = $(nlog[t]), and log = $(log(2, length(pwl[t]) + 1))")
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    @variable(m, y[1 : ndeps], Bin)
    @variable(m, x[1 : ndeps, 1 : ncusts] >= 0)
    @variable(m, w[1 : ndeps] >= 0)
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)
    @variable(m, z[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0)
    @variable(m, v[t in 1 : nperiods, j in 1 : nlog[t]], Bin)

    @constraint(m, [j in 1 : ncusts], sum(x[i, j] for i in 1 : ndeps) >= 1)
    @constraint(m, [i in 1 : ndeps, j in 1 : ncusts], x[i, j] - y[i] <= 0)
    if capcons
        @constraint(m, cap[i in 1 : ndeps], sum(demands[j] * x[i, j] for j in 1 : ncusts) - caps[i] * y[i] <= 0)
    end

    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(z[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - v[t, l] == 0)
        @constraint(m, sum(z[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + v[t, l] == 1)
    end

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * z[t, j] - s[t, j] <= 0)
    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * z[t, j] - s[t, j] >= 0)

    @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)

    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(s[t, j] for j in 1 : length(pwl[t])) == 0)

    @constraint(m, [i in 1 : ndeps], w[i] - sum(demands[j] * x[i, j] for j in 1 : ncusts) >= 0)

    println("cutoff = $cutoff")
    # obj = @expression(m, LinearAlgebra.dot(fixed_costs, y) + LinearAlgebra.dot(asscosts, x) + sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    obj = @expression(m, LinearAlgebra.dot(asscosts, x) + sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
    @objective(m, Min, obj)
    optimize!(m)
    # if termination_status(m) != MOI.OPTIMAL
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        yval = round.(Int64, value.(y))
        xval = value.(x)
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        obj0 = LinearAlgebra.dot(asscosts, xval)
        wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end



function mip_piecewise_uflp_nldep_grb(pwl,
                        fixed_costs,
                        asscosts,
                        caps,
                        demands,
                        capcons,
                        env,
                        model,
                        w,
                        y,
                        x,
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20)

    ndeps = length(caps)
    ncusts = length(demands)

    if !params.warmstartMILP
        model = Gurobi.Optimizer()
    end

    grbparams = gurobi_parameters(0.0, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if params.printon
        MOI.set(model, MOI.Name(), "pbUFLPnldep_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable x_ij
        w=Array{MathOptInterface.VariableIndex, 1}(undef, ndeps)
        for i in 1:ndeps
            w[i],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.set(model, MOI.VariableName(), w[i], "w("*string(i)*")")
        end
        x=Array{MathOptInterface.VariableIndex, 2}(undef, ndeps, ncusts)
        y=Array{MathOptInterface.VariableIndex, 1}(undef, ndeps)
        for i in 1:ndeps
            y[i],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            MOI.set(model, MOI.VariableName(), y[i], "y("*string(i)*")")
            for j in 1:ncusts
                x[i,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                MOI.set(model, MOI.VariableName(), x[i,j], "x("*string(i)*","*string(j)*")")
            end
        end

        for j in 1:ncusts
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i,j]) for i in 1:ndeps], 0.0),
                    #    MOI.GreaterThan(1.0)
                        MOI.EqualTo(1.0)
                   );
            for i in 1:ndeps
                cons = MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i,j]), MOI.ScalarAffineTerm(-1.0, y[i]) ], 0.0),
                           MOI.LessThan(0.0)
                       )
                MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), cons, 1)
            end
        end

        for i in 1:ndeps
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0], [-1.0*demands[j] for j in 1:ncusts]), vcat([w[i]],[x[i,j] for j in 1:ncusts])), 0.0),
                    #    MOI.GreaterThan(0.0)
                        MOI.EqualTo(0.0)
                   );
           if capcons
               MOI.add_constraint(
                          model,
                          MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([-1.0*caps[i]], [demands[j] for j in 1:ncusts]), vcat([y[i]],[x[i,j] for j in 1:ncusts])), 0.0),
                          MOI.LessThan(0.0),
                      );
           end
        end

        println("cutoff = $cutoff")
        # obj = @expression(m, LinearAlgebra.dot(fixed_costs, y) + LinearAlgebra.dot(asscosts, x) + sum(pwl[t][j].a * s[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))
        #start define objective
        MOI.set(
                   model,
                   MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                   MOI.ScalarAffineFunction(vec([MOI.ScalarAffineTerm(asscosts[i,j], x[i,j]) for i in 1:ndeps, j in 1:ncusts]), 0.0),
               );
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    #PWL objective coefs
    # for t in Base.throw_eachindex_mismatch_indices(pwl)
    for t in eachindex(pwl)
        breakx, breaky = CCGetBreakpoints(pwl[t], 1e-7)
        nbreaks = length(breakx)
        old_nbreaks = [copy(Cint(nbreaks))]
        ptr_old_nbreaks = Ptr{Cint}(UInt(pointer(old_nbreaks)))
        error = GRBgetpwlobj(model, t - 1, ptr_old_nbreaks, Ptr{UInt8}(C_NULL), Ptr{UInt8}(C_NULL))
        prevx = zeros(Cdouble, old_nbreaks[1])
        prevy = zeros(Cdouble, old_nbreaks[1])
        if params.warmstartMILP
            error = GRBgetpwlobj(model, t - 1, ptr_old_nbreaks, prevx, prevy)
            @assert(size(prevx, 1) == old_nbreaks[1], "sizes mismatch! $(size(prevx, 1)) != $(old_nbreaks[1])")
        end
        if nbreaks == old_nbreaks[1] && length(breakx) == length(prevx) && length(breaky) == length(prevy) && maximum(abs.(breakx - prevx)) < 1e-5 && maximum(abs.(breaky - prevy)) < 1e-5
            if params.printon
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

    # if print
    #     GRBwrite(model, "pbUFLPnldep_"*string(ndeps)*"_"*string(ncusts)*"_grb.lp")
    # end

    status = GRBoptimize(model)
    if status != 0
        return [], -1, -1, env, model, w, y, x, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        println("writing file")
        # GRBwrite(model, "unbounded.lp")
        return [], -1, -1, env, model, w, y, x, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        yval = round.(Int64, MOI.get(model, MOI.VariablePrimal(), y))
        xval = MOI.get(model, MOI.VariablePrimal(), vec(x))
        for i in 1 : ndeps
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        obj0 = LinearAlgebra.dot(asscosts, xval)
        wval, obj0, obj, env, model, w, y, x, var2D2
    end
end
