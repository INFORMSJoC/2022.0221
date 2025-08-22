mip_piecewise_hevALL = Dict(
                           :jump   => (pwl, Pdem, nbinstants, losscoefSEin, losscoefSEout, p_FCmax, p_FCmin, p_ESmax, p_ESmin, ESmin, ESmax, ESinit, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters::Function=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_hev_jump(pwl,
                                                   Pdem,
                                                   nbinstants,
                                                   losscoefSEin,
                                                   losscoefSEout,
                                                   p_FCmax,
                                                   p_FCmin,
                                                   p_ESmax,
                                                   p_ESmin,
                                                   ESmin,
                                                   ESmax,
                                                   ESinit,
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
                           :grb   => (pwl, Pdem, nbinstants, losscoefSEin, losscoefSEout, p_FCmax, p_FCmin, p_ESmax, p_ESmin, ESmin, ESmax, ESinit, env, model, w, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters::Function=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_hev_grb(pwl,
                                                   Pdem,
                                                   nbinstants,
                                                   losscoefSEin,
                                                   losscoefSEout,
                                                   p_FCmax,
                                                   p_FCmin,
                                                   p_ESmax,
                                                   p_ESmin,
                                                   ESmin,
                                                   ESmax,
                                                   ESinit,
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

function mip_piecewise_hev_jump(pwl,
                        Pdem,
                        nbinstants,
                        losscoefSEin,
                        losscoefSEout,
                        p_FCmax,
                        p_FCmin,
                        p_ESmax,
                        p_ESmin,
                        ESmin,
                        ESmax,
                        ESinit,
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
                        warmstart=false)

    params = with_solver_parameters(currlb, cutoff, tilim)
    m = JuMP.Model(optimizer_with_attributes(with_solver, params...))

    nperiods = length(pwl)
    if nperiods != nbinstants
        println("errrrooorrrr", pwl, " ", nperiods, " ", nbinstants)
    end
    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        #println(t)
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    #println("debut vars")
    @variable(m, 0 <= w[i in 1 : nperiods] <= p_FCmax, base_name="Xfc")
    @variable(m, p_ESmin <= Xesin[i in 1 : nperiods] <=0, base_name="Xesin")
    @variable(m, 0 <= Xesout[i in 1 : nperiods] <= p_ESmax, base_name="Xesout")
    #@variable(m, Prodornot[i in 1 : nperiods], Bin, base_name="Prodornot")

    @variable(m, v[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="v")
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin, base_name="zbinary")
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="scontforbinary")

    #println("fin vars")

    obj = @expression(m, sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t])))

    @objective(m, Min, obj)

    #begin pwl modeling
    @constraint(m, [t in 1 : nperiods], sum(s[t, j] for j in 1 : length(pwl[t])) == 1)
    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * s[t, j] - v[t, j] <= 0)

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * s[t, j] - v[t, j] >= 0)
    # @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)

#    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(v[t, j] for j in 1 : length(pwl[t])) == 0)
    @constraint(m, [t in 1 : nperiods], w[t] - sum(v[t, j] for j in 1 : length(pwl[t])) == 0)

    #end pwl modeling


    #satisfy demand
    @constraint(m, [t in 1 : nperiods], w[t]+Xesout[t]+Xesin[t] == Pdem[t])

    #restore initial level
    @constraint(m, sum( (1.0+losscoefSEout)*Xesout[t]+(1.0+losscoefSEin)*Xesin[t] for t in 1 : nperiods) <= 0)

    #respect minimum energy level
    @constraint(m, [t in 1 : nperiods], sum( (1.0+losscoefSEout)*Xesout[i]+(1.0+losscoefSEin)*Xesin[i] for i in 1:t) <= ESinit - ESmin)

    #respect maximum energy level
    @constraint(m, [t in 1 : nperiods], sum( (1.0+losscoefSEout)*Xesout[i]+(1.0+losscoefSEin)*Xesin[i] for i in 1:t) >= ESinit - ESmax)

    #respect Pmin if Prod
    #@constraint(m, [t in 1 : nperiods], w[t] - p_FCmin*Prodornot[t] >= 0)

    #no prod if no prodornot
    #@constraint(m, [t in 1 : nperiods], w[t] - p_FCmax*Prodornot[t] <= 0)

    #println(m)
    #f = open("model.lp", "w")
    #print(f, m)
    #close(f)
    #lp_file = MathOptFormat.LP.Model()
    #MOI.copy_to(lp_file, m)
    #MOI.copy_to(lp_file, backend(model))
    #MOI.write_to_file(lp_file, "my_model.lp")

    optimize!(m)

    # if termination_status(m) != MOI.OPTIMAL
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        println(termination_status(m))
        return [], -1, -1
    else
        #println("SUNtprint", value.(v))
        wval = vec(value.(w))
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        obj0 = 0.0
        wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end
end


function mip_piecewise_hev_grb(pwl,
                        Pdem,
                        nbinstants,
                        losscoefSEin,
                        losscoefSEout,
                        p_FCmax,
                        p_FCmin,
                        p_ESmax,
                        p_ESmin,
                        ESmin,
                        ESmax,
                        ESinit,
                        env,
                        model,
                        w,
                        var1D2,
                        var2D1,
                        var2D2;
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        tilim = 3600,
                        print = false,
                        warmstart=false)

    nperiods = length(pwl)
    if nperiods != nbinstants
        println("errrrooorrrr", pwl, " ", nperiods, " ", nbinstants)
    end

    if !warmstart
        model = Gurobi.Optimizer()
    end

    params = gurobi_parameters(0.0, cutoff, tilim, warmstart, Integer(print))
    for p in params
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if print
        MOI.set(model, MOI.Name(), "pbHEV_"*string(nperiods)*"_grb.lp")
    end

    if !warmstart
        #println("debut vars")
        w=Array{MathOptInterface.VariableIndex, 1}(undef, nperiods)
        for t in 1:nperiods
            w[t],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0*p_FCmax))
            MOI.set(model, MOI.VariableName(), w[t], "Xfc("*string(t)*")")
        end
        """ SUNhevonoff
        onoff=Array{MathOptInterface.VariableIndex, 1}(undef, nperiods)
        for t in 1:nperiods
            onoff[t],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            MOI.set(model, MOI.VariableName(), onoff[t], "onoff("*string(t)*")")
        end
        """
        Xesin=Array{MathOptInterface.VariableIndex, 1}(undef, nperiods)
        for t in 1:nperiods
            Xesin[t],_ = MOI.add_constrained_variable(model, MOI.Interval(1.0*p_ESmin, 0.0))
            MOI.set(model, MOI.VariableName(), Xesin[t], "Xesin("*string(t)*")")
        end
        Xesout=Array{MathOptInterface.VariableIndex, 1}(undef, nperiods)
        for t in 1:nperiods
            Xesout[t],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0*p_ESmax))
            MOI.set(model, MOI.VariableName(), Xesout[t], "Xesout("*string(t)*")")
        end
        #println("fin vars")

        for t in 1:nperiods
            #satisfy demand
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, w[t]), MOI.ScalarAffineTerm(1.0, Xesout[t]), MOI.ScalarAffineTerm(1.0, Xesin[t])  ], 0.0),
                       MOI.EqualTo(1.0*Pdem[t])
                   );
            """SUNhevonoff
           #respect PFCmin
           MOI.add_constraint(
                      model,
                      MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, w[t]), MOI.ScalarAffineTerm(-1.0*p_FCmin, onoff[t]) ], 0.0),
                      MOI.GreaterThan(0.0)
            );
            """
            #respect minimum energy level
            MOI.add_constraint(
                      model,
                      MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0+losscoefSEout for i in 1:t],[1.0+losscoefSEin for i in 1:t]), vcat([Xesout[i] for i in 1:t],[Xesin[i] for i in 1:t])), 0.0),
                      MOI.LessThan(1.0*(ESinit - ESmin))
                  );
            #respect maximum energy level
            MOI.add_constraint(
                    model,
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0+losscoefSEout for i in 1:t],[1.0+losscoefSEin for i in 1:t]), vcat([Xesout[i] for i in 1:t],[Xesin[i] for i in 1:t])), 0.0),
                    MOI.GreaterThan(1.0*(ESinit - ESmax))
                );

        end

        #restore initial level
        MOI.add_constraint(
                   model,
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0+losscoefSEout for _ in 1:nperiods],[1.0+losscoefSEin for _ in 1:nperiods]), vcat([Xesout[t] for t in 1:nperiods],[Xesin[t] for t in 1:nperiods])), 0.0),
                   MOI.LessThan(0.0)
               );

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
       GRBwrite(model, "pbHEV_"*string(nperiods)*"_grb.lp")
    end

    status = GRBoptimize(model)
    if status != 0
       return [], -1, -1, env, model, w, var1D2, var2D1, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        println(termination_status(model))
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
        wval, obj0, obj, env, model, w, var1D2, var2D1, var2D2
    end

end
