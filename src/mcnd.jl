mip_piecewise_mcnd_nodesALL = Dict(
                           :jump   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_mcnd_jump(pwl,
                                               nodes,
                                               arcs,
                                               commods,
                                               listofinputarcs,
                                               listofoutputarcs,
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
                           :grb   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, w, yy, xx, var2D2, params, currlb=-1e+20, cutoff=1e+20)-> mip_piecewise_mcnd_nodes_grb(pwl,
                                               nodes,
                                               arcs,
                                               commods,
                                               listofinputarcs,
                                               listofoutputarcs,
                                               env,
                                               model,
                                               w,
                                               yy,
                                               xx,
                                               var2D2,
                                               params,
                                               currlb,
                                               cutoff)
)


mip_piecewise_mcnd_arcsALL = Dict(
        :jump   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_mcnd_jump(pwl,
                            nodes,
                            arcs,
                            commods,
                            listofinputarcs,
                            listofoutputarcs,
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
        :grb   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, w, yy, xx, var2D2, params, currlb=-1e+20, cutoff=1e+20)-> mip_piecewise_mcnd_arcs_grb(pwl,
                            nodes,
                            arcs,
                            commods,
                            listofinputarcs,
                            listofoutputarcs,
                            env,
                            model,
                            w,
                            yy,
                            xx,
                            var2D2,
                            params,
                            currlb,
                            cutoff)
)


function add_common_constraints_jump(m, arcs, commods, listofinputarcs, listofoutputarcs,
      nbNode, nbArc, nbCommod, xx, yy, w)

        #println("ccc4")
        for i in 1:nbNode, p in 1:nbCommod
            rhs=0.0
            if(commods[p].from == i)
                rhs=commods[p].qty
            elseif (commods[p].to == i)
                rhs= - commods[p].qty
            end

            #println("ccc5")

            @constraint(m, sum( xx[listofoutputarcs[i][ij], p] for ij in 1:length(listofoutputarcs[i]))
                            - sum( xx[listofinputarcs[i][ij], p] for ij in 1:length(listofinputarcs[i]))
                            == rhs)#, base_name="flowconservation_"*string(i)*"_"*string(p)*"_" )
        end

        #println("ccc6")

        @constraint(m,[ij in 1:nbArc, p in 1:nbCommod], xx[ij, p] - min(commods[p].qty,arcs[ij].capa)*yy[ij] <= 0)
                    #,    base_name="noflowifnoarc_")

                    #println("ccc7")

        @constraint(m,[ij in 1:nbArc], sum(xx[ij, p] for p in 1:nbCommod) - arcs[ij].capa*yy[ij] <= 0)
                    #, base_name="arccapalimit_")

    #println("ccc8")

        @constraint(m,[i in 1:nbNode], sum(xx[listofinputarcs[i][ji], p] for p in 1:nbCommod for ji in 1:length(listofinputarcs[i])) - w[i] == 0)
                    #, base_name="arccapalimit_")

end


function add_common_constraints_grb(model, arcs, commods, listofinputarcs, listofoutputarcs,
      nbNode, nbArc, nbCommod, xx, yy, w, pbtype)

      if pbtype == :mcnd_nodes
        for i in 1:nbNode
            MOI.add_constraint(
                    model,
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat(vec([1.0 for p in 1:nbCommod, ji in 1:length(listofinputarcs[i])]), [-1.0]), vcat(vec([xx[listofinputarcs[i][ji], p] for p in 1:nbCommod, ji in 1:length(listofinputarcs[i])]), [w[i]])), 0.0),
                    MOI.EqualTo(0.0),
                );
        end
    elseif pbtype == :mcnd_arcs
        for i in 1:nbArc
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat(vec([1.0 for p in 1:nbCommod]), [-1.0]), vcat(vec([xx[i, p] for p in 1:nbCommod]), [w[i]])), 0.0),
                       MOI.LessThan(0.0),
                   );
        end
    end
        #println("ccc4")
        for i in 1:nbNode
            for p in 1:nbCommod
                rhs=0.0
                if(commods[p].from == i)
                    rhs=commods[p].qty
                elseif (commods[p].to == i)
                    rhs= - commods[p].qty
                end

                #println("ccc5")
                MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0 for ij in 1:length(listofoutputarcs[i])], [-1.0 for ij in 1:length(listofinputarcs[i])]), vcat([xx[listofoutputarcs[i][ij], p] for ij in 1:length(listofoutputarcs[i])], [xx[listofinputarcs[i][ij], p] for ij in 1:length(listofinputarcs[i])])), 0.0),
                           MOI.EqualTo(1.0*rhs),
                       );
           end

        end

        for ij in 1:nbArc
            let
                MOI.add_constraint(
                        model,
                        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([1.0 for p in 1:nbCommod], [-1.0*arcs[ij].capa]), vcat([xx[ij, p] for p in 1:nbCommod], [yy[ij]])), 0.0),
                        MOI.LessThan(0.0),
                    );#, base_name="arccapalimit_")
                # MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), cons, 3)
            end
            for p in 1:nbCommod
                cons = MOI.add_constraint(
                           model,
                           MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, xx[ij,p]), MOI.ScalarAffineTerm(-1.0*min(commods[p].qty,arcs[ij].capa), yy[ij])], 0.0),
                           MOI.LessThan(0.0),
                       )#,    base_name="noflowifnoarc_")
                MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), cons, 1)
            end
        end
end


function mip_piecewise_mcnd_jump(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
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

    nbNode = length(nodes)
    nbArc = length(arcs)
    nbCommod = length(commods)
    nperiods = length(pwl)

    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        #println(t)
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    #println("debut vars")
    @variable(m, yy[ij in 1 : nbArc], Bin, base_name="yy_("*string(arcs[ij].from)*"_"*string(arcs[ij].to)*")_")
    @variable(m, xx[ij in 1 : nbArc, p in 1 : nbCommod] >= 0, base_name="xx_("*string(arcs[ij].from)*"_"*string(arcs[ij].to)*")_")
    @variable(m, w[i in 1 : nbNode] >= 0, base_name="w_"*string(i)*"_")
    @variable(m, v[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="v_")
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin, base_name="z_binary_")
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="s_contforbinary_")

    #println("fin vars")

    obj = @expression(m, LinearAlgebra.dot([x.fixedcost for x in arcs], yy) + LinearAlgebra.dot([x.varcost for x in arcs], sum(xx,dims=2))    #)
                        + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
                        )
    #println("objjj")

    @objective(m, Min, obj)

    #println("cons")

    #begin pwl modeling
    @constraint(m, [t in 1 : nperiods], sum(s[t, j] for j in 1 : length(pwl[t])) == 1)
    for t in 1 : nperiods, l in 1 : nlog[t]
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if binary[t][j][l]) - z[t, l] == 0)
        @constraint(m, sum(s[t, j] for j in 1 : length(pwl[t]) if !binary[t][j][l]) + z[t, l] == 1)
    end

    #println("ccc1")

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMin * s[t, j] - v[t, j] <= 0)

    #println("ccc2")

    @constraint(m, [t in 1 : nperiods, j in 1 : length(pwl[t])], pwl[t][j].xMax * s[t, j] - v[t, j] >= 0)
    # @constraint(m, [t in 1 : nperiods], sum(z[t, j] for j in 1 : length(pwl[t])) <= 1)

    #println("ccc3")

    @constraint(m, [t in 1 : nperiods], w[CartesianIndices(w)[t]] - sum(v[t, j] for j in 1 : length(pwl[t])) == 0)
    #end pwl modeling

    add_common_constraints_jump(m, arcs, commods, listofinputarcs, listofoutputarcs,
          nbNode, nbArc, nbCommod, xx, yy, w)
#=
    #println("ccc4")
    for i in 1:nbNode, p in 1:nbCommod
        rhs=0.0
        if(commods[p].from == i)
            rhs=commods[p].qty
        elseif (commods[p].to == i)
            rhs= - commods[p].qty
        end

        #println("ccc5")

        @constraint(m, sum( xx[listofoutputarcs[i][ij], p] for ij in 1:length(listofoutputarcs[i]))
                        - sum( xx[listofinputarcs[i][ij], p] for ij in 1:length(listofinputarcs[i]))
                        == rhs)#, base_name="flowconservation_"*string(i)*"_"*string(p)*"_" )
    end

    #println("ccc6")

    @constraint(m,[ij in 1:nbArc, p in 1:nbCommod], xx[ij, p] - min(commods[p].qty,arcs[ij].capa)*yy[ij] <= 0)
                #,    base_name="noflowifnoarc_")

                #println("ccc7")

    @constraint(m,[ij in 1:nbArc], sum(xx[ij, p] for p in 1:nbCommod) - arcs[ij].capa*yy[ij] <= 0)
                #, base_name="arccapalimit_")

#println("ccc8")

    @constraint(m,[i in 1:nbNode], sum(xx[listofinputarcs[i][ji], p] for p in 1:nbCommod for ji in 1:length(listofinputarcs[i])) - w[i] == 0)
                #, base_name="arccapalimit_")
=#
#println("ccc9")


#=    obj = @expression(m, LinearAlgebra.dot([x.fixedcost for x in arcs], yy) + LinearAlgebra.dot([x.fixedcost for x in arcs], sum(xx,dims=2))    #)
                        + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
                        )]
#println("objjj")

    @objective(m, Min, obj)

#println("objjj")

=#
#println(m)

    optimize!(m)
    # if termination_status(m) != MOI.OPTIMAL
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        yval = round.(Int64, value.(yy))
        xval = value.(xx)
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        obj0 = LinearAlgebra.dot([x.fixedcost for x in arcs], yval) + LinearAlgebra.dot([x.varcost for x in arcs], sum(xval,dims=2))
        wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end

end


#=

function minlp_mcnd(func1,
                    func2,
                    pdata,
                    minlpsolver)

    my_f(x,y) = (x - 1)^2 + (y - 2)^2

    m = Model(Ipopt.Optimizer)

    #m = Model(() -> AmplNLWriter.Optimizer(minlpsolver))
    #m=Model(Ipopt.Optimizer)

#=
    m = Model(optimizer_with_attributes(AmplNLWriter.Optimizer("couenne"),
        "Threads" => 1,
        "FeasibilityTol" => 1e-9,
        "OptimalityTol" => 1e-9,
        "BestObjStop" => currlb - 1e-5,
        "IntFeasTol" => 1e-9,
        "Cutoff" => (1 + 1e-5) * cutoff,
        "OutputFlag" => 0,
        "TimeLimit" => tilim + 1)
        )
# =#
    nperiods = pdata.nbNode
    #println("debut vars")
    @variable(m, yy[ij in 1 : pdata.nbArc], Bin, base_name="yy_("*string(pdata.arcs[ij].from)*"_"*string(pdata.arcs[ij].to)*")_")
    @variable(m, xx[ij in 1 : pdata.nbArc, p in 1 : pdata.nbCommod] >= 0, base_name="xx_("*string(pdata.arcs[ij].from)*"_"*string(pdata.arcs[ij].to)*")_")
    @variable(m, w[i in 1 : pdata.nbNode] >= 0, base_name="w_"*string(i)*"_")

    @variable(m, w1[i in 1 : nperiods] >= 0, base_name="w1_"*string(i)*"_")
    @variable(m, w2[i in 1 : nperiods] >= 0, base_name="w2_"*string(i)*"_")
    @variable(m, s[i in 1 : nperiods], Bin, base_name="s_"*string(i)*"_")

    obj1 = @expression(m, LinearAlgebra.dot([x.fixedcost for x in pdata.arcs], yy) + LinearAlgebra.dot([x.varcost for x in pdata.arcs], sum(xx,dims=2))
    )    #)
    #println("fin vars")

    #my_obj_f=sum(func1[t](w1[t]) for t in 1:nperiods)
    #print()
    #register(model, :my_obj, nperiods, my_f, autodiff=true)
    #@NLexpression(m, obj, func1[1](w1[1])
    #            )
    #JuMP.setNLobjective(m, :Min, :(func1[1]($(w1[1]))))
    #=@NLexpression(m, obj,sin($(w1[1])^2))
    for t in 2:nperiods
        obj = obj+sin($(w1[t])^2)
    end
    print(obj)

    #@NLexpression(m, obj, sum(func1[t](w1[t]) for t in 1:nperiods)
    #            )
    #                    + LinearAlgebra.dot([func1[t](w1[t]) for t in 1 : nperiods], 1.0) + LinearAlgebra.dot([func2[t](w2[t]) for t in 1 : nperiods], 1.0)
    #            )
    #println("objjj")
    =#
    #my_f=func1[1]
    #register(m, :my_f, 1, my_f, autodiff=true)

    @NLexpression(m, obj, w1[1]^2)
    @NLexpression(m, obj2, w1[1]^2)
    @NLexpression(m, obj3, obj+obj2)
    #ddd=Algebra.sub(:x=>:($(w1[1])),func1[1])
    #print(ddd)

    #@NLexpression(m, obj4, func1[1])
    #@NLexpression(m, obj4, func1[1](w1[1]^2))

    register(m, :my_f, 2, my_f, autodiff=true)

    #@NLobjective(m, Min, w1[1]^2)
    @NLobjective(m, Min, my_f(w1[1], w2[2]))

    @constraint(m, [t in 1 : nperiods], w[t] - w1[t] - w2[t] == 0)
    @constraint(m, [t in 1 : nperiods], w1[t] + pdata.nodes[t].croisement * s[t] <= pdata.nodes[t].croisement)
    @constraint(m, [t in 1 : nperiods], w2[t] - (pdata.nodes[t].capainit + pdata.nodes[t].capaincrement) * s[t] <= 0.0)
    @constraint(m, [t in 1 : nperiods], w2[t] - pdata.nodes[t].croisement * s[t] >= 0)

    add_common_constraints(m, pdata.arcs, pdata.commods, pdata.listofinputarcs, pdata.listofoutputarcs,
          pdata.nbNode, pdata.nbArc, pdata.nbCommod, xx, yy, w)

    optimize!(m)

print(termination_status(m))

    # if termination_status(m) != MOI.OPTIMAL
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        yval = round.(Int64, value.(yy))
        xval = value.(xx)
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        xval, yval, wval, obj
    end

end

=#



function mip_piecewise_mcnd_nodes_grb(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
                        env,
                        model,
                        w,
                        yy,
                        xx,
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        )

    nbNode = length(nodes)
    nbArc = length(arcs)
    nbCommod = length(commods)
    nperiods = length(pwl)

    if !params.warmstartMILP
        model = Gurobi.Optimizer()
    end

    grbparams = gurobi_parameters(0.0, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if params.printon
        MOI.set(model, MOI.Name(), "pbmcnd_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(nperiods)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable w_i
        w=Array{MathOptInterface.VariableIndex, 1}(undef, nbNode)
        for i in 1:nbNode
            # w[i],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            w[i],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * nodes[i].capainit + nodes[i].capaincrement))
            MOI.set(model, MOI.VariableName(), w[i], "w("*string(i)*")")
        end

        #declare variable yy_i
        yy=Array{MathOptInterface.VariableIndex, 1}(undef, nbArc)
        for ij in 1:nbArc
            yy[ij],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            MOI.set(model, MOI.VariableName(), yy[ij], "yy("*string(ij)*")")
        end    #println("debut vars")


        #declare variable xx_ij
        xx=Array{MathOptInterface.VariableIndex, 2}(undef, nbArc, nbCommod)
        for ij in 1:nbArc
            for j in 1:nbCommod
                # xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * min(commods[j].qty,arcs[ij].capa)))
                MOI.set(model, MOI.VariableName(), xx[ij,j], "xx("*string(ij)*","*string(j)*")")
            end
        end

        #println("fin vars")


        #println("cons")

        add_common_constraints_grb(model, arcs, commods, listofinputarcs, listofoutputarcs,
              nbNode, nbArc, nbCommod, xx, yy, w, :mcnd_nodes)

              #start define objective
              #obj = @expression(m, LinearAlgebra.dot([x.fixedcost for x in arcs], yy) + LinearAlgebra.dot([x.varcost for x in arcs], sum(xx,dims=2))    #)
            #                      + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
            #                      )
        # sumxx = vec(xx)
        MOI.set(
                model,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        #         #  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([x.fixedcost for x in arcs], [x.varcost for x in arcs]), vcat(yy, sum(xx,dims=2))), 0.0),
        #          MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat([x.fixedcost for x in arcs], vec([x.varcost for x in arcs, p in 1:nbCommod])), vcat(yy, vec(xx))), 0.0),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([x.fixedcost for x in arcs], yy), 0.0),
             );
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    #PWL objective coefs
    for t in eachindex(pwl)
        # breakx, breaky = SUNGetBreakpoints(pwl[t])
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

    if params.printon
      GRBwrite(model, "pbMCND_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(eps)*"_grb.lp")
    end

    status = GRBoptimize(model)
    if status != 0
      return [], -1, -1, env, model, w, yy, xx, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, w, yy, xx, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        yval = round.(Int64, MOI.get(model, MOI.VariablePrimal(), yy))
        xval = zeros(nbArc, nbCommod)
        for ij in 1:nbArc, k in 1 : nbCommod
            xval[ij, k] = MOI.get(model, MOI.VariablePrimal(), xx[ij, k])
        end
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        # for t in 1:nperiods
        #     println("pwl[t](w[t]) = $(pwl[t](wval[t]))")
        # end
        obj2 = sum(pwl[t](wval[t]).val for t in eachindex(pwl) if wval[t] >= 1e-9)
        obj = MOI.get(model, MOI.ObjectiveValue())
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        # println.([x.varcost for x in arcs, p in 1 : nbCommod])
        # println()
        # println.(vec(xval))
        obj0 = LinearAlgebra.dot([x.fixedcost for x in arcs], yval)# + LinearAlgebra.dot([x.varcost for x in arcs], sum(xval, dims = 2))
        # @assert(abs(obj2 + obj0 - obj) < 1e-5, "obj = $obj, obj0 = $obj0, obj2 = $obj2")
        wval, obj0, obj, env, model, w, yy, xx, var2D2
    end

end

function mip_piecewise_mcnd_arcs_grb(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
                        env,
                        model,
                        w,
                        yy,
                        xx,
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20)

    nbNode = length(nodes)
    nbArc = length(arcs)
    nbCommod = length(commods)
    nperiods = length(pwl)

    if !params.warmstartMILP
        model = Gurobi.Optimizer()
    end

    grbparams = gurobi_parameters(0.0, cutoff, params.tilim, params.warmstartMILP, Integer(params.printon))
    for p in grbparams
        MOI.set(model, MOI.RawOptimizerAttribute(p[1]), p[2])
    end
    if params.printon
        MOI.set(model, MOI.Name(), "pbmcnd_Arcs_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(nperiods)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable w_i
        w=Array{MathOptInterface.VariableIndex, 1}(undef, nbArc)
        for i in 1:nbArc
            # w[i],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            w[i],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * arcs[i].capa))
            MOI.set(model, MOI.VariableName(), w[i], "w("*string(i)*")")
        end

        #declare variable yy_i
        yy=Array{MathOptInterface.VariableIndex, 1}(undef, nbArc)
        for ij in 1:nbArc
            yy[ij],_ = MOI.add_constrained_variable(model, MOI.ZeroOne())
            MOI.set(model, MOI.VariableName(), yy[ij], "yy("*string(ij)*")")
        end    #println("debut vars")


        #declare variable xx_ij
        xx=Array{MathOptInterface.VariableIndex, 2}(undef, nbArc, nbCommod)
        for ij in 1:nbArc
            for j in 1:nbCommod
                # xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
                xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * min(commods[j].qty,arcs[ij].capa)))
                MOI.set(model, MOI.VariableName(), xx[ij,j], "xx("*string(ij)*","*string(j)*")")
            end
        end

        #println("fin vars")


        #println("cons")

        add_common_constraints_grb(model, arcs, commods, listofinputarcs, listofoutputarcs,
              nbNode, nbArc, nbCommod, xx, yy, w, :mcnd_arcs)

              #start define objective
              #obj = @expression(m, LinearAlgebra.dot([x.fixedcost for x in arcs], yy) + LinearAlgebra.dot([x.varcost for x in arcs], sum(xx,dims=2))    #)
            #                      + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
            #                      )
        sumxx = vec(xx)
        MOI.set(
                model,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([x.fixedcost for x in arcs], yy), 0.0),
             );
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    end

    #PWL objective coefs
    for t in eachindex(pwl)
        # breakx, breaky = SUNGetBreakpoints(pwl[t])
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

    if params.printon
      GRBwrite(model, "pbMCND_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(eps)*"_grb.lp")
    end

    status = GRBoptimize(model)
    if status != 0
      return [], -1, -1, env, model, w, yy, xx, var2D2
    end

    status = MOI.get(model, MOI.TerminationStatus())
    # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, w, yy, xx, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        yval = round.(Int64, MOI.get(model, MOI.VariablePrimal(), yy))
        xval = zeros(nbArc, nbCommod)
        for ij in 1:nbArc, k in 1 : nbCommod
            xval[ij, k] = MOI.get(model, MOI.VariablePrimal(), xx[ij, k])
        end
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        # for t in 1:nperiods
        #     println("pwl[t](w[t]) = $(pwl[t](wval[t]))")
        # end
        obj2 = sum(pwl[t](wval[t]).val for t in eachindex(pwl) if wval[t] >= 1e-9)
        obj = MOI.get(model, MOI.ObjectiveValue())
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        # println.([x.varcost for x in arcs, p in 1 : nbCommod])
        # println()
        # println.(vec(xval))
        obj0 = LinearAlgebra.dot([x.fixedcost for x in arcs], yval)# + LinearAlgebra.dot([x.varcost for x in arcs], sum(xval, dims = 2))
        # @assert(abs(obj2 + obj0 - obj) < 1e-5, "obj = $obj, obj0 = $obj0, obj2 = $obj2")
        wval, obj0, obj, env, model, w, yy, xx, var2D2
    end

end
