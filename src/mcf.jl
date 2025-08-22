mip_piecewise_mcf_nodesALL = Dict(
                           :jump   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_mcf_nodes_jump(pwl,
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
                           :grb   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, w, var1D2, xx, var2D2, params, currlb=-1e+20, cutoff=1e+20)-> mip_piecewise_mcf_nodes_grb(pwl,
                                               nodes,
                                               arcs,
                                               commods,
                                               listofinputarcs,
                                               listofoutputarcs,
                                               env,
                                               model,
                                               w,
                                               var1D2,
                                               xx,
                                               var2D2,
                                               params,
                                               currlb,
                                               cutoff,
                                               )
                    )

mip_piecewise_mcf_arcsALL = Dict(
    :jump   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, var1D1, var1D2, var2D1, var2D2; with_solver=GLPK.Optimizer, with_solver_parameters=empty_parameters, currlb=-1e+20, cutoff=1e+20, tilim=3600, print=false, warmstart=false)-> mip_piecewise_mcf_arcs_jump(pwl,
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
    :grb   => (pwl, nodes, arcs, commods, listofinputarcs, listofoutputarcs, env, model, w, var1D2, xx, var2D2, params, currlb=-1e+20, cutoff=1e+20)-> mip_piecewise_mcf_arcs_grb(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
                        env,
                        model,
                        w,
                        var1D2,
                        xx,
                        var2D2,
                        params,
                        currlb,
                        cutoff,
                        )
)


function mip_piecewise_mcf(pwl,
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
                        var2D2,
                        params,
                        currlb = -1e+20,
                        cutoff = 1e+20)


    grbparams = with_solver_parameters(currlb, cutoff, params.tilim)
    m = JuMP.Model(optimizer_with_attributes(params.with_solver, grbparams...))

    nperiods = length(pwl)
    nlog = [ceil(Int64, log(2, length(pwl[t]) + 1)) for t in 1 : nperiods]
    binary = [[falses(nlog[t]) for i in 1 : length(pwl[t])] for t in 1 : nperiods]
    for t in 1 : nperiods
        #println(t)
        compute_binary_representation!(nlog[t], binary[t], [i for i in 1 : length(pwl[t])], "")
    end

    #println("debut vars")
    #@variable(m, yy[ij in 1 : nbArc], Bin, base_name="yy_("*string(arcs[ij].from)*"_"*string(arcs[ij].to)*")_")
    @variable(m, xx[ij in 1 : nbArc, p in 1 : nbCommod] >= 0, base_name="xx_("*string(arcs[ij].from)*"_"*string(arcs[ij].to)*")_")
    @variable(m, w[i in 1 : nbNode] >= 0, base_name="w_"*string(i)*"_")
    @variable(m, v[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="v_")
    @variable(m, z[t in 1 : nperiods, j in 1 : nlog[t]], Bin, base_name="z_binary_")
    @variable(m, s[t in 1 : nperiods, j in 1 : length(pwl[t])] >= 0, base_name="s_contforbinary_")

    #println("fin vars")

    obj = @expression(m, sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
                        )

    #=
    obj = @expression(m, LinearAlgebra.dot([x.varcost for x in arcs], sum(xx,dims=2))    #)
                        + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * s[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
                    )
    =#

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

    #@constraint(m,[ij in 1:nbArc, p in 1:nbCommod], xx[ij, p] - min(commods[p].qty,arcs[ij].capa)*yy[ij] <= 0)
                #,    base_name="noflowifnoarc_")

                #println("ccc7")

    @constraint(m,[ij in 1:nbArc], sum(xx[ij, p] for p in 1:nbCommod) <= arcs[ij].capa)
                #, base_name="arccapalimit_")

    #println("ccc8")

    @constraint(m,[i in 1:nbNode], sum(xx[listofinputarcs[i][ji], p] for p in 1:nbCommod for ji in 1:length(listofinputarcs[i])) - w[i] == 0)
                #, base_name="arccapalimit_")

                #println("ccc9")


#=    obj = @expression(m, LinearAlgebra.dot([x.fixedcost for x in arcs], yy) + LinearAlgebra.dot([x.fixedcost for x in arcs], sum(xx,dims=2))    #)
                        + sum(pwl[t][j].a * v[t, j] + pwl[t][j].b * z[t, j] for t in 1 : nperiods, j in 1 : length(pwl[t]))
                        )]
println("objjj")

    @objective(m, Min, obj)

println("objjj")

=#
#println(m)

    optimize!(m)
    # if termination_status(m) != MOI.OPTIMAL
    if !in(termination_status(m), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        wval = vec(value.(w))
        #yval = round.(Int64, value.(yy))
        xval = value.(xx)
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = objective_value(m)
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        # obj0 = LinearAlgebra.dot([x.varcost for x in arcs], sum(xval,dims=2))
        obj0 = 0
        wval, obj0, obj, env, model, var1D1, var1D2, var2D1, var2D2
    end

end

#=
function mip_piecewise_mcndtodeleete(pwl,
                        infoGraph,
                        commodity,
                        infoNode,
                        nbNode,
                        nbArc,
                        nbCommod,
                        mg,
                        G;
                        currlb = -1e+20,
                        cutoff = 1e+20,
                        tilim = 3600)
"""
"""
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
        "Threads" => 1,
        "MIPGap" => 0.0,
        "MIPGapAbs" => 0.0,
        "FeasibilityTol" => 1e-9,
        "OptimalityTol" => 1e-9,
        "BestObjStop" => currlb - 1e-5,
        # "Symmetry" => 2,
        "IntFeasTol" => 1e-9,
        "Cutoff" => (1 + 1e-5) * cutoff,
        # "MIPFocus" => 3,
        # "VarBranch" => 3,
        "TimeLimit" => tilim + 1)
        )
"""
"""
    notreModele = Model(optimizer_with_attributes(CPLEX.Optimizer,
        "CPXPARAM_Threads" =>1,
        "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
        "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 0.0,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPRELAX" => 1e-9,
        "CPX_PARAM_EPOPT" => 1e-9,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_CUTUP" => (1 + 1e-5) * cutoff,
        "CPX_PARAM_TILIM" => tilim + 1)
         )


             for i in 1:length(edges(mg))
                 e = Edge(infoGraph[i, :from], infoGraph[i, :to])
                 set_prop!(mg, e, :active, @variable(notreModele,binary = true,base_name="varedgeactive__"*string(infoGraph[i, :from])*"_"*string(infoGraph[i, :to])*"__"*string(i)*"__"))
                 set_prop!(mg, e, :flotCommod, @variable(notreModele, [1:nbCommod],lower_bound = 0,base_name="flotedge__"*string(infoGraph[i, :from])*"_"*string(infoGraph[i, :to])*"__"*string(i)*"__commod"))
                 @constraint(notreModele, sum(get_prop(mg, e, :flotCommod)) - infoGraph[i,:capacity] * get_prop(mg, e, :active) <=0 )

                 #println("for edge e=", e, "i=", i, "edge(mg)(i)", edges(mg), "activevar=", get_prop(mg, e,:active), " flotcommod=", get_prop(mg, e,:flotCommod), "\n\n\n")

             end

             for v in vertices(mg)
                 #println("verticemip=",v)
                 #inEdges = Edge.(inneighbors(G,v),v) #SUN
                 inEdges = Edge[]
                 for ii in 1:length(inneighbors(G,v))
                     push!(inEdges, Edge(inneighbors(G,v)[ii],v))
                 end
                 #println("for vertex", v, "theinedgesare", inEdges, "fromSUN")

                 #outEdges = Edge.(v,outneighbors(G,v)) #SUN
                 outEdges = Edge[]
                 for ii in 1:length(outneighbors(G,v))
                     push!(outEdges, Edge(v,outneighbors(G,v)[ii]))
                 end

                 #totalOut = sum(sum.(get_prop.(mg, outEdges,:flotCommod))) #SUN
                 totalIn = 0
                 for ii in 1:length(inEdges)
                     tmp = get_prop(mg, inEdges[ii],:flotCommod)
                     #println("for vertex", v, "theindge=", inEdges[ii], "tmp=", tmp, "\n\n\n")
                     for vv in 1:length(tmp)
                         totalIn = totalIn+tmp[vv]
                     end
                 end
                 #println("for vertex", v, "totalIn=", totalIn, "fromSUN", "\n\n\n")
                 #totalOut = sum(sum.(get_prop.(mg, outEdges,:flotCommod))) #SUN
                 totalOut = 0
                 for ii in 1:length(outEdges)
                     tmp = get_prop(mg, outEdges[ii],:flotCommod)
                     for vv in 1:length(tmp)
                         totalOut = totalOut+tmp[vv]
                     end
                 end


                 set_prop!(mg, v, :inFlow, totalIn)
                 set_prop!(mg, v, :outFlow, totalOut)
                 #set_prop!(mg, v, :upgraded, @variable(notreModele,binary = true,base_name="nodeupgradedornot__"*string(v)*"__"))

                 #Contraintes de conservation de flow par commodite
                 #inFlow = get_prop.(mg, inEdges,:flotCommod)
                 inFlow = []
                 for ii in 1:length(inEdges)
                     tmp = get_prop(mg, inEdges[ii],:flotCommod)
                     push!(inFlow, tmp)
                 end
                 #outFlow = get_prop.(mg, outEdges,:flotCommod)
                 outFlow = []
                 for ii in 1:length(outEdges)
                     tmp = get_prop(mg, outEdges[ii],:flotCommod)
                     push!(outFlow, tmp)
                 end

                 for c in 1:nbCommod

                     if length(inFlow) >= 1
                         inFlowCommod = sum(map(x->x[c],inFlow))
                     else
                         inFlowCommod = zero(AffExpr)
                     end

                     if length(outFlow) >= 1
                         outFlowCommod = sum(map(x->x[c],outFlow))
                     else
                         outFlowCommod = zero(AffExpr)
                     end


                     #si on est une usine ou un point d arrivee
                     if commodity[c,:from] == v
                    @constraint(notreModele, outFlowCommod - inFlowCommod == commodity[c,:quantity])

                     elseif commodity[c,:to] == v
                    @constraint(notreModele, outFlowCommod - inFlowCommod == -commodity[c,:quantity])
                     else
                        @constraint(notreModele, outFlowCommod - inFlowCommod == 0)
                     end

                 end

                 #Capacite augmentee ou non du noeud

                 @constraint(notreModele, totalIn <= infoNode[v,:InitCapacity] + infoNode[v,:IncremCapacity])
                 #@constraint(notreModele, totalOut - infoNode[v,:InitCapacity] - infoNode[v,:IncremCapacity]* props(mg,v)[:upgraded] <=0)

             end


             for v in vertices(mg)
                 #println("verticemip2=",v)
                 set_prop!(mg, v, :flowPiece, @variable(notreModele,[1:length(v_pwl[v])],base_name="flow__"*string(v)*"__perpiece"))
                 set_prop!(mg, v, :activationPiece, @variable(notreModele,[1:length(v_pwl[v])],binary = true,base_name="active__"*string(v)*"__perpiece"))

                 xMin = getfield.(v_pwl[v],:xMin)
                 xMax = getfield.(v_pwl[v],:xMax)

                 #contrainte pour etre sur un seul morceau de la lineariation

                 @constraint(notreModele, get_prop(mg, v,:flowPiece) .- xMin .* get_prop(mg, v,:activationPiece) .>= 0)
                 @constraint(notreModele, get_prop(mg, v,:flowPiece) .- xMax .* get_prop(mg, v,:activationPiece) .<= 0)
                 @constraint(notreModele, get_prop(mg, v,:inFlow)- sum(get_prop(mg, v,:flowPiece)) == 0)

                 @constraint(notreModele, sum(get_prop(mg, v,:activationPiece)) == 1)
                 cout=sum(getfield.(v_pwl[v],:a) .* get_prop(mg, v,:flowPiece) .+ getfield.(v_pwl[v],:b) .* get_prop(mg, v,:activationPiece))

                 set_prop!(mg, v, :nodeCost, cout)




             end

             activationCostArc=0
             arcFlowCost=0

             for i in 1:length(edges(mg))
                 e = Edge(infoGraph[i, :from], infoGraph[i, :to])
                 activationCostArc = activationCostArc + infoGraph.fixed_cost[i] * get_prop(mg, e, :active)
                 arcFlowCost = arcFlowCost + infoGraph.variable_cost[i] * sum(get_prop(mg, e, :flotCommod))
             end
             #println(activationCostArc)

             #incluant upgrade
             #nodeFlowCost = sum( get_prop.(mg, vertices(mg), :nodeCost))
             nodeFlowCost=0
             for ii in 1:length(vertices(mg))
                 #println("verticemip3=",ii)
                 tmp = get_prop(mg, vertices(mg)[ii],:nodeCost)
                 nodeFlowCost = nodeFlowCost+tmp
             end

             #println("beforesolvemilp")

             @objective(notreModele, Min, activationCostArc + arcFlowCost + nodeFlowCost)




    optimize!(notreModele)
    # if termination_status(m) != MOI.OPTIMAL
    println("status = $(termination_status(notreModele))")
    if !in(termination_status(notreModele), [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1
    else
        """
        wval = vec(value.(w))
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        """

        #SUN to improve BEGIN
        v_val = []
        for v in vertices(mg)
            theval = sum(JuMP.value.(get_prop(mg, v, :flowPiece)))
            push!(v_val, theval)
        end
        wval=v_val
        arcCost = infoGraph.variable_cost .* map(x->value(sum(get_prop(mg,x,:flotCommod))),edges(mg))
        arcFixedCost = infoGraph.fixed_cost .*map(x->value(sum(get_prop(mg,x,:active))),edges(mg))
        #nodeCost = map(x-> value(sum(get_prop(mg,x,:nodeCost))) ,vertices(mg))
        #SUN to improve END


        obj = objective_value(notreModele)
        obj0 = arcFixedCost + arcCost
        return wval, obj0, obj
    end
end
"""
=#


function mip_piecewise_mcf_nodes_grb(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
                        env,
                        model,
                        w,
                        var1D2,
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
        MOI.set(model, MOI.Name(), "pbmcnd_Nodes"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(nperiods)*"_grb.lp")
    end

    if !params.warmstartMILP
        #declare variable w_i
        w=Array{MathOptInterface.VariableIndex, 1}(undef, nbNode)
        for i in 1:nbNode
            w[i],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * nodes[i].capainit + nodes[i].capaincrement))
            MOI.set(model, MOI.VariableName(), w[i], "w("*string(i)*")")
        end

        #declare variable xx_ij
        xx=Array{MathOptInterface.VariableIndex, 2}(undef, nbArc, nbCommod)
        for ij in 1:nbArc
            for j in 1:nbCommod
                xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * min(commods[j].qty, arcs[ij].capa)))
                MOI.set(model, MOI.VariableName(), xx[ij,j], "xx("*string(ij)*","*string(j)*")")
            end
        end



        #println("ccc4")
        for i in 1:nbNode
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat(vec([1.0 for p in 1:nbCommod, ji in 1:length(listofinputarcs[i])]), [-1.0]), vcat(vec([xx[listofinputarcs[i][ji], p] for p in 1:nbCommod, ji in 1:length(listofinputarcs[i])]), [w[i]])), 0.0),
                       MOI.LessThan(0.0),
                   );
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
                       );#, base_name="flowconservation_"*string(i)*"_"*string(p)*"_" )
            end
        end

        for ij in 1:nbArc
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for p in 1:nbCommod], [xx[ij, p] for p in 1:nbCommod]), 0.0),
                       MOI.LessThan(1.0*arcs[ij].capa),
                   );#, base_name="arccapalimit_")
        end
    #     MOI.set(
    #             model,
    #             MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    #             MOI.ScalarAffineFunction(vec([MOI.ScalarAffineTerm(arcs[ij].varcost, xx[ij, p]) for ij in 1:nbArc for p in 1:nbCommod]), 0.0),
    #         );
    #    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
   end

   #PWL objective coefs
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

   if params.printon
     GRBwrite(model, "pbMCND_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(eps)*"_grb.lp")
   end

   status = GRBoptimize(model)
   if status != 0
     return [], -1, -1, env, model, w, var1D2, xx, var2D2
   end

   status = MOI.get(model, MOI.TerminationStatus())
   # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, w, var1D2, xx, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        #yval = round.(Int64, value.(yy))
        xval = Array{Float64, 2}(undef, nbArc, nbCommod)
        for ij in 1:nbArc, k in 1 : nbCommod
            xval[ij, k] = MOI.get(model, MOI.VariablePrimal(), xx[ij, k])
        end
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        # obj0 = LinearAlgebra.dot([x.varcost for x in arcs], sum(xval,dims=2))
        # println("obj0 = $obj0")
        obj0 = 0.0
        wval, obj0, obj, env, model, w, var1D2, xx, var2D2
    end

end

function mip_piecewise_mcf_arcs_grb(pwl,
                        nodes,
                        arcs,
                        commods,
                        listofinputarcs,
                        listofoutputarcs,
                        env,
                        model,
                        w,
                        var1D2,
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
            w[i],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * arcs[i].capa))
            MOI.set(model, MOI.VariableName(), w[i], "w("*string(i)*")")
        end

        #declare variable xx_ij
        xx=Array{MathOptInterface.VariableIndex, 2}(undef, nbArc, nbCommod)
        for ij in 1:nbArc
            for j in 1:nbCommod
                xx[ij,j],_ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0 * min(commods[j].qty, arcs[ij].capa)))
                MOI.set(model, MOI.VariableName(), xx[ij,j], "xx("*string(ij)*","*string(j)*")")
            end
        end



        #println("ccc4")
        for i in 1:nbArc
            # @constraint(model, sum(xx[i, p] for p in 1:nbCommod) - w[i] <= 0)
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat(vec([1.0 for p in 1:nbCommod]), [-1.0]), vcat(vec([xx[i, p] for p in 1:nbCommod]), [w[i]])), 0.0),
                       MOI.LessThan(0.0),
                   );
        end
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
                       );#, base_name="flowconservation_"*string(i)*"_"*string(p)*"_" )
            end
        end

        for ij in 1:nbArc
            MOI.add_constraint(
                       model,
                       MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0 for p in 1:nbCommod], [xx[ij, p] for p in 1:nbCommod]), 0.0),
                       MOI.LessThan(1.0*arcs[ij].capa),
                   );#, base_name="arccapalimit_")
       end
       MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
   end

   #PWL objective coefs
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

   if params.printon
     GRBwrite(model, "pbMCND_"*string(nbNode)*"_"*string(nbArc)*"_"*string(nbCommod)*"_"*string(eps)*"_grb.lp")
   end

   status = GRBoptimize(model)
   if status != 0
     return [], -1, -1, env, model, w, var1D2, xx, var2D2
   end

   status = MOI.get(model, MOI.TerminationStatus())
   # if termination_status(m) != MOI.OPTIMAL
    if !in(status, [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.SOLUTION_LIMIT, MOI.OTHER_LIMIT])
        return [], -1, -1, env, model, w, var1D2, xx, var2D2
    else
        wval = MOI.get(model, MOI.VariablePrimal(), w)
        #yval = round.(Int64, value.(yy))
        # xval = Array{Float64, 2}(undef, nbArc, nbCommod)
        # for ij in 1:nbArc, k in 1 : nbCommod
        #     xval[ij, k] = MOI.get(model, MOI.VariablePrimal(), xx[ij, k])
        # end
        for i in 1 : nperiods
            if abs(wval[i] - round(Int64, wval[i])) < 1e-9
                wval[i] = round(Int64, wval[i])
            end
        end
        obj = MOI.get(model, MOI.ObjectiveValue())
        # obj0 = LinearAlgebra.dot(fixed_costs, yval) + LinearAlgebra.dot(asscosts, xval)
        # obj0 = LinearAlgebra.dot([x.varcost for x in arcs], sum(xval,dims=2))
        obj0 = 0.0
        wval, obj0, obj, env, model, w, var1D2, xx, var2D2
    end

end
