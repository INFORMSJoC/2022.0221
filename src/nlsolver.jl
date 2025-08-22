"""
TODO: add sincost, sinabc1, sinabc2, sinabc costs in gams and nl functions

TODO: create nl file for sc-minlp

"""


# using JuMP, SCIP, Alpine, CPLEX, Ipopt, Juniper

function gams_mcnd(graphfilepathname, commodityfilemathname, nodefilepathname, incremcost, outfile; with_cost_type = :twofuncs, restrict_to_commodities = 0, congestioncoef=1.0)

    pdata, exprs, expr1, expr2, funcs, func1, func2 = load_mcnddata(graphfilepathname, commodityfilemathname,
                                                                        nodefilepathname, incremcost; with_cost_type=with_cost_type,
                                                                        restrict_to_commodities=restrict_to_commodities, printon=false, congestioncoef=congestioncoef)

    out = open(outfile, "w")
    println(out, "\$Title Multicommodity Network Design Problem with Congestion")
    println(out, "\$Ontext")
    println(out, "graph file: $graphfilepathname \n commodity file: $commodityfilemathname \n node file: $nodefilepathname")
    println(out, "cost function: $(string(with_cost_type))")
    println(out, "\n\$Offtext")
    println(out)
    println(out, "\tSets")
    println(out, "\t\tNODES set of nodes / 1 * $(pdata.nbNode) /")
    print(out, "\t\tARCS(NODES, NODES) set of arcs / $(pdata.arcs[1].from).$(pdata.arcs[1].to)")
    for i in 2:pdata.nbArc
        print(out, ", $(pdata.arcs[i].from).$(pdata.arcs[i].to)")
    end
    println(out, " /")
    println(out, "\t\tPRODS set of commodities / 1 * $(pdata.nbCommod) /")
    println(out, ";\n\n\tAlias(NODES, NODESIN, NODESOUT) ;\n")
    println(out)
    println(out, "\tParameters")
    print(out, "\t\tO_arcfixedcost(NODES,NODES) fixed costs of arcs\n\t\t/\t$(pdata.arcs[1].from).$(pdata.arcs[1].to)\t$(pdata.arcs[1].fixedcost)")
    for i in 2:pdata.nbArc
        print(out, "\n\t\t\t$(pdata.arcs[i].from).$(pdata.arcs[i].to)\t$(pdata.arcs[i].fixedcost)")
    end
    println(out, "\t/")
    print(out, "\t\tR_arcroutingcost(NODES,NODES) routing cost of commodity on arc\n\t\t/\t$(pdata.arcs[1].from).$(pdata.arcs[1].to)\t$(pdata.arcs[1].varcost)")
    for i in 2:pdata.nbArc
        print(out, "\n\t\t\t$(pdata.arcs[i].from).$(pdata.arcs[i].to)\t$(pdata.arcs[i].varcost)")
    end
    println(out, "\t/")
    print(out, "\t\tU_arccapacity(NODES,NODES) capacity of arcs\n\t\t/\t$(pdata.arcs[1].from).$(pdata.arcs[1].to)\t$(pdata.arcs[1].capa)")
    for i in 2:pdata.nbArc
        print(out, "\n\t\t\t$(pdata.arcs[i].from).$(pdata.arcs[i].to)\t$(pdata.arcs[i].capa)")
    end
    println(out, "\t/")

    println(out, "\t\tnodecapacity(NODES) capacity of nodes")
    print(out, "\t\t/\t1\t$(pdata.nodes[1].capainit+pdata.nodes[1].capaincrement)")
    for i in 2 : pdata.nbNode
        print(out, "\n\t\t\t$i\t$(pdata.nodes[i].capainit+pdata.nodes[i].capaincrement)")
    end
    println(out, "\t/")

    println(out, "\t\tnodesplitpoint(NODES) split point of the nonlinear congestion function of nodes")
    if with_cost_type == :twofuncs
        print(out, "\t\t/\t1\t$(pdata.nodes[1].croisement)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].croisement)")
        end
        println(out, "\t/")
    else
        print(out, "\t\t/\t1\t0.0")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t0.0")
        end
        println(out, "\t/")
    end


    println(out, "\t\tW_prodtocarry(PRODS) qty of products to transport from origin to destination")
    print(out, "\t\t/\t1\t$(pdata.commods[1].qty)")
    for i in 2 : pdata.nbCommod
        print(out, "\n\t\t\t$i\t$(pdata.commods[i].qty)")
    end
    println(out, "\t/")

    print(out, "\t\tD_proddemand(NODES,PRODS) value equal to W if origin nodes or -W if destination nodes or 0 otherwise")
    for i in 1:pdata.nbNode
        for k in 1:pdata.nbCommod
            dem=0
            if pdata.commods[k].from == i
                dem=pdata.commods[k].qty
            elseif pdata.commods[k].to == i
                dem=-pdata.commods[k].qty
            end
            if ((i==1) && (k ==1))
                print(out, "\n\t\t/\t$i.$k\t$dem")
            else
                print(out, "\n\t\t\t$i.$k\t$dem")
            end
        end
    end
    println(out, "\t/")

    if with_cost_type == :twofuncs
        println(out, "\t\tc_DiFi(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.nodes[1].delaycost*pdata.nodes[1].freeflowdelay)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].delaycost*pdata.nodes[i].freeflowdelay)")
        end
        println(out, "\t/")

        println(out, "\t\tc_DiBeta(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.nodes[1].delaycost)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].delaycost)")
        end
        println(out, "\t/")

        println(out, "\t\tc_Ei(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.nodes[1].fixedcost)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].fixedcost)")
        end
        println(out, "\t/")

        println(out, "\t\tc_Cizero(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.nodes[1].capainit)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].capainit)")
        end
        println(out, "\t/")

        println(out, "\t\tc_CizeroCiplus(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.nodes[1].capainit+pdata.nodes[1].capaincrement)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.nodes[i].capainit+pdata.nodes[i].capaincrement)")
        end
        println(out, "\t/")

    elseif with_cost_type in [:cubic, :exp]
        println(out, "\t\tc_a1(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.cubicexpcoefs[1].a1)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.cubicexpcoefs[i].a1)")
        end
        println(out, "\t/")

        println(out, "\t\tc_a2(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.cubicexpcoefs[1].a2)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.cubicexpcoefs[i].a2)")
        end
        println(out, "\t/")

        println(out, "\t\tc_a3(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.cubicexpcoefs[1].a3)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.cubicexpcoefs[i].a3)")
        end
        println(out, "\t/")

        println(out, "\t\tc_a4(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.cubicexpcoefs[1].a4)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.cubicexpcoefs[i].a4)")
        end
        println(out, "\t/")

        println(out, "\t\tc_a5(NODES) coefficient of nonlinear function")
        print(out, "\t\t/\t1\t$(pdata.cubicexpcoefs[1].a5)")
        for i in 2 : pdata.nbNode
            print(out, "\n\t\t\t$i\t$(pdata.cubicexpcoefs[i].a5)")
        end
        println(out, "\t/")
    else
        println("UNKNOWN cost type")
    end

    println(out, ";\n")


    println(out, "\tScalars")
    println(out, "\t\talpha alphadd \t/ 3 /\n\t\tbeta betadd \t/ 1 /\n;\n")
    println(out, "\tVariables")
    println(out, "\t\tTotalCost asdasd")
    println(out, "\t\tX_IJK(NODES,NODES,PRODS) quantity sent")
    println(out, "\t\tY_IJ(NODES,NODES) active or not")
    println(out, "\t\tV_I(NODES) quantity going through node")
    println(out, "\t\tS_I(NODES) choice to switch to upgrade or not")
    println(out, ";")
    println(out)
    println(out, "\tPositive Variable X_IJK;")
    println(out, "\tBinary Variable Y_IJ;")
    println(out, "\tPositive Variable V_I;")
    println(out, "\tBinary Variable S_I;")
    println(out, "\n\tV_I.up(NODES) = nodecapacity(NODES) ;")
    println(out)
    println(out, "\tEquations")
    println(out, "\t\tComputeCost define objective function")
    println(out, "\t\tFlowConservation flow")
    println(out, "\t\tNoFlowifNoArc ddd")
    println(out, "\t\tCapacityonArc ddd")
    println(out, "\t\tComputeNodeFlow ddd")
    println(out, "\t\tBeforeSplitPoint ddd")
    println(out, "\t\tAfterSplitPoint ddd")
    println(out, ";")
    println(out)
    if with_cost_type == :cubic
        #println(out, "\tComputeCost .. TotalCost =g= sum((ORIGIN, DESTINATION), a1(ORIGIN, DESTINATION) * power(X_IJ(ORIGIN, DESTINATION), 3) + a2(ORIGIN, DESTINATION) * power(X_IJ(ORIGIN, DESTINATION), 2) + a3(ORIGIN, DESTINATION) * X_IJ(ORIGIN, DESTINATION) + a4(ORIGIN, DESTINATION) * ACT_IJ(ORIGIN, DESTINATION));")
        println(out, "\tComputeCost .. TotalCost =e= sum((NODESIN,NODES)\$(ARCS(NODESIN,NODES)) , O_arcfixedcost(NODESIN,NODES) * Y_IJ(NODESIN,NODES))
+   sum((NODESIN,NODES,PRODS)\$(ARCS(NODESIN,NODES)) , R_arcroutingcost(NODESIN,NODES) * X_IJK(NODESIN,NODES,PRODS))
+   sum(NODES, c_a1(NODES)*power(V_I(NODES),3))
+   sum(NODES, c_a2(NODES)*power(V_I(NODES),2))
+   sum(NODES, c_a3(NODES)*V_I(NODES))
+   sum(NODES, c_a4(NODES)*S_I(NODES))
;")
    elseif with_cost_type == :exp
        println(out, "\tComputeCost .. TotalCost =e= sum((NODESIN,NODES)\$(ARCS(NODESIN,NODES)) , O_arcfixedcost(NODESIN,NODES) * Y_IJ(NODESIN,NODES))
+   sum((NODESIN,NODES,PRODS)\$(ARCS(NODESIN,NODES)) , R_arcroutingcost(NODESIN,NODES) * X_IJK(NODESIN,NODES,PRODS))
+   sum(NODES, c_a5(NODES)*S_I(NODES))
+   sum(NODES, c_a3(NODES)*S_I(NODES) / (1+c_a2(NODES)*exp(-c_a1(NODES)*(V_I(NODES)+c_a4(NODES)))) )
;")
    else
        println(out, "\tComputeCost .. TotalCost =e= sum((NODESIN,NODES)\$(ARCS(NODESIN,NODES)) , O_arcfixedcost(NODESIN,NODES) * Y_IJ(NODESIN,NODES))
+   sum((NODESIN,NODES,PRODS)\$(ARCS(NODESIN,NODES)) , R_arcroutingcost(NODESIN,NODES) * X_IJK(NODESIN,NODES,PRODS))
+   sum(NODES, c_DiFi(NODES)*V_I(NODES))
+   sum(NODES, c_DiBeta(NODES)*power(V_I(NODES)/c_Cizero(NODES),alpha)*V_I(NODES)*(1-S_I(NODES)))
+   sum(NODES,c_Ei(NODES)*S_I(NODES))
+   sum(NODES, c_DiBeta(NODES)*power(V_I(NODES)/c_CizeroCiplus(NODES),alpha)*V_I(NODES)*S_I(NODES))
;")
    end
    println(out, "\tFlowConservation(NODES,PRODS) ..
        sum(NODESOUT\$(ARCS(NODES,NODESOUT)), X_IJK(NODES,NODESOUT,PRODS))
        - sum(NODESIN\$(ARCS(NODESIN,NODES)), X_IJK(NODESIN,NODES,PRODS))
        =e= D_proddemand(NODES,PRODS)
;")
    println(out, "\tNoFlowifNoArc(NODESIN,NODES,PRODS)\$(ARCS(NODESIN,NODES)) ..
        X_IJK(NODESIN,NODES,PRODS)
        =l= W_prodtocarry(PRODS) * Y_IJ(NODESIN,NODES)
;")
    println(out, "\tCapacityonArc(NODESIN,NODES)\$(ARCS(NODESIN,NODES)) ..
        sum(PRODS, X_IJK(NODESIN,NODES,PRODS))
        =l=  U_arccapacity(NODESIN,NODES) * Y_IJ(NODESIN,NODES)
;")
    println(out, "\tComputeNodeFlow(NODES) ..
        sum((PRODS,NODESIN)\$(ARCS(NODESIN,NODES)),X_IJK(NODESIN,NODES,PRODS))
        =e= V_I(NODES)
;")
    println(out, "\tBeforeSplitPoint(NODES) ..
        V_I(NODES) =l= nodesplitpoint(NODES)
        + (nodecapacity(NODES) - nodesplitpoint(NODES)) * S_I(NODES)
;")
    println(out, "AfterSplitPoint(NODES) ..
        V_I(NODES) =g= nodesplitpoint(NODES)*S_I(NODES)
;")
    println(out)
    println(out, "\nModel cMCND /all/ ;\n\nSolve cMCND using minlp minimizing ToTalCost ;")
    close(out)
end


function gams_uflp_nldep(filename, cost, outfile = out.gms)
    out = open(outfile, "w")
    if occursin("holmberg", filename)
        pdata = IterativePWLB.read_uflp_holmberg(filename)
    elseif occursin("orlib", filename)
        pdata = IterativePWLB.read_uflp_orlib(filename)
    end
    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    println(out, "\$Title CFLP_nonlinear_warehousing")
    println(out, "\$Ontext")
    println(out, "instance file: $filename")
    println(out, "cost function: $(string(cost))")
    println(out, "\$Offtext")
    println(out)
    println(out, "\tSets")
    println(out, "\tDEPOT set of depots / 1 * $ndeps /")
    println(out, "\tCUST set of customers / 1 * $ncusts /")
    println(out, ";")
    println(out)
    println(out, "\tParameters")
    println(out, "\tD_demand(CUST) demand per customer")
    println(out, "\t/\t1\t$(pdata.demands[1])")
    for i in 2 : ncusts - 1
        println(out, "\t\t$i\t$(pdata.demands[i])")
    end
    println(out, "\t\t$ncusts\t$(pdata.demands[end])\t/")
    println(out)
    println(out, "\tc_IJ(DEPOT, CUST) assignment cost")
    println(out, "\t/\t1.1 $(pdata.asscosts[1, 1])")
    for i in 1 : ndeps, j in 1 : ncusts
        if i == 1 && j == 1 continue
        elseif i == ndeps && j == ncusts continue
        else
            println(out, "\t\t$i.$j $(pdata.asscosts[i, j])")
        end
    end
    println(out, "\t\t$ndeps.$ncusts $(pdata.asscosts[end])\t/")
    println(out)

    println(out, "\tcap(DEPOT) capacity")
    println(out, "\t/\t1 $(pdata.caps[1])")
    for i in 2 : ndeps - 1
        println(out, "\t\t$i $(pdata.caps[i])")
    end
    println(out, "\t\t$ndeps $(pdata.caps[end])\t/")
    println(out)

    costmult = 1
    if cost == :cubic
        cf = zeros(ndeps, 4)
        for i in 1 : ndeps
            maxoffer = min(sum(pdata.demands), pdata.caps[i])
            xf, yf = cubic_function_multipliers(maxoffer, costmult * pdata.fixed_costs[i])
            cf[i, :] = cubic_coefs(xf, yf)
        end
    else
        cf = zeros(ndeps, 5)
        for i in 1 : ndeps
            maxoffer = min(sum(pdata.demands), pdata.caps[i])
            xf, yf = exp_function_multipliers(maxoffer, costmult * pdata.fixed_costs[i])
            cf[i, :] = exp_coefs(xf, yf)
        end
    end

    ncoefs = size(cf, 2)
    for k in 1 : ncoefs
        println(out, "\ta$k(DEPOT) coeff $k in the objective")
        println(out, "\t/\t1\t$(cf[1, k])")
        for i in 2 : ndeps - 1
            println(out, "\t\t$i\t$(cf[i, k])")
        end
        println(out, "\t\t$ndeps\t$(cf[ndeps, k])\t/")
    end

    if cost == :exp
        f0 = zeros(ndeps)
        for i in 1 : ndeps
            f0[i] = cf[i, 5] + cf[i, 3] / (1 + cf[i, 2] * exp(-cf[i, 1] * cf[i, 4]))
        end
        println(out)
        println(out, "\tF0(DEPOT) constant term when s = 0")
        println(out, "\t/\t1\t$(f0[1])")
        for i in 2 : ndeps - 1
            println(out, "\t\t$i\t$(f0[i])")
        end
        println(out, "\t\t$ndeps\t$(f0[ndeps])\t/")
    end

    println(out, ";")
    println(out)
    println(out, "\tVariables")
    println(out, "\tTotalCost asdasd")
    println(out, "\tAssCost asdasd")
    println(out, "\tWHCost asdasd")
    println(out, "\tX_IJ(DEPOT, CUST) asdasd")
    println(out, "\tY_I(DEPOT) asdasd")
    println(out, "\tS_I(DEPOT) asdasd")
    println(out, ";")
    println(out)
    println(out, "\tPositive Variable X_IJ;")
    println(out, "\tPositive Variable S_I;")
    println(out, "\tBinary Variable Y_I;")
    println(out)
    println(out, "\tEquations")
    println(out, "\tComputeCost asdasd")
    println(out, "\tComputeAssCost asdasd")
    println(out, "\tComputeWarehousingCost asdasd")
    println(out, "\tAssignment asdasd")
    println(out, "\tLinking asdasd")
    println(out, "\tWarehousing asdasd")
    println(out, "\tCapacity asdasd")
    println(out, ";")
    println(out)
    println(out, "\tComputeAssCost .. AssCost =g= sum((DEPOT, CUST), c_IJ(DEPOT, CUST) * X_IJ(DEPOT, CUST));")
    if cost == :cubic
        println(out, "\tComputeWarehousingCost .. WHCost =g= sum(DEPOT, a1(DEPOT) * power(S_I(DEPOT), 3) + a2(DEPOT) * power(S_I(DEPOT), 2) + a3(DEPOT) * S_I(DEPOT) + a4(DEPOT) * Y_I(DEPOT));")
    else
        println(out, "\tComputeWarehousingCost .. WHCost =g= sum(DEPOT, F0(DEPOT) * (Y_I(DEPOT) - 1) + a5(DEPOT) + a3(DEPOT) / (1 + a2(DEPOT) * exp(-a1(DEPOT) * (S_I(DEPOT) + a4(DEPOT)))));")
    end
    println(out, "\tComputeCost .. TotalCost =g= AssCost + WHCost;")
    println(out, "\tAssignment(CUST) .. sum(DEPOT, X_IJ(DEPOT, CUST)) =g= 1;")
    println(out, "\tLinking(DEPOT, CUST) .. X_IJ(DEPOT, CUST) =l= Y_I(DEPOT);")
    println(out, "\tWarehousing(DEPOT) .. S_I(DEPOT) =g= sum(CUST, D_demand(CUST) * X_IJ(DEPOT, CUST));")
    println(out, "\tCapacity(DEPOT) .. S_I(DEPOT) =l= cap(DEPOT) * Y_I(DEPOT);")
    println(out)
    println(out, "\nModel FLnldep /all/ ;\n\nSolve FLnldep using minlp minimizing ToTalCost ;")
    #=
    println(out, "\toption resLim = 3600;
    Model FLnldep /all/ ;

    option optCR = 0.001;
    FLnldep.solPrint = 2;

    Solve FLnldep using minlp minimizing ToTalCost ;

    option decimals = 8;
    parameter LBCost  the lower bound of the cost;

    LBCost = FLnldep.objEst;

    display ToTalCost.l, LBCost ;")
    =#
    close(out)
end

function gams_uflp_nlass(filename, cost, outfile = out.gms)
    out = open(outfile, "w")
    if occursin("holmberg", filename)
        pdata = IterativePWLB.read_uflp_holmberg(filename)
    elseif occursin("orlib", filename)
        pdata = IterativePWLB.read_uflp_orlib(filename)
    end
    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    println(out, "\$Title CFLP_nonlinear_assignment")
    println(out, "\$Ontext")
    println(out, "instance file: $filename")
    println(out, "cost function: $(string(cost))")
    println(out, "\$Offtext")
    println(out)
    println(out, "\tSets")
    println(out, "\tDEPOT set of depots / 1 * $ndeps /")
    println(out, "\tCUST set of customers / 1 * $ncusts /")
    println(out, ";")
    println(out)
    println(out, "\tParameters")
    println(out, "\tD_demand(CUST) demand per customer")
    println(out, "\t/\t1\t$(pdata.demands[1])")
    for i in 2 : ncusts - 1
        println(out, "\t\t$i\t$(pdata.demands[i])")
    end
    println(out, "\t\t$ncusts\t$(pdata.demands[end])\t/")
    println(out)
    println(out, "\tf_I(DEPOT) fixed cost")
    println(out, "\t/\t1 $(pdata.fixed_costs[1])")
    for i in 2 : ndeps - 1
        println(out, "\t\t$i $(pdata.fixed_costs[i])")
    end
    println(out, "\t\t$ndeps $(pdata.fixed_costs[end])\t/")
    println(out)

    println(out, "\tcap(DEPOT) capacity")
    println(out, "\t/\t1 $(pdata.caps[1])")
    for i in 2 : ndeps - 1
        println(out, "\t\t$i $(pdata.caps[i])")
    end
    println(out, "\t\t$ndeps $(pdata.caps[end])\t/")
    println(out)

    costmult = 1
    if cost == :cubic
        cf = zeros(ndeps, ncusts, 4)
        for i in 1 : ndeps, j in 1 : ncusts
            maxoffer = 1.0
            xf, yf = cubic_function_multipliers(maxoffer, costmult * pdata.asscosts[i, j])
            cf[i, j, :] = cubic_coefs(xf, yf)
        end
    else
        cf = zeros(ndeps, ncusts, 5)
        for i in 1 : ndeps, j in 1 : ncusts
            maxoffer = 1.0
            xf, yf = exp_function_multipliers(maxoffer, costmult * pdata.asscosts[i, j])
            cf[i, j, :] = exp_coefs(xf, yf)
        end
    end

    ncoefs = size(cf, 3)
    for k in 1 : ncoefs
        println(out, "\ta$k(DEPOT, CUST) coeff $k in the objective")
        println(out, "\t/\t1.1\t$(cf[1, 1, k])")
        for i in 1 : ndeps, j in 1 : ncusts
            if i == 1 && j == 1 continue
            elseif i == ndeps && j == ncusts continue
            else
                println(out, "\t\t$i.$j\t$(cf[i, j, k])")
            end
        end
        println(out, "\t\t$ndeps.$ncusts\t$(cf[ndeps, ncusts, k])\t/")
    end

    if cost == :exp
        f0 = zeros(ndeps, ncusts)
        for i in 1 : ndeps, j in 1 : ncusts
            f0[i, j] = cf[i, j, 5] + cf[i, j, 3] / (1 + cf[i, j, 2] * exp(-cf[i, j, 1] * cf[i, j, 4]))
        end
        println(out)
        println(out, "\tF0(DEPOT, CUST) constant term when s = 0")
        println(out, "\t/\t1.1\t$(f0[1, 1])")
        for i in 1 : ndeps, j in 1 : ncusts
            if i == 1 && j == 1 continue
            elseif i == ndeps && j == ncusts continue
            else
                println(out, "\t\t$i.$j\t$(f0[i, j])")
            end
        end
        println(out, "\t\t$ndeps.$ncusts\t$(f0[ndeps, ncusts])\t/")
    end

    println(out, ";")
    println(out)
    println(out, "\tVariables")
    println(out, "\tTotalCost asdasd")
    println(out, "\tAssCost asdasd")
    println(out, "\tWHCost asdasd")
    println(out, "\tX_IJ(DEPOT, CUST) asdasd")
    println(out, "\tY_I(DEPOT) asdasd")
    println(out, "\tS_I(DEPOT) asdasd")
    println(out, "\tACT_IJ(DEPOT, CUST) activation variable to link with X_IJ")
    println(out, ";")
    println(out)
    println(out, "\tPositive Variable X_IJ;")
    println(out, "\tPositive Variable S_I;")
    println(out, "\tBinary Variable Y_I;")
    println(out, "\tBinary Variable ACT_IJ;")
    println(out)
    println(out, "\tEquations")
    println(out, "\tComputeCost asdasd")
    println(out, "\tComputeAssCost asdasd")
    println(out, "\tComputeWarehousingCost asdasd")
    println(out, "\tAssignment asdasd")
    println(out, "\tLinking asdasd")
    println(out, "\tActivationConstraint asdasd")
    println(out, "\tWarehousing asdasd")
    println(out, "\tCapacity asdasd")
    println(out, ";")
    println(out)
    println(out, "\tComputeWarehousingCost .. WHCost =g= sum(DEPOT, f_I(DEPOT) * Y_I(DEPOT));")
    if cost == :cubic
        println(out, "\tComputeAssCost .. AssCost =g= sum((DEPOT, CUST), a1(DEPOT, CUST) * power(X_IJ(DEPOT, CUST), 3) + a2(DEPOT, CUST) * power(X_IJ(DEPOT, CUST), 2) + a3(DEPOT, CUST) * X_IJ(DEPOT, CUST) + a4(DEPOT, CUST) * ACT_IJ(DEPOT, CUST));")
    else
        println(out, "\tComputeAssCost .. AssCost =g= sum((DEPOT, CUST), F0(DEPOT, CUST) * (ACT_IJ(DEPOT, CUST) - 1) + a5(DEPOT, CUST) + a3(DEPOT, CUST) / (1 + a2(DEPOT, CUST) * exp(-a1(DEPOT, CUST) * (X_IJ(DEPOT, CUST) + a4(DEPOT, CUST)))));")
    end
    println(out, "\tComputeCost .. TotalCost =g= AssCost + WHCost;")
    println(out, "\tAssignment(CUST) .. sum(DEPOT, X_IJ(DEPOT, CUST)) =g= 1;")
    println(out, "\tLinking(DEPOT, CUST) .. X_IJ(DEPOT, CUST) =l= Y_I(DEPOT);")
    println(out, "\tActivationConstraint(DEPOT, CUST) .. X_IJ(DEPOT, CUST) =l= ACT_IJ(DEPOT, CUST);")
    println(out, "\tWarehousing(DEPOT) .. S_I(DEPOT) =g= sum(CUST, D_demand(CUST) * X_IJ(DEPOT, CUST));")
    println(out, "\tCapacity(DEPOT) .. S_I(DEPOT) =l= cap(DEPOT) * Y_I(DEPOT);")
    println(out)
    println(out, "\nModel FLnlass /all/ ;\n\nSolve FLnlass using minlp minimizing ToTalCost ;")
    #=
    println(out, "\toption resLim = 3600;
    Model FLnlass /all/ ;

    option optCR = 0.001;
    FLnlass.solPrint = 2;

    Solve FLnlass using minlp minimizing ToTalCost ;

    option decimals = 8;
    parameter LBCost  the lower bound of the cost;

    LBCost = FLnlass.objEst;

    display ToTalCost.l, LBCost ;")
    =#
    close(out)
end

function gams_transp(filename, cost, outfile = out.gms)
    out = open(outfile, "w")
    pdata = IterativePWLB.read_transportation(filename)
    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    println(out, "\$Title Transportation Problem")
    println(out, "\$Ontext")
    println(out, "instance file: $filename")
    println(out, "cost function: $(string(cost))")
    println(out, "\$Offtext")
    println(out)
    println(out, "\tSets")
    println(out, "\tORIGIN set of origins / 1 * $ndeps /")
    println(out, "\tDESTINATION set of destinations / 1 * $ncusts /")
    println(out, ";")
    println(out)
    println(out, "\tParameters")
    println(out, "\tD_demand(DESTINATION) demand per destination")
    println(out, "\t/\t1\t$(pdata.demands[1])")
    for i in 2 : ncusts - 1
        println(out, "\t\t$i\t$(pdata.demands[i])")
    end
    println(out, "\t\t$ncusts\t$(pdata.demands[end])\t/")
    println(out)
    println(out, "\tO_offer(ORIGIN) offer per origin")
    println(out, "\t/\t1\t$(pdata.caps[1])")
    for i in 2 : ndeps - 1
        println(out, "\t\t$i\t$(pdata.caps[i])")
    end
    println(out, "\t\t$ncusts\t$(pdata.caps[end])\t/")
    println(out)

    costmult = 1
    if cost == :cubic
        cf = zeros(ndeps, ncusts, 4)
        for i in 1 : ndeps, j in 1 : ncusts
            maxoffer = min(pdata.caps[i], pdata.demands[j])
            xf, yf = cubic_function_multipliers(maxoffer, costmult * pdata.asscosts[i, j])
            cf[i, j, :] = cubic_coefs(xf, yf)
        end
    else
        cf = zeros(ndeps, ncusts, 5)
        for i in 1 : ndeps, j in 1 : ncusts
            maxoffer = min(pdata.caps[i], pdata.demands[j])
            xf, yf = exp_function_multipliers(maxoffer, costmult * pdata.asscosts[i, j])
            cf[i, j, :] = exp_coefs(xf, yf)
        end
    end

    ncoefs = size(cf, 3)
    for k in 1 : ncoefs
        println(out, "\ta$k(ORIGIN, DESTINATION) coeff $k in the objective")
        println(out, "\t/\t1.1\t$(cf[1, 1, k])")
        for i in 1 : ndeps, j in 1 : ncusts
            if i == 1 && j == 1 continue
            elseif i == ndeps && j == ncusts continue
            else
                println(out, "\t\t$i.$j\t$(cf[i, j, k])")
            end
        end
        println(out, "\t\t$ndeps.$ncusts\t$(cf[ndeps, ncusts, k])\t/")
    end

    if cost == :exp
        f0 = zeros(ndeps, ncusts)
        for i in 1 : ndeps, j in 1 : ncusts
            f0[i, j] = cf[i, j, 5] + cf[i, j, 3] / (1 + cf[i, j, 2] * exp(-cf[i, j, 1] * cf[i, j, 4]))
        end
        println(out)
        println(out, "\tF0(ORIGIN, DESTINATION) constant term when s = 0")
        println(out, "\t/\t1.1\t$(f0[1, 1])")
        for i in 1 : ndeps, j in 1 : ncusts
            if i == 1 && j == 1 continue
            elseif i == ndeps && j == ncusts continue
            else
                println(out, "\t\t$i.$j\t$(f0[i, j])")
            end
        end
        println(out, "\t\t$ndeps.$ncusts\t$(f0[ndeps, ncusts])\t/")
    end

    println(out, ";")
    println(out)
    println(out, "\tVariables")
    println(out, "\tTotalCost asdasd")
    println(out, "\tX_IJ(ORIGIN, DESTINATION) asdasd")
    println(out, "\tACT_IJ(ORIGIN, DESTINATION) activation variable to link with X_IJ")
    println(out, ";")
    println(out)
    println(out, "\tPositive Variable X_IJ;")
    println(out, "\tBinary Variable ACT_IJ;")
    println(out)
    println(out, "\tEquations")
    println(out, "\tComputeCost asdasd")
    println(out, "\tCapacityConstraint asdasd")
    println(out, "\tDemandConstraint asdasd")
    println(out, "\tActivationConstraint asdasd")
    println(out, ";")
    println(out)
    if cost == :cubic
        println(out, "\tComputeCost .. TotalCost =g= sum((ORIGIN, DESTINATION), a1(ORIGIN, DESTINATION) * power(X_IJ(ORIGIN, DESTINATION), 3) + a2(ORIGIN, DESTINATION) * power(X_IJ(ORIGIN, DESTINATION), 2) + a3(ORIGIN, DESTINATION) * X_IJ(ORIGIN, DESTINATION) + a4(ORIGIN, DESTINATION) * ACT_IJ(ORIGIN, DESTINATION));")
    else
        println(out, "\tComputeCost .. TotalCost =g= sum((ORIGIN, DESTINATION), F0(ORIGIN, DESTINATION) * (ACT_IJ(ORIGIN, DESTINATION) - 1) + a5(ORIGIN, DESTINATION) + a3(ORIGIN, DESTINATION) / (1 + a2(ORIGIN, DESTINATION) * exp(-a1(ORIGIN, DESTINATION) * (X_IJ(ORIGIN, DESTINATION) + a4(ORIGIN, DESTINATION)))));")
    end
    println(out, "\tDemandConstraint(DESTINATION) .. sum(ORIGIN, X_IJ(ORIGIN, DESTINATION)) =g= D_demand(DESTINATION);")
    println(out, "\tCapacityConstraint(ORIGIN) .. sum(DESTINATION, X_IJ(ORIGIN, DESTINATION)) =l= O_offer(ORIGIN);")
    println(out, "\tActivationConstraint(ORIGIN, DESTINATION) .. X_IJ(ORIGIN, DESTINATION) =l= min(D_demand(DESTINATION), O_offer(ORIGIN)) * ACT_IJ(ORIGIN, DESTINATION);")
    println(out)
    println(out, "\nModel TPnl /all/ ;\n\nSolve TPnl using minlp minimizing ToTalCost ;")
    #=
    println(out, "\toption resLim = 3600;
    Model FLnlass /all/ ;

    option optCR = 0.001;
    FLnlass.solPrint = 2;

    Solve FLnlass using minlp minimizing ToTalCost ;

    option decimals = 8;
    parameter LBCost  the lower bound of the cost;

    LBCost = FLnlass.objEst;

    display ToTalCost.l, LBCost ;")
    =#
    close(out)
end

function nl_uflp_nldep(filename, cost, solver, solver_params; solve = true, write_nl = false, with_cap_constraint = false)
    if occursin("holmberg", filename)
        println("reading holmberg data file $filename")
        pdata = read_uflp_holmberg(filename; with_cap_constraint = with_cap_constraint)
    elseif occursin("orlib", filename)
        println("reading orlib data file $filename")
        pdata = read_uflp_orlib(filename; with_cap_constraint = with_cap_constraint)
    elseif occursin("gunluk", filename)
        println("Generating Gunluk and Lee data file $filename")
        pdata = read_gunluk_file(filename)
    else
        println("unknown datafile!")
        return
    end
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))
    
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end


    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    @variable(model, y[1 : ndeps], Bin)
    maxf = [pdata.capcons ? pdata.caps[i] : sum(pdata.demands) for i in 1 : ndeps]
    @variable(model, 0 <= f[i in 1 : ndeps] <= maxf[i])
    @variable(model, 0 <= x[1 : ndeps, j in 1 : ncusts] <= 1)
    @constraint(model, ass[j in 1 : ncusts], sum(x[i, j] for i in 1 : ndeps) == 1)
    if pdata.capcons
        @constraint(model, cap[i in 1 : ndeps], f[i] - pdata.caps[i] * y[i] <= 0)
    end
    @constraint(model, link[i in 1 : ndeps, j in 1 : ncusts], x[i, j] - y[i] <= 0)
    @constraint(model, wh[i in 1 : ndeps], sum(pdata.demands[j] * x[i, j] for j in 1 : ncusts) - f[i] == 0)
    
    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    costmult = 1
    cf = zeros(ndeps, dim)
    for i in 1 : ndeps
        if with_cap_constraint
            maxoffer = pdata.caps[i]
        else maxoffer = sum(pdata.demands)
        end
        scaled_cost = costmult * pdata.fixed_costs[i]
        cf[i, :] = get_nonlinear_coefficients(1.0 * maxoffer, scaled_cost, cost)
    end
    
    @variable(model, -1e+10 <= nlobj[1 : ndeps] <= 1e+10)

    for i in 1 : ndeps
        if solver isa Gurobi.Optimizer
            nlexp, rhs = get_nlgurobi_constraint_expression(model, cf[i, :], f[i], y[i], cost, :MIN)
            @constraint(model, nlexp - nlobj[i] <= rhs)    
        else
            nlexp, rhs = get_nonlinear_constraint_expression(model, cf[i, :], f[i], y[i], cost, :MIN)
            @NLconstraint(model, nlexp - nlobj[i] <= rhs)
        end
    end
    @objective(model, Min, sum(nlobj) + sum(pdata.asscosts[i, j] * x[i, j] for i in 1 : ndeps, j in 1 : ncusts))
    if write_nl
        println("writing file $(filename)-$(cost).nl")
        JuMP.write_to_file(model, "$(filename)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_uflp_nlass(filename, cost, solver, solver_params; solve = true, write_nl = false, with_cap_constraint = false)
    if occursin("holmberg", filename)
        println("reading holmberg data file $filename")
        pdata = read_uflp_holmberg(filename; with_cap_constraint = with_cap_constraint)
    elseif occursin("orlib", filename)
        println("reading orlib data file $filename")
        pdata = read_uflp_orlib(filename; with_cap_constraint = with_cap_constraint)
    elseif occursin("gunluk", filename)
        println("Generating Gunluk and Lee data file $filename")
        pdata = read_gunluk_file(filename)
    else
        println("unknown datafile!")
        return
    end
    
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))
    
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end

    
    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    @variable(model, y[1 : ndeps], Bin)
    @variable(model, 0 <= x[1 : ndeps, 1 : ncusts] <= 1)
    @variable(model, w[1 : ndeps, 1 : ncusts], Bin)
    @constraint(model, ass[j in 1 : ncusts], sum(x[i, j] for i in 1 : ndeps) == 1)
    @constraint(model, link[i in 1 : ndeps, j in 1 : ncusts], w[i, j] - y[i] <= 0)
    @constraint(model, link2[i in 1 : ndeps, j in 1 : ncusts], x[i, j] - w[i, j] <= 0)
    # @constraint(model, link[i in 1 : ndeps, j in 1 : ncusts], x[i, j] - y[i] <= 0)
    
    if pdata.capcons
        @constraint(model, cap[i in 1 : ndeps], sum(pdata.demands[j] * x[i, j] for j in 1 : ncusts) - pdata.caps[i] * y[i] <= 0)
    end
    
    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    costmult = 1
    cf = zeros(ndeps, ncusts, dim)
    for i in 1 : ndeps, j in 1 : ncusts
        scaled_cost = costmult * pdata.asscosts[i, j]
        cf[i, j, :] = get_nonlinear_coefficients(1.0, 1.0 * scaled_cost, cost)
    end
    
    # @variable(model, -1e+30 <= z <= 1e+30)
    @variable(model, -1e+10 <= nlobj[1 : ndeps, 1 : ncusts] <= 1e+10)
    for i in 1 : ndeps, j in 1 : ncusts
        if solver isa Gurobi.Optimizer
            nlexp, rhs = get_nlgurobi_constraint_expression(model, cf[i, j, :], x[i, j], w[i, j], cost, :MIN)
            @constraint(model, nlexp - nlobj[i, j] <= rhs)
        else
            nlexp, rhs = get_nonlinear_constraint_expression(model, cf[i, j, :], x[i, j], w[i, j], cost, :MIN)
            @NLconstraint(model, nlexp - nlobj[i, j] <= rhs)
        end
    end
    @objective(model, Min, sum(nlobj) + dot(pdata.fixed_costs, y))
    if write_nl
        JuMP.write_to_file(model, "$(filename)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_transp(filename, cost, solver, solver_params; solve = true, write_nl = false)
    pdata = IterativePWLB.read_transportation(filename)
    
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))
    
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end


    ndeps = size(pdata.asscosts, 1)
    ncusts = size(pdata.asscosts, 2)
    @variable(model, 0 <= x[i in 1 : ndeps, j in 1 : ncusts] <= min(pdata.caps[i], pdata.demands[j]))
    @variable(model, y[1 : ndeps, 1 : ncusts], Bin)
    # @variable(model, z[1 : ndeps, 1 : ncusts] >= 0)

    # @constraint(model, linkyz[i in 1 : ndeps, j in 1 : ncusts], z[i, j] + y[i, j] == 1)
    @constraint(model, outflow[i in 1 : ndeps], sum(x[i, j] for j in 1 : ncusts) <= pdata.caps[i])
    @constraint(model, inflow[j in 1 : ncusts], sum(x[i, j] for i in 1 : ndeps) >= pdata.demands[j])
    @constraint(model, link[i in 1 : ndeps, j in 1 : ncusts], x[i, j] - min(pdata.demands[j], pdata.caps[i]) * y[i, j] <= 0)
    
    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    cf = zeros(ndeps, ncusts, dim)
    for i in 1 : ndeps, j in 1 : ncusts
        maxoffer = 1.0 * min(pdata.demands[j], pdata.caps[i])
        maxcost = 1.0 * pdata.asscosts[i, j]
        cf[i, j, :] = get_nonlinear_coefficients(maxoffer, maxcost, cost)
    end
    
    @variable(model, -1e+10 <= nlobj[1 : ndeps, 1 : ncusts] <= 1e+10)
    
    for i in 1 : ndeps, j in 1 : ncusts
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, j, :], x[i, j], y[i, j], cost, :MIN)
            @constraint(model, exp - nlobj[i, j] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, j, :], x[i, j], y[i, j], cost, :MIN)
            @NLconstraint(model, exp - nlobj[i, j] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj))
    if write_nl
        JuMP.write_to_file(model, "$(filename)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
        # for i in 1 : ndeps, j in 1 : ncusts
        #     xv = value(x[i, j])
        #     yv = round(Int64, value(y[i, j]))
        #     if abs(xv) > EPS_NZ
        #         f, _ = construct_nonlinear_functions(1.0 * min(pdata.demands[j], pdata.caps[i]), 1.0 * pdata.asscosts[i, j], cost, :MIN)
        #         println("x[$i, $j] = $(xv), f(x) = $(f(xv))")
        #     end
        #     if yv < 1
        #         f, _ = construct_nonlinear_functions(1.0 * min(pdata.demands[j], pdata.caps[i]), 1.0 * pdata.asscosts[i, j], cost, :MIN)
        #         println("x[$i, $j] = $(xv), f(x) = $(f(xv))")
        #     end
        # end
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_scheduling(filename, cost, eps, solver, solver_params; solve = true, write_nl = false)

    pdata = IterativePWLB.read_scheduling(filename)
    # model = JuMP.Model(optimizer_with_attributes(SCIP.Optimizer, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))#, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end


    ntasks = pdata.ntasks
    durations = pdata.durations
    demands = pdata.demands
    tw = pdata.time_windows
    nperiods = maximum([t[2] for t in tw])

    @variable(model, 0 <= x[i in 1 : ntasks, t in 1 : nperiods] <= durations[i])
    @variable(model, 0 <= w[1 : nperiods] <= sum(demands))
    @variable(model, y[1 : nperiods], Bin)

    @constraint(model, [i in 1 : ntasks], sum(x[i, t] for t in 1 : nperiods if tw[i][1] <= t && t < tw[i][2]) >= durations[i])
    @constraint(model, [t in 1 : nperiods], sum(demands[i] * x[i, t] for i in 1 : ntasks) - w[t] <= 0)
    @constraint(model, [t in 1 : nperiods], w[t] - sum(demands) * y[t] <= 0)

    for t in 1 : nperiods - 1
        feast = [i for i in 1 : ntasks if tw[i][1] <= t && t < tw[i][2]]
        feasnext = [i for i in 1 : ntasks if tw[i][1] <= t + 1 && t + 1 < tw[i][2]]
        if feast == feasnext
            println("equivalent periods $t and $(t + 1)")
            @constraint(model, sum(x[i, t] - x[i, t + 1] for i in feast) >= 0)
        end
    end

    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    cf = zeros(nperiods, dim)
    for i in 1 : nperiods
        maxoffer = sum(demands)
        maxcost = 1000
        cf[i, :] = get_nonlinear_coefficients(maxoffer, maxcost, cost)
    end
    
    @variable(model, -1e+10 <= nlobj[1 : nperiods] <= 1e+10)
    for i in 1 : nperiods
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, :], w[i], y[i], cost, :MIN)
            @constraint(model, exp - nlobj[i] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, :], w[i], y[i], cost, :MIN)
            @NLconstraint(model, exp - nlobj[i] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj))
    
    if write_nl
        JuMP.write_to_file(model, "$(filename)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_mcnd(ptype, graph, commodities, nodes, incr, cost, solver, solver_params; solve = true, write_nl = false)
    rawdata = IterativePWLB.read_cMCND(graph, commodities, nodes, incr, 0)
    # return rawdata
    cubicexpcoefs=NamedTuple{(:a1, :a2, :a3, :a4, :a5),Tuple{Float64, Float64,Float64,Float64,Float64}}[]
    pdata = MCNDData(rawdata[1],
                    rawdata[2],
                    rawdata[3],
                    rawdata[4],
                    rawdata[5],
                    length(rawdata[1]),
                    length(rawdata[2]),
                    length(rawdata[3]),
                    0.0,
                    cubicexpcoefs)
    # return pdata
    # model = JuMP.Model(optimizer_with_attributes(SCIP.Optimizer, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))#, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end

    #=
    pdata =
        nodes::Array{NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement),Tuple{Int64, Float64,Int64,Int64,Float64,Float64,Float64}},1}
        arcs::Array{NamedTuple{(:ref, :from, :to, :fixedcost, :varcost, :capa),Tuple{Int64,Int64,Int64,Float64,Float64, Int64}},1}
        commods::Array{NamedTuple{(:ref, :from, :to, :qty),Tuple{Int64,Int64,Int64,Int64}},1}
        listofinputarcs
        listofoutputarcs
        nbNode::Int64
        nbArc::Int64
        nbCommod::Int64
        alpha::Float64
    =#
    narcs = pdata.nbArc
    nnodes = pdata.nbNode
    ncomm = pdata.nbCommod
    commods = pdata.commods
    arcs = pdata.arcs
    nodes = pdata.nodes

    maxoffer = [v.capainit + v.capaincrement for v in pdata.nodes]

    @variable(model, y[1 : narcs], Bin)
    @variable(model, 0 <= x[a in 1 : narcs, p in 1 : ncomm] <= min(commods[p].qty, arcs[a].capa))
    @variable(model, 0 <= w[a in 1 : narcs] <= arcs[a].capa)
    @variable(model, u[1 : nnodes], Bin)
    @variable(model, 0 <= v[i in 1 : nnodes] <= maxoffer[i])

    let
        rhs = zeros(nnodes, ncomm)
        for p in 1 : ncomm
            to = commods[p].to
            from = commods[p].from
            rhs[from, p] += commods[p].qty
            rhs[to, p] -= commods[p].qty
        end
        @constraint(model, [i in 1 : nnodes, p in 1 : ncomm], sum(x[a, p] for a in pdata.listofoutputarcs[i]) - sum(x[a, p] for a in pdata.listofinputarcs[i]) == rhs[i, p])
    end

    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    if ptype == "mcnd_nodes"
        nperiods = nnodes
        maxcost = Float64[nodes[i].fixedcost for i in 1:nnodes]
        maxx = Float64[maxoffer[i] for i in 1:nnodes]
        wvar_vec = v
        fc_vec = u
    elseif ptype == "mcnd_arcs"
        nperiods = narcs
        maxcost = Float64[arcs[i].fixedcost + arcs[i].varcost * arcs[i].capa for i in 1:narcs]
        maxx = Float64[arcs[i].capa for i in 1:narcs]
        wvar_vec = w
        fc_vec = y
    end

    cf = zeros(nperiods, dim)
    for i in 1 : nperiods
        cf[i, :] = get_nonlinear_coefficients(maxx[i], maxcost[i], cost)
    end
    
    @constraint(model, [a in 1 : narcs], sum(x[a, p] for p in 1 : ncomm) - w[a] <= 0)
    @constraint(model, [a in 1 : narcs], w[a] - arcs[a].capa * y[a] <= 0)
    @constraint(model, [a in 1 : narcs, p in 1 : ncomm], x[a, p] - min(commods[p].qty, arcs[a].capa) * y[a] <= 0)
    @constraint(model, [a in 1 : narcs], sum(x[a, p] for p in 1 : ncomm) - arcs[a].capa * y[a] <= 0)
    @constraint(model, [i in 1 : nnodes], sum(x[a, p] for a in pdata.listofinputarcs[i] for p in 1 : ncomm) - v[i] <= 0)
    @constraint(model, [i in 1 : nnodes], v[i] - maxoffer[i] * u[i] <= 0)

    @variable(model, -1e+10 <= nlobj[1 : nperiods] <= 1e+10)
    
    for i in 1 : nperiods
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, :], wvar_vec[i], fc_vec[i], cost, :MIN)
            @constraint(model, exp - nlobj[i] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, :], wvar_vec[i], fc_vec[i], cost, :MIN)
            @NLconstraint(model, exp - nlobj[i] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj) + LinearAlgebra.dot([x.fixedcost for x in arcs], y))
    
    if write_nl
    #    JuMP.write_to_file(model, "$(graph)-$(commodities)-$(nodes)-$(cost).nl")
	#JuMP.write_to_file(model, "baron.nl")
	print(model)
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_mcf(ptype, graph, commodities, nodes, incr, cost, solver, solver_params; restrict_to_commodities = 0, solve = true, write_nl = false)
    rawdata = IterativePWLB.read_cMCND(graph, commodities, nodes, incr, 0)
    cubicexpcoefs=NamedTuple{(:a1, :a2, :a3, :a4, :a5),Tuple{Float64, Float64,Float64,Float64,Float64}}[]
    # return rawdata
    pdata = MCNDData(rawdata[1],
                    rawdata[2],
                    rawdata[3],
                    rawdata[4],
                    rawdata[5],
                    length(rawdata[1]),
                    length(rawdata[2]),
                    length(rawdata[3]),
                    0.0,
                    cubicexpcoefs)
    # return pdata
    # model = JuMP.Model(optimizer_with_attributes(SCIP.Optimizer, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    # model = JuMP.Model(optimizer_with_attributes(solver,
    #                 collect(solver_params)...))#, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end

    #=
    pdata =
        nodes::Array{NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement),Tuple{Int64, Float64,Int64,Int64,Float64,Float64,Float64}},1}
        arcs::Array{NamedTuple{(:ref, :from, :to, :fixedcost, :varcost, :capa),Tuple{Int64,Int64,Int64,Float64,Float64, Int64}},1}
        commods::Array{NamedTuple{(:ref, :from, :to, :qty),Tuple{Int64,Int64,Int64,Int64}},1}
        listofinputarcs
        listofoutputarcs
        nbNode::Int64
        nbArc::Int64
        nbCommod::Int64
        alpha::Float64
    =#
    narcs = pdata.nbArc
    nnodes = pdata.nbNode
    ncomm = pdata.nbCommod
    commods = pdata.commods
    arcs = pdata.arcs
    nodes = pdata.nodes
    if restrict_to_commodities > 0
        ncomm = min(ncomm, restrict_to_commodities)
    end

    maxoffer = [v.capainit + v.capaincrement for v in pdata.nodes]

    @variable(model, 0 <= x[a in 1 : narcs, p in 1 : ncomm] <= min(commods[p].qty, arcs[a].capa))
    @variable(model, 0 <= w[a in 1 : narcs] <= arcs[a].capa)
    @variable(model, 0 <= won[a in 1 : narcs] <= 1)
    @variable(model, u[1 : nnodes], Bin)
    @variable(model, 0 <= v[i in 1 : nnodes] <= maxoffer[i])

    let
        rhs = zeros(nnodes, ncomm)
        for p in 1 : ncomm
            to = commods[p].to
            from = commods[p].from
            rhs[from, p] += commods[p].qty
            rhs[to, p] -= commods[p].qty
        end
        @constraint(model, [i in 1 : nnodes, p in 1 : ncomm], sum(x[a, p] for a in pdata.listofoutputarcs[i]) - sum(x[a, p] for a in pdata.listofinputarcs[i]) == rhs[i, p])
    end

    # @constraint(model, [a in 1 : narcs, p in 1 : ncomm], x[a, p] <= min(commods[p].qty, arcs[a].capa))
    @constraint(model, [a in 1 : narcs], sum(x[a, p] for p in 1 : ncomm) - w[a] == 0)
    @constraint(model, [a in 1 : narcs], w[a] - arcs[a].capa * won[a] <= 0)
    @constraint(model, [i in 1 : nnodes], sum(x[a, p] for a in pdata.listofinputarcs[i] for p in 1 : ncomm) - v[i] <= 0)
    @constraint(model, [i in 1 : nnodes], v[i] - maxoffer[i] * u[i] <= 0)

    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    if ptype == "mcf_nodes"
        nperiods = nnodes
        maxcost = Float64[nodes[i].fixedcost for i in 1:nnodes]
        maxx = Float64[maxoffer[i] for i in 1:nnodes]
        wvar_vec = v
        fc_vec = u
    elseif ptype == "mcf_arcs"
        nperiods = narcs
        maxcost = Float64[arcs[i].fixedcost + arcs[i].varcost * arcs[i].capa for i in 1:narcs]
        maxx = Float64[arcs[i].capa for i in 1:narcs]
        wvar_vec = w
        fc_vec = won
    end

    cf = zeros(nperiods, dim)
    for i in 1 : nperiods
        cf[i, :] = get_nonlinear_coefficients(maxx[i], maxcost[i], cost)
    end
    
    @variable(model, -1e+10 <= nlobj[1 : nperiods] <= 1e+10)
    
    for i in 1 : nperiods
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, :], wvar_vec[i], fc_vec[i], cost, :MIN)
            @constraint(model, exp - nlobj[i] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, :], wvar_vec[i], fc_vec[i], cost, :MIN)
            @NLconstraint(model, exp - nlobj[i] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj))
    
    if write_nl
        JuMP.write_to_file(model, "$(graph)-$(commodities)-$(nodes)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end

function nl_knap(filename, ctype, solver, solver_params; solve = true, write_nl = false)
    pdata, _, _ = read_knapsacknlobj(filename, false, ctype)

    nknap=length(pdata.knapsack_cap)
    nobj=length(pdata.object_U)
    object_coefc=pdata.object_coefc
    object_coefb=pdata.object_coefb
    minuslistea=pdata.minuslistea
    object_coefd=pdata.object_coefd

    # model = JuMP.Model(optimizer_with_attributes(solver, collect(solver_params)...))
    
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    #set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
	if solver_params isa Vector{String}
    		set_attributes(model, solver_params)
	else set_attributes(model, collect(solver_params)...)
	end
    end
    
    @variable(model, 0.0 <= x[i in 1 : nknap, j in 1 : nobj] <= pdata.object_U[j])

    # if cost in [:cubic, :A1, :A2, :R]
        @variable(model, onoff[1 : nknap, p in 1 : nobj], Bin)
        @constraint(model, [i in 1 : nknap, j in 1 : nobj], x[i, j] <= pdata.object_U[j]*onoff[i, j])
    # end

    @constraint(model, [i in 1 : nknap], sum(pdata.object_weight[j] * x[i, j] for j in 1 : nobj) <= pdata.knapsack_cap[i])
    @constraint(model, [j in 1 : nobj], sum(x[i, j] for i in 1 : nknap) <= pdata.object_U[j])

    if ctype == :exp
        dim = 5
    elseif ctype in [:cubic, :divexp]
        dim = 4
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif ctype in [:sincos, :sqrt]
        dim = 2
    elseif ctype == :square
        dim = 1
    end
    
    cf = zeros(nknap, nobj, dim)
    for i in 1 : nknap, j in 1 : nobj
        if ctype == :divexp
            cf[i, j, :] = [minuslistea[j], object_coefb[j], object_coefc[j], object_coefd[j]]
        else
            fkc = x -> ((object_coefc[j]))/(1+(object_coefb[j])*(exp((minuslistea[j])*(x+(object_coefd[j])))))
            maxcost = fkc(pdata.object_U[j])
            maxoffer = pdata.object_U[j]
            cf[i, j, :] = get_nonlinear_coefficients(maxoffer, maxcost, ctype)
        end
    end
    
    @variable(model, -1e+10 <= nlobj[1 : nknap, 1 : nobj] <= 1e+10)
    
    for i in 1 : nknap, j in 1 : nobj
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, j, :], x[i, j], onoff[i, j], ctype, :MAX)
            @constraint(model, exp - nlobj[i, j] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, j, :], x[i, j], onoff[i, j], ctype, :MAX)
            @NLconstraint(model, exp - nlobj[i, j] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj))
    
#    println(model)
    t0 = time_ns()
    if write_nl
        JuMP.write_to_file(model, "$(filename)-$(ctype).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = 0.0#objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    return dualbound, primalbound, elapsed
end


function nl_hev(filename, cost, eps, solver, solver_params; solve = true, write_nl = false)
    #pdata = IterativePWLB.get_hevdata(filename)
    pdata = get_hevdata(filename)

    nperiods = pdata.nbinstants
    pemin = pdata.p_ESmin
    pemax = pdata.p_ESmax
    pfmin = pdata.p_FCmin
    pfmax = pdata.p_FCmax
    pdem = pdata.Pdem
    outloss = pdata.losscoefSEout
    inloss = pdata.losscoefSEin


    # model = JuMP.Model(optimizer_with_attributes(SCIP.Optimizer, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    # model = JuMP.Model(optimizer_with_attributes(solver,
                    # collect(solver_params)...))#, "limits/gap" => eps, "numerics/feastol" => 1e-6))
    if solver isa Gurobi.Optimizer
        model = direct_model(solver)
    else
        model = JuMP.Model(() -> solver)
    end
    # set_attributes(model, collect(solver_params)...)
    if !isempty(solver_params)
        if solver_params isa Vector{String}
                set_attributes(model, solver_params)
        else set_attributes(model, collect(solver_params)...)
	end
    end

    
    @variable(model, pfmin <= x1[1 : nperiods] <= pfmax)
    @variable(model, y1[1 : nperiods], Bin)
    @variable(model, 0 <= x2[1 : nperiods] <= pemax)
    @variable(model, pemin <= x3[1 : nperiods] <= 0)
    # @variable(model, z[1 : nperiods])

    # @constraint(model, [i in 1 : nperiods], z[i] - x1[i] + y1[i] * pfmin == 0)

    @constraint(model, [i in 1 : nperiods], x1[i] - pfmin * y1[i] >= 0)
    @constraint(model, [i in 1 : nperiods], x1[i] - pfmax * y1[i] <= 0)
    @constraint(model, [i in 1 : nperiods], x1[i] + x2[i] + x3[i] >= pdem[i])
    @constraint(model, sum((1 + outloss) * x2[i] + (1 + inloss) * x3[i] for i in 1 : nperiods) <= 0)
    @constraint(model, [t in 1 : nperiods], sum((1 + outloss) * x2[i] + (1 + inloss) * x3[i] for i in 1 : t) <= pdata.ESinit - pdata.ESmin)
    @constraint(model, [t in 1 : nperiods], sum((1 + outloss) * x2[i] + (1 + inloss) * x3[i] for i in 1 : t) >= pdata.ESinit - pdata.ESmax)

    if cost == :exp
        dim = 5
    elseif cost in [:cubic, :pardivexp]
        dim = 4
    elseif cost in [:sinabc1, :sinabc2, :sinabc3, :log]
        dim = 3
    elseif cost in [:sincos, :sqrt]
        dim = 2
    elseif cost == :square
        dim = 1
    end
    
    cf = zeros(nnodes, dim)
    for i in 1 : nperiods
        maxoffer = pfmax
        maxcost = 1000
        cf[i, :] = get_nonlinear_coefficients(maxoffer, maxcost, cost)
    end
    
    @variable(model, -1e+10 <= nlobj[1 : nperiods] <= 1e+10)
    
    for i in 1 : nperiods
        if solver isa Gurobi.Optimizer
            exp, rhs = get_nlgurobi_constraint_expression(model, cf[i, :], x1[i], y1[i], cost, :MIN)
            @constraint(model, exp - nlobj[i] <= rhs)
        else
            exp, rhs = get_nonlinear_constraint_expression(model, cf[i, :], x1[i], y1[i], cost, :MIN)
            @NLconstraint(model, exp - nlobj[i] <= rhs)
        end
    end
    
    @objective(model, Min, sum(nlobj))
    
    if write_nl
        JuMP.write_to_file(model, "$(filename)-$(cost).nl")
    end
    dualbound = primalbound = gap = 0.0
    elapsed = 0.0
    if solve
        t0 = time_ns()
        optimize!(model)
        elapsed = ceil(100.0 * (time_ns() - t0) * 1e-9) / 100.0
        dualbound = objective_bound(model)
        primalbound = objective_value(model)
        gap = abs(primalbound - dualbound) / max(abs(primalbound), abs(dualbound)) * 100
    end
    println("primal_bound = $primalbound, dual_bound = $dualbound, gap = $gap, time = $elapsed")
    dualbound, primalbound, elapsed
end
