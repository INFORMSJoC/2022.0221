function ipwlb(exprs,
                expr1,
                expr2,
                functions,
                func1,
                func2,
                pdata::Data,
                _params::Parameters)
                # with_solver = GLPK.Optimizer,
                # with_solver_parameters = empty_parameters,
                # with_lina_solver = LinA.HeuristicLin,
                # eps0 = 1e-1,
                # tilim = 3600,
                # printon=false,
                # stepdiv=2.0,
                # dynstep=false,
				# activeUpdateLargePwlf=false,
				# warmstartMILP=false,
				# miptype=:grb,
                # with_fixed_charge = false)
    params = deepcopy(_params)
    start_time = time_ns()
    nperiods = max(length(exprs), length(expr1), length(expr2))
    pwl = construct_linear_expressions(functions, 
                                        pdata,
                                        params)#SUN
    
	if params.miptype==:grb
		env = Gurobi.Env()
		model = Gurobi.Optimizer()
		var2D1=Array{MathOptInterface.VariableIndex, 2}(undef, 0, 0)
		var2D2=Array{MathOptInterface.VariableIndex, 2}(undef, 0, 0)
		var1D1=Array{MathOptInterface.VariableIndex, 1}(undef, 0)
		var1D2=Array{MathOptInterface.VariableIndex, 1}(undef, 0)
	else
		env = nothing
		model = nothing
		var2D1 = nothing
		var2D2 = nothing
		var1D1 = nothing
		var1D2 = nothing
	end

	#=
	println("no of pieces = $(length(pwl))")
    for (t, exp) in enumerate(exprs)
		println("expression $(exp)")
    end
    println("data = $(pdata)")
	=#

    lb, ub = -1e+20, 1e+20#0.0, 1e+20 #modified by SUN to handle the knapsack problem
    obj0ub = 0.0
    prev_wvals = Array{Array{Float64, 1}, 1}()
    prev_obj0 = Array{Float64, 1}()
    report = Array{Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Int64, Int64, Float64, Float64}, 1}()
    int_x = true
    nbattempts=0
    mipcpu=0.0


    elapsed = 0.0
    nit_main_loop = 0
    while true
        nit_main_loop += 1
        # if ptype == :mcnd && isassigned(func2, 1)
            # check_linearization_feasible(pwl, pdata.nodes)
        # end
		#println("SUN debug: nbattempts= $nbattempts")
        nbattempts=nbattempts+1
        linearFuncts = [[pwl[t][i].lp for i in 1 : length(pwl[t])] for t in 1 : nperiods]
        cutoff = compute_cutoff(prev_wvals, prev_obj0, pwl, pdata, params.ptype)
        num_pieces = sum(length(pwl[t]) for t in 1 : nperiods)
        elapsed = (time_ns() - start_time) * 1e-9
        timeLimit = params.tilim - elapsed
        println("Iteration $nit_main_loop\tLB = $lb\tUB = $ub\tCPU = $elapsed, #pieces = $num_pieces")
        if timeLimit < 0
            break
        end



		#print("SUN before mip")
        startmipcpu = time_ns()
        params.warmstartMILP = _params.warmstartMILP && (nbattempts >= 2)
        params.tilim = timeLimit
        wval, obj0, newlb, env, model, var1D1, var1D2, var2D1, var2D2 = 
            mip_piecewise(linearFuncts, pdata, env, model, var1D1, var1D2, var2D1, var2D2,
                            params,
                            lb,
                            ub)
                            
        endmipcpu = time_ns()
        mipcpu = mipcpu + (endmipcpu - startmipcpu) * 1e-9
        elapsed = (time_ns() - start_time) * 1e-9
        #print("SUN newlb",newlb)
        #println(obj0, " hh ", newlb-obj0)
        if isempty(wval)
            end_time = time_ns()
            gaptotal = (ub - lb) / abs(ub) * 100.0 #CCtocheck
            #gap = ((ub-obj0ub) - (lb-obj0)) / abs(ub-obj0ub) * 100#CCtocheck
            gap = (ub-lb) / abs(ub-obj0ub) * 100.0 #CCtocheck
            push!(report, (lb, ub, obj0, obj0ub, gap, gaptotal, num_pieces, nbattempts, mipcpu, (end_time - start_time) * 1e-9))#CCtocheck
            break
        else
            otherlb = sum(pwl[t](wval[t]).val for t in 1 : nperiods) + obj0
#            println("recomputed lower bound = $otherlb, original was $newlb")
            newlb = min(otherlb, newlb)
        end
        lb = max(lb, newlb)
        newub = compute_upperbound(params.ptype, pdata, functions, wval, obj0, params.with_fixed_charge)
        if newub < ub
            ub = newub
            obj0ub=obj0
        end
        # println("newlb = $newlb, obj0 = $obj0, newub = $newub, $ub")
        #gap = ((ub-obj0ub) - (lb-obj0)) / abs(ub-obj0ub) * 100
        gap = (newub-newlb) / abs(newub-obj0) * 100.0
        gaptotal = (ub - lb) / abs(ub) * 100.0
        if params.eps0 < 0 || abs(params.eps0 - params.eps) < 1e-12
            end_time = time_ns()
            return [(lb, ub, obj0, obj0ub, gap, gaptotal, num_pieces, nbattempts, mipcpu, (end_time - start_time) * 1e-9)]#CCtocheck
        else
            end_time = time_ns()
            push!(report, (lb, ub, obj0, obj0ub, gap, gaptotal, num_pieces, nbattempts, mipcpu, (end_time - start_time) * 1e-9))#CCtocheck
        end
        push!(prev_wvals, wval)
        push!(prev_obj0, obj0)
        # println("upper bound = $ub")
        if params.printon
            println("gap = $gap")
            println("gaptotal = $gaptotal")
        end
        if gaptotal < params.eps
            break
        end
		if params.ptype == :hev
			abs_errs, rel_errs, indices = compute_all_errors(wval, functions, pwl; minw = pdata.p_FCmin)
                else
                    minw = params.with_fixed_charge ? EPS_NZ : -EPS_NZ
                    abs_errs, rel_errs, indices = compute_all_errors(wval, functions, pwl; minw = minw)
		end
        if sum(indices) <= 0
            println("sum(indices) = $(sum(indices))")
            break
        end
        if ((maximum(rel_errs) < params.eps) && params.printon)
            println("maximum(errs) = $(maximum(rel_errs))")
            println("gaptotal =", gaptotal)
            println(rel_errs)
            println(indices)
            println(wval)
            #=println(pwl[3][4].lp.a*(wval[3])+pwl[3][4].lp.b, pwl[3][4].lp.xMin, pwl[3][4].lp.xMax)
            println(pwl[3][5].lp.a*(wval[3])+pwl[3][5].lp.b, pwl[3][5].lp.xMin, pwl[3][5].lp.xMax)
            println(functions[3](wval[3]))

            println(pwl[2][3].lp.a*(wval[2])+pwl[2][3].lp.b)
            println(functions[2](wval[2]))
            =#
            #println(pwl[4](30.4195876))
            break
        end
        updated = false
        if params.printon
            println("maximum error of $(maximum(rel_errs))")
        end

        multiply(u) = max(abs(compute(functions[u], wval[u], params.with_fixed_charge)), abs(pwl[u](wval[u]).val))
        score_fn(u) = 0.9 * rel_errs[u] * multiply(u) + 0.1 * abs_errs[u]
        
        sorted_periods = [t for t in 1 : nperiods if indices[t] > 0 && rel_errs[t] >= params.eps]
        sort!(sorted_periods, lt = (u, v) -> score_fn(u) > score_fn(v) )
        # if length(sorted_periods) > 1
        #     resize!(sorted_periods, 1)
        # end
        maxscore = 0.0
        let
            p = sorted_periods[begin]
            maxscore = score_fn(p)
        end
        # println.(sorted_periods)
        for t in sorted_periods
            score = score_fn(t)
            if score >= 0.1 * maxscore || !updated
                i = indices[t]
                pwlf = pwl[t][i]
                # println("pwlf.eps = $(pwlf.eps), error $(errs[t]) and score $score, xmin = $(pwlf.lp.xMin), xmax = $(pwlf.lp.xMax), delta = $(pwlf.lp.xMax - pwlf.lp.xMin)")
                # if pwlf.eps < 0.5 * eps continue #SUNtoCC: si ceci arrive c une erreur ?
                # if pwlf.lp.xMax - pwlf.lp.xMin <= EPS_NZ continue # gurobi cannot handle smaller pieces
                # end

                neweps = pwlf.eps / params.stepdiv
                if (params.dynstep)
                    neweps = min( params.eps0 / (params.stepdiv^nbattempts) , neweps)
                end

                while neweps >= 0.75 * rel_errs[t]
                    neweps /= params.stepdiv
                end #while loop uncommented by SUN. Maybe it's the reason my code was not efficient ?

                neweps = max(neweps, 0.75 * params.eps) #modifiied by SUN so that a too small neweps is never requested
                updated = true
                mindelta = max(10 * EPS_NZ, (pwl[t][end].lp.xMax - pwl[t][begin].lp.xMin) * 1e-5)

                if  ((params.ptype == :mcnd) || (params.ptype == :mcf)) && isassigned(expr2, t)
                    if pwl[t][i].lp.xMin >= pdata.nodes[t].croisement
                        # println("updating pwl at w = $(wval[t]) in interval $(pwl[t][i].lp.xMin), $(pwl[t][i].lp.xMax) with croisement $(pdata.nodes[t].croisement) using right side")
                        pwl[t] = update_pwlf(pwl[t], expr2[t], func2[t], i, wval[t], mindelta, neweps, 0.75 * eps, with_lina_solver; printon=printon, activeUpdateLargePwlf=activeUpdateLargePwlf)
                        #=if pwl[t][i].lp.xMin == pdata.nodes[t].croisement
                            println("\n ****** pour info, croisement atteint ******\n")
                        end
                        =#
                    else
                        # println("updating pwl at w = $(wval[t]) in interval $(pwl[t][i].lp.xMin), $(pwl[t][i].lp.xMax) with croisement $(pdata.nodes[t].croisement) using left side")
                        pwl[t] = update_pwlf(pwl[t], expr1[t], func1[t], i, wval[t], mindelta, neweps, 0.75 * eps, with_lina_solver; printon=printon, activeUpdateLargePwlf=activeUpdateLargePwlf)
                    end
                else
                    # println("updating pwl using a single function")
                    # println("expr1 = $expr1")
                    # println("expr2 = $expr2")
                    # println("exprs = $exprs")
                    # @assert(pwl[t] !== nothing, ("pwl is void $(pwl[t])"))
                    # println("pwl = $(pwl[t])")
                    x = update_pwlf(pwl[t], exprs[t], functions[t], i, wval[t], mindelta, neweps, 0.75 * params.eps, params)#params.with_lina_solver; printon=printon, activeUpdateLargePwlf=activeUpdateLargePwlf)
                    # println("x = $x")
                    pwl[t] = x#update_pwlf(pwl[t], exprs[t], functions[t], i, wval[t], mindelta, neweps, eps, with_lina_solver; printon=printon)
                    # check_lowerbound_pwl(pwl[t], functions[t])
                end
            end
        end
        if !updated
            println("not updated")
            break
        end
    end
    report
end

function check_linearization_feasible(pwl, nodes)
    println("entering check linearization feasible")
    nperiods = length(nodes)
    for t in 1 : nperiods
        croisement = nodes[t].croisement
        for pwlf in pwl[t]
            if pwlf.lp.xMin < croisement - 1e-9 && pwlf.lp.xMax > croisement + 1e-9
                println("infeasible linear piece at t = $t with crossing in $croisement and (xmin, xmax) = ($(pwlf.lp.xMin), $(pwlf.lp.xMax))")
            end
        end
    end
end
