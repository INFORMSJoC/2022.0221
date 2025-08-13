"""
TODO: deactivate any output on screen when print==false
TODO: nlsolvers: for all models, output nl file so that it can be solved with sc-minlp for example
TODO: nlsolvers: add sincos, sinabc1, sinabc2, sinabc2 in nlsolvers
TODO: improve the code for generating/reading nl functions
TODO: how to output/save the final solution ?
TODO: implement the naive knapsack pb for the comparison with sc-minlp (no need for absolute delta)
---> maybe this just corresponds to forcing Flim = 0.0 ? oe epsAbs = 0.0 ?
TODO: scheduling: add non linear function FC ??? to match the DAM2016 paper
TODO: write to nl file does write the file in the same folder as the original instance
but need to verify if that nl file can be solved by SC-MINLP
TODO: checker of input parameter
TODO: run mode, that deactivates checker and print

PB:   HEV problem E1400 takes too long to solve because of precision issues
"""


function solve_hev(Pdemprofilename, nlcost, params)
    if init < eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    pdata = get_hevdata(Pdemprofilename)
    nperiods = pdata.nbinstants
    pfmax = pdata.p_FCmax
    pfmin = pdata. p_FCmin

    f, g = construct_nonlinear_functions(pfmax, 1000, nlcost, :MIN)
    functions = [f for t in 1 : nperiods]
    exprs = [g for t in 1 : nperiods]
    report = ipwlb(exprs, [], [], functions, [], [], eps, :hev, pdata, params)
    report
end

function solve_knap(filename, are_vars_integer, params, 
                    Flim = 0.0)
    if params.eps0 < params.eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end

    pdata, exprs, functions = read_knapsacknlobj(filename, are_vars_integer, params.with_cost_type)
    # return functions
    
    report = ipwlb(exprs, [], [], functions, [], [], pdata, params)
    report
end

function solve_mcnd_or_mcf_from_file(graphfilepathname, commodityfilemathname, nodefilepathname, incremcost, 
                                    params,
                                    congestioncoef=1.0,
                                    restrict_to_commodities = 0,
                                    minlpsolver="no")
    if params.eps0 < params.eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    #outfiledat(nodes, arcs, commods, "example.txt", graphfilepathname, commodityfilemathname, nodefilepathname, incremcost)

    pdata, exprs, expr1, expr2, funcs, func1, func2 = load_mcnddata(graphfilepathname, commodityfilemathname,
                                                                            nodefilepathname, incremcost, params.ptype, with_cost_type = params.with_cost_type,
                                                                            restrict_to_commodities = restrict_to_commodities, printon=params.printon, congestioncoef=congestioncoef)

    #println("OKOKO")
    if params.ptype ∈ [:mcnd_nodes, :mcnd_arcs]
        if minlpsolver == "no"
            report = ipwlb(exprs, expr1, expr2, funcs, func1, func2, pdata, params)
        else
            xval, yval, wval, obj = minlp_mcnd(expr1, expr2, pdata, minlpsolver)
        end
    elseif params.ptype ∈ [:mcf_nodes, :mcf_arcs]
        report = ipwlb(exprs, expr1, expr2, funcs, func1, func2, pdata, params)
    else
        print("unknown pb type, should be mncd or mcf")
    end

    report
end


function solve_scheduling_from_file(filename, params)
    if init < eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    pdata = read_scheduling(filename)
    nperiods = maximum([tw[2] for tw in pdata.time_windows])
    ntasks = pdata.ntasks
    sumw = sum(pdata.demands)
    maxw = maximum(pdata.demands)
    avgw = sumw / ntasks
    println("periods = $nperiods")
    println("demands = $(sum(pdata.demands))")
    f, exp = construct_nonlinear_functions(sumw, 1000, with_cost_type, :MIN)
    exprs = [exp for t in 1 : nperiods]
    functions = [f for t in 1 : nperiods]

    report = ipwlb(exprs, [], [], functions, [], [], eps, :scheduling, pdata, params)
    report
end

function solve_uflp_from_file(filename,
                            params::Parameters,
                            with_cap_constraint = false,
                            scale_f = 1.0,
                            scale_ass = 1.0)
    if params.eps0 < params.eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    if occursin("holmberg", filename)
        println("reading holmberg data file $filename")
        pdata = read_uflp_holmberg(filename; with_cap_constraint = with_cap_constraint, scale_f = scale_f, scale_ass = scale_ass)
    elseif occursin("orlib", filename)
        println("reading orlib data file $filename")
        pdata = read_uflp_orlib(filename; with_cap_constraint = with_cap_constraint, scale_f = scale_f, scale_ass = scale_ass)
    elseif occursin("gunluk", filename)
        println("reading gunluk data file $filename")
        pdata = read_gunluk_file(filename)
    else
        println("unknown datafile!")
        return
    end
    ndeps = length(pdata.caps)
    ncusts = length(pdata.demands)
    costmult = 1
    if params.ptype == :uflp_nlass
        funcs = Array{Any, 2}(undef, ndeps, ncusts)
        exprs = Array{Expr, 2}(undef, ndeps, ncusts)
        for i in 1 : ndeps, j in 1 : ncusts
            scaled_cost = costmult * pdata.asscosts[i, j]
            funcs[i, j], exprs[i, j] = construct_nonlinear_functions(1.0, scaled_cost, params.with_cost_type, :MIN)
        end
    else
        funcs = Array{Any, 1}(undef, ndeps)
        exprs = Array{Expr, 1}(undef, ndeps)
        for i in 1 : ndeps
            if with_cap_constraint
                maxoffer = min(sum(pdata.demands), pdata.caps[i])
            else maxoffer = sum(pdata.demands)
            end
            scaled_cost = costmult * pdata.fixed_costs[i]
            funcs[i], exprs[i] = construct_nonlinear_functions(maxoffer, scaled_cost, params.with_cost_type, :MIN)
        end
        for exp in exprs
            println(exp)
        end
        # return
    end

    report = ipwlb(exprs, [], [], funcs, [], [], pdata, params)
end

function solve_transp_random(ndeps, ncusts, params)
    if params.eps0 < params.eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    pdata = generate_random_transp(ndeps, ncusts)
    report = solve_transp(pdata, params)
    pdata, report
end

function solve_transp_from_file(filename, params)
    if params.eps0 < params.eps
        throw(ArgumentError("Your init value $init must be greater than the final eps value $eps !!!!!!"))
    end
    pdata = read_transportation(filename)
    report = solve_transp(pdata, params)
end

function solve_transp(pdata::TranspProblemData, params::Parameters)
    ndeps = length(pdata.caps)
    ncusts = length(pdata.demands)
    funcs = Array{Any, 2}(undef, ndeps, ncusts)
    exprs = Array{Expr, 2}(undef, ndeps, ncusts)
    for i in 1 : ndeps, j in 1 : ncusts
        maxoffer = min(pdata.demands[j], pdata.caps[i])
        cMax = pdata.asscosts[i, j]
        funcs[i, j], exprs[i, j] = construct_nonlinear_functions(maxoffer, cMax, params.with_cost_type, :MIN)
    end
    report = ipwlb(exprs, [], [], funcs, [], [], pdata, params)
end
