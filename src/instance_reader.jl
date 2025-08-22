"""
TODO : add sincos and sinabc1,2,3 cost functions to mcnd

add sincos and sinabc1,2,3 cost functions to knapsacknl

"""


function get_hevdata(Pdemprofilename)
    Pdem=[]
    nbinstants = 0
    open(Pdemprofilename) do f
        tok = split(readline(f))
        nbinstants = parse(Int64, tok[1])
        tok = split(readline(f))
        if (length(tok) != nbinstants)
            print("pb fichier données, incohérence entre le nombre d'instants et le nombre de valeurs\n")
        else
            for i in 1:nbinstants
                #print(i, "\n")
                push!(Pdem, parse(Int64, tok[i]))
            end
        end
    end

    losscoefSEin=-0.070
    losscoefSEout=0.075
    p_FCmax=60
    p_FCmin=1
    p_ESmax=60
    p_ESmin=-60
    ESmin=400
    ESmax=1600
    ESinit=900
    HEVData(Pdem, nbinstants, losscoefSEin, losscoefSEout, p_FCmax, p_FCmin, p_ESmax, p_ESmin, ESmin, ESmax, ESinit)
end

function read_cMCND(graphfilepathname, commodityfilepathname, nodefilepathname, incremcost, congestioncoef)

    infoGraph=CSV.read(graphfilepathname, NamedTuple)
    #println(infoGraph)
    #println(infoGraph[4])

    #sort!(infoGraph)
    commodity=CSV.read(commodityfilepathname, NamedTuple)
    infoNode=CSV.read(nodefilepathname, NamedTuple)
    #println(infoNode)

    #intialisation du graphe
    nbNode=length(infoNode[1]) #infoGraph[end,:from]
    nbArc= length(infoGraph[1])
    nbCommod= length(commodity[1])

    #println("fff", nbNode, "ggg", nbArc, "hhh", nbCommod)

    #alpha = 3.0

    listofinputarcs = [[] for i=1:nbNode]
    listofoutputarcs = [[] for i=1:nbNode]
    #nodes=Array{NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay),Tuple{Int64, Float64,Float64,Float64,Float64,Float64}},1} ()
    nodes=Array{NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement),Tuple{Int64, Float64,Int64,Int64,Float64,Float64,Float64}}, 1}()
    thenames=(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement)
    for i in 1:nbNode
        incremcostvalue=-1
        if incremcost == 6
            incremcostvalue = infoNode.IncremCosts6[i]
        elseif incremcost == 7
            incremcostvalue = infoNode.IncremCosts7[i]
        elseif incremcost == 8
            incremcostvalue = infoNode.IncremCosts7[i]
        else
            println("invalid parameter incremcost, it shoud be 6, 7 ou 8")
        end

        vals=(Int64(i), incremcostvalue, Int64(infoNode.InitCapacity[i]), Int64(infoNode.IncremCapacity[i]), congestioncoef*Float64(infoNode.D[i]), 1.0, Float64(-1.0))
        push!(nodes, NamedTuple{thenames}(vals))
    end

    arcs=Array{NamedTuple{(:ref, :from, :to, :fixedcost, :varcost, :capa),Tuple{Int64,Int64,Int64,Float64,Float64,Int64}}, 1}()
    thenames=(:ref, :from, :to, :fixedcost, :varcost, :capa)
    #println(infoGraph)
    for i in 1:nbArc
        #println(Int64(infoGraph.from[i]))
        #println(Int64(infoGraph[i, :to]), "\n")
        #println(Float64(infoGraph.fixed_cost[i]), "\n")
        #println(Float64(infoGraph.variable_cost[i]), "\n")
        #println(infoGraph[i,:capacity], "\n")
        #println("\n")
        #println("\n")
        #vals=(Int64(i), Int64(infoGraph[i, :from]), Int64(infoGraph[i, :to]), Float64(infoGraph.fixed_cost[i]), Float64(infoGraph.variable_cost[i]), infoGraph[i,:capacity])
        vals=(Int64(i), Int64(infoGraph.from[i]), Int64(infoGraph.to[i]), Float64(infoGraph.fixed_cost[i]), Float64(infoGraph.variable_cost[i]), infoGraph.capacity[i])
        push!(arcs, NamedTuple{thenames}(vals))
        push!(listofinputarcs[infoGraph.to[i]],i)
        push!(listofoutputarcs[infoGraph.from[i]],i)
        #print(vals)
    end

    commods=Array{NamedTuple{(:ref, :from, :to, :qty),Tuple{Int64,Int64,Int64,Int64}}, 1}()
    thenames=(:ref, :from, :to, :qty)
    for i in 1:nbCommod
        vals = (Int64(i), Int64(commodity.from[i]), Int64(commodity.to[i]), Int64(commodity.quantity[i]))
        push!(commods, NamedTuple{thenames}(vals))
    end

    nodes, arcs, commods, listofinputarcs, listofoutputarcs#, nbNode, nbArc, nbCommod#, alpha
    #MCNDData(nodes, arcs, commods, listofinputarcs, listofoutputarcs, nbNode, nbArc, nbCommod, alpha)
end


function load_mcnddata(graphfilepathname, commodityfilemathname, nodefilepathname, incremcost, ptype;
                        with_cost_type = :exp,
                        restrict_to_commodities = 0,
                        printon=true,
                        congestioncoef=1.0)

    nodes, arcs, commods, listofinputarcs, listofoutputarcs = read_cMCND(graphfilepathname, commodityfilemathname, nodefilepathname, incremcost, congestioncoef)

    alpha = 3.0
    beta = 1.0 #not included in the nonlinear function

    nbNode = length(nodes)
    nbArc = length(arcs)
    nbCommod = length(commods)
    if restrict_to_commodities > 0
        nbCommod = min(nbCommod, restrict_to_commodities)
        resize!(commods, nbCommod)
    end

    nodesupdated=NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement),Tuple{Int64, Float64,Int64,Int64,Float64,Float64,Float64}}[]
    thenames=(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement)
    cubicexpcoefs=NamedTuple{(:a1, :a2, :a3, :a4, :a5),Tuple{Float64, Float64,Float64,Float64,Float64}}[]

    if ptype ∈ [:mcnd_nodes, :mcf_nodes]
        funcs = Array{Any, 1}(undef, nbNode)
        func1 = Array{Any, 1}(undef, nbNode)
        func2 = Array{Any, 1}(undef, nbNode)
        exprs = Array{Expr, 1}(undef, nbNode)
        expr1 = Array{Expr, 1}(undef, nbNode)
        expr2 = Array{Expr, 1}(undef, nbNode)
        for i in 1 : nbNode
            maxx = nodes[i].capainit + nodes[i].capaincrement
            f, g = construct_nonlinear_functions(maxx, nodes[i].fixedcost, with_cost_type, :MIN)
            funcs[i] = f
            func1[i] = f
            exprs[i] = g
            expr1[i] = g
        end
    elseif ptype ∈ [:mcnd_arcs, :mcf_arcs]
        funcs = Array{Any, 1}(undef, nbArc)
        func1 = Array{Any, 1}(undef, nbArc)
        func2 = Array{Any, 1}(undef, nbArc)
        exprs = Array{Expr, 1}(undef, nbArc)
        expr1 = Array{Expr, 1}(undef, nbArc)
        expr2 = Array{Expr, 1}(undef, nbArc)
        for i in 1 : nbArc
            maxx = arcs[i].capa
            maxcost = arcs[i].fixedcost + arcs[i].varcost * arcs[i].capa
            f, g = construct_nonlinear_functions(maxx, maxcost, with_cost_type, :MIN)
            funcs[i] = f
            func1[i] = f
            exprs[i] = g
            expr1[i] = g
        end
    end
    for i in 1 : nbNode
        vals=(nodes[i].ref, nodes[i].fixedcost, nodes[i].capainit, nodes[i].capaincrement, nodes[i].delaycost, 1.0, 0.0)
        push!(nodesupdated, NamedTuple{thenames}(vals))
    end
    pdata = MCNDData(nodesupdated, arcs, commods, listofinputarcs, listofoutputarcs, nbNode, nbArc, nbCommod, alpha, cubicexpcoefs)
    pdata, exprs, expr1, expr2, funcs, func1, func2
end

function read_scheduling(filename)
    ntasks = 0
    demands = Int64[]
    durations = Int64[]
    time_windows = Array{Tuple{Int64, Int64}, 1}()
    open(filename) do f
        ntasks = parse(Int64, readline(f))
        for ln in eachline(f)
            tok = split(ln, ' ')
            push!(demands, parse(Int64, tok[1]))
            push!(durations, parse(Int64, tok[2]))
            push!(time_windows, (parse(Int64, tok[3]) + 1, parse(Int64, tok[4]) + 1))
        end
    end
    SchedulingData(ntasks, durations, demands, time_windows)
end

function random_uflp(ndeps, ncusts)
    fixed_costs = rand(5 : 20, ndeps)
    caps = rand(500 : 1000, ndeps)
    demands = rand(50 : 200, ncusts)
    UFLPData(fixed_costs, caps, demands, true)
end

function read_gunluk_file(filename::String)
    ndeps = 0
    ncust = 0
    caps = Int64[]
    fc = Float64[]
    dems = Int64[]
    ass = Matrix{Float64}(undef, 0, 0)
    open(filename) do f
        while !eof(f)
            line = readline(f)
            tok = split(line, [' ', ':']; keepempty = false)
            if tok[1] == "DATASET"
            elseif tok[1] == "FACILITIES"
                ndeps = parse(Int64, tok[2])
            elseif tok[1] == "CUSTOMERS"
                ncust = parse(Int64, tok[2])
            elseif tok[1] == "FIXED_COSTS"
                for i in 1:ndeps
                    line = readline(f)
                    tok = split(line, [' ', ':']; keepempty = false)
                    push!(fc, parse(Float64, tok[2]))
                end
            elseif tok[1] == "CAPACITIES"
                for i in 1:ndeps
                    line = readline(f)
                    tok = split(line, [' ', ':']; keepempty = false)
                    push!(caps, parse(Int64, tok[2]))
                end
            elseif tok[1] == "DEMANDS"
                for i in 1:ncust
                    line = readline(f)
                    tok = split(line, [' ', ':']; keepempty = false)
                    push!(dems, parse(Int64, tok[2]))
                end
            elseif tok[1] == "ASSIGNMENT_COSTS"
                ass = zeros(ndeps, ncust)
                for i in 1:ndeps, j in 1:ncust
                    line = readline(f)
                    tok = split(line, [' ', ':']; keepempty = false)
                    i = parse(Int64, tok[1])
                    j = parse(Int64, tok[2])
                    ass[i, j] = parse(Float64, tok[3])
                end
            end
        end
    end
    return UFLPData(fc, caps, dems, ass, false)
end


function read_uflp_holmberg(filename; with_cap_constraint = true, scale_f = 1.0, scale_ass = 1.0)
    fixed_costs = Float64[]
    caps = Int64[]
    demands = Int64[]
    asscosts = zeros(0, 0)
    ncusts = 0
    ndeps = 0
    open(filename) do f
        let
            first = split(readline(f), [' ', '.']; keepempty = false)
            ndeps = parse(Int64, first[1])
            ncusts = parse(Int64, first[2])
            fixed_costs = zeros(Float64, ndeps)
            caps = zeros(Int64, ndeps)
            demands = zeros(Int64, ncusts)
            asscosts = zeros(Float64, ndeps, ncusts)
        end
        for i in 1 : ndeps
            line = split(readline(f), [' ', '.']; keepempty = false)
            caps[i] = parse(Int64, line[1])
            fixed_costs[i] = parse(Int64, line[2]) * scale_f
        end
        let
            j = 0
            while j < ncusts
                line = split(readline(f), [' ', '.']; keepempty = false)
                for l in line
                    j += 1
                    demands[j] = parse(Int64, l)
                end
            end
        end
        for i in 1 : ndeps
            j = 0
            while j < ncusts
                line = split(readline(f), [' ', '.']; keepempty = false)
                for l in line
                    j += 1
                    asscosts[i, j] = parse(Int64, l) * scale_ass
                end
            end
        end
    end
    UFLPData(fixed_costs, caps, demands, asscosts, with_cap_constraint)
end

function read_uflp_orlib(filename; with_cap_constraint = true, scale_f = 1.0, scale_ass = 1.0)
    fixed_costs = Float64[]
    caps = Int64[]
    demands = Int64[]
    asscosts = zeros(0, 0)
    ncusts = 0
    ndeps = 0
    costmult = 1
    open(filename) do f
        let
            first = split(readline(f), ' '; keepempty = false)
            ndeps = parse(Int64, first[1])
            ncusts = parse(Int64, first[2])
            fixed_costs = zeros(Float64, ndeps)
            caps = zeros(Int64, ndeps)
            demands = zeros(Int64, ncusts)
            asscosts = zeros(Float64, ndeps, ncusts)
        end
        for i in 1 : ndeps
            line = split(readline(f), ' '; keepempty = false)
            caps[i] = parse(Int64, line[1])
            fixed_costs[i] = parse(Float64, line[2]) * scale_f
        end
        for j in 1 : ncusts
            demands[j] = parse(Int64, readline(f))
            i = 0
            while i < ndeps
                line = split(readline(f), ' '; keepempty = false)
                for l in line
                    i += 1
                    asscosts[i, j] = parse(Float64, l) * scale_ass
                end
            end
        end
    end
    UFLPData(fixed_costs, caps, demands, asscosts, with_cap_constraint)
end

function read_transportation(filename)
    caps = Array{Int64, 1}()
    dems = Array{Int64, 1}()
    ass = zeros(Int64, 0, 0)
    open(filename) do f
        tok = split(readline(f))
        ndeps = parse(Int64, tok[1])
        ncusts = parse(Int64, tok[2])
        caps = zeros(Int64, ndeps)
        dems = zeros(Int64, ncusts)
        ass = zeros(Int64, ndeps, ncusts)
        while !eof(f)
            tok = split(readline(f))
            if tok[1] == "cap"
                i = parse(Int64, tok[2])
                c = parse(Int64, tok[3])
                caps[i] = c
            elseif tok[1] == "dem"
                j = parse(Int64, tok[2])
                d = parse(Int64, tok[3])
                dems[j] = d
            elseif tok[1] == "cost"
                i = parse(Int64, tok[2])
                j = parse(Int64, tok[3])
                c = parse(Int64, tok[4])
                ass[i, j] = c
            end
        end
    end
    TranspProblemData(caps, dems, ass)
end


function read_knapsacknlobj(filename, are_vars_integer, ctype)
    accuracy=0.00001::Real
    minxvalueforgain=0.0::Real
    digit=6::Int64

    paramN=0
    paramM=0
    paramUmax=0
    listea = Float64[]
    minuslistea = Float64[]
    object_coefb = Float64[]
    object_coefc = Float64[]
    object_coefd = Float64[]
    knapsack_cap = Float64[]
    object_U = Float64[]
    object_weight = Float64[]

    filed = open(filename)
        while !eof(filed)
            data=readline(filed)
            #println(data)
            data=rstrip(data)
            if cmp(data[1:min(sizeof(data),11)] , String("param N := ") ) == 0
                paramN = parse(Int, data[min(sizeof(data),11):end-1])
                #println(paramN)
            end
            if cmp(data[1:min(sizeof(data),11)] , String("param M := ") ) == 0
                paramM = parse(Int, data[min(sizeof(data),11):end-1])
            elseif cmp(data[1:min(sizeof(data),14)] , String("param Umax := ") ) == 0
                paramUmax = parse(Float64, data[min(sizeof(data),14):end-1])
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: C :=") ) == 0
                for count in 1:paramM
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(knapsack_cap,parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: a :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(listea,parse(Float64,dd[2]))
                    push!(minuslistea,-parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: b :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(object_coefb,parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: c :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(object_coefc,parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: d :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(object_coefd,parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),11)] , String("param: U :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(object_U,parse(Float64,dd[2]))
                end
            elseif cmp(data[1:min(sizeof(data),16)] , String("param: weight :=") ) == 0
                for count in 1:paramN
                    data=readline(filed)
                    data=rstrip(data)
                    dd=split(data,"\t")
                    push!(object_weight,parse(Float64,dd[2]))
                end
            end
            #println(paramN)
        end
    close(filed)

    #println(minuslistea)

    funcs = Matrix{Any}(undef, paramM, paramN)
    minusfuncs = Matrix{Any}(undef, paramM, paramN)
    exprs = Matrix{Expr}(undef, paramM, paramN)
    minusexprs = Matrix{Expr}(undef, paramM, paramN)
    
    for count in 1:paramN
        for k in 1:paramM
            fkc = x -> ((object_coefc[count]))/(1+(object_coefb[count])*(exp((minuslistea[count])*(x+(object_coefd[count])))))
            if ctype == :divexp
                funcs[k, count] = fkc
                exprs[k, count] = :(($(object_coefc[count]))/(1+$(object_coefb[count])*(exp($(minuslistea[count])*(x+$(object_coefd[count]))))))
                println("expr $k, $count = $(string(exprs[k, count]))")
                minusfuncs[k, count] = x -> -1.0*fkc(x)
                minusexprs[k, count] = :(-1.0*($(object_coefc[count]))/(1+$(object_coefb[count])*(exp($(minuslistea[count])*(x+$(object_coefd[count]))))))

            else
                # funcs[k, count], exprs[k, count] = construct_nonlinear_functions(object_U[count], 100.0, ctype, :MIN)
                maxc = fkc(object_U[count])
                funcs[k, count], exprs[k, count] = construct_nonlinear_functions(object_U[count], maxc, ctype, :MIN)
                minusfuncs[k, count] = x -> -funcs[k, count](x)
                minusexprs[k, count] = :(-1.0 * eval(exprs[k, count]))
            end
        end

        #("thefunctionis = cj/(1+bjj*(exp(-aj*(q+dj)))) ;\n")
    end
    # exit(1)


    """
        #using Plots
        thefig = plot(funcs[1], 0, object_U[1], label="function_S", title="ThefunctionS")
        thefig=plot!(funcs[2], 0, object_U[2], label="2");
        plot!([zoLimitZones[1][end]], [funcs[1](zoLimitZones[1][end])], seriestype=:scatter, label="data 1")
        plot!([zoLimitZones[2][end]], [funcs[2](zoLimitZones[2][end])], seriestype=:scatter, label="data 2")
        @show thefig
    """


    #pdata = KnapData(object_U, knapsack_cap, object_weight, elLimitZones, zoLimitZones)
    pdata = KnapData(are_vars_integer, minuslistea, object_coefb, object_coefc, object_coefd, object_U, knapsack_cap, object_weight)
    pdata, minusexprs, minusfuncs

end
