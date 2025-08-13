function generate_gams_files_flp_orlib()
    for i in 4 : 13, j in 1 : 4, cost in [:exp, :cubic]
        if i == 5 && j > 1 continue
        else
            num = 10 * i + j
            IterativePWLB.gams_uflp_nldep("../../instances/uflp/orlib/cap$(num).txt", cost, "../IterativePWLB/neos/python/minlp/Input/flp/cap$(num)_nldep_$(string(cost)).gms")
            IterativePWLB.gams_uflp_nlass("../../instances/uflp/orlib/cap$(num).txt", cost, "../IterativePWLB/neos/python/minlp/Input/flp/cap$(num)_nlass_$(string(cost)).gms")
        end
    end
end

function write_knapsack_trindade(nobjs, filename, seed)
    rng = MersenneTwister(seed)
    w = 1.0 * rand(rng, collect(1 : 100), nobjs)
    C = 50 * sum(w)
    U = 100 * ones(Int64, nobjs)
    a = 0.1 * (1 .+ rand(rng, nobjs))
    b = 100 * rand(rng, nobjs)
    c = 100 * rand(rng, nobjs)
    d = -100 * rand(rng, nobjs)
    open(filename, "w") do f
        @printf(f,"# seed = %d\n", seed)
        println(f)
        @printf(f, "param N := %d;\n", nobjs)
        println(f)
        @printf(f, "param M := %d;\n", 1)
        println(f)
        @printf(f, "param Umax := 100.000000;\n")
        println(f)
        @printf(f, "param: C :=\n")
        @printf(f, "1\t%.6f\n;\n", C)
        
        for (v, nv) in [(a, "a"), (b, "b"), (c, "c"), (d, "d"), (w, "weight"), (U, "U")]
            println(f)
            @printf(f, "param: %s :=\n", nv)
            for n in 1 : nobjs
                @printf(f, "%d\t%.6f\n", n, v[n])
            end
            println(f, ";")
        end
    end
end

function generate_gunluk_uflp(ndeps::Int64, ncusts::Int64, seed::Int64 = 0; with_cap_constraint = false)
    rng = MersenneTwister(seed)
    fixed_costs = rand(rng, 1:100, ndeps)
    dems = rand(rng, 1:100, ncusts)
    caps = with_cap_constraint ? ceil(Int64, sum(dems) / 5) * ones(Int64, ndeps) : sum(dems) * ones(Int64, ndeps)
    pos_deps = rand(rng, ndeps, 2)
    pos_custs = rand(rng, ncusts, 2)
    ass = zeros(ndeps, ncusts)
    for i in 1:ndeps, j in 1:ncusts
        ass[i, j] = 50 * sqrt(LinearAlgebra.dot(pos_deps[i, :] - pos_custs[j, :], pos_deps[i, :] - pos_custs[j, :])) # From Gunluk and Lee 2009
    end
    UFLPData(fixed_costs, caps, dems, ass, false)
end

function write_gunluk_file(data::UFLPData, filename::String)
    ndeps = length(data.caps)
    ncust = length(data.demands)
    open(filename, "w") do f
        println(f, "DATASET: Gunluk and Lee 2009")
        println(f, "FACILITIES: $(length(data.caps))")
        println(f, "CUSTOMERS: $(length(data.demands))")
        println(f, "FIXED_COSTS:")
        for i in 1:ndeps
            println(f, "$i $(data.fixed_costs[i])")
        end
        println(f, "CAPACITIES:")
        for i in 1:ndeps
            println(f, "$i $(data.caps[i])")
        end
        println(f, "DEMANDS:")
        for i in 1:ncust
            println(f, "$i $(data.demands[i])")
        end
        println(f, "ASSIGNMENT_COSTS:")
        for i in 1:ndeps, j in 1:ncust
            println(f, "$i $j $(data.asscosts[i, j])")
        end
    end
end

function generate_random_transp(ndeps, ncusts)
    a = shuffle(collect(1 : 999))
    a = sort(a[1 : ndeps - 1])
    caps = Int64[]
    for i in 1 : ndeps
        c = 0
        if i == 1
            c = a[1]
        elseif i == ndeps
            c = 1000 - a[i - 1]
        else
            c = a[i] - a[i - 1]
        end
        push!(caps, c)
    end
    a = shuffle(collect(1 : 999))
    #a = sort(a[1 : ndeps - 1]) ##modified by SUN
    a = sort(a[1 : ncusts - 1]) ##modified by SUN
    dems = Int64[]
    for i in 1 : ncusts
        c = 0
        if i == 1
            c = a[1]
        elseif i == ncusts
            c = 1000 - a[i - 1]
        else
            c = a[i] - a[i - 1]
        end
        push!(dems, c)
    end
    println("sum caps = $(sum(caps))")
    println("sum dems = $(sum(dems))")
    # assert(sum(dems) == sum(caps))
    ndeps = length(caps)
    ncusts = length(dems)
    dpos = rand(0 : 100, ndeps, 2)
    cpos = rand(0 : 100, ncusts, 2)
    unitcost = zeros(Int64, ndeps, ncusts)
    for i in 1 : ndeps, j in 1 : ncusts
        dist = ceil(Int64, sqrt((dpos[i, 1] - cpos[j, 1])^2 + (dpos[i, 2] - cpos[j, 2])^2))
        unitcost[i, j] = max(1, dist + rand(-3 : 3))
    end
    pdata = TranspProblemData(caps, dems, unitcost)
end

function write_transp(pdata, filename)
    ndeps = length(pdata.caps)
    ncusts = length(pdata.demands)
    open(filename, "w") do f
        println(f, "$ndeps $ncusts")
        for i in 1 : ndeps
            println(f, "cap $i $(pdata.caps[i])")
        end
        for j in 1 : ncusts
            println(f, "dem $j $(pdata.demands[j])")
        end
        for i in 1 : ndeps, j in 1 : ncusts
            println(f, "cost $i $j $(pdata.asscosts[i, j])")
        end
    end
end
