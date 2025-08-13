function SUNGetBreakpoints(pwl, ε = 1e-5)


    bpx = [pwl[1].xMin]
    #bpy = [pwl[1].a*pwl[1].xMin+pwl[1].b]
    bpy = [pwl[1](pwl[1].xMin)]


    for i in 1:length(pwl)-1

        temp = pwl[i].xMax
        push!(bpx,temp)
        #temp = pwl[i].a*pwl[i].xMax+pwl[i].b
        push!(bpy,pwl[i](temp))

        if (abs(pwl[i](temp) - pwl[i+1](temp)) > ε)
            push!(bpx, temp)
            #temp = pwl[i].a*pwl[i].xMax+pwl[i].b
            push!(bpy, pwl[i+1](temp))
        end

    end

    push!(bpx,pwl[end].xMax)
    push!(bpy,pwl[end](pwl[end].xMax))

    return bpx, bpy


end

function CCGetBreakpoints(pwl, ε = 1e-5)
    bpx = Float64[]
    bpy = Float64[]
    for p in pwl
        push!(bpx, p.xMin)
        push!(bpx, p.xMax)
        ymin = p(p.xMin)
        ymax = p(p.xMax)
        if abs(ymin - round(ymin)) < ε 
            ymin = round(ymin)
        end 
        if abs(ymax - round(ymax)) < ε 
            ymax = round(ymax)
        end 
        push!(bpy, ymin)
        push!(bpy, ymax)
    end
    return bpx, bpy
end


function construct_linear_expressions(exprs1, pdata::Data, params::Parameters)
    ptype = params.ptype
    pwl = Array{Array{PWLInfo, 1}, 1}()
    cis = CartesianIndices(exprs1)
    for (t, exp) in enumerate(exprs1)
        minw = 0.0
        midw = 0.0
        maxw = 0.0
        if ptype == :transp
            i, j = cis[t][1], cis[t][2]
            maxw = min(pdata.caps[i], pdata.demands[j])
        elseif ptype == :uflp_nlass
            i, j = cis[t][1], cis[t][2]
            maxw = 1.0
        elseif ptype == :uflp_nldep
            i = cis[t][1]
            maxw = pdata.capcons ? pdata.caps[i] : sum(pdata.demands)
        elseif ((ptype == :mcnd_nodes) || (ptype == :mcf_nodes))
            maxw = pdata.nodes[t].capainit + pdata.nodes[t].capaincrement
        elseif ((ptype == :mcnd_arcs) || (ptype == :mcf_arcs))
            maxw = pdata.arcs[t].capa#pdata.nodes[t].capainit + pdata.nodes[t].capaincrement
            #maxw = capacite max including increment !!!! SUN (max value possible for input value in pwl)##SUNSUN
        elseif ptype == :hev
            maxw = pdata.p_FCmax
            midw = pdata.p_FCmin
        elseif ptype == :scheduling
            maxw = sum(pdata.demands)
        elseif ptype == :knap
            # item = 1+div(t-1,length(pdata.knapsack_cap))
            _, item = cis[t][1], cis[t][2]
            #minw = pdata.zoLimitZones[item][end]
            #maxw = pdata.elLimitZones[item][end] ????
            maxw = pdata.object_U[item]
        else
            maxw = sum(pdata.demands)
        end

        found_before = false
        for s in 1 : t - 1
            if almost_equal_functions(exprs1[s], exp, minw, maxw) && abs(pwl[s][begin].lp.xMin - minw) < 1e-9 && abs(pwl[s][begin].lp.xMax - maxw) < 1e-9
                push!(pwl, pwl[s])
                found_before = true
                # println("found same cost before: $exp")
                break
            end
        end
        if found_before continue 
        else
            # println("expression not found before: $exp")
        end
        # println(exp)
        # println(typeof(exp))
        if abs(exp(minw) - exp(maxw)) < 1e-16
            # println("zero function $exp")
            pwlf = [LinA.LinearPiece(0.0, maxw, 0.0, 0.0, x -> 0.0)]
        else
            pwlf = []
            params.eps0 = min(params.eps0, 50.0)
            minw = 0.0
            pwlf = LinA.Linearize(exp, minw, maxw, Relative(params.eps0), params.with_lina_solver(); bounding = Under())
        end

        if params.with_fixed_charge
            insert!(pwlf, 1, LinA.LinearPiece(0.0, 0.0, 0.0, 0.0, x -> 0.0))
        end
        newpwlf = [PWLInfo(p, params.eps0) for p in pwlf]
        push!(pwl, newpwlf)
        #println("finished LinA for t = $t out of $(length(exprs1))")# with $(length(pwlf)) pieces")
    end
    # println("lower bound for expr 9 = $(pwl[9][end].lp.a) * x + $(pwl[9][end].lp.b), valid in [$(pwl[9][end].lp.xMin), $(pwl[9][end].lp.xMax)]")
    # println.(pwl[9])
    pwl
end

function compute_x_at_absobj(f::Union{Function, Expr}, x0::Real, xf::Real, target::Float64)::Real
    xs = IntervalRootFinding.roots(t -> f(t) - target, interval(x0, xf))
    if isempty(xs) return x0 end
    # @assert(!isempty(xs), "No target found for f between $x0 and $xf")
    x = (xs[begin].region.bareinterval.lo + xs[begin].region.bareinterval.hi)
    @assert(!isnan(x), "Did not find target $target between $x0 and $xf")
    return x
end

function check_lowerbound_pwl(pwl::Vector{PWLInfo}, f::Function)::Bool
    for p in pwl
        lp = p.lp
        xmin, xmax = lp.xMin, lp.xMax
        g(x) = f(x) - lp.fct(x)
        fm = fmin(g, xmin, xmax)
        if fm < -EPS_NZ 
            # println("error, linear piece is not lower bounding between $xmin and $xmax, error of $fm")
            return false
        end
    end
    return true
end

function almost_equal_functions(f, g, l, u)
    step = (u - l) / 1000
    x = collect(l : step : u)
    for u in x
        if abs(f(u) - g(u)) > EPS_NZ return false
        end
    end
    return true
end

function compute_max_error(w, f, pwl; minw = -EPS_NZ)
    nperiods = length(w)
    max_err, argt, argi = 0, 0, 0
    for t in 1 : nperiods
        #if w[t] < 1e-9 continue # deactivated by SUN to handle the knapsack problem
        if w[t] < minw continue # inserted by SUN to handle the knapsack problem
        end
        ub = f[t](w[t])
        lb = ub + 1
        i = -1 #0 modified by SUN
        for (j, pwlf) in enumerate(pwl[t])
            if w[t] >= (pwlf.lp.xMin-1e-9) && w[t] <= (pwlf.lp.xMax+1e-9)
                newlb = pwlf.lp.fct(w[t])
                if newlb < lb
                    lb = newlb
                    i = j
                end
            end
        end
        rel_err = (ub - lb) / abs(ub) * 100.0
        if rel_err > max_err
            max_err = rel_err
            argt = t
            argi = i
        end
    end
    max_err, argt, argi
end

function compute_all_errors(w, f, pwl; minw = -EPS_NZ)
#    println("computing all errors for w = $w, minw = $minw")
    nperiods = length(w)
    rel_err = zeros(Float64, nperiods)
    abs_err = zeros(Float64, nperiods)
    indices = zeros(Int64, nperiods)
    #println("f[3](63.154437002642325)=",f[3](63.154437002642325), "\n")
    for t in 1 : nperiods
        #if w[t] < 1e-9 continue # deactivated by SUN to handle the knapsack problem
        if w[t] < minw continue # inserted by SUN to handle the knapsack problem
        end
        # println("w[t] = $(w[t]) and f(w[t]) = $(f[t](w[t])), xmin = $(pwl[t][begin].lp.xMin), xmax = $(pwl[t][end].lp.xMax)")
        ub = f[t](w[t])
        linveceval = pwl[t](w[t])
        lb = linveceval.val
        if ub - lb < 1e-12 continue
        else
            rel_err[t] = (ub - lb) / max(1e-10, abs(ub), abs(lb)) * 100.0
            abs_err[t] = (ub - lb)
            errs2 = (ub - lb) / abs(lb) * 100.0
            indices[t] = linveceval.piece #0 modified by SUN
#            for (j, pwlf) in enumerate(pwl[t])
                # println("w[t] = $(w[t]) and (min, max) = ($(pwlf.lp.xMin), $(pwlf.lp.xMax))")
#                if w[t] >= (pwlf.lp.xMin - EPS_NZ) && w[t] <= (pwlf.lp.xMax + EPS_NZ) && abs(pwlf.lp(w[t]) - lb) < EPS_NZ
#                        indices[t] = j
#                end
#            end
        end
    end
    # println("done computing all errors")
    abs_err, rel_err, indices
end

function compute_upperbound(ptype, pdata, functions, wval, obj0, with_fixed_charge)
#    minw = ptype == :knap ? -EPS_NZ : EPS_NZ
    if ptype == :scheduling
        nperiods = maximum([tw[2] for tw in pdata.time_windows])
        return obj0 + sum(compute(functions[t], wval[t], with_fixed_charge) for t in 1 : nperiods)
    elseif ((ptype == :mcnd_nodes) || (ptype == :mcf_nodes))
        nperiods = pdata.nbNode
        return obj0 + sum(compute(functions[t], wval[t], with_fixed_charge) for t in 1 : nperiods)
    elseif ((ptype == :mcnd_arcs) || (ptype == :mcf_arcs))
        nperiods = pdata.nbArc
        return obj0 + sum(compute(functions[t], wval[t], with_fixed_charge) for t in 1 : nperiods)
    elseif ptype == :transp
        nperiods = length(functions)
        varobj = sum(compute(functions[t], wval[t], with_fixed_charge) for t in 1 : nperiods)
        return obj0 + varobj
    elseif ptype == :uflp_nlass
        ndeps = length(pdata.caps)
        ncusts = length(pdata.demands)
        w = reshape(wval, ndeps, ncusts)
        used = falses(ndeps)
        cost = 0.0
        for i in 1 : ndeps, j in 1 : ncusts
            cost += compute(functions[i, j], w[i, j], with_fixed_charge)
        end
        return obj0 + cost
    elseif ptype == :hev
        nperiods=pdata.nbinstants
        return obj0 + sum(functions[t](wval[t]) for t in 1 : nperiods if wval[t] >= pdata.p_FCmin-1e-9)
        #return obj0 + sum(functions[t](wval[t]) for t in 1 : nperiods if wval[t] > 1e-9) #SUN

    elseif ptype == :knap
        nperiods = length(functions)
        return obj0 + sum(compute(functions[t], wval[t], with_fixed_charge) for t in 1 : nperiods)
    else
        ndeps = length(pdata.caps)
        ncusts = length(pdata.demands)
        w = reshape(wval, ndeps)
        cost = 0.0
        for i in 1 : ndeps
            cost += compute(functions[i], w[i], with_fixed_charge)
        end
        return obj0 + cost
    end
end

function fmax(f::Function, xmin::Real, xmax::Real)::Real
    res = optimize(x -> -f(x), xmin, xmax)
    return -res.minimum
end

function fmin(f::Function, xmin::Real, xmax::Real)::Real
    res = optimize(f, xmin, xmax)
    return res.minimum
end

function update_null_pwlf(pwlf, exp, f, i, x, delta, neweps, params::Parameters)
    p = pwlf[i].lp
    @assert(is_null_function(f, p.xMin, p.xMax), "component does not seem null!")
    # println("updating [$(p.xMin), $(p.xMax)] at x = $x")
    midp = (p.xMin + p.xMax) / 2
    # while abs(midp - p.xMin) < EPS_NZ midp += EPS_NZ end
    # while abs(midp - p.xMax) < EPS_NZ midp -= EPS_NZ end
    if is_null_function(f, p.xMin, midp)
        minfirst = fmin(f, p.xMin, midp)
        first = [PWLInfo(LinA.LinearPiece(p.xMin, midp, 0.0, minfirst, x -> minfirst), pwlf[i].eps)]
    else
        lps = Linearize(f, p.xMin, midp, Relative(neweps), params.with_lina_solver(); bounding = Under())
        first = [PWLInfo(lp, neweps) for lp in lps]
    end
    if is_null_function(f, midp, p.xMax)
        minsecond = fmin(f, midp, p.xMax)
        second = [PWLInfo(LinA.LinearPiece(midp, p.xMax, 0.0, minsecond, x -> minsecond), pwlf[i].eps)]
    else
        lps = Linearize(f, midp, p.xMax, Relative(neweps), params.with_lina_solver(); bounding = Under())
        second = [PWLInfo(lp, neweps) for lp in lps]
    end
    return vcat(pwlf[begin : i - 1], first, second, pwlf[i + 1 : end])
end

function update_small_pwlf(pwlf, exp, f, i, x, delta, neweps::Float64, params::Parameters)
    # println("running small pwlf")
    first = last = i
    while first > 1 && pwlf[first - 1].lp.fct == f #SUN
        if abs(pwlf[first - 1].eps - neweps) < 1e-9
            first -= 1
        else break
        end
    end
    while last < length(pwlf) && pwlf[last + 1].lp.fct == f #SUN
        if abs(pwlf[last + 1].eps - neweps) < 1e-9
            last += 1
        else break
        end
    end
    interval = pwlf[last].lp.xMax - pwlf[first].lp.xMin
    # delta = 1e-1
    while interval < delta && (first > 1 || last < length(pwlf))
        left = right = -1
        if first > 1 && pwlf[first - 1].lp.fct == f #SUN
            left = pwlf[first - 1].lp.xMax - pwlf[first - 1].lp.xMin
        end
        if last < length(pwlf) && pwlf[last + 1].lp.fct == f #SUN
            right = pwlf[last + 1].lp.xMax - pwlf[last + 1].lp.xMin
        end
        if first > 1 && (right < 0 || left < right)
            first -= 1
        elseif last < length(pwlf)
            last += 1
        end
        neweps = min(min(neweps, pwlf[first].eps), pwlf[last].eps)
        interval = pwlf[last].lp.xMax - pwlf[first].lp.xMin
    end

    # if printon
    #     println("consolidating $(last - first + 1) items")
    #     println([(pwlf[k].lp, pwlf[k].eps) for k in first : last])
    #     println("new eps = $neweps")
    # end
    lhs = first > 1 ? pwlf[1 : first - 1] : []
    rhs = last < length(pwlf) ? pwlf[last + 1 : end] : []
    newlin = LinA.Linearize(f, pwlf[first].lp.xMin, pwlf[last].lp.xMax, Relative(neweps), params.with_lina_solver(); bounding = Under()) #SUN
    newpwlf = [PWLInfo(p, neweps) for p in newlin]
    pwlfvec = vcat(vcat(lhs, newpwlf), rhs)
    # println("small $pwlfvec")
    return pwlfvec
end

function update_large_pwlf(pwlf, exp, f, i, x, delta::Float64, neweps::Float64, params::Parameters)
    # println("running large pwlf")
    xmax::Float64 = min(x + delta / 2, pwlf[i].lp.xMax)
    xmin::Float64 = xmax - delta
    if xmin < pwlf[i].lp.xMin
        d = pwlf[i].lp.xMin - xmin
        xmin = pwlf[i].lp.xMin
        xmax += d
    end

    lhs = i > 1 ? pwlf[1 : i - 1] : []
    rhs = i < length(pwlf) ? pwlf[i + 1 : end] : []
    if pwlf[i].lp.xMin < xmin - 1e-9
        push!(lhs, PWLInfo(LinA.LinearPiece(pwlf[i].lp.xMin, xmin, pwlf[i].lp.a, pwlf[i].lp.b, x-> pwlf[i].lp.a * x + pwlf[i].lp.b), pwlf[i].eps))
        #push!(lhs, PWLInfo(LinA.LinearPiece(pwlf[i].lp.xMin, xmin, pwlf[i].lp.a, pwlf[i].lp.b, f), pwlf[i].eps)) #SUNcorrected f
    end
    if pwlf[i].lp.xMax > xmax + 1e-9
        pushfirst!(rhs, PWLInfo(LinA.LinearPiece(xmax, pwlf[i].lp.xMax, pwlf[i].lp.a, pwlf[i].lp.b, x -> pwlf[i].lp.a * x + pwlf[i].lp.b), pwlf[i].eps))
        #pushfirst!(rhs, PWLInfo(LinA.LinearPiece(xmax, pwlf[i].lp.xMax, pwlf[i].lp.a, pwlf[i].lp.b, f), pwlf[i].eps)) #SUNcorrected f
    end

    newlin = LinA.Linearize(f, xmin, xmax, Relative(neweps), params.with_lina_solver(); bounding = Under()) #SUN
    newpwlf = [PWLInfo(p, neweps) for p in newlin]
    pwlfvec = vcat(vcat(lhs, newpwlf), rhs)
    # println("large $pwlfvec")
    return pwlfvec
end

function update_simple_pwlf(pwlf, exp, f, i, x, delta::Float64, neweps::Float64, params::Parameters)
    # println("running simple pwlf")
    lhs = i > 1 ? pwlf[1 : i - 1] : []
    rhs = i < length(pwlf) ? pwlf[i + 1 : end] : []
    xmin::Float64 = pwlf[i].lp.xMin
    xmax::Float64 = pwlf[i].lp.xMax
    newlin = LinA.Linearize(f, xmin, xmax, Relative(neweps), params.with_lina_solver(); bounding = Under()) #SUN
    newpwlf = [PWLInfo(p, neweps) for p in newlin]
    pwlfvec = vcat(vcat(lhs, newpwlf), rhs)
    # println("simple $pwlfvec")
    return pwlfvec
end

function update_pwlf_burlacu(pwlf, exp, f::Function, i::Integer, x::Real, delta::Float64, neweps::Float64, params::Parameters)
    lhs = i > 1 ? pwlf[1 : i - 1] : []
    rhs = i < length(pwlf) ? pwlf[i + 1 : end] : []
    xmin = pwlf[i].lp.xMin
    xmax = pwlf[i].lp.xMax
    xmid = (xmax + xmin) / 2
    ymin = f(xmin)
    ymax = f(xmax)
    ymid = f(xmid)
    m1 = (ymid - ymin) / (xmid - xmin)
    m2 = (ymax - ymid) / (xmax - xmid)
    # linear function: y = ymin + m * (x - xmin)
    viol_left = 0.0
    viol_right = 0.0
    # function newf(y::Real)
    #     return f(y)
    # end

    let
        mleft = JuMP.Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
        # println("f = $(f(1))")
        register(mleft, :f, 1, f; autodiff = true)
        @variable(mleft, xmin <= xvar <= xmid)
        @NLobjective(mleft, Max, ymin + m1 * (xvar - xmin) - f(xvar))
        # @NLexpression(mleft, w, lambda - m1 * xvar + ff(xvar))
        # @NLconstraint(mleft, w <= ymin - m1 * xmin)
        # @objective(mleft, Max, lambda)
        optimize!(mleft)
        viol_left = objective_value(mleft)
    end
    let
        mright = JuMP.Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
        # println("f = $f")
        register(mright, :f, 1, f; autodiff = true)
        @variable(mright, xmid <= xvar <= xmax)
        @NLobjective(mright, Max, ymid + m2 * (xvar - xmid) - f(xvar))
        # @variable(mright, lambda)
        # @NLexpression(mleft, w, lambda - m2 * xvar + ff(xvar))
        # @NLconstraint(mright, w <= ymid - m2 * xmid)
        # @objective(mright, Max, lambda)
        # @NLobjective(mright, Max, ymid + m2 * (xvar - xmid) - cos(xvar))
        optimize!(mright)
        viol_right = objective_value(mright)
    end

    # shifted linear function: y = ymin - viol + m * (x - xmin) = (ymin - viol - m * xmin) + m * x

    lp_left = LinA.LinearPiece(xmin, xmid, m1, ymin - viol_left - m1 * xmin, x -> m1 * x + ymin - viol_left - m1 * xmin)
    #lp_left = LinA.LinearPiece(xmin, xmid, m1, ymin - viol_left - m1 * xmin, f) #SUNcorrected f
    lp_right = LinA.LinearPiece(xmid, xmax, m2, ymid - viol_right - m2 * xmid, x -> m2 * x + ymid - viol_right - m2 * xmid )
    #lp_right = LinA.LinearPiece(xmid, xmax, m2, ymid - viol_right - m2 * xmid, f) #SUNcorrected f

    pwlfvec = vcat(vcat(lhs, [PWLInfo(lp_left, neweps), PWLInfo(lp_right, neweps)]), rhs)
    return pwlfvec
end

function update_pwlf(pwlf, exp, f::Function, i::Integer, x::Real, delta::Float64, neweps::Float64, mineps::Float64, params::Parameters)
    # return update_pwlf_burlacu(pwlf, exp, f, i, x, delta, neweps; printon = printon)
    # println("iteration $it of update_pwlf with eps = $neweps")
    d = pwlf[i].lp.xMax - pwlf[i].lp.xMin
    # update_fn = is_null_function(f, pwlf[i].lp.xMin, pwlf[i].lp.xMax) ? update_null_pwlf : ((d > 10 * delta) && activeUpdateLargePwlf ? update_large_pwlf : (d < delta ? update_small_pwlf : update_simple_pwlf))
    update_fn = is_null_function(f, pwlf[i].lp.xMin, pwlf[i].lp.xMax) ? update_null_pwlf : (d < delta ? update_small_pwlf : update_simple_pwlf)
    # common_delta = d > 10 * delta ? 10.0 * delta : delta
    lastpwlf = update_fn(pwlf, exp, f, i, x, delta, neweps, params)
    while true
        newpwlf = update_fn(pwlf, exp, f, i, x, delta, neweps, params)
        if length(newpwlf) > length(lastpwlf) || neweps < mineps
            break
        else 
            lastpwlf = newpwlf
            neweps /= 2
        end
    end
    sdif = length(lastpwlf) - length(pwlf)
    # if sdif >= 400
    #     x0 = pwlf[begin].lp.xMin
    #     xf = pwlf[end].lp.xMax
    #     println(stderr, "increase from $(length(pwlf)) to $(length(lastpwlf)) pwl functions for d = $d, neweps = $neweps, x1, x2 = $(pwlf[i].lp.xMin),  $(pwlf[i].lp.xMax), f valid in ($x0, $xf) with C = $(pwlf[end].lp(xf))")
    #     println.(stderr, lastpwlf)
    # end
    return lastpwlf
end

function is_null_piece(lp::LinA.LinearPiece)::Bool
    return min(abs(lp(lp.xMin)), abs(lp(lp.xMax))) < 1e-7
end

function is_null_function(f::Function, xmin::Real, xmax::Real)::Bool
    if min(abs(f(xmin)), abs(f(xmax))) < EPS_NZ return true
    elseif min(abs(fmin(f, xmin, xmax)), abs(fmax(f, xmin, xmax))) < EPS_NZ return true
    else return false
    end 
end

function compute_binary_representation!(nbits, binary, numbers, prefix)
    for n in numbers
        s = string(n; base = 2, pad = nbits)
        # println("n = $n, nnums = $(length(numbers)), string = $s, pad = $nbits")
        for (j, c) in enumerate(s)
            if c == '1'
                binary[n][j] = true
            end
        end
    end
    # nn = length(numbers)
    # if nn <= 1
    #     n = numbers[1]
    #     s = prefix
    #     for (j, c) in enumerate(s)
    #         if c == '1'
    #             binary[n][j] = true
    #         end
    #     end
    # else
    #     nhalf = ceil(Int64, nn / 2)
    #     leftnums = numbers[1 : nhalf]
    #     rightnums = numbers[nhalf + 1 : nn]
    #     compute_binary_representation!(nbits, binary, leftnums, "$(prefix)0")
    #     compute_binary_representation!(nbits, binary, rightnums, "$(prefix)1")
    # end
end

function compute_cutoff(ws, obj0s, pwl, pdata, ptype)
    ##SUNSUN remove all ifs
    cutoff = 1e+20
    #if ptype == :scheduling || ptype == :transp
    for (i, pw) in enumerate(ws)
        obj = obj0s[i]
        for (t, p) in enumerate(pwl)
            if ptype != :knap && pw[t] < 1e-9 
                continue
            end
            minobj = 1e+20
            for pwlf in p
                if pw[t] >= pwlf.lp.xMin && pw[t] <= pwlf.lp.xMax
                    #minobj = min(minobj, pwlf.lp.a * pw[t] + pwlf.lp.b) modified by SUN
                    minobj = min(minobj, pwlf.lp.a * pw[t] + pwlf.lp.b)
                end
            end
            obj += minobj
        end
        cutoff = min(cutoff, obj)
    end
    cutoff
end

function construct_nonlinear_functions(xMax::Real, cMax::Real, ctype::Symbol, sense::Symbol) 
    # println("nonlinear function for xMax = $xMax and cMax = $cMax for cost $(string(ctype))")
    f = zero_functor_function()
    g = zero_expression()
    if abs(cMax) > 1e-7
        if ctype == :A1
            f = FC_A1_functor_function()
            g = FC_A1_expression()
        elseif ctype == :A2
            f = FC_A2_functor_function()
            g = FC_A2_expression()
        elseif ctype == :R
            f = FC_R_functor_function()
            g = FC_R_expression()
        elseif ctype == :cubic
            xf, yf = cubic_function_multipliers(xMax, cMax)
            coefs = cubic_coefs(xf, yf)
            f = cubic_functor_function(coefs)
            g = cubic_expression(coefs) # for now it looks like this, but this is not the same as above!
        elseif ctype == :exp
            xf, yf = exp_function_multipliers(xMax, cMax)
            coefs = exp_coefs(xf, yf)
            f = exp_functor_function(coefs)
            g = exp_expression(coefs)
        elseif ctype == :log
            coefs = logfn_coefs(xMax, cMax)
            f = logfn_functor_function(coefs)
            g = logfn_expression(coefs)
        elseif ctype == :sqrt
            coefs = sqrt_coefs(xMax, cMax)
            f = sqrt_functor_function(coefs)
            g = sqrt_expression(coefs)
        elseif ctype == :square
            coef = square_coef(xMax, cMax)
            f = square_functor_function(coef)
            g = square_expression(coef)
        elseif ctype == :sincos
            f = sincos_functor_function(xMax, cMax)
            g = sincos_expression(xMax, cMax)
        elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
            if ctype == :sinabc1 k  = 1
            elseif ctype == :sinabc2 k = 2
            elseif ctype == :sinabc3 k = 3
            end
            coefs = sinabc_coefs(k, xMax, cMax)
            f = sinabc_functor_function(coefs)
            g = sinabc_expression(coefs)
        elseif ctype == :pardivexp
            coefs = parametric_divexp_coefs(xMax, cMax)
            f = parametric_divexp_functor_function(coefs)
            g = parametric_divexp_expression(coefs)
        else
            throw(ArgumentError("incorrect parameter with_cost_type = $with_cost_type !!!!!!"))
        end
    end
    return sense == :MIN ? (f, g) : (x -> -f(x), :(-1.0 * eval(g)))
end

function string_nonlinear_function(xMax::Real, cMax::Real, ctype::Symbol) 
    # println("nonlinear function for xMax = $xMax and cMax = $cMax for cost $(string(ctype))")
    fstr = ""
    if abs(cMax) > 1e-7
        if ctype == :A1
        elseif ctype == :A2
        elseif ctype == :R
        elseif ctype == :cubic
            xf, yf = cubic_function_multipliers(xMax, cMax)
            coefs = cubic_coefs(xf, yf)
            fstr = cubic_function_string(coefs)
        elseif ctype == :exp
            xf, yf = exp_function_multipliers(xMax, cMax)
            coefs = exp_coefs(xf, yf)
            fstr = exp_function_string(coefs)
        elseif ctype == :log
            coefs = logfn_coefs(xMax, cMax)
            fstr = logfn_function_string(coefs)
        elseif ctype == :sqrt
            coefs = sqrt_coefs(xMax, cMax)
            fstr = sqrt_function_string(coefs)
        elseif ctype == :square
            coef = square_coef(xMax, cMax)
            fstr = square_function_string(coef)
        elseif ctype == :sincos
            fstr = sincos_function_string(xMax, cMax)
        elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
            if ctype == :sinabc1 k  = 1
            elseif ctype == :sinabc2 k = 2
            elseif ctype == :sinabc3 k = 3
            end
            coefs = sinabc_coefs(k, xMax, cMax)
            fstr = sinabc_function_string(coefs)
        elseif ctype == :pardivexp
            coefs = parametric_divexp_coefs(xMax, cMax)
            fstr = parametric_divexp_function_string(coefs)
        else
            throw(ArgumentError("incorrect parameter with_cost_type = $with_cost_type !!!!!!"))
        end
    end
    return fstr
end

function get_nonlinear_coefficients(xMax::Float64, cMax::Float64, ctype::Symbol)
    if ctype == :cubic
        xf, yf = cubic_function_multipliers(xMax, cMax)
        return cubic_coefs(xf, yf)
    elseif ctype == :exp
        xf, yf = exp_function_multipliers(xMax, cMax)
        return exp_coefs(xf, yf)
    elseif ctype == :sqrt
        return sqrt_coefs(xMax, cMax)
    elseif ctype == :square
        return [square_coef(xMax, cMax)]
    elseif ctype == :log
        return logfn_coefs(xMax, cMax)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        if ctype == :sinabc1 k = 1
        elseif ctype == :sinabc2 k = 2
        elseif ctype == :sinabc3 k = 3
        end
        return sinabc_coefs(k, xMax, cMax)
    elseif ctype == :sincos
        return [xMax, cMax]
    elseif ctype == :pardivexp
        return parametric_divexp_coefs(xMax, cMax)
    end
end

function get_nonlinear_constraint_expression(model::JuMP.Model, cf::Vector{Float64}, x::JuMP.VariableRef, y::JuMP.VariableRef, ctype::Symbol, sense::Symbol)
    mult = sense == :MIN ? 1.0 : -1.0
    if ctype == :cubic
        return @NLexpression(model, mult * (cf[1] * x^3 + cf[2] * x^2 + cf[3] * x + cf[4] * y)), 0.0
    elseif ctype == :exp
        f0 = cf[5] + cf[3] / (1 + cf[2] * exp(-cf[1] * cf[4]))
        return @NLexpression(model, mult * (y * f0 + cf[3] / (1 + cf[2] * exp(-cf[1] * (x + cf[4]))))), mult * (f0 - cf[5])
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        return @NLexpression(model, mult * (cf[1] * ((sin(cf[2] * x) + cf[3] * x)^2))), 0.0
    elseif ctype == :sincos
        return @NLexpression(model, mult * (cf[2] * (7.5*sin(pi*(x * 100 / cf[1] - 10.0)/40.0) - 15.0*cos(pi*(x * 100 / cf[1] -10.0)/80.0) + 19.5) / 38.66)), 0.0
    elseif ctype in [:divexp, :pardivexp] #TODO: valider que je n'ai pas mis une parenthese au mauvais endroit
        f0 =  cf[3] / (1 + cf[2] * exp(cf[1] * cf[4]))
        println("f0 = $f0")
        return @NLexpression(model, mult * (y * f0 + cf[3] / (1 + cf[2] * exp(cf[1] * (x + cf[4]))))), mult * f0
    elseif ctype == :log
        return @NLexpression(model, mult * cf[3] * log(cf[1] + cf[2] * x)), 0.0
    elseif ctype == :sqrt
        return @NLexpression(model, mult * cf[2] * ((cf[1] * x)^(0.5))), 0.0
    elseif ctype == :square
        return @NLexpression(model, mult * cf[1] * (x^2)), 0.0
    end
end
