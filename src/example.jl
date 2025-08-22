using Gurobi
using JuMP
using LinA
using IterativePWLB

column(x::VariableRef) = Gurobi.c_column(JuMP.backend(owner_model(x)), index(x))
function ff(a, b, c, d)
    f(x) = sin(a * pi * x) + b * x + (c * x)^2 + d
    return f
end

function tikz_string(pwl)
    return "\\addplot[color = blue, thick, domain = $(pwl.xMin) : $(pwl.xMax)]{$(pwl.a) * x + $(pwl.b)};"
end

F = [ff(4, -0.4, 2.4, 5), ff(4, -0.4, 2.9, 4)]

function solve_example(pwl, model, x)::Vector{Float64}
    for t in 1 : 2
        breakx, breaky = IterativePWLB.CCGetBreakpoints(pwl[t], 1e-7)
        nbreaks = length(breakx)
        error = GRBsetpwlobj(JuMP.backend(model), column(x[t]), nbreaks, breakx, breaky)
    end
    optimize!(model)
    xval = value.(x)
    return xval
end

function tighten_pwl(f, pwl, x, eps)
    pieces = Int64[]
    for (i, p) in enumerate(pwl)
        if x >= p.xMin - 1e-9 && x <= p.xMax + 1e-9
            push!(pieces, i)
        end
    end
    first = pieces[begin]
    last = pieces[end]
    xmin = pwl[first].xMin
    xmax = pwl[last].xMax
    epsf = minimum(eps[pieces]) / 2
    newpwl = Linearize(f, xmin, xmax, LinA.Relative(epsf), LinA.ExactLin(); bounding = Under())
    neweps = [epsf for i in eachindex(newpwl)]
    pwlret = vcat(pwl[1 : first - 1], newpwl, pwl[last + 1 : end])
    epsret = vcat(eps[1 : first - 1], neweps, eps[last + 1 : end])
    return pwlret, epsret
end

function solve_example(maxits = 4)
    model = direct_model(Gurobi.Optimizer())
    @variable(model, 0 <= x[1:2] <= 1)
    @constraint(model, 2* x[1] + x[2] >= 1)
    @constraint(model, 2 * x[1] + 5 * x[2] <= 4)
    @constraint(model, x[2] >= 0.3)
    
    pwl = [Linearize(F[i], 0, 1, LinA.Relative(20.0), LinA.ExactLin(); bounding = Under()) for i in 1 : 2]
    eps = [[20.0 for i in eachindex(p)] for p in pwl]
    for nit in 1 : maxits
        println("########### iteration $nit ##############")
        println("pwl[1]:")
        println.(tikz_string.(pwl[1]))
        println("pwl[2]:")
        println.(tikz_string.(pwl[2]))
        xopt = solve_example(pwl, model, x)
        println("opt = $xopt")
        for t in 1 : 2
            if F[t](xopt[t]) - pwl[t](xopt[t]) > 0.0001 * F[t](xopt[t])
                pwl[t], eps[t] = tighten_pwl(F[t], pwl[t], xopt[t], eps[t])
            end
        end
    end
end