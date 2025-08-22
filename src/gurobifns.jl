
# function get_nonlinear_constraint_expression(model::JuMP.Model, cf::Vector{Float64}, x::JuMP.VariableRef, y::JuMP.VariableRef, ctype::Symbol, sense::Symbol)
#     mult = sense == :MIN ? 1.0 : -1.0
#     if ctype == :cubic
#         return @NLexpression(model, mult * (cf[1] * x^3 + cf[2] * x^2 + cf[3] * x + cf[4] * y)), 0.0
#     elseif ctype == :exp
#         f0 = cf[5] + cf[3] / (1 + cf[2] * exp(-cf[1] * cf[4]))
#         return @NLexpression(model, mult * (y * f0 + cf[3] / (1 + cf[2] * exp(-cf[1] * (x + cf[4]))))), mult * (f0 - cf[5])
#     elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
#         return @NLexpression(model, mult * (cf[1] * ((sin(cf[2] * x) + cf[3] * x)^2))), 0.0
#     elseif ctype == :sincos
#         return @NLexpression(model, mult * (cf[2] * (7.5*sin(pi*(x * 100 / cf[1] - 10.0)/40.0) - 15.0*cos(pi*(x * 100 / cf[1] -10.0)/80.0) + 19.5) / 38.66)), 0.0
#     elseif ctype in [:divexp, :pardivexp] #TODO: valider que je n'ai pas mis une parenthese au mauvais endroit
#         f0 = (cf[3])/(1 + (cf[2]) * (exp(cf[1] * cf[4])))
#         println("f0 = $f0")
#         return @NLexpression(model, mult * (y * f0 + cf[3]/(1 + (cf[2]) * (exp(cf[1] * (x + cf[4])))))), mult * f0
#     elseif ctype == :log
#         return @NLexpression(model, mult * cf[3] * log(cf[1] + cf[2] * x)), 0.0
#     elseif ctype == :sqrt
#         return @NLexpression(model, mult * cf[2] * ((cf[1] * x)^(0.5))), 0.0
#     elseif ctype == :square
#         return @NLexpression(model, mult * cf[1] * (x^2)), 0.0
#     end
# end


grbcolumn(x::VariableRef) = Gurobi.c_column(backend(owner_model(x)), index(x))

function get_nlgurobi_constraint_expression(model::JuMP.Model, cf::Vector{Float64}, x::JuMP.VariableRef, y::JuMP.VariableRef, ctype::Symbol, sense::Symbol)
    grbx = grbcolumn(x)
    grby = grbcolumn(y)
    mult = sense == :MIN ? 1.0 : -1.0
    if ctype == :cubic
        p = [cf[1], cf[2], cf[3], 0.0]
        z = @variable(model)
        grbz = grbcolumn(z)
        GRBaddgenconstrPoly(backend(model), "poly", grbx, grbz, 4, p, "")
        return mult * (z + cf[4] * y), 0.0 
    elseif ctype == :exp
        f0 = cf[5] + cf[3] / (1 + cf[2] * exp(-cf[1] * cf[4]))
        z1 = @variable(model)
        @constraint(model, z1 + cf[1] * x == -cf[1] * cf[4])
        z2 = @variable(model)
        GRBaddgenconstrExp(backend(model), "exp", grbcolumn(z1), grbcolumn(z2), "")
        z3 = @variable(model, lower_bound = 0.0)
        @constraint(model, z3 - cf[2] * z2 == 1)
        z4 = @variable(model)
        GRBaddgenconstrPow(backend(model), "pow", grbcolumn(z3), grbcolumn(z4), -1.0, "")
        return mult * (f0 * y + cf[3] * z4), mult * (f0 - cf[5])
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        z1 = @variable(model)
        @constraint(model, z1 - cf[2] * x == 0)
        z2 = @variable(model, lower_bound = -1.0, upper_bound = 1.0)
        GRBaddgenconstrSin(backend(model), "sin", grbcolumn(z1), grbcolumn(z2), "")
        z3 = @variable(model)
        @constraint(model, z3 - z2 - cf[3] * x == 0)
        z4 = @variable(model, lower_bound = 0.0)
        GRBaddgenconstrPow(backend(model), "sq", grbcolumn(z3), grbcolumn(z4), 2.0, "")
        return mult * (cf[1] * z4), 0.0
    elseif ctype == :sincos
        z1 = @variable(model)
        @constraint(model, 40 * z1 - (100 * pi / cf[1]) * x == -10 * pi)
        z2 = @variable(model)
        @constraint(model, 80 * z2 - (100 * pi / cf[1]) * x == -10 * pi)
        z3 = @variable(model, lower_bound = -1.0, upper_bound = 1.0)
        z4 = @variable(model, lower_bound = -1.0, upper_bound = 1.0)
        GRBaddgenconstrSin(backend(model), "sin", grbcolumn(z1), grbcolumn(z3), "")
        GRBaddgenconstrCos(backend(model), "cos", grbcolumn(z2), grbcolumn(z4), "")
        return mult * cf[2] * (7.5 * z3 - 15.0 * z4 + 19.5) / 38.66, 0.0
    elseif ctype in [:divexp, :pardivexp]
        f0 = cf[3] / (1 + cf[2] * exp(cf[1] * cf[4]))
        z1 = @variable(model)
        @constraint(model, z1 - cf[1] * x == cf[1] * cf[4])
        z2 = @variable(model, lower_bound = 0.0)
        GRBaddgenconstrExp(backend(model), "exp", grbcolumn(z1), grbcolumn(z2), "")
        z3 = @variable(model)
        @constraint(model, z3 - cf[2] * z2 == 1.0)
        z4 = @variable(model)
        GRBaddgenconstrPow(backend(model), "inv", grbcolumn(z3), grbcolumn(z4), -1.0, "")
        return mult * (f0 * y + cf[3] * z4), mult * f0
    elseif ctype == :log
        z1 = @variable(model, lower_bound = 0.0)
        @constraint(model, z1 - cf[2] * x == cf[1])
        z2 = @variable(model)
        GRBaddgenconstrLog(backend(model), "log", grbcolumn(z1), grbcolumn(z2), "")
        return mult * cf[3] * z2, 0.0
    elseif ctype == :sqrt
        z1 = @variable(model, lower_bound = 0.0)
        @constraint(model, z1 - cf[1] * x == 0)
        z2 = @variable(model)
        GRBaddgenconstrPow(backend(model), "log", grbcolumn(z1), grbcolumn(z2), 0.5, "")
        return mult * cf[2] * z2, 0.0
    elseif ctype == :square
        z = @variable(model)
        GRBaddgenconstrPow(backend(model), "log", grbcolumn(x), grbcolumn(z), 2.0, "")
        return mult * cf[1] * z, 0.0
    end
end