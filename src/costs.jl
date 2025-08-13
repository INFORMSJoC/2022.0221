function ngueveu2019_x2_function(x)
    return x^2
end

function ngueveu2019_x2_expression()
    return :(x^2)
end

function log_expression()
    return :(log(x + 1) + 1)
end
function log_function(x)
    return log(x + 1) + 1
end

############################################################

function FC_R_expression()
    return :(0.0000002*x^5 - 0.0000274*x^4 + 0.00151450*x^3 - 0.02453270*x^2 +
1.92434870*x + 5.90568630)
end
function FC_R_function(x)
    return 0.0000002*x^5 - 0.0000274*x^4 + 0.00151450*x^3 - 0.02453270*x^2 +
1.92434870*x + 5.90568630
end
function FC_R_functor_function()::Function
    return FC_R_function
end

############################################################

function sincos_expression(Q, C)
    return :(C * (7.5*sin(pi*(x * 100 / Q -10.0)/40.0) - 15.0*cos(pi*(x * 100 / Q -10.0)/80.0) + 19.5) / 38.66)
end
function sincos_function(Q, C, x)
    return C * (7.5*sin(pi*(x * 100 / Q -10.0)/40.0) - 15.0*cos(pi*(x * 100 / Q -10.0)/80.0) + 19.5) / 38.66
end

function sincos_functor_function(Q, C)::Function
    fn(q) = sincos_function(Q, C, q)
    return fn
end

function sincos_function_string(Q, C)::String
    return "$C * (7.5*sin(pi*(x * 100 / $Q - 10.0) / 40.0) - 15.0 * cos(pi * (x * 100 / $Q -10.0) / 80.0) + 19.5) / 38.66"
end

############################################################
function sinabc_coefs(k, Q = 100, C = 100)
    if k == 1
        a = 15.0 / 54.68125 ; b = 2.0; c = 1.0
    elseif k == 2
        a = 25.0 / 408.25733; b = 5.0; c  = 5.0
    elseif k == 3
        a = 25.0 / 638.779 ; b = 10.0; c  = 5.0
    else
        throw(ArgumentError("the parameter k of sinabc function must be in [1, 2, 3], so it can not be equal to $k !!!!!!"))
    end
    [a * C, b / Q, c / Q]
end

function sinabc_function(x, q)
    return x[1]*(( sin(x[2]*q) + x[3]*q )^2 )
end
function sinabc_expression(y)
    return :($(y[1])*( ( sin($(y[2])*x) + $(y[3])*x )^2 ))
end

function sinabc_functor_function(x)::Function
    sinabc_fn(q) = sinabc_function(x, q)
    return sinabc_fn
end

function sinabc_function_string(x)::String
    return "$(x[1]) * ((sin($(x[2]) * x) + $(x[3]) * x)^2)"
end

############################################################

function FC_A1_expression()
    return :(0.001*x^3 - 0.024*x^2 + 1.92*x + 5.91)
end
function FC_A1_function(x)
    return 0.001*x^3 - 0.024*x^2 + 1.92*x + 5.91
end

function FC_A1_functor_function()::Function
#    FC_A1_fn(x) = FC_A1_function(x)
#   return FC_A1_fn
    return FC_A1_function
end

function FC_A2_expression()
    return :(-0.005*x^3 + 0.5*x^2 - 0.8*x + 10.0)
end
function FC_A2_function(x)
    return -0.005*x^3 + 0.5*x^2 - 0.8*x + 10.0
end
function FC_A2_functor_function()::Function
    return FC_A2_function
end


function cubic_coefs(x::Array{Float64, 1}, y::Array{Float64, 1})
    Q = maximum(x)
    C = maximum(y)
    if C > EPS_NZ
        function f!(F, w)
            for i in 1 : 4
                F[i] = cubic_function(w, x[i] / Q) - (y[i] / C)
            end
        end
        # println("finding cubic coefs for x = $x, y = $y")
        zeros = nlsolve(f!, [1.0, 1.0, 1.0, 1.0]).zero
    else
        zeros = [0, 0, 0, 0]
    end
    # println("found zeros $zeros")
    return C * [zeros[1] / (Q^3), zeros[2] / (Q^2), zeros[3] / Q, zeros[4]]
    # push!(zeros, y[1])
    # zeros
    # zeros = round.(Int64, zeros * 1000) / 1000
    # p = CubicPoly(zeros)
    # zeros
end

function zero_function()
    return 0.0
end

function zero_expression()
    return :(0.0 * x)
end

function zero_functor_function()::Function
    zero_fn(q) = zero_function()
    return zero_fn
end

function cubic_function(x, q)
    return x[1] * q^3 + x[2] * q^2 + x[3] * q + x[4]
end
function cubic_expression(y)
    return :($(y[1]) * x^3 + $(y[2]) * x^2 + $(y[3]) * x^1 + $(y[4]))
end

function cubic_functor_function(x)::Function
    cubic_fn(q) = cubic_function(x, q)
    return cubic_fn
end

function cubic_function_string(x):String
    return "$(x[1]) * x^3 + $(x[2]) * x^2 + $(x[3]) * x + $(x[4])"
end

function exp_coefs(x::Array{Float64, 1}, y::Array{Float64, 1})
    k = y[1]
    c = y[2] - y[1]
    a = 10 / x[2]
    b = x[2] / 2
    d = -x[2] / 10
    [a, b, c, d, k]
end

function exp_function(x, q)
    return x[5] + x[3] / (1 + x[2] * exp(-x[1] * (q + x[4])))
end
function exp_expression(y)
    return :($(y[5]) + $(y[3]) / (1 + $(y[2]) * exp(-$(y[1]) * (x + $(y[4])))))
end

function exp_functor_function(x)::Function
    exp_fn(q) = exp_function(x, q)
    return exp_fn
end

function exp_function_string(x)::String
    return "$(x[5]) + $(x[3]) / (1 + $(x[2]) * e^(-$(x[1]) * (x + $(x[4]))))"
end

function (p::CubicPoly)(q::Float64)
    p.coefs[1] * (q - p.coefs[2])^3 + p.coefs[3] * (q - p.coefs[2]) + p.coefs[4]
end

# -(p::CubicPoly) =

function real_function(q)
    return ((484360120612810239*q^16)/212676479325586539664609129644855132160 -
            (55975000056945233*q^15)/51922968585348276285304963292200960 +
            (374578569900980301*q^14)/1622592768292133633915780102881280 -
            (374961211368334637*q^13)/12676506002282294014967032053760 +
            (100147218409756427*q^12)/39614081257132168796771975168 -
            (235467440160020499*q^11)/1547425049106725343623905280 +
            (642510758331249011*q^10)/96714065569170333976494080 -
            (161525572054452217*q^9)/755578637259143234191360 +
            (482916429491142721*q^8)/94447329657392904273920 -
            (269581733069483441*q^7)/2951479051793528258560 +
            (5646351665002773*q^6)/4611686018427387904 -
            (71461922790103043*q^5)/5764607523034234880 +
            (137606181927263737*q^4)/1441151880758558720 -
            (101741333933636769*q^3)/180143985094819840 +
            (471783432168476983*q^2)/180143985094819840 +
            (298634054576571277*q)/23058430092136939520 - 6/5)
end

function real_expression()
    return :((484360120612810239*q^16)/212676479325586539664609129644855132160 -
            (55975000056945233*q^15)/51922968585348276285304963292200960 +
            (374578569900980301*q^14)/1622592768292133633915780102881280 -
            (374961211368334637*q^13)/12676506002282294014967032053760 +
            (100147218409756427*q^12)/39614081257132168796771975168 -
            (235467440160020499*q^11)/1547425049106725343623905280 +
            (642510758331249011*q^10)/96714065569170333976494080 -
            (161525572054452217*q^9)/755578637259143234191360 +
            (482916429491142721*q^8)/94447329657392904273920 -
            (269581733069483441*q^7)/2951479051793528258560 +
            (5646351665002773*q^6)/4611686018427387904 -
            (71461922790103043*q^5)/5764607523034234880 +
            (137606181927263737*q^4)/1441151880758558720 -
            (101741333933636769*q^3)/180143985094819840 +
            (471783432168476983*q^2)/180143985094819840 +
            (298634054576571277*q)/23058430092136939520 - 6/5)
end

function trigonometric_function(x)
    return 0.01 * (x  + x^2 * (1 + sin(3 * x / 150)) + 5000)
end

function trigonometric_expression()
    # println(:(x  + x^2 * (1 + sin(3 * x / D)) + 1000))
    return :(0.01 * (x  + x^2 * (1 + sin(3 * x / 150)) + 5000))
end


function congestioncost_functions(x, q)
    # return min( x[2]*(q+q^(x[1]+1.0)/(x[3])^x[1]) , x[5]+x[2]*(q+(q^(x[1]+1.0)/(x[3] + x[4])^x[1])) )
    return min(congestioncost_function1(x, q), congestioncost_function2(x, q))
    #minimum( D*(x+x^(alpha+1.0)/(initCapa)^alpha) , incremCosts6+D*(x+(x^(alpha+1.0)/(initCapa + incremCapa)^alpha)) )
    """
    alpha = coefs[1]
    D = coefs[2]
    initCapa = coefs[3]
    incremCapa = coefs[4]
    incremCost6  = coefs[5]
    capaMax = initCapacity + incremCapacity
    """
    #return minimum( D*(x+x^(alpha+1.0)/(initCapa)^alpha) , incremCosts6+D*(x+(x^(alpha+1.0)/(initCapa + incremCapa)^alpha)) )
    #v_g = infoNode[v,:D]*(x+x^(alpha+1.0)/(infoNode[v,:InitCapacity])^alpha)
    #v_gu = infoNode[v,:IncremCosts6]+  infoNode[v,:D]*(x+(x^(alpha+1.0)/v_capaMax[v]^alpha))
end

function congestioncost_function1(x, q)
    return x[2]*(q+q^(x[1]+1.0)/(x[3])^x[1])
    #return min( x[2]*(q+q^(x[1]+1.0)/(x[3])^x[1]) , x[5]+x[2]*(q+(q^(x[1]+1.0)/(x[3] + x[4])^x[1])) )
    """
    alpha = coefs[1]
    D = coefs[2]
    initCapa = coefs[3]
    incremCapa = coefs[4]
    incremCost6  = coefs[5]
    capaMax = initCapacity + incremCapacity
    """
end

function congestioncost_function2(x, q)
    return x[5]+x[2]*(q+(q^(x[1]+1.0)/(x[3] + x[4])^x[1]))
    #return min( x[2]*(q+q^(x[1]+1.0)/(x[3])^x[1]) , x[5]+x[2]*(q+(q^(x[1]+1.0)/(x[3] + x[4])^x[1])) )
    """
    alpha = coefs[1]
    D = coefs[2]
    initCapa = coefs[3]
    incremCapa = coefs[4]
    incremCost6  = coefs[5]
    capaMax = initCapacity + incremCapacity
    """
end

function congestioncost_functor_functions(x)::Function
    congestioncost_fns(q) = congestioncost_functions(x,q)
    return congestioncost_fns
end
function congestioncost_functor_function1(x)::Function
    congestioncost_fn1(q) = congestioncost_function1(x,q)
    return congestioncost_fn1
end
function congestioncost_functor_function2(x)::Function
    congestioncost_fn2(q) = congestioncost_function2(x,q)
    return congestioncost_fn2
end

function congestioncost_expressions(y)
    return  :(min( $(y[2])*(x+x^($(y[1])+1.0)/($(y[3]))^$(y[1])) , $(y[5])+$(y[2])*(x+(x^($(y[1])+1.0)/($(y[3]) + $(y[4]))^$(y[1]))) )),
            :( $(y[2])*(x+x^($(y[1])+1.0)/($(y[3]))^$(y[1])) ) ,
            :( $(y[5])+$(y[2])*(x+(x^($(y[1])+1.0)/($(y[3]) + $(y[4]))^$(y[1]))) )
    #return :(min( $(y[2])*(x+x^($(y[1])+1.0)/($(y[3]))^$(y[1])) , $(y[5])+$(y[2])*(x+(x^($(y[1])+1.0)/($(y[3]) + $(y[4]))^$(y[1]))) ))
end

function cubic_function_x_multipliers(maxx)
    return maxx * [0, 1/3, 2/3, 1]
end

function cubic_function_y_multipliers(nomcost)
    return nomcost * [1/2, 3/4, 7/8, 2] / 2
end

function cubic_function_multipliers(maxx, nomcost)
    cubic_function_x_multipliers(maxx), cubic_function_y_multipliers(nomcost)
end

function exp_function_x_multipliers(maxx)
    return [0.0, maxx]
end

function exp_function_y_multipliers(nomcost)
    return nomcost * [1/2, 2]
end

function exp_function_multipliers(maxx, nomcost)
    exp_function_x_multipliers(maxx), exp_function_y_multipliers(nomcost)
end


#=
function knapcost_functions(x, q)
    return x[2]*(q+q^(x[1]+1.0)/(x[3])^x[1])
end

function knapcost_functor_functions(x)::Function
    knapcost_fns(q) = knapcost_functions(x,q)
    return knapcost_fns
end

=#

############### Parametrized divex function ##################
function parametric_divexp_coefs(Q = 100, C = 100)
    a = 0.152216 * 100 / Q
    b = 44.742557
    c = 72.273751 * C / 72.25368996908502
    d = 21.229801 * Q / 100
    [-a, b, c, -d]
end

function parametric_divexp_function(x, q)
    return x[3] / (1 + (x[2] * exp(x[1] * (q + x[4]))))
end
function parametric_divexp_expression(y)
    return :($(y[3]) / (1 + ($(y[2]) * exp($(y[1]) * (x + $(y[4]))))))
end

function parametric_divexp_functor_function(x)::Function
    fn(q) = parametric_divexp_function(x, q)
    return fn
end

function parametric_divexp_function_string(x)::String
    return "$(x[3]) / (1 + ($(x[2]) * e^($(x[1]) * (x + $(x[4])))))"
end

##############################################################
############### Logarithmic function ##################
function logfn_coefs(Q = 100, C = 100)
    a = 1.0
    b = (MathConstants.e^3 - 1) / Q
    c = C * 1.0 / 3
    [a, b, c]
end

function logfn_function(x, q)
    return x[3] * log(x[1] + x[2] * q)
end
function logfn_expression(y)
    return :($(y[3]) * log($(y[1]) + $(y[2]) * x))
end

function logfn_functor_function(x)::Function
    fn(q) = logfn_function(x, q)
    return fn
end

function logfn_function_string(x)::String
    return "$(x[3]) * log($(x[1]) + $(x[2]) * x)"
end 
##############################################################
############### Square-root function ##################
function sqrt_coefs(Q = 100, C = 100)
    a = 1 / Q
    b = C
    [a, b]
end

function sqrt_function(x, q)
    return x[2] * sqrt(x[1] * q)
end
function sqrt_expression(y)
    return :($(y[2]) * sqrt($(y[1]) * x))
end

function sqrt_functor_function(x)::Function
    fn(q) = sqrt_function(x, q)
    return fn
end

function sqrt_function_string(x)::String
    return "$(x[2]) * sqrt($(x[1]) * x)"
end

##############################################################
############### x^2 function ##################
function square_coef(Q = 100, C = 100)
    return (C / Q) / Q
end

function square_function(x, q)
    return x * q * q
end
function square_expression(y)
    return :($(y) * x * x)
end

function square_functor_function(x)::Function
    fn(q) = square_function(x, q)
    return fn
end

function square_function_string(x)::String
    return "$x * x^2"
end
##############################################################