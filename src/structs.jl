import Base.-

# struct CommonParameters
#     ptype::Symbol
#     eps::Float64
#     with_solver::Any
#     with_solver_parameters::Dict{String, Any}
#     with_lina_solver::LinA.AbstractAlgorithm
#     eps0::Float64
#     tilim::Float64
#     printon::Bool
#     stepdiv::Float64
#     dynstep::Bool
#     activeUpdateLargePwlf::Bool
#     warmstartMILP::Bool
#     miptype::Symbol
#     with_fixed_charge::Bool
# end

struct PWLInfo
    lp::LinA.LinearPiece
    eps::Float64
end

abstract type Data end

mutable struct Parameters
    with_solver::DataType# = GLPK.Optimizer,
    with_solver_parameters::Any# = empty_parameters,
    with_lina_solver::DataType# = LinA.HeuristicLin,
    with_cost_type::Symbol
    ptype::Symbol
    eps0::Float64# = 1e-1,
    eps::Float64
    tilim::Float64 # = 3600,
    printon::Bool # =false,
    stepdiv::Float64# =2.0,
    dynstep::Bool # =false,
    activeUpdateLargePwlf::Bool # =false,
    warmstartMILP::Bool # =false,
    miptype::Symbol # =:grb,
    with_fixed_charge::Bool # = false
end

struct HEVData <: Data
    Pdem::Array{Int64, 1}
    nbinstants::Int64
    losscoefSEin::Float64
    losscoefSEout::Float64
    p_FCmax::Float64
    p_FCmin::Float64
    p_ESmax::Float64
    p_ESmin::Float64
    ESmin::Float64
    ESmax::Float64
    ESinit::Float64
end

struct KnapData <: Data
    are_vars_integer::Bool # input parameter that specifies if we should consider
    # integer variables or not (integer or real variables) 
    minuslistea::Array{Float64, 1}
    object_coefb::Array{Float64, 1}
    object_coefc::Array{Float64, 1}
    object_coefd::Array{Float64, 1}
    object_U::Array{Float64, 1}
    knapsack_cap::Array{Float64, 1}
    object_weight::Array{Float64, 1}
end


struct MCNDData <: Data
    nodes::Array{NamedTuple{(:ref, :fixedcost, :capainit, :capaincrement, :delaycost, :freeflowdelay, :croisement),Tuple{Int64, Float64,Int64,Int64,Float64,Float64,Float64}},1}
    arcs::Array{NamedTuple{(:ref, :from, :to, :fixedcost, :varcost, :capa),Tuple{Int64,Int64,Int64,Float64,Float64, Int64}},1}
    commods::Array{NamedTuple{(:ref, :from, :to, :qty),Tuple{Int64,Int64,Int64,Int64}},1}
    listofinputarcs
    listofoutputarcs
    nbNode::Int64
    nbArc::Int64
    nbCommod::Int64
    alpha::Float64
    cubicexpcoefs::Array{NamedTuple{(:a1, :a2, :a3, :a4, :a5),Tuple{Float64, Float64,Float64,Float64,Float64}},1}
end

struct UFLPData <: Data
    fixed_costs::Array{Float64, 1}
    caps::Array{Int64, 1}
    demands::Array{Int64, 1}
    asscosts::Array{Float64, 2}
    capcons::Bool
end

struct SchedulingData <: Data
    ntasks::Int64
    durations::Array{Int64, 1}
    demands::Array{Int64, 1}
    time_windows::Array{Tuple{Int64, Int64}, 1}
end

struct TranspProblemData <: Data
    caps::Array{Int64, 1}
    demands::Array{Int64, 1}
    asscosts::Array{Int64, 2}
end

struct CubicPoly
    coefs::Array{Float64, 1}
end

const EPS_NZ = 1e-6
# -(f::Function) = x -> -f(x)

Parameters() = Parameters(GLPK.Optimizer,
                            empty_parameters,
                            LinA.HeuristicLin,
                            :cubic,
                            :transp,
                            1e+1,
                            1e-1,
                            3600.0,
                            true,
                            2.0,
                            false,
                            false,
                            false,
                            :glpk,
                            false)

(pwl::Vector{LinA.LinearPiece})(x::Float64) = pwl(x, LinA.Under(), EPS_NZ)
(pwl::Vector{PWLInfo})(x::Float64) = [p.lp for p in pwl](x)
compute(f::Function, x::Float64, with_fixed_charge::Bool) = with_fixed_charge && x < EPS_NZ ? 0.0 : f(x)
