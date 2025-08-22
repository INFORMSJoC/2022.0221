module IterativePWLB
    using NLsolve
    using Optim
    using JuMP
    using MathOptInterface
    using GLPK
    using CSV
    using ForwardDiff
    using Pkg
    using Printf
    using LinA
    using Dates
    using SCIP
    using Alpine
    using CPLEX
    using Ipopt
    using Juniper
    using Gurobi
    using IntervalArithmetic
    using IntervalRootFinding
    using Random
    using AmplNLWriter
    using Couenne_jll
    using Dates
    using LinearAlgebra
    using PrecompileTools

    const MOI = MathOptInterface

    include("structs.jl")
    include("amplwriter.jl")
    include("amplutils.jl")
    include("scminlputils.jl")
    include("costs.jl")
    include("hev.jl")
    include("instance_gen.jl")
    include("instance_reader.jl")
    include("gurobifns.jl")
    include("minlpsolvers.jl")
    include("nlsolver.jl")
    include("optimizer.jl")
    include("util.jl")
    include("knap.jl")
    include("uflp.jl")
    include("transp.jl")
    include("mcf.jl")
    include("mcnd.jl")
    include("scheduling.jl")
    include("solve.jl")
    include("mip.jl")
    include("ipwlb.jl")
    
    export solve_hev, solve_knap, solve_mcnd_or_mcf_from_file, solve_scheduling_from_file, solve_uflp_from_file, solve_transp_random, solve_transp_from_file, solve_transp
    export setLinaSolver
    export nl_knap, nl_hev, nl_mcf, nl_mcnd, nl_scheduling, nl_transp, nl_uflp_nlass, nl_uflp_nldep ### SUNSUN update here

    @setup_workload begin
        home = pkgdir(IterativePWLB)
        datadir = joinpath([home, "data"])
        @compile_workload begin
            params = Parameters(Gurobi.Optimizer,
                                gurobi_parameters,
                                LinA.ExactLin,
                                :cubic,
                                :transp,
                                1e+1,
                                1e+0,
                                60.0,
                                false,
                                2.0,
                                false,
                                false,
                                true,
                                :grb,
                                false)
            for ctype in [:cubic]#:sincos, :sinabc2, :sqrt, :square, :log, :exp, :cubic, :divexp, :pardivexp, :sinabc1, :sinabc3]
                params.with_cost_type = ctype
                with_fixed_charge = ctype in [:exp, :cubic]
                println(stdout, "precompiling for cost of type $(string(ctype))")
                if ctype != :divexp    
                    for ptype in ["mcnd_nodes", "mcf_nodes", "mcnd_arcs", "mcf_arcs"]
                        params.ptype = Symbol(ptype)
                        println("doing $ptype")
                        redirect_stdout(devnull) do
                            solve_mcnd_or_mcf_from_file(joinpath([datadir, "test-mcnd.csv"]),
                                                    joinpath([datadir, "test-mcndCommod.csv"]),
                                                    joinpath([datadir, "nodes_test-mcnd_03.csv"]),
                                                    6, params)
                        end
                    end
                end
                 if ctype != :divexp && ctype != :twofuncs
                    for ptype in [:uflp_nldep, :uflp_nlass]
                        params.ptype = ptype
                        println("doing $ptype")
                        redirect_stdout(devnull) do
                            solve_uflp_from_file(joinpath([datadir, "holmberg-p1-5"]),
                            params)
                        end
                    end
                end
                if ctype != :twofuncs
                    println("doing knapsack")
                    params.ptype = :knap
                    redirect_stdout(devnull) do
                        solve_knap(joinpath([datadir, "knap.dat"]),
                                    false,
                                    params)
                    end
                end
                if ctype != :divexp && ctype != :twofuncs
                    params.ptype = :transp
                    println("doing transp")
                    redirect_stdout(devnull) do
                        solve_transp_from_file(joinpath([datadir, "p_2_1"]), params)
                    end
                end    
            end
        end
    end

end # module
