function get_juniper()
    return optimizer_with_attributes(Juniper.Optimizer,
                                    MOI.Silent() => true,
					"processors" => 1,
                                    "mip_solver" => get_gurobi_alpine(),
                                    "nl_solver" => get_ipopt())
end

function get_cplex()
    return optimizer_with_attributes(CPLEX.Optimizer,
                                    MOI.Silent() => true,
                                    "CPX_PARAM_PREIND" => 1,
                                    "CPX_PARAM_THREADS" => 1)
end

function get_ipopt()
    return optimizer_with_attributes(Ipopt.Optimizer,
                                    MOI.Silent() => true,
                                    "sb" => "yes",
                                    "max_iter" => Int(1E4))
end

function get_gurobi_alpine()
	return optimizer_with_attributes(Gurobi.Optimizer,
                                    MOI.Silent() => true,
                                    "Presolve" => 1)
                                    # "max_iter" => Int(1E4))
end

function get_scip()
    return SCIP.Optimizer
end

function get_alpine()
    return Alpine.Optimizer
end

function minlp_scip_params(reps = 1e-1, tilim = 3600)
    return Dict("limits/gap" => reps,
                "limits/time" => tilim)
end

function minlp_alpine_params(reps = 1e-1, tilim = 3600)
    return Dict("nlp_solver" => get_ipopt(),
                "mip_solver" => get_gurobi_alpine(),
                "rel_gap" => reps,
                "minlp_solver" => get_juniper(),
		        "presolve_bt" => false,
                "time_limit" => tilim)
end

function get_amplnlwriter(solver)
    return AmplNLWriter.Optimizer(solver)
end

function minlp_bonmin_params(reps = 1e-1, tilim = 3600)
    return Dict("bonmin.time_limit" => tilim,
                "bonmin.allowable_fraction_gap" => reps)
end

function minlp_couenne_params(reps = 1e-1, tilim = 3600)
    return Dict()#"time_limit" => tilim,
                # "feas_tolerance" => reps)
end

function minlp_baron_params(reps = 1e-1, tilim = 3600, sumfile = "sumfile.lst")
    #return Dict("epsr" => reps / 100, "maxtime" => tilim, "sumfile" => sumfile, "lpsolver" => "cplex")
	return ["epsr = $(reps/100)", "maxtime = $tilim", "sumfile = $sumfile", "lpsolver = cplex", "barstats"]
end

function minlp_scminlp_params(tilim = 3600)
    return Dict("time" => tilim)
end

function get_baron(reps = 1e-1, 
                    tilim = 3600, 
                    licpath = "/opt/baron/23.11.13/baronlice.txt", 
                    cpxpath = "/opt/ibm/cplex/22.1.1/cplex/bin/x86-64_linux/libcplex2211.so")
    return BARON.Optimizer(EpsR = reps, 
                        EpsA = 1e-12,
                        AbsConFeasTol = 1e-12,
                        RelConFeasTol = 1e-12,
                        LicName = licpath, 
                        MaxTime = tilim,
                        CplexLibName = cpxpath)
end
