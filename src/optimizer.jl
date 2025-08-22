function cplex_parameters(currlb, cutoff, tilim)::Dict{String, Any}
    return  Dict(
        "CPXPARAM_Threads" => 1,
        "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
        "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 0.0,
        # "CPXPARAM_MIPEmphasis" => 3,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPRELAX" => 1e-9,
        "CPX_PARAM_EPOPT" => 1e-9,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_CUTUP" => (1 + 1e-5) * cutoff,
        "CPX_PARAM_TILIM" => tilim + 1
    )
end

function gurobi_parameters(currlb, cutoff, tilim, warmstart, printint)::Dict{String, Any}
    if warmstart
        return  Dict(
            # "BestObjStop" => currlb - 1e-5,
            # "Cutoff" => cutoff > 0 ? (1 + 1e-5) * cutoff : (1 - 1e-5) * cutoff,
            "TimeLimit" => tilim + 1
        )
    else
        return  Dict(
            "OutputFlag" => printint,
            "Threads" => 1,
            "MIPGap" => 1e-9,
            "MIPGapAbs" => 1e-9, 
            "FeasibilityTol" => 1e-9,
            "OptimalityTol" => 1e-9,
            # "BestObjStop" => currlb - 1e-5,
            "IntFeasTol" => 1e-9,
            # "Cutoff" => cutoff > 0 ? (1 + 1e-5) * cutoff : (1 - 1e-5) * cutoff,
            "TimeLimit" => tilim + 1
        )
    end
end

function empty_parameters(currlb, cutoff, tilim)::Dict{String, Any}
    return Dict{String, Any}()
end
