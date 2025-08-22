function write_ampl(pdata::Data, costs, filename)
    if typeof(pdata) == IterativePWLB.TranspProblemData
        write_ampl_transportation(pdata, costs, filename)
    end
end

function write_ampl_transportation(pdata::TranspProblemData, costs, filename)
    open(filename, "w") do f
        println(f, "set ORIGINS;")
        println(f, "set DESTINATIONS;")
        println(f)
        println(f, "param cubic_coef {ORIGINS, DESTINATIONS};")
        println(f, "param quadr_coef {ORIGINS, DESTINATIONS};")
        println(f, "param linear_coef {ORIGINS, DESTINATIONS};")
        println(f, "param const_coef {ORIGINS, DESTINATIONS};")
        println(f)
        println(f, "param offer {ORIGINS} >= 0;")
        println(f, "param demand {DESTINATIONS} >= 0;")
        println(f, "param maxflow {i in ORIGINS, j in DESTINATIONS} = min{offer[i], demand[j]};")
        println(f)
        println(f, "var x {ORIGINS, DESTINATIONS} >= 0;")
        println(f, "var y {ORIGINS, DESTINATIONS} binary;")
        println(f)
        println(f, "minimize cost: sum {i in ORIGINS, j in DESTINATIONS}(cubic_coef[i, j] * x[i, j]^3 + quadr_coef[i, j] * x[i, j]^2 + linear_coef[i, j] * x[i, j] + const_coef[i, j] * y[i, j]);")
        println(f)
        println(f, "subject to Offer {i in ORIGINS}:")
        println(f, "\tsum {j in DESTINATIONS} x[i, j] = offer[i];")
        println(f)
        println(f, "subject to Demand {j in DESTINATIONS}:")
        println(f, "\tsum {i in ORIGINS} x[i, j] = demand[j];")
        println(f)
        println(f, "subject to LINKING {i in ORIGINS, j in DESTINATIONS}:")
        println(f, "\tx[i, j] - maxflow[i, j] * y[i, j] <= 0;")
    end
end