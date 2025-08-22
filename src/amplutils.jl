# FUNCTIONS TO WRITE AMPL MOD/DAT FILES

#------------------------------ COMMON FUNCTIONS -------------------------------#

function get_nonlinear_coefficients_mod(xMax::Real, cMax::Real, ctype::Symbol)::Vector{Float64}
    coefs = get_nonlinear_coefficients(xMax, cMax, ctype)
    if ctype == :cubic
        coefs[1] = coefs[1]^(1/3)
    end
    #return [abs(t) < 1e-15 ? 1e-15 : t for t in coefs]
    return coefs
end

function copyright_header(f::IOStream, filename::String)
    println(f, "#########################################################")
    println(f, "#  Claudio Contardo, Concordia U")
    println(f, "#  Sandra U. Ngueveu, INP Toulouse")
    println(f, "#  20 September 2023")
    println(f, "#  $filename")
    println(f, "#########################################################")
end

function nlobj_string_cubic(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)((COEF1$(istr) * $(xstr)$(istr))^3 + (COEF2$(istr) * $(xstr)$(istr)) * $(xstr)$(istr) + COEF3$(istr) * $(xstr)$(istr))"
    # return "nlobj$(istr) >= $(sign)($(xstr)$(istr) * (COEF3$(istr) + $(xstr)$(istr) * (COEF2$(istr) + $(xstr)$(istr) * COEF1$(istr))))"
end
function nlobj_string_exp(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF3$(istr) / (1 + COEF2$(istr) * exp(-COEF1$(istr) * ($(xstr)$(istr) + COEF4$(istr)))))"
end
function nlobj_string_sinabc(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF1$(istr) * ((sin(COEF2$(istr) * $(xstr)$(istr)) + COEF3$(istr) * $(xstr)$(istr))^2))"
end
function nlobj_string_sincos(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF2$(istr) * (7.5 * sin(PI * ($(xstr)$(istr) * 100 / COEF1$(istr) - 10) / 40) - 15 * cos(PI * ($(xstr)$(istr) * 100 / COEF1$(istr) - 10) / 80) + 19.5) / 38.66)"
end
function nlobj_string_divexp(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF3$(istr) / (1 + COEF2$(istr) * (exp(COEF1$(istr) * ($(xstr)$(istr) + COEF4$(istr))))))"
end
function nlobj_string_log(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF3$(istr) * log(COEF1$(istr) + COEF2$(istr) * $(xstr)$(istr)))"
end
function nlobj_string_sqrt(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF2$(istr) * ((COEF1$(istr) * $(xstr)$(istr))^(0.5)))"
end
function nlobj_string_square(istr::String, xstr::String, sign::String = "")::String
    return "nlobj$(istr) >= $(sign)(COEF1$(istr) * ($(xstr)$(istr)^2))"
end

function mincost_string_cubic(istr::String, ystr::String, sign::String = "+")::String
    return "(nlobj$(istr) $sign COEF4$(istr) * $(ystr)$(istr))"
end

function mincost_string_exp(istr::String, ystr::String, sign::String = "+")::String
    C6 = "COEF5$(istr) + COEF3$(istr) / (1 + COEF2$(istr) * exp(-COEF1$(istr) * COEF4$(istr)))"
    return "(nlobj$(istr) $sign ($(C6)) * ($(ystr)$(istr) - 1) $sign COEF5$(istr))"
end

function mincost_string_divexp(istr::String, ystr::String, sign::String = "+")::String
    C5 = "COEF3$(istr) / (1 + COEF2$(istr) * (exp(COEF1$(istr) * COEF4$(istr))))"
    return "(nlobj$(istr) $sign ($(C5)) * ($(ystr)$(istr) - 1))"
end

function mincost_string_other(istr::String, ystr::String, sign::String = "+")::String
    return "nlobj$(istr)"
end

function get_nl_condition(istr::String, ctype::Symbol)::String
	eps = "1e-15"
	if ctype == :cubic
		return " : abs(COEF1$(istr)) > $eps or abs(COEF2$(istr)) > $eps or abs(COEF3$(istr)) > $eps or abs(COEF4$(istr)) > $eps"
	elseif ctype == :exp
		return " : abs(COEF3$(istr)) > $eps"
	elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
		return " : abs(COEF1$(istr)) > $eps "
	elseif ctype == :sincos
		return " : abs(COEF2$(istr)) > $eps"
	elseif ctype in [:divexp, :pardivexp]
		return " : abs(COEF3$(istr)) > $eps"
	elseif ctype == :log
		return " : abs(COEF3$(istr)) > $eps"
	elseif ctype == :sqrt
		return " : abs(COEF2$(istr)) > $eps"
	elseif ctype == :square
		return " : abs(COEF1$(istr)) > $eps" 
	else return ""
	end
end

#--------------------------------------- END OF COMMON FUNCTIONS ---------------------------------------#

#--------------------------------------- TRANSPORTATION ----------------------------------#

function write_amplmod_transp(ctype, outfile)
    if ctype == :cubic
        write_amplmod_transp(outfile, 4, nlobj_string_cubic, mincost_string_cubic)
    elseif ctype == :exp
        write_amplmod_transp(outfile, 5, nlobj_string_exp, mincost_string_exp)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_transp(outfile, 3, nlobj_string_sinabc, mincost_string_other)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_transp(outfile, 4, nlobj_string_divexp, mincost_string_divexp)
    elseif ctype == :sincos
        write_amplmod_transp(outfile, 2, nlobj_string_sincos, mincost_string_other)
    elseif ctype == :log
        write_amplmod_transp(outfile, 3, nlobj_string_log, mincost_string_other)
    elseif ctype == :sqrt
        write_amplmod_transp(outfile, 2, nlobj_string_sqrt, mincost_string_other)
    elseif ctype == :square    
        write_amplmod_transp(outfile, 1, nlobj_string_square, mincost_string_other)
    end
end

function write_amplmod_transp_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nDeps;")
    println(f, "param nCust;")
    println(f)
    println(f, "param SUPPLY {1..nDeps};")
    println(f, "param DEMAND {1..nCust};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var x {i in 1..nDeps, j in 1..nCust} >= 0, <= min(SUPPLY[i], DEMAND[j]);")
    println(f, "var y {i in 1..nDeps, j in 1..nCust} binary;")
    println(f, "var nlobj {i in 1..nDeps, j in 1..nCust};")
    println(f)
    println(f, "subject to outflow {i in 1..nDeps}:")
    println(f, "\tsum {j in 1..nCust} x[i, j] <= SUPPLY[i];")
    println(f)
    println(f, "subject to inflow {j in 1..nCust}:")
    println(f, "\tsum {i in 1..nDeps} x[i, j] >= DEMAND[j];")
    println(f)
    println(f, "subject to link {i in 1..nDeps, j in 1..nCust}:")
    println(f, "\tx[i, j] - min(SUPPLY[i], DEMAND[j]) * y[i, j] <= 0;")
end


function write_amplmod_transp(outfile, ncoefs, nlobj_fn, mincost_fn)
    istr = "[i, j]"
    xstr = "x"
    ystr = "y"
    open(outfile, "w") do f
        write_amplmod_transp_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nDeps, 1..nCust};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nDeps, j in 1..nCust}:")
        println(f, "\t$(nlobj_fn(istr, xstr));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nDeps, j in 1..nCust} $(mincost_fn(istr, ystr));")    
    end
end

function write_ampldat_transp(filename::String, ctype::Symbol, outfile::String)
    pdata = read_transportation(filename)
    
    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_transp(pdata, outfile, ctype, K)
end

function write_ampldat_transp_common(pdata, f)
    ndeps = size(pdata.asscosts, 1)
    ncust = size(pdata.asscosts, 2)
    
    println(f, "param nDeps := $(ndeps);")
    println(f, "param nCust := $(ncust);")
    
    println(f)

    print(f, "param SUPPLY :=")
    for i in 1 : ndeps
        print(f, " $i $(pdata.caps[i])")
    end
    println(f, ";")

    print(f, "param DEMAND :=")
    for j in 1 : ncust
        print(f, " $j $(pdata.demands[j])")
    end
    println(f, ";")
end

function write_ampldat_transp(pdata::TranspProblemData, outfile::String, ctype::Symbol, K::Int64)
    ndeps = size(pdata.asscosts, 1)
    ncust = size(pdata.asscosts, 2)

    cf = zeros(ndeps, ncust, K)
    
    for i in 1 : ndeps, j in 1 : ncust
        maxoffer = 1.0 * min(pdata.demands[j], pdata.caps[i])
        maxcost = 1.0 * pdata.asscosts[i, j]
        cf[i, j, :] = get_nonlinear_coefficients_mod(maxoffer, maxcost, ctype)
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_transp_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k):=")
            for i in 1 : ndeps
                print(f, "\n\t[$i, *]")
                for j in 1 : ncust
                    print(f, "\t$j $(cf[i, j, k])")
                end
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF TRANSPORTATION ------------------------------#

#---------------------------- UFLP NLDEP ---------------------------#

function write_amplmod_uflp_nldep(ctype, outfile)
    nlcond = get_nl_condition("[i]", ctype)
    if ctype == :cubic
        write_amplmod_uflp_nldep(outfile, 4, nlobj_string_cubic, mincost_string_cubic, nlcond)
    elseif ctype == :exp
        write_amplmod_uflp_nldep(outfile, 5, nlobj_string_exp, mincost_string_exp, nlcond)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_uflp_nldep(outfile, 3, nlobj_string_sinabc, mincost_string_other, nlcond)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_uflp_nldep(outfile, 4, nlobj_string_divexp, mincost_string_divexp, nlcond)
    elseif ctype == :sincos
        write_amplmod_uflp_nldep(outfile, 2, nlobj_string_sincos, mincost_string_other, nlcond)
    elseif ctype == :log
        write_amplmod_uflp_nldep(outfile, 3, nlobj_string_log, mincost_string_other, nlcond)
    elseif ctype == :sqrt
        write_amplmod_uflp_nldep(outfile, 2, nlobj_string_sqrt, mincost_string_other, nlcond)
    elseif ctype == :square    
        write_amplmod_uflp_nldep(outfile, 1, nlobj_string_square, mincost_string_other, nlcond)
    end
end

function write_amplmod_uflp_nldep_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nDeps;")
    println(f, "param nCust;")
    println(f)
    println(f, "param DEMAND {1..nCust};")
    println(f, "param ASSIGNMENT_COST {1..nDeps, 1..nCust};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var f {i in 1..nDeps, j in 1..nCust} >= 0, <= 1;")
    println(f, "var x {i in 1..nDeps} >= 0, <= sum {j in 1..nCust} DEMAND[j];")
    println(f, "var y {i in 1..nDeps} binary;")
    println(f, "var nlobj {i in 1..nDeps};")
    println(f)
    println(f, "subject to degree {j in 1..nCust}:")
    println(f, "\tsum {i in 1..nDeps} f[i, j] = 1;")
    println(f)
    println(f, "subject to link {i in 1..nDeps, j in 1..nCust}:")
    println(f, "\tf[i, j] - y[i] <= 0;")
    println(f)
    println(f, "subject to wh {i in 1..nDeps}:")
    println(f, "\tsum {j in 1..nCust} (DEMAND[j] * f[i, j]) - x[i] = 0;")
end

function write_amplmod_uflp_nldep(outfile, ncoefs, nlobj_fn, mincost_fn, nlcond)
    istr = "[i]"
    xstr = "x"
    ystr = "y"
    open(outfile, "w") do f
        write_amplmod_uflp_nldep_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nDeps};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nDeps$(nlcond)}:")
        println(f, "\t$(nlobj_fn(istr, xstr));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nDeps$(nlcond)} ($(mincost_fn(istr, ystr))) + sum {i in 1..nDeps, j in 1..nCust} (ASSIGNMENT_COST[i, j] * f[i, j]);")    
    end
end

function write_ampldat_uflp_nldep(filename::String, ctype::Symbol, outfile::String)
    if occursin("holmberg", filename)
        pdata = read_uflp_holmberg(filename; with_cap_constraint = false)
    else
        pdata = read_uflp_orlib(filename; with_cap_constraint = false)
    end

    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_uflp_nldep(pdata, outfile, ctype, K)
end

function write_ampldat_uflp_nldep_common(pdata, f)
    ndeps = size(pdata.asscosts, 1)
    ncust = size(pdata.asscosts, 2)
    
    println(f, "param nDeps := $(ndeps);")
    println(f, "param nCust := $(ncust);")
    
    println(f)

    print(f, "param DEMAND :=")
    for j in 1 : ncust
        print(f, "\n\t$j $(pdata.demands[j])")
    end
    println(f, " ;")
    print(f, "param ASSIGNMENT_COST :=")
    for i in 1 : ndeps
        print(f, "\n\t[$i, *]")
        for j in 1 : ncust
            print(f, "\t$j $(pdata.asscosts[i, j])")
        end
    end
    println(f, " ;")
end

function write_ampldat_uflp_nldep(pdata::UFLPData, outfile::String, ctype::Symbol, K::Int64)
    ndeps = size(pdata.asscosts, 1)
    
    cf = zeros(ndeps, K)
    
    for i in 1 : ndeps
        maxoffer = min(pdata.caps[i], sum(pdata.demands))
        maxcost = 1.0 * pdata.fixed_costs[i]
        cf[i, :] = get_nonlinear_coefficients_mod(1.0 * maxoffer, maxcost, ctype)
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_uflp_nldep_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k) :=")
            for i in 1 : ndeps
                print(f, "\n\t$i $(cf[i, k])")
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF UFLP NLDEP ---------------------------#

#---------------------------- UFLP NLASS ---------------------------#

function write_amplmod_uflp_nlass(ctype, outfile)
    if ctype == :cubic
        write_amplmod_uflp_nlass(outfile, 4, nlobj_string_cubic, mincost_string_cubic)
    elseif ctype == :exp
        write_amplmod_uflp_nlass(outfile, 5, nlobj_string_exp, mincost_string_exp)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_uflp_nlass(outfile, 3, nlobj_string_sinabc, mincost_string_other)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_uflp_nlass(outfile, 4, nlobj_string_divexp, mincost_string_divexp)
    elseif ctype == :sincos
        write_amplmod_uflp_nlass(outfile, 2, nlobj_string_sincos, mincost_string_other)
    elseif ctype == :log
        write_amplmod_uflp_nlass(outfile, 3, nlobj_string_log, mincost_string_other)
    elseif ctype == :sqrt
        write_amplmod_uflp_nlass(outfile, 2, nlobj_string_sqrt, mincost_string_other)
    elseif ctype == :square    
        write_amplmod_uflp_nlass(outfile, 1, nlobj_string_square, mincost_string_other)
    end
end

function write_amplmod_uflp_nlass_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nDeps;")
    println(f, "param nCust;")
    println(f)
    println(f, "param DEMAND {1..nCust};")
    println(f, "param FIXED_COST {1..nDeps};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var x {i in 1..nDeps, j in 1..nCust} >= 0, <= 1;")
    println(f, "var y {i in 1..nDeps} binary;")
    println(f, "var w {i in 1..nDeps, j in 1..nCust} binary;")
    println(f, "var nlobj {i in 1..nDeps, j in 1..nCust};")
    println(f)
    println(f, "subject to assignment {j in 1..nCust}:")
    println(f, "\tsum {i in 1..nDeps} x[i, j] = 1;")
    println(f)
    println(f, "subject to link {i in 1..nDeps, j in 1..nCust}:")
    println(f, "\tw[i, j] - y[i] <= 0;")
    println(f)
    println(f, "subject to link2 {i in 1..nDeps, j in 1..nCust}:")
    println(f, "\tx[i, j] - w[i, j] <= 0;")
end

function write_amplmod_uflp_nlass(outfile, ncoefs, nlobj_fn, mincost_fn)
    istr = "[i, j]"
    xstr = "x"
    ystr = "w"
    open(outfile, "w") do f
        write_amplmod_uflp_nlass_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nDeps, 1..nCust};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nDeps, j in 1..nCust}:")
        println(f, "\t$(nlobj_fn(istr, xstr));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nDeps, j in 1..nCust} ($(mincost_fn(istr, ystr))) + sum {i in 1..nDeps} (FIXED_COST[i] * y[i]);")    
    end
end

function write_ampldat_uflp_nlass(filename::String, ctype::Symbol, outfile::String)
    if occursin("holmberg", filename)
        pdata = read_uflp_holmberg(filename; with_cap_constraint = false)
    else
        pdata = read_uflp_orlib(filename; with_cap_constraint = false)
    end

    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_uflp_nlass(pdata, outfile, ctype, K)
end

function write_ampldat_uflp_nlass_common(pdata, f)
    ndeps = size(pdata.asscosts, 1)
    ncust = size(pdata.asscosts, 2)
    
    println(f, "param nDeps := $(ndeps);")
    println(f, "param nCust := $(ncust);")
    
    println(f)

    print(f, "param DEMAND :=")
    for j in 1 : ncust
        print(f, "\n\t$j $(pdata.demands[j])")
    end
    println(f, " ;")
    print(f, "param FIXED_COST :=")
    for i in 1 : ndeps
        print(f, "\n\t$i $(pdata.fixed_costs[i])")
    end
    println(f, " ;")
end

function write_ampldat_uflp_nlass(pdata::UFLPData, outfile::String, ctype::Symbol, K::Int64)
    ndeps = size(pdata.asscosts, 1)
    ncust = size(pdata.asscosts, 2)
    
    cf = zeros(ndeps, ncust, K)
    
    for i in 1 : ndeps, j in 1 : ncust
        maxoffer = 1.0
        maxcost = 1.0 * pdata.asscosts[i, j]
        cf[i, j, :] = get_nonlinear_coefficients_mod(maxoffer, maxcost, ctype)
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_uflp_nlass_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k) :=")
            for i in 1 : ndeps
                print(f, "\t[$i, *]")
                for j in 1 : ncust
                    print(f, "\t$j $(cf[i, j, k])")
                end
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF UFLP NLASS ---------------------------#

#---------------------------- MCF ---------------------------#

function write_amplmod_mcf(ctype, outfile)
    if ctype == :cubic
        write_amplmod_mcf(outfile, 4, nlobj_string_cubic, mincost_string_cubic)
    elseif ctype == :exp
        write_amplmod_mcf(outfile, 5, nlobj_string_exp, mincost_string_exp)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_mcf(outfile, 3, nlobj_string_sinabc, mincost_string_other)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_mcf(outfile, 4, nlobj_string_divexp, mincost_string_divexp)
    elseif ctype == :sincos
        write_amplmod_mcf(outfile, 2, nlobj_string_sincos, mincost_string_other)
    elseif ctype == :log
        write_amplmod_mcf(outfile, 3, nlobj_string_log, mincost_string_other)
    elseif ctype == :sqrt
        write_amplmod_mcf(outfile, 2, nlobj_string_sqrt, mincost_string_other)
    elseif ctype == :square    
        write_amplmod_mcf(outfile, 1, nlobj_string_square, mincost_string_other)
    end
end

function write_amplmod_mcf_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nArcs;")
    println(f, "param nComm;")
    println(f, "param nNodes;")
    println(f)
    println(f, "param CAPACITY {1..nArcs};")
    println(f, "param FROM {1..nArcs};")
    println(f, "param TO {1..nArcs};")
    println(f, "param MAXOFFER {1..nNodes};")
    println(f, "param RHS {1..nNodes, 1..nComm};")
    println(f, "param QTY {1..nComm};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var x {a in 1..nArcs, p in 1..nComm} >= 0, <= min(QTY[p], CAPACITY[a]);")
    println(f, "var u {i in 1..nNodes} binary;")
    println(f, "var v {i in 1..nNodes} >= 0, <= MAXOFFER[i];")
    println(f, "var nlobj {i in 1..nNodes};")
    println(f)
    println(f, "subject to flowcons {i in 1..nNodes, p in 1..nComm}:")
    println(f, "\tsum {a in 1..nArcs : FROM[a] = i} x[a, p] - sum {a in 1..nArcs : TO[a] = i} x[a, p] = RHS[i, p];")
    println(f)
    println(f, "subject to capacity {a in 1..nArcs}:")
    println(f, "\tsum {p in 1..nComm} x[a, p] <= CAPACITY[a];")
    println(f)
    println(f, "subject to cutflow {i in 1..nNodes}:")
    println(f, "\tsum {a in 1..nArcs, p in 1..nComm : TO[a] = i} x[a, p] - v[i] = 0;")
    println(f)
    println(f, "subject to nodecap {i in 1..nNodes}:")
    println(f, "\tv[i] - MAXOFFER[i] * u[i] <= 0;")
end

function write_amplmod_mcf(outfile, ncoefs, nlobj_fn, mincost_fn)
    istr = "[i]"
    xstr = "v"
    ystr = "u"
    open(outfile, "w") do f
        write_amplmod_mcf_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nNodes};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nNodes}:")
        println(f, "\t$(nlobj_fn(istr, xstr));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nNodes} $(mincost_fn(istr, ystr));")    
    end
end

function write_ampldat_mcf(pathname::String, filename::String, incr::Int64, cong::Int64, ctype::Symbol, outfile::String)
    gfilename = joinpath(pathname, "$filename.csv")
    cfilename = joinpath(pathname, "$(filename)Commod.csv")
    nfilename = joinpath(pathname, "infoNodes", "nodes_$(filename)_0$(cong).csv")
    pdata = load_mcnddata(gfilename,
                        cfilename,
                        nfilename,
                        incr;
                        congestioncoef = cong)[1]
    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_mcf(pdata, outfile, ctype, K)
end

function write_ampldat_mcf_common(pdata, f)
    narcs = length(pdata.arcs)
    ncomm = length(pdata.commods)
    nnodes = length(pdata.nodes)

    println(f, "param nArcs := $(narcs);")
    println(f, "param nComm := $(ncomm);")
    println(f, "param nNodes := $(nnodes);")
    
    println(f)

    print(f, "param CAPACITY :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:capa])")
    end
    println(f, " ;")
    println(f)
    print(f, "param FROM :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:from])")
    end
    println(f, " ;")
    println(f)
    print(f, "param TO :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:to])")
    end
    println(f, " ;")
    println(f)
    print(f, "param MAXOFFER :=")
    for i in 1 : nnodes
        print(f, "\n\t$i $(pdata.nodes[i][:capainit] + pdata.nodes[i][:capaincrement])")
    end
    println(f, " ;")
    println(f)

    rhs = zeros(nnodes, ncomm)
    for p in 1 : ncomm
        to = pdata.commods[p].to
        from = pdata.commods[p].from
        rhs[from, p] += pdata.commods[p].qty
        rhs[to, p] -= pdata.commods[p].qty
    end

    print(f, "param RHS :=")
    for i in 1 : nnodes
        print(f, "\n\t[$i, *]")
        for p in 1:ncomm
            print(f, " $p $(rhs[i, p])")
        end
    end
    println(f, ";")
    
    println(f)
    print(f, "param QTY :=")
    for p in 1 : ncomm
        print(f, "\n\t$p $(pdata.commods[p][:qty])")
    end
    println(f, " ;")
end

function write_ampldat_mcf(pdata::MCNDData, outfile::String, ctype::Symbol, K::Int64)
    nnodes = size(pdata.nodes, 1)
    
    cf = zeros(nnodes, K)
    
    for i in 1 : nnodes
        maxoffer = pdata.nodes[i][:capainit] + pdata.nodes[i][:capaincrement]
        maxcost = pdata.nodes[i].fixedcost
        cf[i, :] = get_nonlinear_coefficients_mod(1.0 * maxoffer, maxcost, ctype)
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_mcf_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k) :=")
            for i in 1 : nnodes
                print(f, "\n\t$i $(cf[i, k])")
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF MCF ---------------------------#

#---------------------------- MCND ---------------------------#

function write_amplmod_mcnd(ctype, outfile)
    if ctype == :cubic
        write_amplmod_mcnd(outfile, 4, nlobj_string_cubic, mincost_string_cubic)
    elseif ctype == :exp
        write_amplmod_mcnd(outfile, 5, nlobj_string_exp, mincost_string_exp)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_mcnd(outfile, 3, nlobj_string_sinabc, mincost_string_other)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_mcnd(outfile, 4, nlobj_string_divexp, mincost_string_divexp)
    elseif ctype == :sincos
        write_amplmod_mcnd(outfile, 2, nlobj_string_sincos, mincost_string_other)
    elseif ctype == :log
        write_amplmod_mcnd(outfile, 3, nlobj_string_log, mincost_string_other)
    elseif ctype == :sqrt
        write_amplmod_mcnd(outfile, 2, nlobj_string_sqrt, mincost_string_other)
    elseif ctype == :square    
        write_amplmod_mcnd(outfile, 1, nlobj_string_square, mincost_string_other)
    end
end

function write_amplmod_mcnd_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nArcs;")
    println(f, "param nComm;")
    println(f, "param nNodes;")
    println(f)
    println(f, "param CAPACITY {1..nArcs};")
    println(f, "param FROM {1..nArcs};")
    println(f, "param TO {1..nArcs};")
    println(f, "param MAXOFFER {1..nNodes};")
    println(f, "param RHS {1..nNodes, 1..nComm};")
    println(f, "param QTY {1..nComm};")
    println(f, "param FIXED_COST {1..nArcs};")
    println(f, "param VAR_COST {1..nArcs};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var y {a in 1..nArcs} binary;")
    println(f, "var x {a in 1..nArcs, p in 1..nComm} >= 0, <= min(QTY[p], CAPACITY[a]);")
    println(f, "var u {i in 1..nNodes} binary;")
    println(f, "var v {i in 1..nNodes} >= 0, <= MAXOFFER[i];")
    println(f, "var nlobj {i in 1..nNodes};")
    println(f)
    println(f, "subject to flowcons {i in 1..nNodes, p in 1..nComm}:")
    println(f, "\tsum {a in 1..nArcs : FROM[a] = i} x[a, p] - sum {a in 1..nArcs : TO[a] = i} x[a, p] = RHS[i, p];")
    println(f)
    println(f, "subject to xy {a in 1..nArcs, p in 1..nComm}:")
    println(f, "\tx[a, p] - min(QTY[p], CAPACITY[a]) * y[a] <= 0;")
    println(f)
    println(f, "subject to capacity {a in 1..nArcs}:")
    println(f, "\tsum {p in 1..nComm} x[a, p] - CAPACITY[a] * y[a] <= 0;")
    println(f)
    println(f, "subject to cutflow {i in 1..nNodes}:")
    println(f, "\tsum {a in 1..nArcs, p in 1..nComm : TO[a] = i} x[a, p] - v[i] = 0;")
    println(f)
    println(f, "subject to nodecap {i in 1..nNodes}:")
    println(f, "\tv[i] - MAXOFFER[i] * u[i] <= 0;")
end

function write_amplmod_mcnd(outfile, ncoefs, nlobj_fn, mincost_fn)
    istr = "[i]"
    xstr = "v"
    ystr = "u"
    open(outfile, "w") do f
        write_amplmod_mcnd_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nNodes};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nNodes}:")
        println(f, "\t$(nlobj_fn(istr, xstr));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nNodes} ($(mincost_fn(istr, ystr))) + sum {a in 1..nArcs} (FIXED_COST[a] * y[a]) + sum {a in 1..nArcs, p in 1..nComm} (VAR_COST[a] * x[a, p]);")    
    end
end

function write_ampldat_mcnd(pathname::String, filename::String, incr::Int64, cong::Int64, ctype::Symbol, outfile::String)
    gfilename = joinpath(pathname, "$filename.csv")
    cfilename = joinpath(pathname, "$(filename)Commod.csv")
    nfilename = joinpath(pathname, "infoNodes", "nodes_$(filename)_0$(cong).csv")
    pdata = load_mcnddata(gfilename,
                        cfilename,
                        nfilename,
                        incr;
                        congestioncoef = cong)[1]
    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_mcnd(pdata, outfile, ctype, K)
end

function write_ampldat_mcnd_common(pdata, f)
    narcs = length(pdata.arcs)
    ncomm = length(pdata.commods)
    nnodes = length(pdata.nodes)

    println(f, "param nArcs := $(narcs);")
    println(f, "param nComm := $(ncomm);")
    println(f, "param nNodes := $(nnodes);")
    
    println(f)

    print(f, "param CAPACITY :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:capa])")
    end
    println(f, " ;")
    println(f)
    print(f, "param FROM :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:from])")
    end
    println(f, " ;")
    println(f)
    print(f, "param TO :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:to])")
    end
    println(f, " ;")
    println(f)
    print(f, "param MAXOFFER :=")
    for i in 1 : nnodes
        print(f, "\n\t$i $(pdata.nodes[i][:capainit] + pdata.nodes[i][:capaincrement])")
    end
    println(f, " ;")
    println(f)

    rhs = zeros(nnodes, ncomm)
    for p in 1 : ncomm
        to = pdata.commods[p].to
        from = pdata.commods[p].from
        rhs[from, p] += pdata.commods[p].qty
        rhs[to, p] -= pdata.commods[p].qty
    end

    print(f, "param RHS :=")
    for i in 1 : nnodes
        print(f, "\n\t[$i, *]")
        for p in 1:ncomm
            print(f, " $p $(rhs[i, p])")
        end
    end
    println(f, ";")
    
    println(f)
    print(f, "param QTY :=")
    for p in 1 : ncomm
        print(f, "\n\t$p $(pdata.commods[p][:qty])")
    end
    println(f, " ;")
    print(f, "param FIXED_COST :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:fixedcost])")
    end
    println(f, " ;")
    print(f, "param VAR_COST :=")
    for a in 1 : narcs
        print(f, "\n\t$a $(pdata.arcs[a][:varcost])")
    end
    println(f, " ;")
end

function write_ampldat_mcnd(pdata::MCNDData, outfile::String, ctype::Symbol, K::Int64)
    nnodes = size(pdata.nodes, 1)
    
    cf = zeros(nnodes, K)
    
    for i in 1 : nnodes
        maxoffer = pdata.nodes[i][:capainit] + pdata.nodes[i][:capaincrement]
        maxcost = pdata.nodes[i].fixedcost
        cf[i, :] = get_nonlinear_coefficients_mod(1.0 * maxoffer, maxcost, ctype)
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_mcnd_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k) :=")
            for i in 1 : nnodes
                print(f, "\n\t$i $(cf[i, k])")
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF MCND ---------------------------#

#---------------------------- KNAPSACK ---------------------------#

function write_amplmod_knapsack(ctype, outfile)
    if ctype == :cubic
        write_amplmod_knapsack(outfile, 4, nlobj_string_cubic, mincost_string_cubic)
    elseif ctype == :exp
        write_amplmod_knapsack(outfile, 5, nlobj_string_exp, mincost_string_exp)
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3]
        write_amplmod_knapsack(outfile, 3, nlobj_string_sinabc, mincost_string_other)
    elseif ctype in [:pardivexp, :divexp]
        write_amplmod_knapsack(outfile, 4, nlobj_string_divexp, mincost_string_divexp)
    elseif ctype == :sincos
        write_amplmod_knapsack(outfile, 2, nlobj_string_sincos, mincost_string_other)
    elseif ctype == :log
        write_amplmod_knapsack(outfile, 3, nlobj_string_log, mincost_string_other)
    elseif ctype == :sqrt
        write_amplmod_knapsack(outfile, 2, nlobj_string_sqrt, mincost_string_other)
    elseif ctype == :square    
        write_amplmod_knapsack(outfile, 1, nlobj_string_square, mincost_string_other)
    end
end

function write_amplmod_knapsack_common(f, outfile)
    copyright_header(f, outfile)
    println(f)
    println(f, "param nKnap;")
    println(f, "param nObj;")
    println(f)
    println(f, "param OBJECT_U {1..nObj};")
    println(f, "param OBJECT_WEIGHT {1..nObj};")
    println(f, "param KNAPSACK_CAPACITY {1..nKnap};")
    println(f, "param PI := 4 * atan(1);")
    println(f)
    println(f, "var x {i in 1..nKnap, j in 1..nObj} >= 0, <= OBJECT_U[j];")
    println(f, "var onoff {i in 1..nKnap, j in 1..nObj} binary;")
    println(f, "var nlobj {i in 1..nKnap, j in 1..nObj} >= -1E+7, <= 1E+7;")
    println(f)
    println(f, "subject to link {i in 1..nKnap, j in 1..nObj}:")
    println(f, "\tx[i, j] - OBJECT_U[j] * onoff[i, j] <= 0;")
    println(f)
    println(f, "subject to cap {i in 1..nKnap}:")
    println(f, "\tsum {j in 1..nObj} (OBJECT_WEIGHT[j] * x[i, j]) <= KNAPSACK_CAPACITY[i];")
    println(f)
    println(f, "subject to objlen {j in 1..nObj}:")
    println(f, "\tsum {i in 1..nKnap} x[i, j] <= OBJECT_U[j];")
end

function write_amplmod_knapsack(outfile, ncoefs, nlobj_fn, mincost_fn)
    istr = "[i, j]"
    xstr = "x"
    ystr = "onoff"
    sign = "-"
    open(outfile, "w") do f
        write_amplmod_knapsack_common(f, outfile)
        println(f)
        for k in 1:ncoefs
            println(f, "param COEF$k {1..nKnap, 1..nObj};")
        end
        println(f)
        println(f, "subject to nllink {i in 1..nKnap, j in 1..nObj}:")
        println(f, "\t$(nlobj_fn(istr, xstr, sign));")
        println(f)
        println(f, "minimize Cost:")
        println(f, "\tsum {i in 1..nKnap, j in 1..nObj} $(mincost_fn(istr, ystr, sign));")    
    end
end

function write_ampldat_knapsack(filename::String, ctype::Symbol, outfile::String)
    pdata = read_knapsacknlobj(filename, false, ctype)[1]
    if ctype in [:cubic, :divexp, :pardivexp]
        K = 4
    elseif ctype == :exp
        K = 5
    elseif ctype in [:sinabc1, :sinabc2, :sinabc3, :log]
        K = 3
    elseif ctype in [:sincos, :sqrt]
        K = 2
    elseif ctype == :square
        K = 1
    end
    write_ampldat_knapsack(pdata, outfile, ctype, K)
end

function write_ampldat_knapsack_common(pdata, f)
    nobjs = length(pdata.object_U)
    nknap = length(pdata.knapsack_cap)
    
    println(f, "param nObj := $(nobjs);")
    println(f, "param nKnap := $(nknap);")
    
    println(f)

    print(f, "param OBJECT_U :=")
    for j in 1 : nobjs
        print(f, "\n\t$j $(pdata.object_U[j])")
    end
    println(f, " ;")
    println(f)
    print(f, "param OBJECT_WEIGHT :=")
    for j in 1 : nobjs
        print(f, "\n\t$j $(pdata.object_weight[j])")
    end
    println(f, " ;")
    println(f)
    print(f, "param KNAPSACK_CAPACITY :=")
    for i in 1 : nknap
        print(f, "\n\t$i $(pdata.knapsack_cap[i])")
    end
    println(f, " ;")
end

function write_ampldat_knapsack(pdata::KnapData, outfile::String, ctype::Symbol, K::Int64)
    nobjs = size(pdata.object_U, 1)
    nknap = size(pdata.knapsack_cap, 1)
    
    cf = zeros(nknap, nobjs, K)
    
    for i in 1 : nknap, j in 1 : nobjs
        if ctype == :divexp
            cf[i, j, :] = [pdata.minuslistea[j], pdata.object_coefb[j], pdata.object_coefc[j], pdata.object_coefd[j]]
        else
            fkc = x -> ((pdata.object_coefc[j]))/(1+(pdata.object_coefb[j])*(exp((pdata.minuslistea[j])*(x+(pdata.object_coefd[j])))))
            maxoffer = pdata.object_U[j]
            maxcost = fkc(pdata.object_U[j])
            cf[i, j, :] = get_nonlinear_coefficients_mod(1.0 * maxoffer, maxcost, ctype)
        end
    end

    open(outfile, "w") do f
        copyright_header(f, outfile)
        
        println(f)
        
        write_ampldat_knapsack_common(pdata, f)
    
        println(f)

        for k in 1 : K
            print(f, "param COEF$(k) :=")
            for i in 1 : nknap
                print(f, "\t[$i, *]")
                for j in 1 : nobjs
                    print(f, "\t$j $(cf[i, j, k])")
                end
            end
            println(f, " ;")
        end
    end
end

#---------------------------- END OF KNAPSACK ---------------------------#
