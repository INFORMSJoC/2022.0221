machine=gethostname()

if (machine != "bodo.laas.fr")
        println("setenv GRB_LICENSE_FILE /home/sungueve/Licences/Gurobi/gurobi-",machine,".lic")
        #run(`setenv GRB_LICENSE_FILE /home/sungueve/Licences/Gurobi/gurobi-$machine.lic`)
        #run(`GRB_LICENSE_FILE=/home/sungueve/Licences/Gurobi/gurobi-$machine.lic`)
        #run(`export GRB_LICENCE_FILE`)
        ENV["GRB_LICENSE_FILE"] = "/home/sungueve/Licences/Gurobi/gurobi-$machine.lic"
else
        println("machine %s: the license file should be at its place", gethostname())
end

if (ARGS[1] == "hev")
        target = parse(Float64,ARGS[4])
        init=parse(Float64,ARGS[5])
        stepdiv=parse(Float64,ARGS[6])
        dynstep=!iszero(parse(Int,ARGS[7]))
        report = solve_hev("mini", ARGS[2], target; init = init, printon=false, stepdiv=stepdiv, dynstep=dynstep)
        #first execution to precompile all functions and second to solve ?
        report = solve_hev(ARGS[3], ARGS[2], target; init = init, printon=false, stepdiv=stepdiv, dynstep=dynstep)
        t=now()
        gaptotal=report[end][6]
        finished = false
        if (gaptotal <= target)
                finished=true
        end
        println(report)
        println("pb=", ARGS[1], "; FCname=", ARGS[2],"; Pdemprofilename=", ARGS[3],
                "; init=", init,"; target=", target,
                "; stepdiv=", stepdiv, "; dynstep=", dynstep,
                "; lb=", report[end][1],
                "; ub=", report[end][2], "; obj0=", report[end][3], "; obj0ub=", report[end][4],
                "; gap=", report[end][5], "; gaptotal=", report[end][6], "; num_pieces=", report[end][7],
                "; nbattempts=", report[end][8], "; mipcpu=", report[end][9],
                "; cpu=", report[end][10], "; finished=", finished, "; machine=", gethostname(),
                "; datetime=", Dates.format(t, "yyyy-mm-dd_HH:MM:SS"))

elseif (ARGS[1] == "mcf")
        println(ARGS)
        inst=ARGS[2]
        version=ARGS[3]
        #3 ou 8
        incremcost=parse(Int,ARGS[4])
        congestioncoef=parse(Float64,ARGS[5])
        #6,7 ou 8
        filepathname=ARGS[6]

        graphfilepathname = string(filepathname,inst,".csv")
        commodityfilepathname = string(filepathname,inst,"Commod.csv")
        nodefilepathname = string(filepathname,"infoNodes/nodes_",inst,"_0",version,".csv")

        testgraphfilepathname = string(filepathname,"test.csv")
        testcommodityfilepathname = string(filepathname,"testCommod.csv")
        testnodefilepathname = string(filepathname,"infoNodes/nodes_test_03.csv")
        target = parse(Float64,ARGS[7])
        init=parse(Float64,ARGS[8])
        stepdiv=parse(Float64,ARGS[9])
        dynstep=!iszero(parse(Int,ARGS[10]))
        report = solve_mcnd_or_mcf_from_file("mcf", testgraphfilepathname, testcommodityfilepathname, testnodefilepathname, incremcost, target; init=init, printon=false, congestioncoef=congestioncoef, stepdiv=stepdiv)
        #first execution to precompile all functions and second to solve ?
        report = solve_mcnd_or_mcf_from_file("mcf", graphfilepathname, commodityfilepathname, nodefilepathname, incremcost, target; init=init, printon=false, congestioncoef=congestioncoef, stepdiv=stepdiv)
        t=now()
        gaptotal=report[end][6]
        finished = false
        if (gaptotal <= target)
                finished=true
        end
        println(report)
        println("pb=mcf", "; inst=",inst,"; version=",version, "; incremcost=",incremcost,
                "; congestioncoef=",congestioncoef, "; init=",init,"; target=",target,
                "; stepdiv=", stepdiv, "; dynstep=", dynstep,
                "; lb=", report[end][1],
                "; ub=", report[end][2], "; obj0=", report[end][3], "; obj0ub=", report[end][4],
                "; gap=", report[end][5], "; gaptotal=", report[end][6], "; num_pieces=", report[end][7],
                "; nbattempts=", report[end][8], "; mipcpu=", report[end][9],
                "; cpu=", report[end][10], "; finished=", finished, "; machine=", gethostname(),
                "; datetime=", Dates.format(t, "yyyy-mm-dd_HH:MM:SS"))

elseif (ARGS[1] == "mcnd")
        println(ARGS)
        inst=ARGS[2]
        version=ARGS[3]
        #3 ou 8
        incremcost=parse(Int,ARGS[4])
        congestioncoef=parse(Float64,ARGS[5])
        #6,7 ou 8
        filepathname=ARGS[6]

        graphfilepathname = string(filepathname,inst,".csv")
        commodityfilepathname = string(filepathname,inst,"Commod.csv")
        nodefilepathname = string(filepathname,"infoNodes/nodes_",inst,"_0",version,".csv")

        testgraphfilepathname = string(filepathname,"test.csv")
        testcommodityfilepathname = string(filepathname,"testCommod.csv")
        testnodefilepathname = string(filepathname,"infoNodes/nodes_test_03.csv")

        target = parse(Float64,ARGS[7])
        init=parse(Float64,ARGS[8])
        stepdiv=parse(Float64,ARGS[9])
        dynstep=!iszero(parse(Int,ARGS[10]))
        report = solve_mcnd_or_mcf_from_file("mcnd", testgraphfilepathname, testcommodityfilepathname, testnodefilepathname, incremcost, target; init=init, printon=false, congestioncoef=congestioncoef, stepdiv=stepdiv, dynstep=dynstep)
        #first execution to precompile all functions and second to solve ?
        report = solve_mcnd_or_mcf_from_file("mcnd", graphfilepathname, commodityfilepathname, nodefilepathname, incremcost, target; init=init, printon=false, congestioncoef=congestioncoef, stepdiv=stepdiv, dynstep=dynstep)
        t=now()
        gaptotal=report[end][6]
        finished = false
        if (gaptotal <= target)
                finished=true
        end
        println(report)
        println("pbtype=mcnd ", "; inst_version_incremcost_congestioncoef=",inst,"+",version, "+",incremcost,
                "+",congestioncoef, "; init=", init, "; target=", target,
                "; stepdiv=", stepdiv, "; dynstep=", dynstep,
                "; lb=", report[end][1],
                "; ub=", report[end][2], "; obj0=", report[end][3], "; obj0ub=", report[end][4],
                "; gap=", report[end][5], "; gaptotal=", report[end][6], "; num_pieces=", report[end][7],
                "; nbattempts=", report[end][8], "; mipcpu=", report[end][9],
                "; cpu=", report[end][10], "; finished=", finished, "; machine=", gethostname(),
                "; datetime=", Dates.format(t, "yyyy-mm-dd_HH:MM:SS"))
else
        println("incorrect problem type")
end
