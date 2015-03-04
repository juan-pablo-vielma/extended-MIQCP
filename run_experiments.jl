require("conic.jl")

CplexQCP  = 1
CplexLP   = 2
GurobiQCP = 0
GurobiLP  = 1

function LiftedLpBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = liftedlpsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[implementSOCPStandard],
							CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_PREIND=0),
	 						CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0)
	 						)
	solveInstance(misocp, results, "LiftedLpBN", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function CplexSepLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end

	solver() = lazycutsolve(misocp,
						    [addTower(BNGLIThreeDimCone,0.01), implementSOCPSeparableBase],
						    [implementSOCPSeparable],
						    CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 1
	 						)
	solveInstance(misocp, results, "CplexSepLazyBN", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiSepLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01), implementSOCPSeparableBase],
							[implementSOCPSeparable],
							GurobiSolver(TimeLimit=timelimit, Threads=1),
	 						GurobiSolver(TimeLimit=timelimit, Threads=1, OutputFlag=0);
	 						separable = 1
	 						)
	solveInstance(misocp, results, "GurobiSepLazyBN", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexTowerLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "CplexTowerLazyBN", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiTowerLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							GurobiSolver(TimeLimit=timelimit, Threads=1),
	 						GurobiSolver(TimeLimit=timelimit, Threads=1, OutputFlag=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "GurobiTowerLazyBN", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexTowerSepLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = lazycutsolve(misocp,
							[addTower(BNGLISepThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "CplexTowerSepLazyBN", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiTowerSepLazyBN(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = lazycutsolve(misocp,
							[addTower(BNGLISepThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							GurobiSolver(TimeLimit=timelimit, Threads=1),
	 						GurobiSolver(TimeLimit=timelimit, Threads=1, OutputFlag=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "GurobiTowerSepLazyBN", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexTowerLp(misocp,results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
			  			  [addTower(StandardThreeDimCone)],
			  			  CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
			  			  )
	solveInstance(misocp, results, "CplexTowerLp", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiTowerLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [addTower(StandardThreeDimCone)],
						  GurobiSolver(TimeLimit=timelimit, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiTowerLp", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexTowerSepLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [addTower(SeparableThreeDimCone)],
						  CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexTowerSepLp", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiTowerSepLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [addTower(SeparableThreeDimCone)],
						  GurobiSolver(TimeLimit=timelimit, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiTowerSepLp", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexSepLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPSeparable],
						  CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexSepLp", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiSepLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPSeparable],
						  GurobiSolver(TimeLimit=timelimit, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiSepLp", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexLp", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiLp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  GurobiSolver(TimeLimit=timelimit, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiLp", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function CplexQcp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  CplexSolver(CPX_PARAM_TILIM=timelimit, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexQCP)
						  )
	solveInstance(misocp, results, "CplexQCP", solver, cplexbasesolver=true;writefile=writefile,filename=filename)
end

function GurobiQcp(misocp, results; writefile=false,filename="")
	if writefile
		timelimit=0
	else
		timelimit=3600
	end
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  GurobiSolver(TimeLimit=timelimit, Threads=1, MIQCPMethod=GurobiQCP)
						  )
	solveInstance(misocp, results, "GurobiQCP", solver, cplexbasesolver=false;writefile=writefile,filename=filename)
end

function test(instancesize, resultfilename, solvers, mark=true, short=true, robust=true; writefile=false)
	if !writefile
		resultfile = open(resultfilename, "a")
	end
	resultsMark = {}
	resultsShort = {}
	resultsRobust = {}

	files = filter!(Regex(string("portfolio_",instancesize,"_.*\.por")), readdir("portfolios"))
	if mark
		for file in files
			misocp = buildMarkowitzPor(joinpath("portfolios", file))
			for solver in solvers
				println(solver)
				solver(misocp, resultsMark; writefile=writefile, filename=string("mpsfiles/Mark_",file,"_",solver,".mps"))
				if !writefile
					println(resultfile, "Mark,", instancesize, ",", file, ",",
							join(resultsMark[length(resultsMark)], ","))
					flush(resultfile)
				end
			end
		end
	end
	if short
		for file in files
			misocp = buildShortfallPor(joinpath("portfolios", file))
			for solver in solvers
				solver(misocp, resultsShort; writefile=writefile, filename=string("mpsfiles/Short_",file,"_",solver,".mps"))
				if !writefile
					println(resultfile, "Short,", instancesize, ",", file, ",",
							join(resultsShort[length(resultsShort)], ","))
					flush(resultfile)
				end
			end
		end
	end

	if robust
		files = filter!(Regex(string("robustportfolio_", instancesize, "_.*\.por")), readdir("robustportfolios"))
		for file in files
			misocp = buildRobustMarkowitzPor(joinpath("robustportfolios", file))
			for solver in solvers
				solver(misocp, resultsRobust; writefile=writefile, filename=string("mpsfiles/Robust_",file,"_",solver,".mps"))
				if !writefile
					println(resultfile, "Robust,", instancesize, ",", file, ",",
							join(resultsRobust[length(resultsRobust)],","))
					flush(resultfile)
				end
			end
		end
	end
end

if length(ARGS) > 1
	test(int(ARGS[2]),ARGS[1],[eval(parse(ARGS[i])) for i in 3:length(ARGS)])
elseif length(ARGS) > 0

	for n in 20:10:60
		test(n, ARGS[1], [CplexSepLp, GurobiSepLp, CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp, CplexQcp, GurobiQcp])
	end

	for n in 20:10:60
		test(n, ARGS[1], [CplexLp, GurobiLp])
	end

	test(100, ARGS[1], [CplexSepLp, GurobiSepLp])
	test(200, ARGS[1], [CplexSepLp, GurobiSepLp],false,false)
	test(300, ARGS[1], [CplexSepLp, GurobiSepLp],false,false)

	test(100, ARGS[1], [CplexLp, GurobiLp, CplexQcp, GurobiQcp],false,false)
	test(200, ARGS[1], [CplexLp, GurobiLp, CplexQcp, GurobiQcp],false,false)

	test(100, ARGS[1], [CplexLp, GurobiLp],false,false)
	test(200, ARGS[1], [CplexQcp, GurobiQcp],false,false)

	test(100, ARGS[1], [CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp],false,false)
	test(200, ARGS[1], [CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp],false,false)

	for n in 20:10:60
		test(n, ARGS[1], [LiftedLpBN, CplexSepLazyBN, GurobiSepLazyBN])
	end
	test(100, ARGS[1], [LiftedLpBN, CplexSepLazyBN, GurobiSepLazyBN])
	test(200, ARGS[1], [LiftedLpBN, CplexSepLazyBN, GurobiSepLazyBN],false,false)
	test(300, ARGS[1], [LiftedLpBN, CplexSepLazyBN, GurobiSepLazyBN],false,false)
else
	for n in 20:10:60
		test(n, "", [CplexSepLp, GurobiSepLp, CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp, CplexQcp, GurobiQcp];writefile=true)
	end
	test(100, "", [CplexSepLp, GurobiSepLp, CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp, CplexQcp, GurobiQcp];writefile=true)
	test(200, "", [CplexSepLp, GurobiSepLp, CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp, CplexQcp, GurobiQcp],false,false;writefile=true)
	test(300, "", [CplexSepLp, GurobiSepLp, CplexTowerLp, GurobiTowerLp, CplexTowerSepLp, GurobiTowerSepLp, CplexQcp, GurobiQcp],false,false;writefile=true)
end

