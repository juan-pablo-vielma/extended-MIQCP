require("conic.jl")

CplexQCP  = 1
CplexLP   = 2
GurobiQCP = 0
GurobiLP  = 1

function LiftedLpBN(misocp, results)
	solver() = liftedlpsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[implementSOCPStandard],
							CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_PREIND=0),
	 						CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0)
	 						)
	solveInstance(misocp, results, "LiftedLpBN", solver, cplexbasesolver=true)
end

function CplexSepLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
						    [addTower(BNGLIThreeDimCone,0.01), implementSOCPSeparableBase],
						    [implementSOCPSeparable],
						    CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 1
	 						)
	solveInstance(misocp, results, "CplexSepLazyBN", solver, cplexbasesolver=true)
end

function GurobiSepLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01), implementSOCPSeparableBase],
							[implementSOCPSeparable],
							GurobiSolver(TimeLimit=3600, Threads=1),
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0);
	 						separable = 1
	 						)
	solveInstance(misocp, results, "GurobiSepLazyBN", solver, cplexbasesolver=false)
end

function CplexTowerLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "CplexTowerLazyBN", solver, cplexbasesolver=true)
end

function GurobiTowerLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
							[addTower(BNGLIThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							GurobiSolver(TimeLimit=3600, Threads=1),
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "GurobiTowerLazyBN", solver, cplexbasesolver=false)
end

function CplexTowerSepLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
							[addTower(BNGLISepThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1),
	 						CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "CplexTowerSepLazyBN", solver, cplexbasesolver=true)
end

function GurobiTowerSepLazyBN(misocp, results)
	solver() = lazycutsolve(misocp,
							[addTower(BNGLISepThreeDimCone,0.01)],
							[addTower(StandardThreeDimCone)],
							GurobiSolver(TimeLimit=3600, Threads=1),
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0);
	 						separable = 2
	 						)
	solveInstance(misocp, results, "GurobiTowerSepLazyBN", solver, cplexbasesolver=false)
end

function CplexTowerLp(misocp,results)
	solver() = buildModel(misocp,
			  			  [addTower(StandardThreeDimCone)],
			  			  CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
			  			  )
	solveInstance(misocp, results, "CplexTowerLp", solver, cplexbasesolver=true)
end

function GurobiTowerLp(misocp, results)
	solver() = buildModel(misocp,
						  [addTower(StandardThreeDimCone)],
						  GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiTowerLp", solver, cplexbasesolver=false)
end

function CplexTowerSepLp(misocp, results)
	solver() = buildModel(misocp,
						  [addTower(SeparableThreeDimCone)],
						  CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexTowerSepLp", solver, cplexbasesolver=true)
end

function GurobiTowerSepLp(misocp, results)
	solver() = buildModel(misocp,
						  [addTower(SeparableThreeDimCone)],
						  GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiTowerSepLp", solver, cplexbasesolver=false)
end

function CplexSepLp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPSeparable],
						  CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexSepLp", solver, cplexbasesolver=true)
end

function GurobiSepLp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPSeparable],
						  GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiSepLp", solver, cplexbasesolver=false)
end

function CplexLp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexLP)
						  )
	solveInstance(misocp, results, "CplexLp", solver, cplexbasesolver=true)
end

function GurobiLp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP)
						  )
	solveInstance(misocp, results, "GurobiLp", solver, cplexbasesolver=false)
end

function CplexQcp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2, CPX_PARAM_MIQCPSTRAT=CplexQCP)
						  )
	solveInstance(misocp, results, "CplexQCP", solver, cplexbasesolver=true)
end

function GurobiQcp(misocp, results)
	solver() = buildModel(misocp,
						  [implementSOCPStandard],
						  GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiQCP)
						  )
	solveInstance(misocp, results, "GurobiQCP", solver, cplexbasesolver=false)
end

function test(instancesize, resultfilename, solvers, mark=true, short=true, robust=true)

	resultfile = open(resultfilename, "a")
	resultsMark = {}
	resultsShort = {}
	resultsRobust = {}

	files = filter!(Regex(string("portfolio_",instancesize,"_.*\.por")), readdir("portfolios"))

	if mark
		for file in files
			misocp = buildMarkowitzPor(joinpath("portfolios", file))
			for solver in solvers
				solver(misocp, resultsMark)
				println(resultfile, "Mark,", instancesize, ",", file, ",",
						join(resultsMark[length(resultsMark)], ","))
				flush(resultfile)
			end
		end
	end
	if short
		for file in files
			misocp = buildShortfallPor(joinpath("portfolios", file))
			for solver in solvers
				solver(misocp, resultsShort)
				println(resultfile, "Short,", instancesize, ",", file, ",",
						join(resultsShort[length(resultsShort)], ","))
				flush(resultfile)
			end
		end
	end

	if robust
		files = filter!(Regex(string("robustportfolio_", instancesize, "_.*\.por")), readdir("robustportfolios"))
		for file in files
			misocp = buildRobustMarkowitzPor(joinpath("robustportfolios", file))
			for solver in solvers
				solver(misocp, resultsRobust)
				println(resultfile, "Robust,", instancesize, ",", file, ",",
						join(resultsRobust[length(resultsRobust)],","))
				flush(resultfile)
			end
		end
	end
end


test(20, ARGS[1], [CplexSepLazyBN])

# test(20,"results.csv","stats.csv",[ LiftedLpBN,
# 									CplexSepLazyBN, GurobiSepLazyBN, CplexTowerLazyBN, GurobiTowerLazyBN, CplexTowerSepLazyBN, GurobiTowerSepLazyBN,
# 									CplexSepLp,	    GurobiSepLp,     CplexTowerLp,     GurobiTowerLp,     CplexTowerSepLp,     GurobiTowerSepLp,
# 									CplexLp,     GurobiLp,
# 									CplexQcp,       GurobiQcp    ])
# test(30,"results.csv","stats.csv",[ LiftedLpBN,
# 									CplexSepLazyBN, GurobiSepLazyBN, CplexTowerLazyBN, GurobiTowerLazyBN, CplexTowerSepLazyBN, GurobiTowerSepLazyBN,
# 									CplexSepLp,	    GurobiSepLp,     CplexTowerLp,     GurobiTowerLp,     CplexTowerSepLp,     GurobiTowerSepLp,
# 									CplexLp,     GurobiLp,
# 									CplexQcp,       GurobiQcp    ])

