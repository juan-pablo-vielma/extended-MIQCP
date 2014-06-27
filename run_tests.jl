require("conic.jl")

CplexQCP=1
CplexLP=2
GurobiQCP=0
GurobiLP=1

function LiftedLpBN(misocp,results)
	solver() = liftedlpsolve(misocp, [addTower(BNGLIThreeDimCone,0.01)], [implementSOCPStandard] , CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_PREIND=0), 
	 						 CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0) )
	solveInstance(misocp,results,"LiftedLpBN",solver,cplexbasesolver=true)
end
function CplexSepLazyBN(misocp,results)
	solver() = lazycutsolve(misocp, [addTower(BNGLIThreeDimCone,0.01),implementSOCPSeparableBase], [implementSOCPSeparable] , CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1), 
	 						 CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0) ;  separable = 1)
	solveInstance(misocp,results,"CplexSepLazyBN",solver,cplexbasesolver=true)
end
function GurobiSepLazyBN(misocp,results)
	solver() =  lazycutsolve(misocp, [addTower(BNGLIThreeDimCone,0.01),implementSOCPSeparableBase], [implementSOCPSeparable] , GurobiSolver(TimeLimit=3600, Threads=1), 
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0) ;  separable = 1)
	solveInstance(misocp,results,"GurobiSepLazyBN",solver,cplexbasesolver=false)
end
function CplexTowerLazyBN(misocp,results)
	solver() = lazycutsolve(misocp, [addTower(BNGLIThreeDimCone,0.01)], [addTower(StandardThreeDimCone)] , CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1), 
	 						 CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0) ;  separable = 2)
	solveInstance(misocp,results,"CplexTowerLazyBN",solver,cplexbasesolver=true)
end
function GurobiTowerLazyBN(misocp,results)
	solver() =  lazycutsolve(misocp, [addTower(BNGLIThreeDimCone,0.01)], [addTower(StandardThreeDimCone)] , GurobiSolver(TimeLimit=3600, Threads=1), 
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0) ;  separable = 2)
	solveInstance(misocp,results,"GurobiTowerLazyBN",solver,cplexbasesolver=false)
end
function CplexTowerSepLazyBN(misocp,results)
	solver() = lazycutsolve(misocp, [addTower(BNGLISepThreeDimCone,0.01)], [addTower(StandardThreeDimCone)] , CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1), 
	 						 CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0) ;  separable = 2)
	solveInstance(misocp,results,"CplexTowerSepLazyBN",solver,cplexbasesolver=true)
end
function GurobiTowerSepLazyBN(misocp,results)
	solver() =  lazycutsolve(misocp, [addTower(BNGLISepThreeDimCone,0.01)], [addTower(StandardThreeDimCone)] , GurobiSolver(TimeLimit=3600, Threads=1), 
	 						GurobiSolver(TimeLimit=3600, Threads=1, OutputFlag=0) ;  separable = 2)
	solveInstance(misocp,results,"GurobiTowerSepLazyBN",solver,cplexbasesolver=false)
end
function CplexTowerLp(misocp,results)
	solver()=buildModel(misocp,[addTower(StandardThreeDimCone)],CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_MIQCPSTRAT=CplexLP))
	solveInstance(misocp,results,"CplexTowerLp",solver,cplexbasesolver=true)
end
function GurobiTowerLp(misocp,results)
	solver()=buildModel(misocp,[addTower(StandardThreeDimCone)],GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP))
	solveInstance(misocp,results,"GurobiTowerLp",solver,cplexbasesolver=false)
end
function CplexTowerSepLp(misocp,results)
	solver()=buildModel(misocp,[addTower(SeparableThreeDimCone)],CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_MIQCPSTRAT=CplexLP))
	solveInstance(misocp,results,"CplexTowerSepLp",solver,cplexbasesolver=true)
end
function GurobiTowerSepLp(misocp,results)
	solver()=buildModel(misocp,[addTower(SeparableThreeDimCone)],GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP))
	solveInstance(misocp,results,"GurobiTowerSepLp",solver,cplexbasesolver=false)
end
function CplexSepLp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPSeparable],CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_MIQCPSTRAT=CplexLP))
	solveInstance(misocp,results,"CplexSepLp",solver,cplexbasesolver=true)
end
function GurobiSepLp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPSeparable],GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP))
	solveInstance(misocp,results,"GurobiSepLp",solver,cplexbasesolver=false)
end
function CplexLp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPStandard],CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_MIQCPSTRAT=CplexLP))
	solveInstance(misocp,results,"CplexLp",solver,cplexbasesolver=true)
end
function GurobiLp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPStandard],GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiLP))
	solveInstance(misocp,results,"GurobiLp",solver,cplexbasesolver=false)
end
function CplexQcp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPStandard],CplexSolver(CPX_PARAM_TILIM=3600, CPX_PARAM_THREADS=1, CPX_PARAM_MIPDISPLAY=2,CPX_PARAM_MIQCPSTRAT=CplexQCP))
	solveInstance(misocp,results,"CplexQCP",solver,cplexbasesolver=true)
end
function GurobiQcp(misocp,results)
	solver()=buildModel(misocp,[implementSOCPStandard],GurobiSolver(TimeLimit=3600, Threads=1, MIQCPMethod=GurobiQCP))
	solveInstance(misocp,results,"GurobiQCP",solver,cplexbasesolver=false)
end

function printstats(results,statsfile)
	resultData = Dict()
	for result in results
		if !haskey(resultData,result[1])
				resultData[result[1]]=Array(Array{Float64},3)
			fill!(resultData[result[1]],[])
		end
		push!(resultData[result[1]][1],result[4])
		push!(resultData[result[1]][2],result[6])
		push!(resultData[result[1]][3],maximum([result[7],result[8],result[9]]))
	end
	println(statsfile,"Time")
	for solver in keys(resultData)
		methodData = resultData[solver][1]
		for i in [solver,minimum(methodData),  mean(methodData), maximum(methodData),std(methodData)]
			print(statsfile,i,repeat(" ",25-length(string(i))))
		end
		println(statsfile,"")
	end
	println(statsfile,"Nodes")
	for solver in keys(resultData)
		methodData = resultData[solver][2]
		for i in [solver,minimum(methodData),  mean(methodData), maximum(methodData),std(methodData)]
			print(statsfile,i,repeat(" ",25-length(string(i))))
		end
		println(statsfile,"")
	end
	println(statsfile,"Quality")
	for solver in keys(resultData)
		methodData = resultData[solver][3]
		for i in [solver,minimum(methodData),  mean(methodData), maximum(methodData),std(methodData)]
			print(statsfile,i,repeat(" ",25-length(string(i))))
		end
		println(statsfile,"")
	end
	flush(statsfile)
end

function test(instancesize,resultfilename,statsfilename,solvers)

	resultfile = open(resultfilename,"a")
	statsfile = open(statsfilename,"a")
	resultsMark = {}
	resultsShort = {}
	resultsRobust = {}
	
	files = filter!(Regex(string("portfolio_",instancesize,"_.*\.por")), readdir("portfolios"))
	
	for file in files
		misocp = buildMarkowitzPor(joinpath("portfolios", file))
		for solver in solvers
			solver(misocp,resultsMark)
			println(resultfile, "Mark,",instancesize,",",file,",",join(resultsMark[length(resultsMark)],",")); flush(resultfile)
		end
	end
	println(statsfile,"Mark ",instancesize)
	printstats(resultsMark,statsfile)
	for file in files
		misocp = buildShortfallPor(joinpath("portfolios", file))
		for solver in solvers
			solver(misocp,resultsShort)
			println(resultfile, "Short,",instancesize,",",file,",",join(resultsShort[length(resultsShort)],",")); flush(resultfile)
		end
	end
	println(statsfile,"Short ",instancesize)
	printstats(resultsShort,statsfile)
	files = filter!(Regex(string("robust_portfolio_",instancesize,"_.*\.por")), readdir("robust_portfolios"))
	for file in files
		misocp = buildRobustMarkowitzPor(joinpath("robust_portfolios", file))
		for solver in solvers
			solver(misocp,resultsRobust)
			println(resultfile, "Robust,",instancesize,",",file,",",join(resultsRobust[length(resultsRobust)],",")); flush(resultfile)
		end
	end
	println(statsfile,"Robust ",instancesize)
	printstats(resultsRobust,statsfile)
end


test(20,"results.csv","stats.csv",[ LiftedLpBN, 
									CplexSepLazyBN, GurobiSepLazyBN, CplexTowerLazyBN, GurobiTowerLazyBN, CplexTowerSepLazyBN, GurobiTowerSepLazyBN,
									CplexSepLp,	    GurobiSepLp,     CplexTowerLp,     GurobiTowerLp,     CplexTowerSepLp,     GurobiTowerSepLp,
									CplexLp,     GurobiLp,
									CplexQcp,       GurobiQcp    ])
test(30,"results.csv","stats.csv",[ LiftedLpBN, 
									CplexSepLazyBN, GurobiSepLazyBN, CplexTowerLazyBN, GurobiTowerLazyBN, CplexTowerSepLazyBN, GurobiTowerSepLazyBN,
									CplexSepLp,	    GurobiSepLp,     CplexTowerLp,     GurobiTowerLp,     CplexTowerSepLp,     GurobiTowerSepLp,
									CplexLp,     GurobiLp,
									CplexQcp,       GurobiQcp    ])

