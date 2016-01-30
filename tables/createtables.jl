


function printStats(data, filename, problemtype, methods, size, datacolumn, methodcolumn, instancecolumn, labels,formatstring, wins = true)
	tempdata = data[(data[:,1] .== problemtype) & (data[:,2] .== size) & (data[:,4] .== methods[1]),[methodcolumn,instancecolumn,datacolumn ]]
	for i  in 2:length(methods)
		tempdata = [tempdata ; data[(data[:,1] .== problemtype) & (data[:,2] .== size) & (data[:,4] .== methods[i]),[methodcolumn,instancecolumn,datacolumn ]]]
	end
	instances = Set(tempdata[:,2]);
	stream = open(filename,"w")
	println(stream,"\\begin{tabular}{lrrrrrrr}")
	println(stream,"Algorithm & min & avg & max & std")
	if wins	
		println(stream,"& wins & 1\\% win & 10\\% win ")
	end
	println(stream,"\\\\")
	println(stream,"\\hline")
	#println(problemtype)
	for i  in 1:length(methods)
		numwins = 0
		numwins10 = 0
		numwins1 = 0
	#	println(methods[i])
		if wins 
			for instance in instances
	#			println(instance)
				if minimum(tempdata[tempdata[:,2] .== instance,3]) >= tempdata[(tempdata[:,2] .== instance) & (tempdata[:,1] .== methods[i]),3][1]
					numwins += 1
				end
				if 1.01*minimum(tempdata[tempdata[:,2] .== instance,3]) >= tempdata[(tempdata[:,2] .== instance) & (tempdata[:,1] .== methods[i]),3][1]
					numwins1 += 1
				end
				if 1.10*minimum(tempdata[tempdata[:,2] .== instance,3]) >= tempdata[(tempdata[:,2] .== instance) & (tempdata[:,1] .== methods[i]),3][1]
					numwins10 += 1
				end
			end
		end
		temp = convert(Array{Float64,1},tempdata[tempdata[:,1] .== methods[i],3]);
		print(stream,labels[i],"&")
		@eval @printf($stream,$formatstring,minimum($temp))
		print(stream,"& \\bf")
		@eval @printf($stream,$formatstring,mean($temp))
		print(stream,"& ")
		@eval @printf($stream,$formatstring,maximum($temp))
		print(stream,"& ")
		@eval @printf($stream,$formatstring,std($temp))
		if wins 
			print(stream,"&",numwins,"&",numwins1,"&",numwins10)
		end
		if i == length(methods)
			println(stream,"")
		else
			println(stream,"\\\\")
		end
	end
	println(stream,"\\end{tabular}")
	close(stream)
end

data = readdlm(ARGS[1],',');

const floatformat = "%.2f"
const scientificformat = "%.2e"

if length(ARGS) > 1
	printStats(data,string("test_time.tex"),ARGS[3],ARGS[4:end],int(ARGS[2]),6,4,3,ARGS[4:end],floatformat);
	printStats(data,string("test_quality.tex"),ARGS[3],ARGS[4:end],int(ARGS[2]),8,4,3,ARGS[4:end],floatformat);
else

	for k in [20,30,40,50,60]
		println(k)
		printStats(data,string("mark",k,".tex"),"Mark",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,6,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],floatformat);
		println(k)
		printStats(data,string("short",k,".tex"),"Short",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,6,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],floatformat);
		println(k)
		printStats(data,string("robust",k,".tex"),"Robust",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,6,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],floatformat);
	end

	for k in [100]
		println(k)
		printStats(data,string("mark",k,".tex"),"Mark",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,6,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],floatformat);
		println(k)
		printStats(data,string("short",k,".tex"),"Short",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,6,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],floatformat);
		println(k)
		printStats(data,string("robust",k,".tex"),"Robust",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,6,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],floatformat);
	end

	for k in [200]

		println(k)
		printStats(data,string("robust",k,".tex"),"Robust",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,6,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],floatformat);
	end

	for k in [300]
		println(k)
		printStats(data,string("robust",k,".tex"),"Robust",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,6,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],floatformat);
	end

	for k in [20,30,40,50,60]
		println(k)
		printStats(data,string("qualitymark",k,".tex"),"Mark",[ "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
		println(k)
		printStats(data,string("qualityshort",k,".tex"),"Short",[ "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
		println(k)
		printStats(data,string("qualityrobust",k,".tex"),"Robust",[ "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
	end

	for k in [100]
		println(k)
		printStats(data,string("qualitymark",k,".tex"),"Mark",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
		println(k)
		printStats(data,string("qualityshort",k,".tex"),"Short",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
		println(k)
		printStats(data,string("qualityrobust",k,".tex"),"Robust",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,8,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
	end

	for k in [200]

		println(k)
		printStats(data,string("qualityrobust",k,".tex"),"Robust",["CplexQCP", "GurobiQCP", "CplexLp", "GurobiLp", "LiftedLpBN","CplexSepLp", "CplexTowerLp", "CplexTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp","CplexSepLazyBN","GurobiSepLazyBN" ],k,8,4,3,
								      ["CPLEXCP", "GurobiCP", "CPLEXLP", "GurobiLP","LiftedLP","CPLEXSepLp", "CPLEXTowerLp", "CPLEXTowerSepLp","GurobiSepLp", "GurobiTowerLp", "GurobiTowerSepLp", "CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
	end

	for k in [300]
		println(k)
		printStats(data,string("qualityrobust",k,".tex"),"Robust",[ "LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN"],k,8,4,3,
								      ["LiftedLP","CPLEXSepLp","GurobiSepLp","CPLEXSepLazy","GurobiSepLazy"],scientificformat,false);
	end

end
