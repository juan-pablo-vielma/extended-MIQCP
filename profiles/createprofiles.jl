function createFileForProfile(data,filename,method,name,instances,datacolumn)
	stream = open(filename,"w")

	println(stream,"---")
	println(stream,"algname: ",name)
	println(stream,"success: solved")
	println(stream,"free_format: True")
	println(stream,"---")
	for i in 1:size(data,1)
		if data[i,4] == method
			instance = join(data[i,1:2],"-")
			if in(instance,instances)
				println(stream,join(data[i,1:3],"-")," solved ",data[i,datacolumn])
			end
		end
	end
end



data = readdlm("../results.csv",',');


solvers = ["LiftedLpBN","CplexSepLp","GurobiSepLp","CplexSepLazyBN","GurobiSepLazyBN","CplexQcp", "GurobiQcp","CplexLp", "GurobiLp","CplexTowerLp", "CplexTowerSepLp","GurobiTowerLp", "GurobiTowerSepLp"]
names = ["LiftedLP","CPLEXSepLP","GurobiSepLP","CPLEXSepLazy","GurobiSepLazy","CPLEXCP", "GurobiCP","CPLEXLP", "GurobiLP","CPLEXTowerLP", "CPLEXTowerSepLP", "GurobiTowerLP", "GurobiTowerSepLP"]


for i in 1:length(solvers)
	createFileForProfile(data,string(names[i],"-prof"),solvers[i],names[i], 
		[[string("Mark-",i) for i in [10.0:10:60]];[string("Short-",i) for i in [10.0:10:60]];[string("Robust-",i) for i in [[10.0:10:60];100.0;200.0]]] 
		, 6)
end





