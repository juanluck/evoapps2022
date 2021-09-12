//currentDate = getdate()
//rand('seed',getdate('s')+currentDate(10))

args=sciargs();

exec("./scilab-scripts/benchmarksRIM.sce");
exec("./scilab-scripts/modelsRIM.sce");
exec("./scilab-scripts/optimBenchmarkRIM.sce");

iter=args(6)
nmodel=args(7)
nbench=args(8)

//disp(iter)
//disp(nmodel)
//disp(nbench)

// Defining the model to test
model=modelRIM_K_L;
constraints=constRIM_K_L;
if strcmp(nmodel,"RIM_K_L")==0 then 
	model=modelRIM_K_L;
	constraints=constRIM_K_L;
elseif strcmp(nmodel,"RIM_Kir_K_L")==0 then 
	model=modelRIM_Kir_K_L;
	constraints=constRIM_Kir_K_L;
end

// Defining the benchmark
benchmark=benchmarkRIM_K_L;
if strcmp(nbench,"RIM_K_L")==0 then 
	benchmark=benchmarkRIM_K_L;
elseif strcmp(nbench,"RIM_Kir_K_L")==0 then 
	benchmark=benchmarkRIM_Kir_K_L;
end

[bM,valBest]=optimization(model,constraints,benchmark);

csvWrite(valBest,"./Results/result.csv")
csvWrite(bM,"./Results/bM.csv")
quit();

