from subprocess import run, Popen, PIPE
import multiprocessing
import time
import csv

def nbContainers():
	result=run("docker ps -aq | wc -l", shell=True, stdout=PIPE).stdout.decode('utf-8')
	return int(result)-1

models = ["RIM_K_L","RIM_Kir_K_L","RIM_Cap_K_L","RIM_Cap_Kir_K_L"]
benchmarks = ["RIM_K_L","RIM_Kir_K_L","RIM_Cap_K_L","RIM_Cap_Kir_K_L"]

#models = ["RIM_K_L"]
#benchmarks = ["RIM_K_L"]

nbCpus = multiprocessing.cpu_count()
maxExp = 30
simId = 0

for model in models:
	for benchmark in benchmarks:
		noExp = 0
		while noExp < maxExp:
			resultPath = "Results_Model"+model+"_Benchmark"+benchmark+"_"+str(noExp)
			run("mkdir Results/"+resultPath, shell = True)
			command = ["docker run --name scilab"+str(simId)+" --rm  -v $(pwd)/scilab-scripts:/scilab-scripts -v $(pwd)/Results/"+resultPath+":/Results scilab scilab-cli -nb -f main.sce -args "+str(noExp)+" "+model+" "+benchmark]
			Popen(command, shell = True)
			noExp += 1
			simId += 1
			time.sleep(2)
			
			# Loop if the current number of running containers is greater or equal than the number of available cpus
			while nbCpus <= nbContainers():
				time.sleep(1)
		

