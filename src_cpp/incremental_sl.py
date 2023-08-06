import os
import numpy as np


def read_log(res_log):
	with open(res_log, 'r') as f:
		line = f.readline().split("\n")[0]		
		if (line == "Infeasible"):
			return False, None, None

		line = f.readline().split(" ")[:-1]	# 0 1 1 0
		shelter_assignment = [int(i) for i in line]

		line = f.readline()	# skip a line
		line = f.readline().split(" ")[:-1]	# 1,2,3
		shelter_index = [int(i) for i in line]
	
	return True, shelter_assignment, shelter_index
	

def find_best_shelter():
	best_shelters_sofar = []
	# prepare params
	graph_path = "./graphs/graph1.txt"
	T = 2
	M = 3
	source = "0,1,2"

	temp_out_log = "temp.log"
	
	# incremental qlist, run experiments with timeout
	for i in range(11):
		out_dir = f"LOG/incremental_{i}"
		qlist = f"{i},{i},{i}"
		os.system(f'timeout 1h ./SMC -graph {graph_path} -T {T} -M {M} -source {source} -qlist {qlist} -output {out_dir} > {temp_out_log}')
	
		res_log = out_dir + "/result.log"
		if (os.path.isfile(res_log)):
			is_feasible, shelter_assignment, shelter_index = read_log(res_log)

			if(is_feasible):
				print("Found shelter plan: ", shelter_index)
				best_shelters_sofar.append([shelter_assignment, shelter_index])
			else:	# Infesible at i
				print(f"Found best shelter plan when i = {i-1}")
				break

		else: # Timeout at i
			print(f"Found best shelter plan when i = {i-1}")
			break
	
	if len(best_shelters_sofar) > 0:
		print(best_shelters_sofar[-1][-1])
	return

if __name__ == "__main__":
    find_best_shelter()
