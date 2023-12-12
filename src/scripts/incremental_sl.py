import os
import numpy as np
import argparse


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
	


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Incremental Shelter Finding')
	parser.add_argument(
		'-g', '--graph', default="./graphs/graph1.txt", help='Graph file required')
	parser.add_argument(
		'-t', '--T', default=2, type=int, help='Number of repeat experiments')
	parser.add_argument(
		'-m', '--M', default=3, type=int, help='Maximum number of shelters')
	parser.add_argument(
		'-s','--source', nargs='+', help='Source nodes (residential areas)', required=True)
	parser.add_argument(
		'-q','--qlist', nargs='+', help='Requirement of the number of paths', required=True)
	parser.add_argument(
		'-o','--output', default="./LOG", help='Output base path')
	parser.add_argument(
		'-l','--timelimit', default="1h", help='Set time limit')

	args = parser.parse_args()

	# params
	graph = args.graph
	T = args.T
	M = args.M
	source = args.source
	qlist = args.qlist
	output = args.output
	timelimit = args.timelimit

	print("Parameters")
	print(f"graph: {graph}")
	print(f"T: {T}")
	print(f"M: {M}")
	print(f"source: {source}")
	print(f"qlist: {qlist}")
	print(f"output: {output}")
	print(f"timelimit: {timelimit}")


	best_shelters_sofar = []
	
	# incremental qlist, run experiments with timeout
	for i in range(10):
		inc_out_dir = output + f"/incremental_{i}"
		temp_cpp_out = output + f"/cpp_out_{i}"

		source_str = ",".join([str(s) for s in source])
		qlist_str = ",".join([str(int(q) + i) for q in qlist])

		os.system(f'timeout {timelimit} ./SMC -graph {graph} -T {T} -M {M} -source {source_str} -qlist {qlist_str} -output {inc_out_dir} > {temp_cpp_out}')
	
		res_log = inc_out_dir + "/result.log"
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


