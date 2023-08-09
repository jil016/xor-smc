import numpy as np
from graph import Graph
from utils_data import *
from pysat.formula import CNF
from pysat.solvers import Solver
from cnf2uai import cnf_to_uai
from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling
import os
import sys

def calc_heuristic(cnf_file, n_free, n_samples, sampler='gibbs'):
    # Sampling based counter.
    # generate K samples, and ...

    pass



def draw_from_gibbs_sampling(num_samples, input_file, cnf_instance, prob):
    cnf_to_uai(cnf_instance, prob, input_file + ".uai")
    returned_samples = Gibbs_Sampling(input_file + ".uai", num_samples, burnin_time=20)
    return returned_samples


def draw_from_quicksampler(num_samples, input_file):
    sampled_assignments = []
    formula = CNF(input_file)
    solver = Solver(bootstrap_with=formula.clauses)
    iter = 0
    max_iter = 1000
    while len(sampled_assignments) < num_samples:
        iter += 1
        if(iter > max_iter):
            break
        cmd = """./scripts/sampler/uniformSATSampler/quicksampler -n {} -t 180.0 {} >/tmp/tmp.log""".format(num_samples, input_file)

        os.system(cmd)
        print(iter, len(sampled_assignments), end="\r")
        sys.stdout.flush()
        lines = []
        with open(input_file + '.samples', 'r') as fr:
            for li in fr.read().split("\n"):
                one_sol = li.strip().split(" ")[-1].strip()
                if len(one_sol) <= 1:
                    continue
                lines.append(one_sol)
        os.remove(input_file + '.samples')
        cnt = 0
        for k in lines:
            assignment = []
            for idx, x in enumerate(k):
                if x == '0':
                    assignment.append(-idx - 1)
                else:
                    assignment.append(idx + 1)
            if solver.solve(assignment):
                sampled_assignments.append([1 if int(x) > 0 else 0 for x in k])
            else:
                cnt += 1

    return np.array(sampled_assignments)


def draw_from_unigen(num_samples, input_file):
    tmpfile = '/tmp/unigen.txt'
    cmd = """./scripts/sampler/uniformSATSampler/unigen --input {} --samples {} --sampleout {} > /tmp/tmp.txt""".format(input_file,
                                                                                                                        num_samples,
                                                                                                                        tmpfile)

    os.system(cmd)
    sampled_assignments = []
    with open(tmpfile, 'r') as fr:
        for li in fr.read().split("\n"):
            one_sol = li.strip().split(" ")[:-1]
            # print(one_sol)
            if len(li) <= 1:
                continue
            sampled_assig = [1 if int(x) > 0 else 0 for x in one_sol]
            # print(sampled_assig)
            sampled_assignments.append(sampled_assig)
    return np.array(sampled_assignments)


def shelter_location_local_search(graph: Graph, sources, sampler = 'gibbs'):
    N = graph.N
    ctx = bx.Context()
    


    start = np.random.randint(0, N+1)
    # start = 50 #94 #113


    print(f"Random start: {start}")

    best_everseen = [start]

    max_steps = 100
    max_neighbor = 10
    n_samples = 20
    n_sampling_bits = 14

    neighbors = [start]

    for step in range(max_steps):
        print(">> Step: ", step)
        # Evaluate start node

        for i in range(N):  # collect all neighnors
            if (graph.Adj[start][i] == 1):
                neighbors.append(i)
                if(len(neighbors) > max_neighbor):
                    break
        
        evaluation = []
        
        print(f"Step: {step}. Looking at neighbors of {start}")
        for nb in neighbors: # evaluate each neighbors
            nb_eval = []
            for s in sources:
                flow = [[ctx.get_var(f'x_{i}_{j}') for j in range(N)] for i in range(N)]
                assign_values = []  # 0 -- 0; 1 -- 1; -1 -- keep
                acc_prob = 1.0

                for nsb in range(n_sampling_bits):
                    print(f"Sampling the {nsb} round")
                    cnf, _ = flow2CNF_nbf(graph, flow, s, nb)

                    if(cnf == bx.ZERO):
                        print("Infeasible!")
                        break

                    if not cnf.is_cnf():
                        cnf.to_cnf()

                    if not cnf.sat()[0]:
                        print("Infeasible!")
                        break


                    vars = exportCNF(cnf, "flow_cnf")

                    # sample from cnf
                    st = time.time()

                    if(sampler == 'quick'):
                        returned_samples = draw_from_quicksampler(n_samples, "flow_cnf.cnf")

                    if(sampler == 'unigen'):
                        returned_samples = draw_from_unigen(n_samples, "flow_cnf.cnf")


                    if(sampler == 'gibbs'):
                        cnf_inst = CNF("flow_cnf.cnf")
                        prob = np.ones(cnf_inst.nv)
                        returned_samples = draw_from_gibbs_sampling(n_samples, "gibbs.uniform", cnf_inst, prob)


                    sampled_time = time.time() - st
                    print(f"Sampling takes: {sampled_time}")

                    
                    # generate a sub-counting problem
                    
                    samples_prob = np.sum(returned_samples, axis=0) / returned_samples.shape[0]
                    prob = samples_prob[0]
                    
                    if(prob > 0.5):
                        assign_values.append(1)
                        acc_prob *= prob
                    else:
                        assign_values.append(0)
                        acc_prob *= (1 - prob)

                    # vars[0] should be assigned with assign_values[-1]
                    for i in range(N):
                        for j in range(N):
                            if(str(flow[i][j]) == vars[0]):
                                flow[i][j] = bx.ZERO if (assign_values[-1] == 0) else bx.ONE


                cnf, _ = flow2CNF_nbf(graph, flow, s, nb)
                cnt = 0
                for sat in cnf.iter_sat():
                    # print(sat)
                    cnt += 1
                
                nb_eval.append([cnt, acc_prob])

                            
        
            evaluation.append(nb_eval)
        
        # All neighborss are evaluated
        # neignbor[] <---> evaluation[]
        # find best next
        best_cnt = 0.0
        best_neighbor = start
        for i, nb in enumerate(neighbors):
            eval_nb = []
            for s in range(len(sources)):
                eval_nb.append((evaluation[i][s][0] / evaluation[i][s][1]))

            min_eval_nb = min(eval_nb)
            
            if(min_eval_nb > best_cnt):
                best_cnt = min_eval_nb
                best_neighbor = nb
    
        if (best_neighbor == start):
            print("Finished! No better neighbor!")
            break
        if best_neighbor in best_everseen:
            print("Finished! No better neighbor, trying to go back!")
            break


        best_everseen.append(best_neighbor)
        start = best_neighbor  # move to the next
        neighbors = []  # initialize

    print(best_everseen)
    return best_everseen[-1]



if __name__ == '__main__':
    graph = Graph()
    graph.readFromFile("graphs/graph_hawaii_200.txt")
    np.random.seed(1086)
    best_everseen = []
    run_time_list = []
    for i in range(10):
        st = time.time()
        res = shelter_location_local_search(graph, [0, 10, 20], 'gibbs')
        best_everseen.append(res)
        run_time = time.time() - st
        run_time_list.append(run_time)
    
    print(best_everseen)
    print(run_time_list)