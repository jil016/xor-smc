import numpy as np
from graph import Graph
from utils_data import *
from pysat.formula import CNF
from cnf2uai import cnf_to_uai
from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling


def calc_heuristic(cnf_file, n_free, n_samples, sampler='gibbs'):
    # Sampling based counter.
    # generate K samples, and ...
    

    pass



def local_search(graph: Graph, sources, sampler = 'gibbs'):
    np.random.seed(10086)
    N = graph.N
    ctx = bx.Context()
    flow = [[ctx.get_var(f'x_{i}_{j}') for j in range(N)]
                            for i in range(N)]


    start = np.random.randint(0, N+1)
    print(f"Random start: {start}")

    best_everseen = [start]

    max_steps = 100
    max_neignbor = 10
    n_samples = 50

    
    for step in range(max_steps):
        print(">> Step: ", step)
        # Evaluate start node

        neighbor = []
        for i in range(N):
            if (graph.Adj[start][i] == 1):
                neighbor.append(i)
                if(len(neighbor) > max_neignbor):
                    break
        
        evaluation = []
        
        for nb in neighbor:
            nb_eval = []
            for s in sources:
                cnf, _ = flow2CNF_nbf(graph, flow, s, nb)
                vars = exportCNF(cnf, "temp_cnf")

                cnf_inst = CNF("temp_cnf.cnf")

                # sample from uai
                st = time.time()
                prob = np.ones(cnf_inst.nv)
                cnf_to_uai(cnf_inst, prob, "temp.uniform.gibbs.uai")
                returned_samples = Gibbs_Sampling("temp.uniform.gibbs.uai", n_samples)

                sampled_time = time.time() - st
                print(f"Gibbs sampling takes {sampled_time}")

                flow_ = flow.copy()

                break
            break
            evaluation.append(nb_eval)
        
        break
    return 



def exact_solver():
    pass


if __name__ == '__main__':
    graph = Graph()
    graph.readFromFile("graphs/graph_hawaii_200.txt")

    local_search(graph, [0, 10, 20])