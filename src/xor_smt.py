import numpy as np
import shelter_location as sl
import supply_chain as sc
from sympy import *
from graph import Graph
from utils import *
import time


def gen_XOR_constraints(var_list: list, K: int):
    return True


def export_CNF():
    pass


def extract_variable_list(graph: Graph, x: list):
    var_list = []
    for i in range(graph.N):
        for j in range(graph.N):
            if(graph.Adj[i][j] == 1):
                var_list.append(x[i][j])

    return var_list


def shelter_design(phi, q_list, eta, c, N):
    T = 5
    m = 3

    var_list = []
    tgt = [Symbol(f't{i}') for i in range(N)]
    graph_list = []
    for i in range(N):
        graph = Graph(N, 'full')
        graph.forceBreakLoopsBFS(i)
        graph_list.append(graph)

    x_e_list = [[[Symbol(f'x{i}_{j}_{k}') 
                  for k in range(N)]
                  for j in range(N)]
                  for i in range(N)]
    
    var_list = []
    for i in range(N):
        vars = extract_variable_list(graph_list[i], x_e_list[i])
        var_list.append(vars)

    psi_t_list = []
    for t in range(T):
        psi_i_list = []
        for i in range(N):
            psi_i = sl.pathIdentifierUltra(graph_list[i], x_e_list[i], i, tgt, m)
            const_xor = gen_XOR_constraints(var_list[i], q_list[i])
            psi_i = And(psi_i, const_xor)
            psi_i_list.append(psi_i)
        psi_t = And(*psi_i_list)
        psi_t_list.append(psi_t)
    
    psi_star = majority(psi_t_list)
    print(psi_star)
    return psi_star, var_list, tgt


def run_shelter_design():
    tic = time.perf_counter()
    N = 4
    psi_star = q_list=[i for i in range(N)]
    psi_star, var_list, tgt = shelter_design(True, q_list, 0, 0, N)
    psi_star_cnf = to_cnf(psi_star)
    toc = time.perf_counter()
    print(f"N = {N}. Converted to SAT in {toc - tic:0.4f} seconds")

    with open(f'psi_star_N{N}.txt', 'w') as f:
        f.write(str(psi_star))
    
    with open(f'psi_star_cnf_N{N}.txt', 'w') as f:
        f.write(str(psi_star_cnf))
    
    with open(f'variables_N{N}.txt', 'w') as f:
        f.write(str(var_list))
        f.write('\n')
        f.write(str(tgt))

    return


if __name__ == '__main__':
    run_shelter_design()