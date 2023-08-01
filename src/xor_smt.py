import numpy as np
import shelter_location as sl
from graph import Graph
from utils import *
import time
from datetime import datetime


def gen_XOR_constraints(var_list: list, K: int):
    n_var = len(var_list)
    xor_assignments = np.random.randint(2, size=(K, n_var + 1))

    list_of_xor_terms = []
    for k in range(K):
        temp = [xor_assignments[k,-1]] # add 0 or 1 for the negation
        for i in range(n_var):
            if xor_assignments[k, i] == 1:
                temp.append(var_list[i])
        list_of_xor_terms.append(Xor(*temp))
    
    res = And(*list_of_xor_terms)
    return res


def extract_variable_list(graph: Graph, x: list):
    var_list = []
    for i in range(graph.N):
        for j in range(graph.N):
            if(graph.Adj[i][j] == 1):
                var_list.append(x[i][j])

    return var_list


def shelter_design_for_test(graph, phi, q_list, eta, c, N, T, m):
    var_list = []
    tgt = [Symbol(f't{i}') for i in range(N)]

    x_e_list = [[[Symbol(f'x{i}_{j}_{k}') 
                  for k in range(N)]
                  for j in range(N)]
                  for i in range(N)]
    
    var_list = []
    for i in range(N):
        vars = extract_variable_list(graph, x_e_list[i])
        var_list.append(vars)

    # number of shelter targets <= m
    m_binlist = int2binlist(m)
    tgt_sum = [0]
    for i in range(N):
        tgt_sum = bin_add_int(tgt_sum, [tgt[i]])

    const_sum = bin_leq_int(tgt_sum, m_binlist)

    psi_t_list = []
    for t in range(T):
        psi_i_list = []
        for i in range(N):
            psi_i = sl.pathIdentifierUltra(graph, x_e_list[i], i, tgt)
            const_xor = gen_XOR_constraints(var_list[i], q_list[i])
            psi_i = And(psi_i, const_xor)
            psi_i_list.append(psi_i)
        psi_t = And(*psi_i_list)
        psi_t_list.append(psi_t)
    
    psi_star = And(phi, const_sum, majority(psi_t_list))
    print(psi_star)
    return psi_star, var_list, tgt



def run_shelter_design_test(prefix = 'full', suffix = ''):
    time0 = time.perf_counter()
    N = 8
    T = 1
    m = 3
    graph = Graph(N, prefix, 'false') # graph without loop
    phi = True
    eta = 0
    c = 0
    q_list=[0] * N

    print(graph.Adj)

    psi_star, var_list, tgt = shelter_design_for_test(graph, phi, q_list, eta, c, N, T, m)
    time1 = time.perf_counter()
    time_sat = time1 - time0
    print(f"N = {N}. Converted to SAT problem in {time_sat:0.4f} seconds")

    with open(f'{prefix}_psi_star_N_{N}_{suffix}.txt', 'w') as f:
        f.write(str(psi_star))

    
    with open(f'{prefix}_variables_N_{N}_{suffix}.txt', 'w') as f:
        f.write(str(var_list))
        f.write('\n')
        f.write(str(tgt))
    
    with open(f'{prefix}_params_N_{N}_{suffix}.txt', 'w') as f:
        f.write("Graph: " + str(graph.Adj))
        f.write('\n')
        f.write(f"N = {N}, m = {m}, T = {T} \n")
        f.write(f"phi = {str(phi)}, eta = {eta}, c = {c} \n")
        f.write(f"q_list = {str(q_list)}")

    sat = None
    sat = satisfiable(psi_star)
    time2 = time.perf_counter()
    time_solve = time2 - time1
    print(f"N = {N}. SAT solved in {time_solve:0.4f} seconds")
    
    with open(f'{prefix}_satisfiability_N_{N}_{suffix}.txt', 'w') as f:
        f.write(str(sat))
    

    psi_star_cnf = None
    psi_star_cnf = to_cnf(psi_star)
    time3 = time.perf_counter()
    time_cnf = time3 - time2
    print(f"N = {N}. Converted to CNF in {time_cnf:0.4f} seconds")

    with open(f'{prefix}_psi_star_cnf_N_{N}_{suffix}.txt', 'w') as f:
        f.write(str(psi_star_cnf))


    with open(f'{prefix}_timer_N_{N}_{suffix}.txt', 'w') as f:
        f.write(f"N = {N}. Converted to SAT problem in {time_sat:0.4f} seconds \n")
        f.write(f"N = {N}. Converted to CNF in {time_cnf:0.4f} seconds \n")
        f.write(f"N = {N}. Solve the SAT in {time_solve:0.4f} seconds")
    

    return


if __name__ == '__main__':
    now = datetime.now()
    date_time = now.strftime("%m-%d-%Y-%H_%M_%S")
    run_shelter_design_test(suffix=date_time)
