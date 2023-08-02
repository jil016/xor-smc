import os
import time
import numpy as np
import boolexpr as bl
import shelter_location as sl
from utils_bx import *
from graph import Graph
from datetime import datetime

def gen_XOR_constraints(var_list: list, K: int):
    # generate XOR constraints that are not self-conflict
    n_var = len(var_list)
    xor_assignments = np.random.randint(2, size=(K, n_var + 1))

    list_of_xor_terms = []
    for k in range(K):
        temp = [xor_assignments[k,-1]] # add 0 or 1 for the negation
        for i in range(n_var):
            if xor_assignments[k, i] == 1:
                temp.append(var_list[i])
        list_of_xor_terms.append(bx.xor_s(*temp))
    
    res = bx.and_s(*list_of_xor_terms)
    return res


def extract_variable_list(graph: Graph, x: list):
    var_list = []
    for i in range(graph.N):
        for j in range(graph.N):
            if(graph.Adj[i][j] == 1):
                var_list.append(x[i][j])

    return var_list


def shelter_design(graph, src, q_list, N, T, M):
    ctx = bx.Context()

    var_list = []
    shelter_assign = [ctx.get_var(f'a{i}') for i in range(N)]

    sat_problems_T = []

    for t in range(T):
        flow = []
        for i_f in range(len(src)):
            flow.append([[ctx.get_var(f'x{t}_s{i_f}_{i}_{j}') for j in range(N)]
                          for i in range(N)])

        reached_shelter = []
        for i_s in range(len(src)):
            reached_shelter.append([ctx.get_var(f's{t}_{i_s}_{i}') for i in range(N)])

        # registers that count to M
        registers = [[ctx.get_var(f'reg{t}_{i}_{j}') for j in range(M)] for i in range(N)]

        sat_problems_t = []

        # number of shelter targets <= m
        const_shelter_number = sl.shelterNumberIdentifier(shelter_assign, registers, M)
        sat_problems_t.append(const_shelter_number)

        for i, s in enumerate(src):
            const_valid_flow, _ = sl.pathIdentifier(graph, flow[i], s, reached_shelter[i])
            const_valid_shelter = sl.sinkIdentifier(shelter_assign, reached_shelter[i])
            # add XOR here (raw XOR without encoding for now) 
            var_list = extract_variable_list(graph, flow[i])
            var_list = var_list + reached_shelter[i]
            const_xor = gen_XOR_constraints(var_list, q_list[i])
            
            sat_problems_t.append(const_valid_shelter)
            sat_problems_t.append(const_valid_flow)
            sat_problems_t.append(const_xor)

        
        sat_problems_t = bx.and_s(*sat_problems_t)
        sat_problems_T.append(sat_problems_t)


    # add majority here among sat_problems_T
    all_sat_problems = bx.and_(*sat_problems_T)
    if_sat = all_sat_problems.sat()

    # print one valid assignment
    if if_sat[0]:
        print("Shelter locations:")
        for i, s in enumerate(shelter_assign):
            if (s != bx.ZERO) and (s != bx.ONE):
                print(f"Node {i}: ", if_sat[1][s])
    else:
        print("Not Satisfiable!")
    
    return if_sat, sat_problems_T

def shelter_location_exact():
    pass

def run_shelter_design(save_to = ''):
    time0 = time.perf_counter()
    N = 100
    T = 2
    M = 3
    map_type = '2'  # max in-degree = max out-degree = 2
    allow_loop = 'false' # graph without loop
    graph = Graph(N, map_type, allow_loop) 

    src = [0,1,2]

    eta = 0
    c = 0
    q_list=[2] * len(src)
    print(graph.Adj)
    sat_problems_T = shelter_design(graph, src, q_list, N, T, M)


    os.mkdir(save_to)
    with open(f'{save_to}//graph_N_{N}_map_{map_type}.txt', 'w') as f:
        f.write(str(graph.Adj))

    with open(f'{save_to}//params_N_{N}_map_{map_type}.txt', 'w') as f:
        f.write("Source nodes: " + str(src))
        f.write('\n')
        f.write(f"N = {N}, M = {M}, T = {T} \n")
        f.write(f"eta = {eta}, c = {c} \n")
        f.write(f"q_list = {str(q_list)}")

    return


if __name__ == '__main__':
    now = datetime.now()
    date_time = now.strftime("%m-%d-%Y-%H_%M_%S")
    result_folder = f"result-{date_time}"
    run_shelter_design(result_folder)
