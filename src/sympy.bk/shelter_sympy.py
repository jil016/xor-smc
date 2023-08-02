from sympy import *
from graph import Graph
from sympy.logic.boolalg import is_cnf, is_dnf
from utils_spy import *
from sympy.logic.inference import satisfiable
import itertools

import time


if __name__ == '__main__':
    N = 8
    M = 1
    # g = Graph(N, 'full', 'false')
    g = Graph(8, '2', 'false')

    x_e = [[Symbol(f'x{i}_{j}') for j in range(N)] for i in range(N)]

    print(g.Adj)

    # x_e = [[0,0,0,0],
    #         [0,0,1,0],
    #         [0,0,0,1],
    #         [0,0,0,0]] # flow from 0->1->2->3

    snk = [Symbol(f't{i}') for i in range(N)]
    src = 1

    # preprocess
    edges_in = []
    edges_out = []
    for i in range(N):
        edges_in_i = []
        edges_out_i = []
        for j in range(N):
            if (g.Adj[j][i] == 1):
                edges_in_i.append(x_e[j][i])
            if (g.Adj[i][j] == 1):
                edges_out_i.append(x_e[i][j])
        edges_in.append(edges_in_i)
        edges_out.append(edges_out_i)    


    # src input sums to zero, in form of CNF/DNF: 
    # ~x1 & ~x2 & ~x3 & ...
    const_src_in = True
    for e in edges_in[src]:
        const_src_in = And(const_src_in, Not(e))


    # src output sums to one, in form of DNF
    # (x1 & ~x2 & ~x3) | (~x1 & x2 & ~x3) | (~x1 & ~x2 & x3) | (~x1 & ~x2 & ~x3 & t)
    # or src is one of those snk
    clauses = []
    for i, e in enumerate(edges_out[src]):
        clause = [Not(x) for j, x in enumerate(edges_out[src]) if j!=i]
        clause.append(e)
        clauses.append(And(*clause))
    
    all_zero_case = [Not(x) for x in edges_out[src]] + [snk[src]]
    all_zero_case = And(*all_zero_case)

    const_src_out = Or(*clauses, all_zero_case)


    # middle node, flow in == flow out
    # Each const_mid[i] is a DNF
    # "AND" between each line
    const_mid = []
    for i in range(N):
        if (i == src):
            continue
        const_mid_i = []
        # all zero
        for e in edges_in[i]:
            const_mid_i.append(Not(e))

        for e in edges_out[i]:
            const_mid_i.append(Not(e))
        
        const_mid_i = [And(*const_mid_i)]

        for j, e_i in enumerate(edges_in[i]):
            for k, e_o in enumerate(edges_out[i]):
                temp_in = [Not(x) for l, x in enumerate(edges_in[i]) if j!=l]
                temp_out = [Not(x) for l, x in enumerate(edges_out[i]) if k!=l]
                clauses = temp_in + temp_out + [e_i,e_o]
                const_mid_i.append(And(*clauses))

        const_mid_i = Or(Or(*const_mid_i), snk[i])
        const_mid.append(const_mid_i)

    
    # sink node
    # flow in could be either 1 or zero
    # "OR" between each line
    const_snk = []
    for i in range(N):
        const_snk_i = [snk[i]]
        # and 
        const_snk_i_in = []
        const_snk_i_out = []
            
        for e in edges_out[i]:
            const_snk_i_out.append(Not(e))

        const_snk_i_out = And(*const_snk_i_out)
        const_snk_i.append(const_snk_i_out) # out must be all zero

        # in must sum to one 
        if(i == src):
            clauses = [Not(x) for x in edges_in[i]]
            const_snk_i_in.append(And(*clauses))
        else:
            for j, e_i in enumerate(edges_in[i]):
                temp_in = [Not(x) for l, x in enumerate(edges_in[i]) if j!=l]
                clauses = temp_in + [e_i]
                const_snk_i_in.append(And(*clauses))
        
        const_snk_i_in = Or(*const_snk_i_in)
        const_snk_i.append(const_snk_i_in)

        const_snk.append(And(*const_snk_i))

    const_snk_io = Or(*const_snk)

    # limit the number of shelters -- DNF
    const_snk_number = []
    for selected_snk in itertools.combinations(snk, len(snk) - M):
        temp = [Not(x) for x in selected_snk]
        const_snk_number.append(And(*temp))
    
    const_snk_number = Or(*const_snk_number)
    print(const_snk_number)
    
    sat_problem = And(const_src_in, const_src_out)
    sat_problem = And(sat_problem, *const_mid)
    sat_problem = And(sat_problem, const_snk_io)
    sat_problem = And(sat_problem, const_snk_number)

    time_list = []
    time_list.append(time.perf_counter())
    
    to_cnf(const_src_in) 
    time_list.append(time.perf_counter())
    print(f"const_src_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    to_cnf(const_src_out) 
    time_list.append(time.perf_counter())
    print(f"const_src_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")


    mid_cnf = to_cnf(And(*const_mid)) 
    time_list.append(time.perf_counter())
    print(f"All const_mid to CNF:{time_list[-1] - time_list[-2]:0.4f}")
    
    print("number of clauses in And(*const_mid)")
    print(len(clausesCNF(mid_cnf)))

    snk_io_dnf = to_dnf(const_snk_io) 
    time_list.append(time.perf_counter())
    print(f"(easy) const_snk_io to DNF:{time_list[-1] - time_list[-2]:0.4f}")
    print(snk_io_dnf)

    to_cnf(const_snk_io) 
    time_list.append(time.perf_counter())
    print(f"const_snk_io to CNF:{time_list[-1] - time_list[-2]:0.4f}")


    # to_cnf(const_snk_number) 
    # time_list.append(time.perf_counter())
    # print(f"const_snk_number to CNF:{time_list[-1] - time_list[-2]:0.4f}")



    print(satisfiable(sat_problem))
