import boolexpr as bx

import shelter_location_bx as sl
from graph import Graph
from utils_bx import majority
import time


def gen_XOR_constraints(var_list: list, K: int):
    return True


def export_CNF():
    pass


def extract_variable_list(graph: Graph, x: list):
    var_list = []
    for i in range(graph.N):
        for j in range(graph.N):
            if (graph.Adj[i][j] == 1):
                var_list.append(x[i][j])

    return var_list


def shelter_design(phi, q_list, eta, c, N):
    ctx = bx.Context()
    T = 1
    m = 3

    var_list = []
    tgt = [ctx.get_var(f't{i}') for i in range(N)]
    graph_list = []
    for i in range(N):
        graph = Graph(N, 'full')
        graph.forceBreakLoopsBFS(i)
        graph_list.append(graph)

    x_e_list = [[[ctx.get_var(f'x{i}_{j}_{k}')
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
            psi_i = psi_i & const_xor
            psi_i_list.append(psi_i)
        psi_t = None
        for x in psi_i_list:
            if not psi_t:
                psi_t=x
            else:
                psi_t=psi_t & x
        psi_t_list.append(psi_t)

    # psi_star = majority(psi_t_list)
    print(psi_t_list[0])
    print("psi_star cnf")
    psi_star_cnf = psi_t_list[0].to_cnf()
    print(psi_star_cnf)
    print("done")
    return psi_t_list, var_list, tgt


def run_shelter_design():
    # tic = time.time()
    N = 4
    q_list = [i for i in range(N)]
    psi_star, var_list, tgt = shelter_design(True, q_list, 0, 0, N)


    # toc = time.time()
    # print(f"N = {N}. Converted to SAT in {toc - tic:0.4f} seconds")

    # with open(f'psi_star_N{N}.txt', 'w') as f:
    #
    #     # f.write(psi_star[0])
    #     print("psi_star")
    #     print(psi_star[0])
    # print('done1')

    with open(f'psi_star_cnf_N{N}.txt', 'w') as f:
        psi_star_cnf=psi_star[0].to_cnf()
        print("psi_star cnf")
        print(psi_star_cnf)
        # f.write(psi_star_cnf)
        # print('done2')

    with open(f'variables_N{N}.txt', 'w') as f:
        f.write(str(var_list))
        f.write('\n')
        f.write(str(tgt))
        print('done3')

    return


if __name__ == '__main__':
    run_shelter_design()
