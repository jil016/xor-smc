# Evaluate the shelter plan

# (1) Put the flow model with given source and sink into a CNF

# (2) Count the number of paths between sources and shelters

import numpy as np
import boolexpr as bx
import time
from graph import Graph
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx 



def flow2CNF(graph: Graph, flow: list, src: int, sink: list):
    # preprocess
    edges_in = []
    edges_out = []
    N = graph.N

    for i in range(N):
        edges_in_i = []
        edges_out_i = []
        for j in range(N):
            if (graph.Adj[j][i] == 1):
                edges_in_i.append(flow[j][i])
            if (graph.Adj[i][j] == 1):
                edges_out_i.append(flow[i][j])
        edges_in.append(edges_in_i)
        edges_out.append(edges_out_i)    

    # src input sums to zero, in form of CNF/DNF: 
    # ~x1 & ~x2 & ~x3 & ...
    const_src_in = bx.ONE
    for e in edges_in[src]:
        const_src_in = bx.and_s(const_src_in, bx.not_(e))

    # src output sums to one, in form of DNF
    # (x1 & ~x2 & ~x3) | (~x1 & x2 & ~x3) | (~x1 & ~x2 & x3) | (~x1 & ~x2 & ~x3 & t)
    # or src is one of those sink
    clauses = []
    for i, e in enumerate(edges_out[src]):
        clause = [bx.not_(x) for j, x in enumerate(edges_out[src]) if j!=i]
        clause.append(e)
        clauses.append(bx.and_s(*clause))
    
    if (src == sink):
        all_zero_case = [bx.not_(x) for x in edges_out[src]]
        all_zero_case = bx.and_s(*all_zero_case)
        const_src_out = bx.or_(*clauses, all_zero_case)
    else:
        const_src_out = bx.or_(*clauses)

    # middle node, flow in == flow out
    # Each const_mid[i] is a DNF
    # "AND" between each line
    const_mid = []
    for i in range(N):
        if (i == src or i == sink):
            continue
        const_mid_i = []
        # all zero
        for e in edges_in[i]:
            const_mid_i.append(bx.not_(e))

        for e in edges_out[i]:
            const_mid_i.append(bx.not_(e))
        
        const_mid_i = [bx.and_s(*const_mid_i)]

        for j, e_i in enumerate(edges_in[i]):
            for k, e_o in enumerate(edges_out[i]):
                temp_in = [bx.not_(x) for l, x in enumerate(edges_in[i]) if j!=l]
                temp_out = [bx.not_(x) for l, x in enumerate(edges_out[i]) if k!=l]
                clauses = temp_in + temp_out + [e_i,e_o]
                const_mid_i.append(bx.and_s(*clauses))

        const_mid_i = bx.or_(*const_mid_i)
        const_mid.append(const_mid_i)

    
    # sink output sums to zero: 
    # ~x1 & ~x2 & ~x3 & ...
    const_sink_out = bx.ONE
    for e in edges_out[sink]:
        const_sink_out = bx.and_s(const_sink_out, bx.not_(e))

    # sink input sums to one, in form of DNF
    # or sink == src
    clauses = []
    for i, e in enumerate(edges_in[sink]):
        clause = [bx.not_(x) for j, x in enumerate(edges_in[sink]) if j!=i]
        clause.append(e)
        clauses.append(bx.and_s(*clause))
    
    if (src == sink):
        all_zero_case = [bx.not_(x) for x in edges_in[sink]]
        all_zero_case = bx.and_s(*all_zero_case)
        const_sink_in = bx.or_(*clauses, all_zero_case)
    else:
        const_sink_in = bx.or_(*clauses)


    ################# CNF Conversion ###################
    time_list = []
    time_list.append(time.perf_counter())
    
    const_src_in_cnf = const_src_in.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_src_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_src_out_cnf = const_src_out.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_src_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    all_mid_cnf = bx.and_s(*const_mid)
    all_mid_cnf = all_mid_cnf.to_cnf()
    time_list.append(time.perf_counter())
    print(f"All const_mid to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    # sink_io_cnf = const_sink_io.to_cnf() 
    # time_list.append(time.perf_counter())
    # print(f"const_sink_io to CNF:{time_list[-1] - time_list[-2]:0.4f}")
    const_sink_in_cnf = const_sink_in.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_sink_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_sink_out_cnf = const_sink_out.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_sink_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    ####################################################
    sub_cnfs = [const_src_in_cnf, 
                const_src_out_cnf, 
                all_mid_cnf, 
                const_sink_in_cnf, 
                const_sink_out_cnf,]
    
    sat_problem = bx.and_s(*sub_cnfs)
    return sat_problem, sub_cnfs


def flow2CNF_nbf(graph: Graph, flow: list, src: int, sink: list):
    # preprocess
    edges_in = []
    edges_out = []
    N = graph.N

    for i in range(N):
        edges_in_i = []
        edges_out_i = []
        for j in range(N):
            if (graph.Adj[j][i] == 1):
                edges_in_i.append(flow[j][i])
            if (graph.Adj[i][j] == 1):
                edges_out_i.append(flow[i][j])
        edges_in.append(edges_in_i)
        edges_out.append(edges_out_i)    

    # src input sums to zero, in form of CNF/DNF: 
    # ~x1 & ~x2 & ~x3 & ...
    const_src_in = bx.ONE
    for e in edges_in[src]:
        const_src_in = bx.and_s(const_src_in, bx.not_(e))

    # src output sums to one, in form of DNF
    # (x1 & ~x2 & ~x3) | (~x1 & x2 & ~x3) | (~x1 & ~x2 & x3) | (~x1 & ~x2 & ~x3 & t)
    # or src is one of those sink
    clauses = []
    for i, e in enumerate(edges_out[src]):
        clause = [bx.not_(x) for j, x in enumerate(edges_out[src]) if j!=i]
        clause.append(e)
        clauses.append(bx.and_s(*clause))
    
    if (src == sink):
        all_zero_case = [bx.not_(x) for x in edges_out[src]]
        all_zero_case = bx.and_s(*all_zero_case)
        const_src_out = bx.or_(*clauses, all_zero_case)
    else:
        const_src_out = bx.or_(*clauses)

    # middle node, flow in == flow out
    # Each const_mid[i] is a DNF
    # "AND" between each line
    const_mid = []
    for i in range(N):
        if (i == src or i == sink):
            continue
        const_mid_i = []
        # all zero
        for e in edges_in[i]:
            const_mid_i.append(bx.not_(e))

        for e in edges_out[i]:
            const_mid_i.append(bx.not_(e))
        
        const_mid_i = [bx.and_s(*const_mid_i)]

        for j, e_i in enumerate(edges_in[i]):
            for k, e_o in enumerate(edges_out[i]):
                temp_in = [bx.not_(x) for l, x in enumerate(edges_in[i]) if j!=l]
                temp_out = [bx.not_(x) for l, x in enumerate(edges_out[i]) if k!=l]
                clauses = temp_in + temp_out + [e_i,e_o]
                const_mid_i.append(bx.and_s(*clauses))

        const_mid_i = bx.or_(*const_mid_i)
        const_mid.append(const_mid_i)
    
    bf_edge_pair = []
    const_nbf = []
    for i in range(N):
        for j in range(i, N):
            if (graph.Adj[j][i] == 1) and (graph.Adj[i][j] == 1):
                bf_edge_pair.append([i,j])
                temp = bx.not_(bx.and_(flow[i][j], flow[j][i]))
                const_nbf.append(temp)
    

    
    
    # sink output sums to zero: 
    # ~x1 & ~x2 & ~x3 & ...
    const_sink_out = bx.ONE
    for e in edges_out[sink]:
        const_sink_out = bx.and_s(const_sink_out, bx.not_(e))

    # sink input sums to one, in form of DNF
    # or sink == src
    clauses = []
    for i, e in enumerate(edges_in[sink]):
        clause = [bx.not_(x) for j, x in enumerate(edges_in[sink]) if j!=i]
        clause.append(e)
        clauses.append(bx.and_s(*clause))
    
    if (src == sink):
        all_zero_case = [bx.not_(x) for x in edges_in[sink]]
        all_zero_case = bx.and_s(*all_zero_case)
        const_sink_in = bx.or_(*clauses, all_zero_case)
    else:
        const_sink_in = bx.or_(*clauses)


    ################# CNF Conversion ###################
    time_list = []
    time_list.append(time.perf_counter())
    
    const_src_in_cnf = const_src_in.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_src_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_src_out_cnf = const_src_out.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_src_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    all_mid_cnf = bx.and_s(*const_mid)
    all_mid_cnf = all_mid_cnf.to_cnf()
    time_list.append(time.perf_counter())
    print(f"All const_mid to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_nbf_cnf = bx.and_s(*const_nbf)
    const_nbf_cnf = const_nbf_cnf.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_nbf to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    # sink_io_cnf = const_sink_io.to_cnf() 
    # time_list.append(time.perf_counter())
    # print(f"const_sink_io to CNF:{time_list[-1] - time_list[-2]:0.4f}")
    const_sink_in_cnf = const_sink_in.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_sink_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_sink_out_cnf = const_sink_out.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_sink_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    ####################################################
    sub_cnfs = [const_src_in_cnf, 
                const_src_out_cnf, 
                all_mid_cnf, 
                const_nbf_cnf,
                const_sink_in_cnf, 
                const_sink_out_cnf,]
    
    sat_problem = bx.and_s(*sub_cnfs)
    return sat_problem, sub_cnfs


def extractVarsFromCNF(expr: bx.And):
    if not(expr.is_cnf()):
        expr.to_cnf()
    
    vars = []

    for clause in expr.args:
        if isinstance(clause, bx.Or):
            for literal in clause.args:
                vars.append(str(literal).replace('~', ''))

        else:
            vars.append(str(bx.not_(clause)).replace('~', ''))
    
    vars = np.unique(vars).tolist()

    return vars


def exportCNF(expr: bx.And, file_name):
    if not(expr.is_cnf()):
        expr.to_cnf()
    
    vars = extractVarsFromCNF(expr)

    with open(file_name + '.vars', 'w') as f:
        f.write(str(vars))
    
    clauses = []

    for clause in expr.args:
        literals = []

        if isinstance(clause, bx.Or):
            for literal in clause.args:
                str_literal = str(literal)
                if (str_literal[0] == '~'):
                    idx = vars.index(str_literal[1:]) + 1
                    literals.append(-idx)
                else:
                    idx = vars.index(str_literal) + 1
                    literals.append(idx)

        else:
            str_literal = str(clause)
            if (str_literal[0] == '~'):
                idx = vars.index(str_literal[1:]) + 1
                literals.append(-idx)
            else:
                idx = vars.index(str_literal) + 1
                literals.append(idx)
        
        clauses.append(literals)

    # write clauses to a file in a standard form
    with open(file_name + '.cnf', 'w') as f:
        f.write(f"p cnf {len(vars)} {len(clauses)}\n")
        for clause in clauses:
            [f.write(f"{l} ") for l in clause]
            f.write("0\n")

    return vars


def checkGraphDetails(graph: Graph):
    in_degrees = []
    out_degrees = []

    count_backandforth = 0
    for i in range(graph.N):
        in_degree = 0
        out_degree = 0

        for j in range(graph.N):
            if(graph.Adj[i][j] == 1):
                out_degree += 1

            if(graph.Adj[j][i] == 1):
                in_degree += 1

            if(graph.Adj[j][i] == 1) and (graph.Adj[i][j] == 1):
                count_backandforth += 1
            
        
        in_degrees.append(in_degree)
        out_degrees.append(out_degree)

    print("Max in-degree: ", max(in_degrees))
    print("Max out-degree: ", max(out_degrees))
    print("Back-and-forth edges: ", count_backandforth)
    print("Total edges: ", graph.M)

    # plt.figure(figsize=(4.5,2.5))
    plt.figure(0)
    sns.histplot(data=in_degrees)
    plt.show()

    plt.figure(1)
    Gd = nx.DiGraph()
    for i in range(graph.N): 
        for j in range(graph.N): 
            if graph.Adj[i][j] == 1: 
                Gd.add_edge(i,j) 

    nx.draw(Gd, node_size=1) 
    plt.show() 

    # plt.figure(1)
    # sns.histplot(data=out_degrees)
    # plt.show()
    return



def extractLargestComponent(graph: Graph, outfile=""):
    Gd = nx.DiGraph()
    for i in range(graph.N): 
        for j in range(graph.N): 
            if graph.Adj[i][j] == 1: 
                Gd.add_edge(i,j) 

    # nx.draw(Gd, node_size=1) 
    # plt.show() 
    sccs = nx.weakly_connected_components(Gd) 
    # sccs = nx.strongly_connected_components(Gd) 
    largest_component = []

    for scc in sccs: 
        if len(scc) > len(largest_component): 
            largest_component = list(scc)


    graph_out = Graph(N = len(largest_component))
    for i in range(len(largest_component)):
        for j in range(len(largest_component)):
            graph_out.Adj[i][j] = graph.Adj[largest_component[i]][largest_component[j]]
    
    graph_out.exportToFile(outfile)
    return

########################################################################################

def test_flow2CNF():
    graph = Graph()
    src = 20
    sink = 90

    graph.readFromFile("graphs/graph_hawaii_200.txt")


    N = graph.N
    
    ctx = bx.Context()

    flow = [[ctx.get_var(f'x_{i}_{j}') for j in range(N)]
                          for i in range(N)]

    cnf, sub_cnfs = flow2CNF_nbf(graph, flow, src, sink)

    exportCNF(cnf, f"src{src}_sink{sink}")

    # sink_list = [bx.ZERO] * N
    # sink_list[sink] = bx.ONE
    # cnf2, sub_cnfs2 = pathIdentifier(graph, flow, src, sink_list)
    
    # if not cnf.is_cnf():
    #     cnf.to_cnf()

    # # print(cnf.sat())
    cnt = 0
    for sat in cnf.iter_sat():
        # print(sat)
        cnt += 1
    
    print(cnt)
    return


def test_checkGraphDetails():
    graph = Graph()
    graph.readFromTSPFile("graphs/tsphcp/SNm_500.hcp")
    # graph.readFromTSPFile("graphs/tsphcp/COL_1000.hcp")
    # graph.readFromTSPFile("graphs/tsphcp/FLS3_408.hcp")
    # graph.readFromHawaiiFile("graphs/hawaii/hawaii_400.txt")
    # graph.exportToFile("graphs/graph_hawaii_400.txt")
    checkGraphDetails(graph)
    return 


def test_extractLargestComponent():
    graph = Graph()
    graph.readFromTSPFile("graphs/tsphcp/SNm_500.hcp")
    extractLargestComponent(graph, "graphs/SNm_500_connected.txt")


if __name__ == "__main__":
    # test_checkGraphDetails()
    test_flow2CNF()
    # 