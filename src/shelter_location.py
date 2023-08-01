from graph import Graph
import itertools

import boolexpr as bx

from utils_bx import *
import time


# def expandCNF(cnf_expr, clause):
#     new_expr = []
#     clause = clause.to_cnf()
#     if isinstance(cnf_expr, bx.Or) or isinstance(clause, bx.Or):
#         return bx.and_s(cnf_expr, clause)
#     elif isinstance(clause, bx.And):
#         for e in cnf_expr.args:
#             for c in clause.args:
#                 new_expr.append(bx.or_s(e, c))
    
#     res = bx.and_s(*new_expr)
#     return res

# def expandCNFwSimplify(cnf_expr, clause):
#     new_expr = []
#     clause = clause.to_cnf()
#     if isinstance(cnf_expr, bx.Or) or isinstance(clause, bx.Or):
#         return bx.and_s(cnf_expr, clause)
#     elif isinstance(clause, bx.And):
#         for e in cnf_expr.args:
#             for c in clause.args:
#                 new_expr.append(bx.or_s(e, c))
    
#     res = bx.and_s(*new_expr)
#     return res.to_cnf()


def shelterNumberIdentifier(sink_assign: list, registers: list, M: int):
    # Alan M. Frisch and Paul A. Giannaros. 
    # SAT Encodings of the At-Most-k Constraint - ModRef 2010
    eq1 = []
    N = len(sink_assign)
    for i in range(N - 1):
        eq1.append(bx.or_s(bx.not_(sink_assign[i]), registers[i][0]))

    eq1 = bx.and_s(*eq1)

    eq2 = []
    for j in range(1, M):
        eq2.append(bx.not_(registers[0][j]))

    eq2 = bx.and_s(*eq2)

    eq3 = []
    for i in range(1, N-1):
        for j in range(M):
            eq3.append(bx.or_s(bx.not_(registers[i-1][j]), registers[i][j]))
    
    eq3 = bx.and_s(*eq3)

    eq4 = []
    for i in range(1, N-1):
        for j in range(1, M):
            eq4.append(bx.or_s(bx.not_(sink_assign[i]),
                               bx.not_(registers[i-1][j-1]), 
                               registers[i][j]))
    
    eq4 = bx.and_s(*eq4)

    eq5 = []
    for i in range(N):
        eq5.append(bx.or_s(bx.not_(sink_assign[i]), registers[i-1][-1]))

    eq5 = bx.and_s(*eq5)

    return bx.and_(eq1, eq2, eq3, eq4, eq5).to_cnf()


# def shelterNumberIdentifier(sink_assign: list, M: int):
#     # limit the number of assigned shelters -- DNF
#     # sum of sink_assign <= 2^M

#     const_list = []
#     sum_to_i = [sink_assign[0]]

#     for i in range(50):
#         sum_to_i = bin_add_int([sink_assign[i]], sum_to_i)
#         if (len(sum_to_i) > M):
#             sum_to_i = sum_to_i[:M+1]
#             const_list.append(bx.not_(sum_to_i[M]).to_cnf())
        
#         if(i % 10 == 0):
#             print(i)

    # const_sink_number = []
    # for selected_sink in itertools.combinations(sink_assign, len(sink_assign) - M):
    #     temp = [bx.not_(x) for x in selected_sink]
    #     const_sink_number.append(bx.and_s(*temp))
    
    # const_sink_number = bx.or_(*const_sink_number)
    # return const_sink_number.to_cnf()


def sinkIdentifier(sink_assign: list, sink: list):
    # only 
    if len(sink_assign) != len(sink):
        return bx.ZERO

    # sink_assign[i] >= sink[i]
    # i.e., sink[i] -> sink_assign[i]
    # NOT(sink[i]) OR sink_assign[i]
    ele_and = []
    for i in range(len(sink_assign)):
        ele_and.append(bx.or_s(sink_assign[i], bx.not_(sink[i])))
    
    const_be_sink = bx.and_s(*ele_and)
    return const_be_sink.to_cnf()

    # const_one_sink = []
    # for i, s in enumerate(sink):
    #     # only one of sink[i] can be one
    #     temp = [bx.not_(x) for j, x in enumerate(sink) if i!=j]
    #     temp.append(s)
    #     const_one_sink.append(bx.and_s(*temp))
    
    # const_one_sink = bx.or_(*const_one_sink)
    # const_valid_sink = bx.and_s(const_be_sink, const_one_sink)
    # return const_valid_sink.to_cnf()

def pathIdentifier(graph: Graph, flow: list, src: int, sink: list):
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
    
    all_zero_case = [bx.not_(x) for x in edges_out[src]] + [sink[src]]
    all_zero_case = bx.and_s(*all_zero_case)

    const_src_out = bx.or_(*clauses, all_zero_case)

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

        const_mid_i = bx.or_(bx.or_(*const_mid_i), sink[i])
        const_mid.append(const_mid_i)

    
    # sink node
    # flow in could be either 1 or zero
    # "OR" between each line
    const_sink = []
    # actually only one of these needs to be considered
    # how to select the sink represented by [0,0,1,0,0,0]?
    for i in range(N):
        const_sink_i = [bx.not_(sink[i])]
        # and 
        const_sink_i_in = []
        const_sink_i_out = []
            
        for e in edges_out[i]:
            const_sink_i_out.append(bx.not_(e))

        const_sink_i_out = bx.and_s(*const_sink_i_out)
        const_sink_i.append(const_sink_i_out) # out must be all zero

        # in must sum to one 
        if(i == src):
            clauses = [bx.not_(x) for x in edges_in[i]]
            const_sink_i_in.append(bx.and_s(*clauses))
        else:
            for j, e_i in enumerate(edges_in[i]):
                temp_in = [bx.not_(x) for l, x in enumerate(edges_in[i]) if j!=l]
                clauses = temp_in + [e_i]
                const_sink_i_in.append(bx.and_s(*clauses))
        
        const_sink_i_in = bx.or_(*const_sink_i_in)
        const_sink_i.append(const_sink_i_in)

        const_sink.append(bx.or_(*const_sink_i))

    const_sink_io = bx.and_s(*const_sink)

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

    # generate CNF by purely expanding
    # sink_io_cnf_exp = const_sink[0].to_cnf()
    # for i in range(1, len(const_sink)):
    #     # sink_io_cnf_exp = expandCNFwSimplify(sink_io_cnf_exp, const_sink[i])
    #     sink_io_cnf_exp = expandCNF(sink_io_cnf_exp, const_sink[i])
    # time_list.append(time.perf_counter())
    # print(f"const_sink_io to CNF by expanding:{time_list[-1] - time_list[-2]:0.4f}")

    sink_io_cnf = const_sink_io.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_sink_io to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    ####################################################
    sub_cnfs = [const_src_in_cnf, 
                const_src_out_cnf, 
                all_mid_cnf, 
                sink_io_cnf]
    
    sat_problem = bx.and_s(*sub_cnfs)
    return sat_problem, sub_cnfs





def test_pathIdentifier():
    N = 500
    M = 1
    T = 1 # repeat T times

    ctx = bx.Context()

    # g = Graph(N, 'full', 'false')
    g = Graph(N, '2', 'false')
    print(g.Adj)


    sink_assign = [ctx.get_var(f't{i}') for i in range(N)]
    # sink_assign = [bx.ZERO for i in range(N - 1)]
    # sink_assign.append(ctx.get_var(f't{0}'))
    src = [0,1,2,99]

    flow = []
    for i_f in range(len(src)):
        flow.append([[ctx.get_var(f'x{i_f}_{i}_{j}') for j in range(N)] for i in range(N)])

    sink = []
    for i_s in range(len(src)):
        sink.append([ctx.get_var(f's{i_s}_{i}') for i in range(N)])

    # registers that count to M
    registers = [[ctx.get_var(f'r{i}_{j}') for j in range(M)] for i in range(N)]

    sat_problems = []

    const_sink_number = shelterNumberIdentifier(sink_assign, registers, M)
    sat_problems.append(const_sink_number)

    for i, s in enumerate(src):
        const_valid_flow, _ = pathIdentifier(g, flow[i], s, sink[i])
        const_valid_sink = sinkIdentifier(sink_assign, sink[i])
        sat_problems.append(const_valid_sink)
        sat_problems.append(const_valid_flow)

    
    all_sat_problems = bx.and_s(*sat_problems)

    c_length = [0] * 10
    for c in all_sat_problems.args:
        if isinstance(c, bx.Or):
            c_length[len(c.args)] += 1
        else:
            c_length[1] += 1
    if_sat = all_sat_problems.sat()
    print(if_sat)

    if if_sat[0]:
        print("Shelter locations:")
        for i, s in enumerate(sink_assign):
            if (s != bx.ZERO) and (s != bx.ONE):
                print(f"Node {i}: ", if_sat[1][s])

    return 

if __name__ == '__main__':
    test_pathIdentifier()
    #######################################################
