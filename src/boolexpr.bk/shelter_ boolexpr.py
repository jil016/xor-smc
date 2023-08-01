from graph import Graph
import itertools

import boolexpr as bx
import time


def expandCNF(cnf_expr, clause):
    new_expr = []
    clause = clause.to_cnf()
    if isinstance(cnf_expr, bx.Or) or isinstance(clause, bx.Or):
        return bx.and_s(cnf_expr, clause)
    elif isinstance(clause, bx.And):
        for e in cnf_expr.args:
            for c in clause.args:
                new_expr.append(bx.or_s(e, c))
    
    res = bx.and_s(*new_expr)
    return res

def expandCNFwSimplify(cnf_expr, clause):
    new_expr = []
    clause = clause.to_cnf()
    if isinstance(cnf_expr, bx.Or) or isinstance(clause, bx.Or):
        return bx.and_s(cnf_expr, clause)
    elif isinstance(clause, bx.And):
        for e in cnf_expr.args:
            for c in clause.args:
                new_expr.append(bx.or_s(e, c))
    
    res = bx.and_s(*new_expr)
    return res.to_cnf()


def pathIdentifier(graph: Graph, flow: list, src: int, sink: list, M: int):
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
        const_src_in = bx.and_(const_src_in, bx.not_(e))


    # src output sums to one, in form of DNF
    # (x1 & ~x2 & ~x3) | (~x1 & x2 & ~x3) | (~x1 & ~x2 & x3) | (~x1 & ~x2 & ~x3 & t)
    # or src is one of those sink
    clauses = []
    for i, e in enumerate(edges_out[src]):
        clause = [bx.not_(x) for j, x in enumerate(edges_out[src]) if j!=i]
        clause.append(e)
        clauses.append(bx.and_(*clause))
    
    all_zero_case = [bx.not_(x) for x in edges_out[src]] + [sink[src]]
    all_zero_case = bx.and_(*all_zero_case)

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
        
        const_mid_i = [bx.and_(*const_mid_i)]

        for j, e_i in enumerate(edges_in[i]):
            for k, e_o in enumerate(edges_out[i]):
                temp_in = [bx.not_(x) for l, x in enumerate(edges_in[i]) if j!=l]
                temp_out = [bx.not_(x) for l, x in enumerate(edges_out[i]) if k!=l]
                clauses = temp_in + temp_out + [e_i,e_o]
                const_mid_i.append(bx.and_(*clauses))

        const_mid_i = bx.or_(bx.or_(*const_mid_i), sink[i])
        const_mid.append(const_mid_i)

    
    # sink node
    # flow in could be either 1 or zero
    # "OR" between each line
    const_sink = []
    for i in range(N):
        const_sink_i = [sink[i]]
        # and 
        const_sink_i_in = []
        const_sink_i_out = []
            
        for e in edges_out[i]:
            const_sink_i_out.append(bx.not_(e))

        const_sink_i_out = bx.and_(*const_sink_i_out)
        const_sink_i.append(const_sink_i_out) # out must be all zero

        # in must sum to one 
        if(i == src):
            clauses = [bx.not_(x) for x in edges_in[i]]
            const_sink_i_in.append(bx.and_(*clauses))
        else:
            for j, e_i in enumerate(edges_in[i]):
                temp_in = [bx.not_(x) for l, x in enumerate(edges_in[i]) if j!=l]
                clauses = temp_in + [e_i]
                const_sink_i_in.append(bx.and_(*clauses))
        
        const_sink_i_in = bx.or_(*const_sink_i_in)
        const_sink_i.append(const_sink_i_in)

        const_sink.append(bx.and_(*const_sink_i))

    const_sink_io = bx.or_(*const_sink)

    # limit the number of shelters -- DNF
    const_sink_number = []
    for selected_sink in itertools.combinations(sink, len(sink) - M):
        temp = [bx.not_(x) for x in selected_sink]
        const_sink_number.append(bx.and_(*temp))
    
    const_sink_number = bx.or_(*const_sink_number)

    ################# CNF Conversion ###################
    time_list = []
    time_list.append(time.perf_counter())
    
    const_src_in_cnf = const_src_in.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_src_in to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_src_out_cnf = const_src_out.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_src_out to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    all_mid_cnf = bx.and_(*const_mid)
    all_mid_cnf = all_mid_cnf.to_cnf()
    time_list.append(time.perf_counter())
    print(f"All const_mid to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    # sink_io_dnf = const_sink_io.to_dnf() 
    # sink_io_cnf = sink_io_dnf.to_cnf() 
    # time_list.append(time.perf_counter())
    # print(f"const_sink_io DNF to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    # generate CNF by purely expanding
    sink_io_cnf_exp = const_sink[0].to_cnf()
    for i in range(1, len(const_sink)):
        # sink_io_cnf_exp = expandCNFwSimplify(sink_io_cnf_exp, const_sink[i])
        sink_io_cnf_exp = expandCNF(sink_io_cnf_exp, const_sink[i])

    time_list.append(time.perf_counter())
    print(f"const_sink_io to CNF by expanding:{time_list[-1] - time_list[-2]:0.4f}")

    sink_io_cnf = const_sink_io.to_cnf() 
    time_list.append(time.perf_counter())
    print(f"const_sink_io to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    const_sink_number_cnf = const_sink_number.to_cnf()
    time_list.append(time.perf_counter())
    print(f"const_sink_number to CNF:{time_list[-1] - time_list[-2]:0.4f}")

    ####################################################
    sub_cnfs = [const_src_in_cnf, 
                const_src_out_cnf, 
                all_mid_cnf, 
                sink_io_cnf, 
                const_sink_number_cnf]
    
    sat_problem = bx.and_(*sub_cnfs)
    return sat_problem, sub_cnfs



def test_pathIdentifier():
    N = 100
    M = 1
    ctx = bx.Context()

    # g = Graph(N, 'full', 'false')
    g = Graph(N, '2', 'false')
    print(g.Adj)


    sink = [ctx.get_var(f't{i}') for i in range(N)]
    # for i in range(1,6):
    #     sink[i] = bx.ZERO
    src = [0,1,4]

    flow = []
    for i_f in range(len(src)):
        flow.append([[ctx.get_var(f'x{i_f}_{i}_{j}') for j in range(N)] for i in range(N)])


    sat_problems = []
    for i, s in enumerate(src):
        sat, _ = pathIdentifier(g, flow[i], s, sink, M)
        sat_problems.append(sat)

    all_sat_problems = bx.and_(*sat_problems)
    if_sat = all_sat_problems.sat()
    print(if_sat)

    if if_sat[0]:
        print("Shelter locations:")
        for i, s in enumerate(sink):
            if (s != bx.ZERO) and (s != bx.ONE):
                print(f"Node {i}: ", if_sat[1][s])

    return 

if __name__ == '__main__':
    test_pathIdentifier()
    #######################################################
