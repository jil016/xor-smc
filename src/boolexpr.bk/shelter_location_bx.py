import numpy as np
import boolexpr as bx
from boolexpr import *
from utils_bx import *
from graph import Graph


def pathIdentifier(graph: Graph, flow: list, src: int, tgt: int):
    # Return a logic formula that determine whether the flow represents a valid path from node src to tgt

    # Assume the graph has **NO LOOP**! 
    # Otherwise flow and path are not one-to-one corresponded

    constraints = []
    # include those constraints
    # src & tgt should also be encoded
    one = [bx.ONE]
    zero = [bx.ZERO]

    # what if src == tgt 
    # then simply examine whether the rest are zero
    const_s_eq_t = False
    if (src == tgt):
        for i in range(graph.N):
            for j in range(graph.N):
                if (graph.Adj[i][j] == 1):
                    const_s_eq_t = const_s_eq_t | flow[i][j]

        return ~const_s_eq_t

    src_out = [bx.ZERO]
    src_in = [bx.ZERO]
    for i in range(graph.N):
        if (graph.Adj[src][i] == 1):  # there is an edge from src to i
            src_out = bin_add_int([flow[src][i]], src_out)
        if (graph.Adj[i][src] == 1):  # there is an edge from i to src
            src_in = bin_add_int([flow[i][src]], src_in)

    # src node has zero in-flow, unit out-flow
    const_src = bin_eq_int(src_out, one) & bin_eq_int(src_in, zero)
    constraints.append(const_src)

    tgt_out = [bx.ZERO]
    tgt_in = [bx.ONE]
    for i in range(graph.N):
        if (graph.Adj[tgt][i] == 1):  # there is an edge from tgt to i
            tgt_out = bin_add_int([flow[tgt][i]], tgt_out)
        if (graph.Adj[i][tgt] == 1):  # there is an edge from i to tgt
            tgt_in = bin_add_int([flow[i][tgt]], tgt_in)

    # tgt node has zero out-flow, unit in-flow
    const_tgt = bin_eq_int(tgt_out, zero) & bin_eq_int(tgt_in, one)
    constraints.append(const_tgt)

    for i in range(graph.N):
        if (i == src) or (i == tgt):  # other nodes
            continue
        node_in = [bx.ZERO]
        node_out = [bx.ZERO]
        for j in range(graph.N):
            if (graph.Adj[i][j] == 1):
                node_out = bin_add_int([flow[i][j]], node_out)
            if (graph.Adj[j][i] == 1):
                node_in = bin_add_int([flow[j][i]], node_in)
        # Each middle node should have equal in-flow and out-flow
        const_node_io = bin_eq_int(node_out, node_in)
        # The maximum in-flow or out-flow of a node is 1, to guarantee single visit.
        # If removed there will be finite loops, i.e., each node can be visited finite many times in a path
        const_node_once = bin_leq_int(node_out, one)

        constraints.append(const_node_io)
        constraints.append(const_node_once)

    # how to recover the path?

    res = None
    for x in constraints:
        if not res:
            res=x
        else:
            res=bx.and_(res, x)
    return res


def pathIdentifierUltra(graph: Graph, flow: list, src: int, tgt: list, m_tgt: int):
    # tgt is also encoded as {0,1}^N symbol vector
    # path can terminate at any one of those target node in tgt

    # sum over all tgt <= m_tgt
    m_binlist = int2binlist(m_tgt)
    const_sum = bin_leq_int(m_binlist, m_binlist)

    res_list = None
    for t_idx, t in enumerate(tgt):
        res_t = t & pathIdentifier(graph, flow, src, t_idx)
        if not res_list:
            res_list = res_t
        else:
            res_list = res_t

    res = const_sum & res_list
    return res


def recoverPath(graph: Graph, flow: list, src: int, tgt: int):
    # recover a path from src to tgt

    assert (NotImplementedError)


def detectLoop(graph: Graph, flow: list, src: int, tgt: int):
    # Given flow, detect loops.
    # Assume the flow already represents a valid path, 
    # but there could be unnecessary loops that are disjoint from path
    assert (NotImplementedError)


def hasPath(graph: Graph, flow: list, src: int, tgt: int):
    assert (NotImplementedError)


def pathIdentifierTester():
    N = 4
    x_e = [[Symbol(f'x{i}_{j}') for j in range(N)] for i in range(N)]

    graph = Graph(N, 'empty')
    for i in range(N - 1):
        graph.addEdge(i, i + 1)
    graph.addEdge(0, 3)

    x_e = [[0, 1, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1],
           [0, 0, 0, 0]]

    src = 1
    tgt = N - 1
    print(pathIdentifier(graph, x_e, src, tgt))

    x_e = [[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]]
    src = 0
    tgt = 0
    print(pathIdentifier(graph, x_e, src, tgt))


def pathIdentifierUltraTester():
    N = 4
    x_e = [[Symbol(f'x{i}_{j}') for j in range(N)] for i in range(N)]
    tgt = [Symbol(f't{i}') for i in range(N)]

    graph = Graph(N, 'empty')
    for i in range(N - 1):
        graph.addEdge(i, i + 1)
    graph.addEdge(0, 3)

    x_e = [[0, 1, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1],
           [0, 0, 0, 0]]  # flow from 0->1->2->3

    src = 0
    tgt = [0, 0, 1, 1]
    m = 3

    print(pathIdentifierUltra(graph, x_e, src, tgt, m))

    x_e = [[0, 0, 0, 1],
           [0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]]  # flow from 0->3

    src = 0
    tgt = [0, 0, 1, 1]
    m = 3
    print(pathIdentifierUltra(graph, x_e, src, tgt, m))

    x_e = [[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]]  # flow from 0->3

    src = 0
    tgt = [1, 0, 1, 1]
    m = 3
    print(pathIdentifierUltra(graph, x_e, src, tgt, m))


# Driver Code
if __name__ == '__main__':
    # pathIdentifierTester()
    pathIdentifierUltraTester()
