import numpy as np
from sympy import *
from utils import *

# np.random.seed(10086)


class Graph:
    def __init__(self, N, mode='empty', allow_loop='true'):
        self.N = N
        self.M = 0
        
        if(mode == 'empty'):
            self.Adj = [[0 for j in range(self.N)] for i in range(self.N)]
        elif(allow_loop == 'true'):
            self.normalInit(mode)
        else:
            self.loopFreeInit(mode)

    def normalInit(self, mode):
        if(mode == 'full'):
            self.Adj = [[1 for j in range(self.N)] for i in range(self.N)]
            for i in range(self.N):
                self.Adj[i][i] = 0 # no self loop

            self.M = np.sum(self.Adj)
        
        if(mode == 'random'):
            self.Adj = np.random.randint(low=0, high=2, 
                                         size=(self.N,self.N),
                                         dtype=int).tolist()
            self.M = np.sum(self.Adj)

    def loopFreeInit(self, mode):
        if(mode == 'full'):
            self.Adj = [[1 for j in range(self.N)] for i in range(self.N)]
            for i in range(self.N):
                for j in range(i + 1):
                    self.Adj[i][j] = 0 # no self loop

            self.M = np.sum(self.Adj)
        if(mode == 'random'):
            self.Adj = [[0 for j in range(self.N)] for i in range(self.N)]

            idx = [i for i in range(self.N)]
            n_hierarchy = np.random.randint(low=0, high=self.N)

            perm = np.random.permutation(idx[:self.N-1])[:n_hierarchy]
            perm = np.sort(perm)

            hierarchies = []

            l = 0
            for p in perm:
                hierarchies.append(idx[l:p+1])
                l = p+1
            hierarchies.append(idx[l:])

            for hierarchy in hierarchies:
                for u in hierarchy:
                    if (hierarchy[-1] != self.N - 1):
                        for v in idx[hierarchy[-1] + 1:]:
                            self.Adj[u][v] = np.random.randint(low=0, high=2)

            self.M =np.sum(self.Adj)
        return

    def addEdge(self, u, v):
        # Add v to u’s list.
        self.M += 1
        self.Adj[u][v] = 1
        return
 
   
    def countEdgeOncePaths(self, ):
        # Count paths that pass each edge at most once.

        return
    
    def countLoopFreePathsDFS(self, s, d):
        # Find all paths by DFS, loop allowed
        visited = [False] * self.N
 
        pathCount = 0
        queue = [s]
        visited_que = [visited]
        while queue:
            u = queue.pop()
            visited = visited_que.pop()

            visited[u] = True
            if (u == d):
                pathCount += 1
            else: 
                for v, value in enumerate(self.Adj[u]):
                    if (value == 1) & (not visited[v]):
                        queue.append(v)
                        visited_que.append([j for j in visited])

        return pathCount



def pathIdentifier(graph: Graph, flow: list, src: int, tgt: int):
    # return a logic formula that determine whether it is a valid path
    # Assume the graph has **NO LOOP**!
    # from node src to tgt

    constraints = []
    # include those constraints
    # src & tgt should also be encoded
    one = [1]
    zero = [0]

    src_out = [0]
    src_in = [0]
    for i in range(graph.N):
        if(graph.Adj[src][i] == 1): # there is an edge from src to i
            src_out = bin_add_int([flow[src][i]], src_out)
        if(graph.Adj[i][src] == 1): # there is an edge from i to src
            src_in = bin_add_int([flow[i][src]], src_in)
    
    const_src = And(bin_eq_int(src_out, one), bin_eq_int(src_in, zero))
    constraints.append(const_src)

    tgt_out = [0]
    tgt_in = [0]
    for i in range(graph.N):
        if(graph.Adj[tgt][i] == 1): # there is an edge from tgt to i
            tgt_out = bin_add_int([flow[tgt][i]], tgt_out)
        if(graph.Adj[i][tgt] == 1): # there is an edge from i to tgt
            tgt_in = bin_add_int([flow[i][tgt]], tgt_in)
    
    const_tgt = And(bin_eq_int(tgt_out, zero), bin_eq_int(tgt_in, one))
    constraints.append(const_tgt)
    
    for i in range(graph.N):
        if(i == src) or (i == tgt): # other nodes
            continue
        node_in = [0]
        node_out = [0]
        for j in range(graph.N):
            if(graph.Adj[i][j] == 1):
                node_out = bin_add_int([flow[i][j]], node_out)
            if(graph.Adj[j][i] == 1):
                node_in = bin_add_int([flow[j][i]], node_in)
                
        const_node_io = bin_eq_int(node_out, node_in)
        const_node_once = bin_leq_int(node_out, one)

        constraints.append(const_node_io)
        constraints.append(const_node_once)

    # how to recover the path?

    res = And(*constraints)
    return res


def pathIdentifierUltra(graph: Graph, flow: list, src: int, tgt: list, m_tgt: int):
    # tgt is also encoded as {0,1}^m symbol vector
    # path can terminate at any one of those tgt
    
    # sum over all tgt <= m_tgt
    m_binlist = int2binlist(m_tgt)
    const_sum = bin_leq_int(m_binlist, m_binlist)
    
    res_list = []
    for t_idx, t in enumerate(tgt):
        res_t = And(t, pathIdentifier(graph, flow, src, t_idx))
        res_list.append(res_t)
    
    res = And(const_sum, Or(*res_list))
    return res



def graphTester():
    # Create a graph given in the
    # above diagram
    N = 4
    g = Graph(N, 'random', 'false')
    # g.addEdge(0, 1)
    # g.addEdge(1, 2)
    # g.addEdge(2, 3)
    # g.addEdge(3, 1)
    # g.addEdge(1, 4)
    s = 0
    d = N - 1
     
    # Function call
    print(g.countLoopFreePathsDFS(s, d))
    print("N: ", g.N, "; M: ", g.M)
    print(g.Adj)
    return


def pathIdentifierTester():
    N = 4
    x_e = [[Symbol(f'x{i}_{j}') for j in range(N)] for i in range(N)]

    graph = Graph(N, 'empty')
    for i in range(N - 1):
        graph.addEdge(i,i+1)
    graph.addEdge(0,3)

    x_e = [[0,1,0,1],
           [0,0,1,0],
           [0,0,0,1],
           [0,0,0,0]]

    src = 0
    tgt = N-1
    print(pathIdentifier(graph, x_e, src, tgt))



def pathIdentifierUltraTester():
    N = 4
    x_e = [[Symbol(f'x{i}_{j}') for j in range(N)] for i in range(N)]
    tgt = [Symbol(f't{i}') for i in range(N)]

    graph = Graph(N, 'empty')
    for i in range(N - 1):
        graph.addEdge(i,i+1)
    graph.addEdge(0,3)

    x_e = [[0,1,0,0],
           [0,0,1,0],
           [0,0,0,1],
           [0,0,0,0]]

    src = 0
    tgt = [0, 0, 1, 1]
    m = 3

    print(pathIdentifierUltra(graph, x_e, src, tgt, m))

    x_e = [[0,0,0,1],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0]]
    

    print(pathIdentifierUltra(graph, x_e, src, tgt, m))


# Driver Code
if __name__ == '__main__':
    pathIdentifierUltraTester()