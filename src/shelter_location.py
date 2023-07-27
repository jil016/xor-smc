import numpy as np
from sympy import *

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
        # this is fully at random
        if(mode == 'full'):
            self.Adj = [[1 for j in range(self.N)] for i in range(self.N)]
            for i in range(self.N):
                for j in range(i + 1):
                    self.Adj[i][j] = 0 # no self loop

            self.M = np.sum(self.Adj)
        if(mode == 'random'):
            n_hierarchy = np.random.randint(low=0, high=self.N)
            
        
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
        # Mark all the vertices
        # as not visited
        visited = [False] * self.N
 
        # Call the recursive helper
        # function to print all paths
        pathCount = 0
        queue = [s]
        visited_que = [visited]
        while queue:
            u = queue.pop()
            visited = visited_que.pop()

            visited[u] = True
            # If current vertex is same as
            # destination, then increment count
            if (u == d):
                pathCount += 1
            # If current vertex is not destination
            else: 
                # Recur for all the vertices
                # adjacent to current vertex
                for v, value in enumerate(self.Adj[u]):
                    if (value == 1) & (not visited[v]):
                        queue.append(v)
                        visited_que.append([j for j in visited])

        return pathCount



def pathIdentifier(graph: Graph, flow: list, hasloop = 'false'):
    # return a logic formula that determine whether it is a valid path
    return


def graphTester():
    # Create a graph given in the
    # above diagram
    g = Graph(4, 'full', 'false')
    # g.addEdge(0, 1)
    # g.addEdge(1, 2)
    # g.addEdge(2, 3)
    # g.addEdge(3, 1)
    # g.addEdge(1, 4)
    s = 0
    d = 3
     
    # Function call
    print(g.countLoopFreePathsDFS(s, d))
    print(g.N, g.M)
    print(g.Adj)
    return

# Driver Code
if __name__ == '__main__':
    M = 10
    x_e = [Symbol(f'x_e{i}') for i in range(M)]
    graphTester()
