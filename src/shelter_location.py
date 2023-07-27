import numpy as np
from sympy import *


class Graph:
    def __init__(self, N, mode='empty'):
        self.N = N
        self.M = 0
        self.adj = [[] for i in range(N)]
        self.Adj = [[0 for j in range(N)] for i in range(N)]

        if(mode == 'full'):
            for i in range(N):
                for j in range(N):
                    if (j != i):
                        self.adj[i].append(j)
                        self.Adj[i][j] = 1
                        self.M += 1


    def randomInit(self):
        pass
 
    def addEdge(self, u, v):
        # Add v to u’s list.
        self.M += 1
        self.adj[u].append(v)
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
                for v in self.adj[u]:
                    if (not visited[v]):
                        queue.append(v)
                        visited_que.append([j for j in visited])

        return pathCount



def pathIdentifier(graph: Graph, flow: list, hasloop = 'false'):
    return



def graphTester():
    # Create a graph given in the
    # above diagram
    g = Graph(4, )
    g.addEdge(0, 1)
    g.addEdge(1, 2)
    g.addEdge(2, 3)
    g.addEdge(3, 1)
    # g.addEdge(1, 4)
    s = 0
    d = 3
     
    # Function call
    print(g.countLoopFreePaths(s, d))
    print(g.countLoopFreePaths2(s, d))
    print(g.N, g.M)
    return

# Driver Code
if __name__ == '__main__':
    
    x_e = [Symbol(f'x{i}_{j}') for i in range(n)]
