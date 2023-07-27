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
 
    def addEdge(self, u, v):
        # Add v to u’s list.
        self.M += 1
        self.adj[u].append(v)
        self.Adj[u][v] = 1
 
    # Returns count of paths from 's' to 'd'

 
    # A recursive function to print all paths
    # from 'u' to 'd'. visited[] keeps track
    # of vertices in current path. path[]
    # stores actual vertices and path_index
    # is current index in path[]
    def countPathsUtil(self, u, d,
                       visited, pathCount):
        visited[u] = True
 
        # If current vertex is same as
        # destination, then increment count
        if (u == d):
            pathCount[0] += 1
 
        # If current vertex is not destination
        else: 
            # Recur for all the vertices
            # adjacent to current vertex
            i = 0
            while i < len(self.adj[u]):
                if (not visited[self.adj[u][i]]):
                    self.countPathsUtil(self.adj[u][i], d,
                                        visited, pathCount)
                i += 1
 
        visited[u] = False
        return
    
    def countEdgeOncePaths(self, ):
        # Count paths that pass each edge at most once.

        return
    
    def countLoopFreePaths(self, s, d):
        # Mark all the vertices
        # as not visited
        visited = [False] * self.N
 
        # Call the recursive helper
        # function to print all paths
        pathCount = [0]
        self.countPathsUtil(s, d, visited, pathCount)
        return pathCount[0]



def pathIdentifier(flow: list, hasloop = 'false'):
    return



# Driver Code
if __name__ == '__main__':
 
    # Create a graph given in the
    # above diagram
    g = Graph(4, 'full')
    # g.addEdge(0, 1)
    # g.addEdge(1, 2)
    # g.addEdge(2, 3)
    # g.addEdge(3, 1)
    # g.addEdge(1, 4)
    s = 0
    d = 3
     
    # Function call
    print(g.countLoopFreePaths(s, d))
    print(g.N, g.M)
