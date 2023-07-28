import numpy as np

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
        assert(NotImplementedError)
    
    def countLoopFreePathsDFS(self, s, d):
        # Find all paths by DFS, loop allowed, but each node visited at most once.
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
    
    def breakLoops(self, s):
        assert(NotImplementedError)
    

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


if __name__ == '__main__':
    graphTester()