
import numpy as np
import boolexpr as bx
import time
from graph import Graph
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx 

from utils_data import *





def reduceCycle(graph: Graph):
    checkGraphDetails(graph)

def reduceBackForthEdges(graph: Graph):
    checkGraphDetails(graph)


def modifyGraph():
    graph = Graph()
    graph.readFromFile("graphs/graph_hawaii_1000.txt")
    
    n_to_rm = 100
    max_in = 3
    max_out = 3

    for _ in range(n_to_rm):
        rand_remove = np.random.randint(low=0, high=graph.N)
        print(rand_remove)
        graph.removeNode(rand_remove)
    
    graph.randomReduceDegree(max_in, max_out)

    # print(graph.Adj)
    extractLargestComponent(graph, "graphs/graph_hawaii_1000_r100_d3.txt")
    graph2 = Graph()
    graph2.readFromFile("graphs/graph_hawaii_1000_r100_d3.txt")
    checkGraphDetails(graph2)
    pass



if __name__ == '__main__':
    modifyGraph()