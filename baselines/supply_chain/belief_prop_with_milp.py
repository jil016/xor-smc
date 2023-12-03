import os.path
import networkx as nx

import numpy as np
from pgmpy.models import BayesianNetwork
from pgmpy.readwrite import UAIWriter
from pgmpy.factors.discrete.CPD import TabularCPD
from pgmpy.inference import BeliefPropagation

def gen(num_of_layers: int, num_of_nodes: list, basepath: str, iter):
    total_nodes = np.sum(num_of_nodes)
    rand_budget = np.random.randint(10, 100, size=total_nodes)
    node_str = "_".join([str(x) for x in num_of_nodes])
    folder_name = f"Layers_{num_of_layers}/Nodes_{node_str}"
    os.makedirs(os.path.join(basepath, folder_name), exist_ok=True)
    budget_filename = os.path.join(basepath, folder_name, f"prog_{iter}_budget.txt")
    with open(budget_filename, 'w') as fw:
        fw.write(" ".join([str(x) for x in rand_budget]))

    adjacency_matrix = np.random.randint(1, 100, size=(total_nodes, total_nodes))
    G = nx.bipartite.random_graph(num_of_nodes[0], num_of_nodes[1], p=0.6)

    Adj = nx.to_numpy_array(G).astype(int)
    print(Adj)

    cap_filename = os.path.join(basepath, folder_name, f"prog_{iter}_capacity.txt")
    all_edges = []
    with open(cap_filename, 'w') as fw:
        fw.write(f"{total_nodes}\n")
        for i in range(total_nodes):
            temp = []
            for j in range(total_nodes):
                if j <= i:
                    temp.append(0)
                    continue
                if Adj[i][j] == 1:
                    temp.append(adjacency_matrix[i][j])
                    all_edges.append("e({},{})".format(i, j))
                else:
                    temp.append(0)
            fw.write(" ".join([str(x) for x in temp]) + '\n')
    cost_matrix = np.zeros((total_nodes, total_nodes))
    for i in range(total_nodes):
        for j in range(total_nodes):
            cost_matrix[i, j] = np.random.randint(1, rand_budget[j])
    cost_filename = os.path.join(basepath, folder_name, f"prog_{iter}_cost.txt")
    with open(cost_filename, 'w') as fw:
        for i in range(total_nodes):
            fw.write(" ".join([str(int(x)) for x in cost_matrix[i]]) + '\n')
    #
    left = np.sum(num_of_nodes[:-1])
    demand_nodes = []
    while len(demand_nodes) == 0:
        rand_bits = np.random.randint(2, size=total_nodes)
        demand_nodes = [i for i in range(left, total_nodes) if rand_bits[i]]
    demand_filename = os.path.join(basepath, folder_name, f"prog_{iter}_demand.txt")
    with open(demand_filename, 'w') as fw:
        fw.write("{}\n".format(len(demand_nodes)))
        fw.write(" ".join([str(x) for x in demand_nodes]) + '\n')
        fw.write("{}".format(sum([int(rand_budget[x] * 0.7) for x in demand_nodes])))

    #
    disaster_filename = os.path.join(basepath, folder_name, f"prog_{iter}_disaster.txt")

    pairs_of_edges = [tuple(np.random.choice(all_edges, size=2)) for _ in range(4)]
    pairs_of_edges = list(set(pairs_of_edges))
    print(pairs_of_edges)
    disaster_graph = BayesianNetwork(pairs_of_edges)

    unique_edges = set([ej for ei, ej in pairs_of_edges])
    base_unit = 32
    for ei in unique_edges:
        rand_pick = np.random.randint(1, base_unit)
        cpd_ei = TabularCPD(ei, 2, [[rand_pick / base_unit], [1 - rand_pick / base_unit]])
        disaster_graph.add_cpds(cpd_ei)

    for ei, ej in pairs_of_edges:
        rand_pick1 = np.random.randint(1, base_unit)
        rand_pick2 = np.random.randint(1, base_unit)
        cond_prob = TabularCPD(ei, 2,
                               [[rand_pick1 / base_unit, 1 - rand_pick1 / base_unit], [rand_pick2 / base_unit, 1 - rand_pick2 / base_unit]],
                               evidence=[ej],
                               evidence_card=[2])
        disaster_graph.add_cpds(cond_prob)

    writer = UAIWriter(disaster_graph)

    writer.write_uai(disaster_filename)


if __name__ == '__main__':
    iter = 0
    while iter < 100:
        try:
            gen(2, [3, 4], "./", iter)
            iter += 1
        except Exception as e:
            print(e)
