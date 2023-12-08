import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


from pgmpy.models import BayesianNetwork
from pgmpy.readwrite import UAIWriter
from pgmpy.factors.discrete.CPD import TabularCPD


def generate_random_dag(nodes, edges, max_parents):
    G = nx.DiGraph()
    for i in range(nodes):
        G.add_node(i)
    while edges > 0:
        a = np.random.randint(0, nodes)
        b = a
        while b == a or len(list(G.predecessors(b))) >= max_parents:
            b = np.random.randint(0, nodes)
        if not G.has_edge(a, b) and not nx.has_path(G, b, a):
            G.add_edge(a, b)
            edges -= 1
    return G


def generate_limited_precision_probs(precision):
    num_values = 2 ** precision
    return [i / (num_values - 1) for i in range(num_values)]


def assign_random_cpts(G, precision):
    cpts = {}
    allowed_probs = generate_limited_precision_probs(precision)
    for node in G.nodes():
        num_parents = len(list(G.predecessors(node)))

        cpt_size = 2 ** num_parents
        cpt = np.random.choice(allowed_probs, (cpt_size, 2))
        cpt[:,-1] = 1 - cpt[:,0]  # Normalize

        cpts[node] = cpt
    return cpts


def export_to_uai(G, cpts, filename):
    with open(filename, 'w') as file:
        file.write("BAYES\n")
        file.write(f"{len(G.nodes())}\n")

        file.write(" ".join(["2"] * len(G.nodes())) + "\n")

        file.write(f"{len(G.nodes())}\n")
        for node in G.nodes():
            parents = list(G.predecessors(node))
            scope = sorted([node] + parents)  # Sort the scope in ascending order
            file.write(f"{len(scope)} " + " ".join(map(str, scope)) + "\n")

        file.write("\n")
        for node, cpt in cpts.items():
            file.write(f"{cpt.size}\n")
            file.write(" ".join(map(str, cpt.flatten())) + "\n\n")


def generate_disaster_model(num_nodes, num_edges, precision, max_parents, out_path):

    os.makedirs(out_path, exist_ok=True)
    num_edges = np.min([num_edges, (num_nodes * (num_nodes - 1)) // 2]) # within maximum possible number
    G = generate_random_dag(num_nodes, num_edges, max_parents)
    cpts = assign_random_cpts(G, precision)
    export_to_uai(G, cpts, os.path.join(out_path, "disaster.uai"))

    # Visualize the graph
    nx.draw(G, with_labels=True)
    plt.savefig(os.path.join(out_path, "disaster_bn.png"))


def generate_budget(net_struct, min_bgt, max_bgt, out_path):
    total_nodes = np.sum(net_struct)
    rand_budget = np.random.randint(min_bgt, max_bgt+1, size=total_nodes)

    os.makedirs(out_path, exist_ok=True)
    budget_filename = os.path.join(out_path, f"budget.txt")
    with open(budget_filename, 'w') as fw:
        fw.write(" ".join([str(x) for x in rand_budget]))
    pass


def generate_capacity(net_struct, min_cap, max_cap, p_edge, out_path):
    total_nodes = sum(net_struct)
    adjacency_matrix = np.zeros((total_nodes, total_nodes), dtype=int)

    start_index = 0
    for i, num_nodes in enumerate(net_struct[:-1]):
        next_layer_start = start_index + num_nodes
        next_layer_nodes = net_struct[i + 1]

        for j in range(start_index, next_layer_start):
            for k in range(next_layer_start, next_layer_start + next_layer_nodes):
                if np.random.uniform(0,1) < p_edge:
                    adjacency_matrix[j, k] = np.random.randint(min_cap, max_cap+1)

        start_index += num_nodes

    os.makedirs(out_path, exist_ok=True)
    capacity_filename = os.path.join(out_path, f"capacity.txt")

    with open(capacity_filename, 'w') as file:
        file.write(f"{total_nodes}\n")
        for row in adjacency_matrix:
            row_str = ' '.join(map(str, row))
            file.write(f"{row_str}\n")

    return adjacency_matrix


def generate_cost(adj_matrix, min_cst, max_cst, out_path):
    random_matrix = np.random.randint(min_cst, max_cst + 1, adj_matrix.shape)

    # Set elements to zero where adj_matrix has zeros
    cost_matrix = np.where(adj_matrix != 0, random_matrix, 0)

    os.makedirs(out_path, exist_ok=True)
    capacity_filename = os.path.join(out_path, f"cost.txt")

    with open(capacity_filename, 'w') as file:
        for row in cost_matrix:
            row_str = ' '.join(map(str, row))
            file.write(f"{row_str}\n")


def generate_demand(net_struct, n_nodes, q, out_path):
    total_nodes = sum(net_struct)
    end_nodes = net_struct[-1]

    demand_nodes = np.random.choice(range(total_nodes - end_nodes, total_nodes), n_nodes, replace=False)

    os.makedirs(out_path, exist_ok=True)
    capacity_filename = os.path.join(out_path, f"demand.txt")

    with open(capacity_filename, 'w') as file:
        file.write(f"{n_nodes}\n")
        nodes_str = ' '.join(map(str, demand_nodes))
        file.write(f"{nodes_str}\n")
        file.write(f"{q}\n")


def generate_disaster_edges(adj_matrix, n_edges, out_path):
    disaster_edges = []
    # Find the indices of non-zero elements
    non_zero_indices = np.argwhere(adj_matrix != 0)

    for i in range(np.min([len(non_zero_indices), n_edges])):
        # Uniformly sample one index pair
        sampled_index = non_zero_indices[np.random.choice(len(non_zero_indices))]
        # Get the corresponding value
        assert (adj_matrix[sampled_index[0], sampled_index[1]] != 0)
        disaster_edges.append(sampled_index)

    os.makedirs(out_path, exist_ok=True)
    disaster_edges_filename = os.path.join(out_path, f"disaster.uai.edges")

    with open(disaster_edges_filename, 'w') as file:
        file.write(f"{len(disaster_edges)}\n")
        for i in range(2):
            for d in disaster_edges:
                file.write(f"{d[i]} ")
            file.write(f"\n")


if __name__ == "__main__":
    # supply network
    net_struct = [10, 10, 10]
    out_path = "./network"

    # Generate Budgets
    min_bgt = 100
    max_bgt = 100
    generate_budget(net_struct, min_bgt, max_bgt, out_path)

    # Generate Capacities
    min_cap = 10
    max_cap = 10
    p_edge = 1
    adjacency_matrix = generate_capacity(net_struct, min_cap, max_cap, p_edge, out_path)

    # Generate Costs
    min_cst = 20
    max_cst = 20
    generate_cost(adjacency_matrix, min_cst, max_cst, out_path)

    # Generate Demands
    n_demands = 3
    q = 4
    generate_demand(net_struct, n_demands, q, out_path)

    # Generate Disasters
    num_dedges = 5
    num_bayes_edges = 9
    precision = 4  # k-digit precision
    max_parents = 1000 # maximum of parents allowed

    generate_disaster_edges(adjacency_matrix, num_dedges, out_path)
    generate_disaster_model(num_dedges, num_bayes_edges, precision, max_parents, out_path)
