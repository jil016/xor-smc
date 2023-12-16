import os
import random
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


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
        cpt[:, -1] = 1 - cpt[:, 0]  # Normalize

        cpts[node] = cpt
    return cpts


def export_to_uai(G, cpts, filename):
    with open(filename, 'w') as file:
        file.write("BAYES\n")
        file.write(f"{len(G.nodes())}\n")

        # Assuming binary variables for simplicity
        file.write(" ".join(["2"] * len(G.nodes())) + "\n")

        file.write(f"{len(G.nodes())}\n")
        for node in G.nodes():
            parents = list(G.predecessors(node))
            scope = parents + [node]  # the order matters! the last one X is P(X|others)
            file.write(f"{len(scope)} " + " ".join(map(str, scope)) + "\n")

        file.write("\n")
        for node, cpt in cpts.items():
            file.write(f"{cpt.size}\n")
            # for li in range(cpt.shape[0]):
            file.write(" ".join(f"{num:.2f}" for num in cpt.flatten()) + "\n\n")  # precision
            # file.write("\n")


def generate_disaster_model(num_nodes, num_edges, precision, max_parents, out_path):
    num_edges = np.min([num_edges, (num_nodes * (num_nodes - 1)) // 2])  # within maximum possible number
    G = generate_random_dag(num_nodes, num_edges, max_parents)
    cpts = assign_random_cpts(G, precision)

    os.makedirs(out_path, exist_ok=True)
    export_to_uai(G, cpts, os.path.join(out_path, "disaster.uai"))
    # Visualize the graph
    nx.draw(G, with_labels=True)
    plt.savefig(os.path.join(out_path, "disaster_bn.png"))


def generate_budget(net_struct, min_bgt, max_bgt, out_path):
    total_nodes = np.sum(net_struct)
    rand_budget = np.random.randint(min_bgt, max_bgt + 1, size=total_nodes)

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
                if np.random.uniform(0, 1) < p_edge:
                    adjacency_matrix[j, k] = np.random.randint(min_cap, max_cap + 1)

        start_index += num_nodes

    os.makedirs(out_path, exist_ok=True)
    capacity_filename = os.path.join(out_path, f"capacity.txt")

    with open(capacity_filename, 'w') as file:
        file.write(f"{total_nodes}\n")
        for row in adjacency_matrix:
            row_str = ' '.join(map(str, row))
            file.write(f"{row_str}\n")

    return adjacency_matrix


def generate_extreme_capacity(net_struct, cap_list, p_edge, out_path):
    total_nodes = sum(net_struct)
    adjacency_matrix = np.zeros((total_nodes, total_nodes), dtype=int)

    start_index = 0
    for i, num_nodes in enumerate(net_struct[:-1]):
        next_layer_start = start_index + num_nodes
        next_layer_nodes = net_struct[i + 1]

        for j in range(start_index, next_layer_start):
            for k in range(next_layer_start, next_layer_start + next_layer_nodes):
                if np.random.uniform(0, 1) < p_edge:
                    adjacency_matrix[j, k] = np.random.choice(cap_list, 1)

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

    return demand_nodes

def generate_disaster_edges(adj_matrix, n_edges, out_path):
    non_zero_indices = np.nonzero(adj_matrix)
    non_zero_indices_list = list(zip(non_zero_indices[0], non_zero_indices[1]))

    num_samples = n_edges

    sampled_indices = random.sample(non_zero_indices_list, num_samples)

    os.makedirs(out_path, exist_ok=True)
    disaster_edges_filename = os.path.join(out_path, f"disaster.uai.edges")

    with open(disaster_edges_filename, 'w') as file:
        file.write(f"{len(sampled_indices)}\n")
        for i in range(2):
            for d in sampled_indices:
                file.write(f"{d[i]} ")
            file.write(f"\n")


def sample_edges_only_tonodes(adj_matrix, n_edges, node_list, out_path):
    N = adj_matrix.shape[0]

    candidate_edges = []
    for node in node_list:
        for i in range(N):
            if adj_matrix[i][node] > 0:
                candidate_edges.append([i, node])
    

    n_alloc = np.min([len(candidate_edges), n_edges])
    sampled_indices = random.sample(candidate_edges, n_alloc)
        
    return sampled_indices, n_alloc


def generate_disaster_edges_only_end(adj_matrix, n_edges, out_path):
    end_nodes = []
    other_nodes = []
    for i in range(adj_matrix.shape[0]):
        out_degree = 0
        for j in range(adj_matrix.shape[1]):
            if adj_matrix[i,j] > 0:
                out_degree += 1
                break
        
        if out_degree == 0:
            end_nodes.append(i)
        else:
            other_nodes.append(i)
    
    sampled_edges, n_alloc = sample_edges_only_tonodes(adj_matrix, n_edges, end_nodes, out_path)
    if(n_alloc < n_edges):
        sampled_edges_add,_ = sample_edges_only_tonodes(adj_matrix, n_edges - n_alloc, other_nodes, out_path)
        sampled_edges.extend(sampled_edges_add)

    os.makedirs(out_path, exist_ok=True)
    disaster_edges_filename = os.path.join(out_path, f"disaster.uai.edges")

    with open(disaster_edges_filename, 'w') as file:
        file.write(f"{len(sampled_edges)}\n")
        for i in range(2):
            for d in sampled_edges:
                file.write(f"{d[i]} ")
            file.write(f"\n")


def generate_disaster_edges_only_tonodes(adj_matrix, n_edges, node_list, out_path):
    other_nodes = []
    for i in range(adj_matrix.shape[0]):
        if i not in node_list:
            other_nodes.append(i)
    
    sampled_edges, n_alloc = sample_edges_only_tonodes(adj_matrix, n_edges, node_list, out_path)
    if(n_alloc < n_edges):
        sampled_edges_add,_ = sample_edges_only_tonodes(adj_matrix, n_edges - n_alloc, other_nodes, out_path)
        sampled_edges.extend(sampled_edges_add)

    os.makedirs(out_path, exist_ok=True)
    disaster_edges_filename = os.path.join(out_path, f"disaster.uai.edges")

    with open(disaster_edges_filename, 'w') as file:
        file.write(f"{len(sampled_edges)}\n")
        for i in range(2):
            for d in sampled_edges:
                file.write(f"{d[i]} ")
            file.write(f"\n")


if __name__ == "__main__":
    # supply network
    net_struct = [5, 5, 5, 5]

    # Generate Budgets
    min_bgt = 20
    max_bgt = 50

    # Generate Capacities
    min_cap = 8
    max_cap = 12
    p_edge = 1

    # Generate Costs
    min_cst = 5
    max_cst = 15

    # Generate Demands
    n_demands = 4
    q = 6

    # Generate Disasters
    num_dedges = 35  #
    num_bayes_edges = 242
    precision = 4  # k-digit precision 2^n-1 0~15
    max_parents = 8  # maximum of parents allowed

    out_path = f"./small_sized_testnet"

    # generates everything
    generate_budget(net_struct, min_bgt, max_bgt, out_path)
    adjacency_matrix = generate_capacity(net_struct, min_cap, max_cap, p_edge, out_path)
    generate_cost(adjacency_matrix, min_cst, max_cst, out_path)
    generate_demand(net_struct, n_demands, q, out_path)
    generate_disaster_edges(adjacency_matrix, num_dedges, out_path)
    generate_disaster_model(num_dedges, num_bayes_edges, precision, max_parents, out_path)

    # save parameters
    with open(os.path.join(out_path, "params.txt"), "w") as fp:
        fp.write("net_struct: ")
        [fp.write(f"{x} ") for x in net_struct]
        fp.write("\n")

        # Generate Budgets
        fp.write(f"min_bgt: {min_bgt}, max_bgt: {max_bgt}\n")

        # Generate Capacities
        fp.write(f"min_cap: {min_cap}, max_cap: {max_cap}, p_edge: {p_edge}\n")

        # Generate Costs
        fp.write(f"min_cst: {min_cst}, max_cst: {max_cst}\n")

        # Generate Demands
        fp.write(f"n_demands: {n_demands}, q: {q}\n")

        # Generate Disasters
        fp.write(
            f"num_dedges: {num_dedges}, num_bayes_edges: {np.min([num_bayes_edges, (num_dedges * (num_dedges - 1)) // 2])} / {(num_dedges * (num_dedges - 1)) // 2}\n")
        fp.write(f"precision: {precision}, max_parents: {max_parents}\n")
