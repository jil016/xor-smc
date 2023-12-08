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


def generate_disaster_model():
    # Parameters for the disaster distribution
    num_nodes = 5
    num_edges = 12
    num_edges_max = num_nodes * (num_nodes - 1) // 2
    num_edges = np.min([num_edges, num_edges_max])

    precision = 4  # k-digit precision
    max_parents = 1000 # maximum of parents allowed

    G = generate_random_dag(num_nodes, num_edges, max_parents)
    cpts = assign_random_cpts(G, precision)
    export_to_uai(G, cpts, "disaster.uai")

    # Visualize the graph
    nx.draw(G, with_labels=True)
    plt.savefig("disaster.png")


if __name__ == "__main__":
    generate_disaster_model()
