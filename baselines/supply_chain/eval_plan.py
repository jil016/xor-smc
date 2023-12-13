import os
import numpy as np
from supply_net_with_disaster import SupplyNet


def find_reachable_nodes(start_nodes, adjacency_matrix):
    num_nodes = len(start_nodes)
    reachable = np.array(start_nodes, dtype=bool)
    matrix_power = np.array(adjacency_matrix, dtype=bool)

    # Iterate until no new nodes are reached
    while True:
        new_reachable = reachable | matrix_power.any(axis=0)
        if np.array_equal(new_reachable, reachable):
            break

        reachable = new_reachable
        matrix_power = np.dot(matrix_power, adjacency_matrix) > 0

    return reachable.astype(int)


def calc_actual_production(supply_net, trade_plan, disaster_sample):
    # trade_plan can be num_edges x 1 or N x N matrix
    # don't verify budget
    # only check connectivity!

    # # Tester: for maximum inference
    # trade_plan = [1] * len(trade_plan)
    # disaster_sample = [0] * len(disaster_sample)

    if trade_plan is None:
        return 0

    capacity_matrix = supply_net.capacity.copy()
    connection_list = [0] * supply_net.num_nodes

    if isinstance(trade_plan, list):
        for i, e in enumerate(supply_net.edges):
            if(trade_plan[i] == 0):
                capacity_matrix[e[0], e[1]] = 0
    else:
        # N x N numpy array
        capacity_matrix[trade_plan == 0] = 0

    for i, de in enumerate(supply_net.disaster_edges):
        if(disaster_sample[i] == 1):
            capacity_matrix[de[0], de[1]] = 0

    for i in range(supply_net.num_nodes):
        in_degree = 0
        for j in range(supply_net.num_nodes):
            in_edge_idx = supply_net.edge_map[j, i]
            if (in_edge_idx > -1):
                in_degree += 1

        if (in_degree == 0):
            connection_list[i] = 1

    reachable_nodes = find_reachable_nodes(connection_list, capacity_matrix)

    total_production = 0
    for node in supply_net.demand_node:
        for i in range(supply_net.num_nodes):
            if(capacity_matrix[i, node] > 0 and reachable_nodes[i] == 1):
                total_production += capacity_matrix[i, node]

    return total_production