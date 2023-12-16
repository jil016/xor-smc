import os
import numpy as np
from supply_net_with_disaster import SupplyNet
from docplex.cp.model import CpoModel
from docplex.cp.config import context

context.solver.agent = 'local'
# context.solver.local.execfile = '/Applications/CPLEX_Studio2211/cpoptimizer/bin/arm64_osx/cpoptimizer'
context.solver.local.execfile = '/home/jinzhao/.local/ibm/ILOG/CPLEX_Studio2211/cpoptimizer/bin/x86-64_linux/cpoptimizer'
context.log_output = None   # sys.stdout
context.verbose = 0         # 0~9


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


def calc_no_disaster_fc_production(supply_net, timelimit=99):
    # no disaster, fully connected
    mdl = CpoModel()

    N = supply_net.num_nodes
    total_capacity = 0
    const_all = []

    candidate_edges = []
    for node in supply_net.demand_node:
        in_cost = 0
        for i in range(N):
            if supply_net.capacity[i][node] > 0:
                var_i_node = mdl.binary_var(f"n_{i}_{node}")
                in_cost += var_i_node * supply_net.cost[i][node]
                total_capacity += var_i_node * supply_net.capacity[i][node]
                candidate_edges.append([i, node])

        const_all.append(in_cost <= supply_net.budget[node])


    # add constraints
    mdl.add(const_all)
    mdl.maximize(total_capacity)

    print("\nSolving model....")
    # mdl.export_model("baseline_sampled.lp")
    msol = mdl.solve(TimeLimit = timelimit)
        
    return msol.get_objective_value()