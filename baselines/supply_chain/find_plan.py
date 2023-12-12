import os
import numpy as np
import argparse
import time
import random
from supply_net_with_disaster import SupplyNet, process_samples

from docplex.cp.model import CpoModel
from docplex.cp.config import context

context.solver.agent = 'local'
context.solver.local.execfile = '/Applications/CPLEX_Studio2211/cpoptimizer/bin/arm64_osx/cpoptimizer'
# context.solver.local.execfile =
# '/home/jinzhao/.local/ibm/ILOG/CPLEX_Studio2211/cpoptimizer/bin/x86-64_linux/cpoptimizer'


def mip_find_plan(supply_net, disaster_sample, find_best=False):
    mdl = CpoModel()

    # trade plan variables
    trade_plan = [mdl.binary_var(f"s_{eg[0]}_{eg[1]}") for eg in supply_net.edges]

    # connection variables
    x_nodes = [mdl.binary_var(f"n_{i}") for i in range(supply_net.num_nodes)]
    x_edges = [mdl.binary_var(f"e_{eg[0]}_{eg[1]}") for eg in supply_net.edges]

    # connection constraints
    edge_connection_const = []
    node_connection_const = []
    for i in range(supply_net.num_nodes):
        in_degree = 0
        for j in range(supply_net.num_nodes):
            out_edge_idx = supply_net.edge_map[i, j]
            dedge_idx = supply_net.dedge_map[i, j]

            no_dis = 1
            if (dedge_idx != -1 and disaster_sample[dedge_idx] == 1):
                no_dis = 0

            if (out_edge_idx > -1):
                edge_connection_const.append(x_edges[out_edge_idx] <= x_nodes[i])
                edge_connection_const.append(x_edges[out_edge_idx] <= trade_plan[out_edge_idx])
                edge_connection_const.append(x_edges[out_edge_idx] <= no_dis)  # can be pruned ...
                edge_connection_const.append(
                    x_edges[out_edge_idx] >= x_nodes[i] + trade_plan[out_edge_idx] + no_dis - 2)

            in_edge_idx = supply_net.edge_map[j, i]
            if (in_edge_idx > -1):
                node_connection_const.append(x_nodes[i] >= x_edges[in_edge_idx])
                in_degree += 1

        if (in_degree == 0):
            node_connection_const.append(x_nodes[i] == 1)

    # budget constraint
    budget_const = []
    for i in range(supply_net.num_nodes):
        total_cost = 0
        for j in range(supply_net.num_nodes):
            in_edge_idx = supply_net.edge_map[j, i]

            if (in_edge_idx > -1):
                total_cost += supply_net.cost[j, i] * trade_plan[in_edge_idx]

        budget_const.append((total_cost <= supply_net.budget[i]))

    # objective
    total_production = 0
    for edge in supply_net.edges:
        if (edge[1] in supply_net.demand_node):
            total_production += supply_net.capacity[edge[0], edge[1]] * x_edges[supply_net.edge_map[edge[0], edge[1]]]

    production_const = (total_production >= 2 ** supply_net.total_demand)
    production_const = (total_production >= 128)

    # optional!! maximize ...
    if (find_best):
        mdl.maximize(total_production)

    # add constraints
    for c in edge_connection_const:
        mdl.add(c)
    for c in node_connection_const:
        mdl.add(c)
    for c in budget_const:
        mdl.add(c)

    mdl.add(production_const)

    print("\nSolving model....")
    mdl.export_model("baseline_sampled.lp")
    msol = mdl.solve() #(TimeLimit=10)

    # msol.stop_cause = 'SearchStoppedByLimit'
    # msol.solve_status = 'Feasible'

    if msol:
        print("Solution status: " + msol.get_solve_status())
        for edge in trade_plan:
            print("   " + edge.get_name() + f": {msol[edge.get_name()]}")
    else:
        print("No solution found")
    return [msol[edge.get_name()] for edge in trade_plan]


def saa_find_plan(supply_net, disaster_samples, find_best=False, time_limit="2h"):
    mdl = CpoModel()

    # trade plan variables, for all samples
    trade_plan = [mdl.binary_var(f"s_{eg[0]}_{eg[1]}") for eg in supply_net.edges]


    x_nodes = []

    for idx_sample, sample in enumerate(disaster_samples):
        # connection variables
        x_nodes = [mdl.binary_var(f"n_{i}") for i in range(supply_net.num_nodes)]
        x_edges = [mdl.binary_var(f"e_{eg[0]}_{eg[1]}") for eg in supply_net.edges]

        # connection constraints
        edge_connection_const = []
        node_connection_const = []
        for i in range(supply_net.num_nodes):
            in_degree = 0
            for j in range(supply_net.num_nodes):
                out_edge_idx = supply_net.edge_map[i, j]
                dedge_idx = supply_net.dedge_map[i, j]

                no_dis = 1
                if (dedge_idx != -1 and sample[dedge_idx] == 1):
                    no_dis = 0

                if (out_edge_idx > -1):
                    edge_connection_const.append(x_edges[out_edge_idx] <= x_nodes[i])
                    edge_connection_const.append(x_edges[out_edge_idx] <= trade_plan[out_edge_idx])
                    edge_connection_const.append(x_edges[out_edge_idx] <= no_dis)  # can be pruned ...
                    edge_connection_const.append(
                        x_edges[out_edge_idx] >= x_nodes[i] + trade_plan[out_edge_idx] + no_dis - 2)

                in_edge_idx = supply_net.edge_map[j, i]
                if (in_edge_idx > -1):
                    node_connection_const.append(x_nodes[i] >= x_edges[in_edge_idx])
                    in_degree += 1

            if (in_degree == 0):
                node_connection_const.append(x_nodes[i] == 1)

    # budget constraint, unique
    budget_const = []
    for i in range(supply_net.num_nodes):
        total_cost = 0
        for j in range(supply_net.num_nodes):
            in_edge_idx = supply_net.edge_map[j, i]

            if (in_edge_idx > -1):
                total_cost += supply_net.cost[j, i] * trade_plan[in_edge_idx]

        budget_const.append((total_cost <= supply_net.budget[i]))

    # objective, for all samples -- on average
    total_production = 0
    for edge in supply_net.edges:
        if (edge[1] in supply_net.demand_node):
            total_production += supply_net.capacity[edge[0], edge[1]] * x_edges[supply_net.edge_map[edge[0], edge[1]]]

    production_const = (total_production >= 2 ** supply_net.total_demand)
    production_const = (total_production >= 128)

    # optional!! maximize ...
    if (find_best):
        mdl.maximize(total_production)

    # add constraints
    for c in edge_connection_const:
        mdl.add(c)
    for c in node_connection_const:
        mdl.add(c)
    for c in budget_const:
        mdl.add(c)

    mdl.add(production_const)

    print("\nSolving model....")
    mdl.export_model("baseline_sampled.lp")
    msol = mdl.solve()  # (TimeLimit=10)

    # msol.stop_cause = 'SearchStoppedByLimit'
    # msol.solve_status = 'Feasible'

    if msol:
        print("Solution status: " + msol.get_solve_status())
        for edge in trade_plan:
            print("   " + edge.get_name() + f": {msol[edge.get_name()]}")
    else:
        print("No solution found")
    return [msol[edge.get_name()] for edge in trade_plan]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath",
                        default="/home/jinzhao/jinzhao/xor_smt/data/supply_chain/large_sized_network_medium_distribution",
                        help="the filename.")

    args = parser.parse_args()

    seed = int(time.perf_counter() * 10000) % 1000007
    random.seed(seed)
    print('random seed=', seed)

    seed = int(time.perf_counter() * 10000) % 1000007
    np.random.seed(seed)
    print('np.random seed=', seed)
    print(args)

    # load network
    sn = SupplyNet(args.filepath)
    n_samples = 10
    disaster_samples = sn.sample_disaster_via_loopy_gibbs_sampling(n_samples, "./sampled_output_")
    disaster_samples = process_samples(disaster_samples, sn.num_dedges)

    # generate a MIP plan
    n_samples = 10
    baseline_trade_plans = []
    for i in range(n_samples):
        baseline_trade_plans.append(mip_find_plan(sn, disaster_samples[i], False))
