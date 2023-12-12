import os
import numpy as np
import argparse
import time
import random
from pgmpy import sampling
from pgmpy.inference import ApproxInference
from supply_net_with_disaster import SupplyNet, process_samples


from docplex.cp.model import CpoModel, CpoIntVar
from docplex.cp.config import context

context.solver.agent = 'local'
context.solver.local.execfile = '/home/jinzhao/.local/ibm/ILOG/CPLEX_Studio2211/cpoptimizer/bin/x86-64_linux/cpoptimizer'

from sampler.bayes_net_sampler.bayesian import Bayesian_Sampling
from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling

def mip_find_plan(supply_net, disaster_sample, find_best = False):
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
            out_edge_idx = supply_net.edge_map[i,j]
            dedge_idx = supply_net.dedge_map[i,j]

            no_dis = 1
            if(dedge_idx != -1 and disaster_sample[dedge_idx] == 1):
                no_dis = 0

            if(out_edge_idx > -1):
                edge_connection_const.append(x_edges[out_edge_idx] <= x_nodes[i])
                edge_connection_const.append(x_edges[out_edge_idx] <= trade_plan[out_edge_idx])
                edge_connection_const.append(x_edges[out_edge_idx] <= no_dis)   # can be pruned ...
                edge_connection_const.append(x_edges[out_edge_idx] >= x_nodes[i] + trade_plan[out_edge_idx] + no_dis - 2)

            in_edge_idx = supply_net.edge_map[j, i]
            if(in_edge_idx > -1):
                node_connection_const.append(x_nodes[i] >= x_edges[in_edge_idx])
                in_degree += 1

        if(in_degree == 0):
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
        if(edge[1] in supply_net.demand_node):
            total_production += supply_net.capacity[edge[0], edge[1]] * x_edges[supply_net.edge_map[edge[0],edge[1]]]

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

    if msol:
        print("Solution status: " + msol.get_solve_status())
        for edge in trade_plan:
            print("   " + edge.get_name() + f": {msol[edge.get_name()]}")
    else:
        print("No solution found")
    return [msol[edge.get_name()] for edge in trade_plan]


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


def parse_SPC_results(filepath):
    matrix = []
    with open(filepath, "r") as fp:
        line1 = fp.readline()[:-1]
        if(line1 != "Optimal"):
            print("SPC No solution!")
            exit(0)
        for line in fp:
            matrix.append([int(num) for num in line.split()])
    return np.array(matrix)

def run_SPC_program(smc_binary = "", network_folder="", out_path=""):
    random_seed = 20
    T = 1
    os.system(f"{smc_binary} {network_folder} -seed {random_seed} -T {T} -output {out_path}")
    print("XMC_SPC finished!")


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

    disaster_samples = Gibbs_Sampling(os.path.join(args.filepath, "disaster.uai"), n_samples, None, 50)
    
    # disaster_samples = sn.disaster_exact_LazyPropagation(n_samples, './sampled_output_')
    # disaster_samples = sn.sample_disaster_via_importance_sampling(n_samples, "./sampled_output_")
    # disaster_samples = sn.sample_disaster_via_weighted_sampling(n_samples, "./sampled_output_")
    # disaster_samples = sn.sample_disaster_via_loopy_weighted_sampling(n_samples, "./sampled_output_")
    # disaster_samples = sn.sample_disaster_via_loopy_gibbs_sampling(n_samples, "./sampled_output_")
    # disaster_samples = sn.sample_disaster_via_loopy_belief_propagation(n_samples, "./sampled_output_")
    # disaster_samples = sn.sample_disaster_via_loopy_importance_sampling(n_samples, "./sampled_output_")

    # disaster_samples = sn.sample_disaster_via_loopy_belief_propagation(n_samples, "./sampled_output_")
    disaster_samples = process_samples(disaster_samples, sn.num_dedges)


    # generate a MIP plan
    n_samples = 1
    baseline_trade_plans = []
    for i in range(n_samples):
        baseline_trade_plans.append(mip_find_plan(sn, disaster_samples[i], False))


    # generate a SMC plan
    smc_binary = "/home/jinzhao/jinzhao/xor_smt/xor_smc/supply_chain/SPC"
    smc_outpath = "/home/jinzhao/jinzhao/xor_smt/xor_smc/supply_chain/LOG-SPC"
    run_SPC_program(smc_binary, args.filepath, smc_outpath)
    smc_trade_plan = parse_SPC_results(os.path.join(smc_outpath, "result.log"))


    n_eval_samples = 100
    gt_disaster_samples = Bayesian_Sampling(os.path.join(args.filepath, "disaster.uai"),
                                          n_eval_samples)
    
    smc_total_prodcution = []
    for i in range(n_eval_samples):
        smc_total_prodcution.append(calc_actual_production(sn, smc_trade_plan, gt_disaster_samples[i]))
    print(f"smc_total_prodcution: {sum(smc_total_prodcution) / n_eval_samples}")


    baseline_total_prodcution_list = []
    for i_plan in range(n_samples):
        baseline_total_prodcution = []
        for i in range(n_eval_samples):
            baseline_total_prodcution.append(calc_actual_production(sn, baseline_trade_plans[i_plan], gt_disaster_samples[i]))

        baseline_total_prodcution_list.append(sum(baseline_total_prodcution) / n_eval_samples)
        print(f"baseline_total_prodcution {i_plan}: {sum(baseline_total_prodcution) / n_eval_samples}")
    

    print(f"baseline_entire_average: {sum(baseline_total_prodcution_list) / n_samples}")

    pass
