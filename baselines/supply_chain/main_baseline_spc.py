import os
import numpy as np
import argparse
import time
import random
from pgmpy import sampling
from pgmpy.inference import ApproxInference
from supply_net import SupplyNet


from docplex.cp.model import CpoModel, CpoIntVar
from docplex.cp.config import context

context.solver.agent = 'local'
context.solver.local.execfile = '/Applications/CPLEX_Studio2211/cpoptimizer/bin/arm64_osx/cpoptimizer'



# def find_best_plan(supply_net, disaster_sample):
#     from mip import Model, xsum, maximize, BINARY
#     from docplex.cp.model import CpoModel
#
#     supply_net = SupplyNet("/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network")
#     mip_model = Model()
#
#     # trade plan variables
#     trade_plan = []
#     for e in supply_net.edges:
#         trade_plan.append(mip_model.add_var(name=f"s_{e[0]}_{e[1]}", var_type=BINARY))
#
#     # node connection
#     x_nodes = [mip_model.add_var(name=f"n_{i}", var_type=BINARY) for i in range(supply_net.num_nodes)]
#     x_edges = [mip_model.add_var(name=f"e_{eg[0]}_{eg[1]}", var_type=BINARY) for eg in supply_net.edges]
#
#     # production >= threshold
#     return

def cplex_color_example():
    # Create CPO model
    mdl = CpoModel()

    # Create model variables containing colors of the countries
    Belgium = mdl.integer_var(0, 3, "Belgium")
    Denmark = mdl.integer_var(0, 3, "Denmark")
    France = mdl.integer_var(0, 3, "France")
    Germany = mdl.integer_var(0, 3, "Germany")
    Luxembourg = mdl.integer_var(0, 3, "Luxembourg")
    Netherlands = mdl.integer_var(0, 3, "Netherlands")
    ALL_COUNTRIES = (Belgium, Denmark, France, Germany, Luxembourg, Netherlands)

    # Create constraints
    mdl.add(Belgium != France)
    mdl.add(Belgium != Germany)
    mdl.add(Belgium != Netherlands)
    mdl.add(Belgium != Luxembourg)
    mdl.add(Denmark != Germany)
    mdl.add(France != Germany)
    mdl.add(France != Luxembourg)
    mdl.add(Germany != Luxembourg)
    mdl.add(Germany != Netherlands)

    # Solve model
    print("\nSolving model....")
    msol = mdl.solve()

    if msol:
        print("Solution status: " + msol.get_solve_status())
        colors = ("Yellow", "Red", "Green", "Blue")
        for country in ALL_COUNTRIES:
            print("   " + country.get_name() + ": " + colors[msol[country]])
    else:
        print("No solution found")


def find_best_plan(supply_net, disaster_sample):
    supply_net = SupplyNet("/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network")
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
            total_production += supply_net.capacity[edge[0], edge[1]] * trade_plan[supply_net.edge_map[edge[0],edge[1]]]

    production_const = (total_production >= 2 ** supply_net.total_demand)
    # add constraints


    # load all constraints
    return



def calc_actual_production(supply_net, trade_plan, disaster_sample):

    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath",
                        default="/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network",
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

    disaster_sample = np.random.randint(0,2,size=sn.num_dedges)

    # generate a MIP instance
    find_best_plan(sn, disaster_sample)



# def evaluateXORres():
#     N = 100
#     file_path = "./LOG-SPC-test/result.log"
#     network_folder = "./test_net"
#
#     # read XOR result
#     res = []
#     with open(file_path, "r") as fp:
#         fp.readline()
#         for i in range(N):
#             line = fp.readline()
#             line = line[:-2].split(' ')
#             row = [int(l) for l in line]
#             res.append(row)
#
#     # read GT
#     supply_net = SupplyNet(network_folder, 10)
#     res_dict, cost_dict = supply_net.ground_truth_production()
#
#     # calc XOR quality
#     xor_res = []
#     for e in supply_net.end_users:
#         cap_acc = 0
#         for i in range(N):
#             if res[i][e] == 1:
#                 cap_acc += res_dict[e][i]
#         xor_res.append(cap_acc)
#
#     print(supply_net.end_users)
#     print(xor_res)
#
#     gt_res = []
#     for e in supply_net.end_users:
#         # temp = mixedIntegerProgramming(res_dict[e], cost_dict[e])
#         # gt_res.append(-temp)
#         best = 0
#         # simply pick best
#         for i, c in enumerate(cost_dict[e]):
#             if c <= 1000:
#                 if res_dict[e][i] > best:
#                     best = res_dict[e][i]
#
#         gt_res.append(best)
#
#     print(gt_res)
#     return