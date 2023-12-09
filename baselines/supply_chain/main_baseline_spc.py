import os
import numpy as np
import argparse
import time
import random
from pgmpy import sampling
from pgmpy.inference import ApproxInference
from supply_net import SupplyNet

from mip import Model, xsum, maximize, BINARY

def mip_knapsack_example():
    p = [10, 13, 18, 31, 7, 15]
    w = [11, 15, 20, 35, 10, 33]
    c, I = 47, range(len(w))

    m = Model("knapsack")

    x = [m.add_var(var_type=BINARY) for i in I]

    m.objective = maximize(xsum(p[i] * x[i] for i in I))

    m += xsum(w[i] * x[i] for i in I) <= c

    m.optimize()

    selected = [i for i in I if x[i].x >= 0.99]
    print("selected items: {}".format(selected))

def find_best_plan(supply_net, disaster_sample):
    supply_net = SupplyNet("/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network")

    mip_model = Model()

    # trade plan variables
    trade_plan = []
    for e in supply_net.edges:
        trade_plan.append(mip_model.add_var(name=f"s_{e[0]}_{e[1]}", var_type=BINARY))

    # node connection
    x_nodes = [mip_model.add_var(name=f"n_{i}", var_type=BINARY) for i in range(supply_net.num_nodes)]
    x_edges = [mip_model.add_var(name=f"e_{eg[0]}_{eg[1]}", var_type=BINARY) for eg in supply_net.edges]

    # for i in range(supply_net.num_nodes):
    #     for j in range(supply_net.num_nodes):

    a = 1


    # production >= threshold



    return


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

    # generate a MIP instance
    find_best_plan(sn, None)



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