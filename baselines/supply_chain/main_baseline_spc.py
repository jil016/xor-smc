import os
import numpy as np
import argparse
from pgmpy.models import BayesianNetwork
from pgmpy.readwrite import UAIReader
import time
import random
from pgmpy.factors.discrete.CPD import TabularCPD
from pgmpy.inference import BeliefPropagation
from pgmpy import sampling
from pgmpy.inference import ApproxInference


class SupplyNet(object):
    def __init__(self, network_folder):
        self.folder = network_folder
        self.initialize_from_file()

    def initialize_from_file(self):
        with open(os.path.join(self.folder, "demand.txt"), "r") as fp:
            # n_lines = int([:-1])
            fp.readline()
            line = fp.readline().strip().split(" ")
            self.demand_node = [int(x) for x in line]
            self.total_demand = int(fp.readline().strip())

        # cost
        self.cost = []
        with open(os.path.join(self.folder, "cost.txt"), "r") as fp:
            for line in fp:
                sp_line = line.strip().split(" ")
                line_costs = [int(s) for s in sp_line]
                self.cost.append(line_costs)
        self.cost = np.asarray(self.cost).astype(int)

        # capacity
        self.capacity = []
        with open(os.path.join(self.folder, "capacity.txt"), "r") as fp:
            fp.readline()
            for line in fp:
                sp_line = line.strip().split(" ")
                line_caps = [int(s) for s in sp_line]
                self.capacity.append(line_caps)
        self.capacity = np.asarray(self.capacity).astype(int)

        # budget
        with open(os.path.join(self.folder, "capacity.txt"), "r") as fp:
            line = fp.readline()
            self.budget = [int(c) for c in line.strip().split(" ")]

        # disaster edges
        self.dedges = []
        with open(os.path.join(self.folder, "disaster.uai.edges"), "r") as fp:
            fp.readline()
            for i in range(2):
                line = fp.readline()
                sp_line = line.strip().split(" ")
                v_idx = [int(s) for s in sp_line]
                self.dedges.append(v_idx)


    def sample_disaster_via_bayesian_sampling(self):
        inference = sampling.BayesianModelSampling(self.disasters)
        return inference.forward_sample(size=1)

    def sample_disaster_via_gibbs_sampling(self):

        gibbs_chain = sampling.GibbsSampling(self.disasters)
        return gibbs_chain.sample(size=1)

    def sample_disaster_via_approx_infer(self):
        infer = ApproxInference(self.disasters)


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


def test_gibbs_sampler():
    from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling
    n_samples = 200
    samples = Gibbs_Sampling("/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network/disaster.uai",
                   n_samples,
                   None,
                   50)

    probs = samples.sum(axis=0)
    probs = probs / n_samples

    return probs





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

    # uai_path = os.path.join(args.filepath, "disaster.uai")
    probs = test_gibbs_sampler()

    print(probs)



    # # load network
    # sn = SupplyNet(args.filepath)
    #
    # # sample a disaster
    # print(sn.sample_disaster_via_gibbs_sampling())
    # print(sn.sample_disaster_via_bayesian_sampling())
