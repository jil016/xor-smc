import os
import numpy as np
import argparse
import time
import random
from pgmpy import sampling
from pgmpy.inference import ApproxInference


class SupplyNet(object):
    def __init__(self, network_folder):
        self.folder = network_folder
        # Parameters
        self.num_nodes = 0
        self.demand_node = None
        self.total_demand = None
        self.cost = None
        self.capacity = None
        self.budget = None
        self.dedges = None
        self.disaster_uai = None
        #
        self.initialize_from_file()
        self.generate_advance_info()

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
        self.num_nodes = self.cost.shape[0]

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
        with open(os.path.join(self.folder, "budget.txt"), "r") as fp:
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

        self.dedges = [list(row) for row in zip(* self.dedges)]
        self.disaster_uai = os.path.join(self.folder, "disaster.uai")
        return

    def generate_advance_info(self):
        self.edges = []
        self.edge_map = - np.ones_like(self.capacity)

        for i, row in enumerate(self.capacity):
            for j, cap in enumerate(row):
                if(cap > 0):
                    self.edges.append([i,j])
                    self.edge_map[i,j] = len(self.edges) - 1

        self.num_edges = len(self.edges)

        self.dedge_map = - np.ones_like(self.capacity)
        for i, e in enumerate(self.dedges):
            self.dedge_map[e[0], e[1]] = i

        self.num_dedges = len(self.dedges)
        return

    def sample_disaster_via_bayesian_sampling(self):
        # inference = sampling.BayesianModelSampling(self.disasters)
        # return inference.forward_sample(size=1)
        pass

    def sample_disaster_via_gibbs_sampling(self):
        # gibbs_chain = sampling.GibbsSampling(self.disasters)
        # return gibbs_chain.sample(size=1)
        pass

    def sample_disaster_via_approx_infer(self):
        # infer = ApproxInference(self.disasters)
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
