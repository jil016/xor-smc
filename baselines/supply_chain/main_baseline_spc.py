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

        with open(self.folder + "_demand.txt", "r") as fp:
            # n_lines = int([:-1])
            fp.readline()
            line = fp.readline().strip().split(" ")
            self.demand_node = [int(x) for x in line]
            self.total_demand = int(fp.readline().strip())

        # cost
        self.cost = []
        with open(self.folder + "_cost.txt", "r") as fp:
            for line in fp:
                sp_line = line.strip().split(" ")
                line_costs = [int(s) for s in sp_line]
                self.cost.append(line_costs)
        self.cost = np.asarray(self.cost).astype(int)

        # capacity
        self.capacity = []
        with open(self.folder + "_capacity.txt", "r") as fp:
            fp.readline()
            for line in fp:
                sp_line = line.strip().split(" ")
                line_caps = [int(s) for s in sp_line]
                self.capacity.append(line_caps)
        self.capacity = np.asarray(self.capacity).astype(int)

        # budget
        with open(self.folder + "_budget.txt", "r") as fp:
            line = fp.readline()
            self.budget = [int(c) for c in line.strip().split(" ")]

        # disasters
        self.disasters = []
        model=BayesianNetwork.load("/home/jiangnan/PycharmProjects/xor_smt/data/supply_chain/Layers_2/Nodes_3_4/prog_0_disaster.bif", filetype='bif')
        print(model)
        self.disasters = BayesianNetwork.load(self.folder + "_disaster.bif", filetype='bif')

    def sample_disaster_via_bayesian_sampling(self):
        inference = sampling.BayesianModelSampling(self.disasters)
        return inference.forward_sample(size=2)

    def sample_disaster_via_gibbs_sampling(self):

        gibbs_chain = sampling.GibbsSampling(self.disasters)
        return gibbs_chain.sample(size=3)

    def sample_disaster_via_approx_infer(self):
        infer = ApproxInference(self.disasters)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath", default="/home/jiangnan/PycharmProjects/xor_smt/data/supply_chain/Layers_2/Nodes_3_4",
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

    # sample a disaster
    print(sn.sample_disaster_via_gibbs_sampling())
    print(sn.sample_disaster_via_bayesian_sampling())
