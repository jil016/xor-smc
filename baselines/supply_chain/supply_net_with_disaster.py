import os
import numpy as np
import argparse
import time
import random

import pyAgrum as gum
from pyAgrum import BayesNet
import os

import matplotlib.pyplot as plt
import timeit


class Timer:
    def __enter__(self):
        self.start = timeit.default_timer()
        return self

    def __exit__(self, *args):
        self.end = timeit.default_timer()
        self.duration = self.end - self.start


def execute(bn, ie):
    with Timer() as t:
        ie.makeInference()
        for i in bn.nodes():
            a = ie.posterior(i)
    return "duration : {:3.3f}s".format(t.duration)


def vals(bn, ie):
    exact = []

    for node in bn.nodes():
        # potentials as numpy array
        exact += ie.posterior(node).tolist()

    return exact


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
        self.disaster_edges = None
        self.disaster_uai_filepath = None

        self.edges = None
        self.edge_map = None
        self.num_edges = 0
        self.dedge_map = None
        self.num_dedges = 0
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
        self.disaster_edges = []
        with open(os.path.join(self.folder, "disaster.uai.edges"), "r") as fp:
            fp.readline()
            for i in range(2):
                line = fp.readline()
                sp_line = line.strip().split(" ")
                v_idx = [int(s) for s in sp_line]
                self.disaster_edges.append(v_idx)

        self.disaster_edges = [list(row) for row in zip(*self.disaster_edges)]
        disaster_uai_filepath = os.path.join(self.folder, "disaster.uai")

        self.disasters = BayesNet()
        self.disasters.loadUAI(name=disaster_uai_filepath)
        print(self.disasters)
        # parameters used for baye net sampling
        self.evs = None
        self.maxtime = 10
        self.epsilon = 0.01

    def generate_advance_info(self):
        self.edges = []
        self.edge_map = - np.ones_like(self.capacity)

        for i, row in enumerate(self.capacity):
            for j, cap in enumerate(row):
                if cap > 0:
                    self.edges.append([i, j])
                    self.edge_map[i, j] = len(self.edges) - 1

        self.num_edges = len(self.edges)

        self.dedge_map = - np.ones_like(self.capacity)
        for i, e in enumerate(self.disaster_edges):
            self.dedge_map[e[0], e[1]] = i

        self.num_dedges = len(self.disaster_edges)

    def sample_disaster_via_gibbs_sampling(self, sample_size, output_filepath):
        """
        return the marginal distribution of nodes
        """
        ie2 = gum.GibbsSampling(self.disasters)
        if self.evs is not None:
            ie2.setEvidence(self.evs)
        ie2.setMaxTime(self.maxtime)
        ie2.setEpsilon(self.epsilon)
        txt = "Gibbs : " + execute(self.disasters, ie2) + "\n" + ie2.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie2)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'gibbs_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def sample_disaster_via_weighted_sampling(self, sample_size, output_filepath):
        ie4 = gum.WeightedSampling(self.disasters)
        if self.evs is not None:
            ie4.setEvidence(self.evs)
        ie4.setMaxTime(self.maxtime)
        ie4.setEpsilon(self.epsilon)
        txt = "Weighted : " + execute(self.disasters, ie4) + "\n" + ie4.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie4)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'weighted_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def sample_disaster_via_importance_sampling(self, sample_size, output_filepath):
        ie5 = gum.ImportanceSampling(self.disasters)
        if self.evs is not None:
            ie5.setEvidence(self.evs)
        ie5.setMaxTime(self.maxtime)
        ie5.setEpsilon(self.epsilon)
        txt = "Importance: " + execute(self.disasters, ie5) + "\n" + ie5.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie5)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'importance_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled


    def sample_disaster_via_loopy_belief_propagation(self, sample_size, output_filepath):
        ie6 = gum.LoopyBeliefPropagation(self.disasters)
        if self.evs is not None:
            ie6.setEvidence(self.evs)
        ie6.setMaxTime(self.maxtime)
        ie6.setEpsilon(self.epsilon)
        txt = "LBP: " + execute(self.disasters, ie6) + "\n" + ie6.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie6)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'loopy_belief_propagation.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def sample_disaster_via_loopy_weighted_sampling(self, sample_size, output_filepath):
        ie7 = gum.LoopyWeightedSampling(self.disasters)
        if self.evs is not None:
            ie7.setEvidence(self.evs)
        ie7.setMaxTime(self.maxtime)
        ie7.setEpsilon(self.epsilon)
        txt = "LoopyWeighted: " + execute(self.disasters, ie7) + "\n" + ie7.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie7)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'loopy_weighted_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def sample_disaster_via_loopy_gibbs_sampling(self, sample_size, output_filepath):
        ie8 = gum.LoopyGibbsSampling(self.disasters)
        if self.evs is not None:
            ie8.setEvidence(self.evs)
        ie8.setMaxTime(self.maxtime)
        ie8.setEpsilon(self.epsilon)
        txt = "LoopyGibbs: " + execute(self.disasters, ie8) + "\n" + ie8.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie8)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'loopy_gibbs_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def sample_disaster_via_loopy_importance_sampling(self, sample_size, output_filepath):
        ie9 = gum.LoopyImportanceSampling(self.disasters)
        if self.evs is not None:
            ie9.setEvidence(self.evs)
        ie9.setMaxTime(self.maxtime)
        ie9.setEpsilon(self.epsilon)
        txt = "LoopyImportance: " + execute(self.disasters, ie9) + "\n" + ie9.messageApproximationScheme()
        marginal_prob = vals(self.disasters, ie9)
        print(marginal_prob)
        print(txt)
        total_sampled = []
        with open(output_filepath + 'loopy_importance_sampling.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled

    def disaster_exact_LazyPropagation(self, sample_size, output_filepath):
        ie = gum.LazyPropagation(self.disasters)
        if self.evs is not None:
            ie.setEvidence(self.evs)
        marginal_prob = vals(self.disasters, ie)
        total_sampled = []
        with open(output_filepath + 'exact_lazy_propagation.txt', 'w') as fw:
            for si in range(sample_size):
                str_sampled, sampled = sample_(marginal_prob)
                fw.write(str_sampled + "\n")
                total_sampled.append(sampled)
        return total_sampled


def sample_(marg_prob: list):
    """return the list of sampled nodes"""
    sampled_nodes = []
    for i in range(int(len(marg_prob))):
        if i % 2 == 0:
            val = np.random.choice([True, False], p=[marg_prob[i], 1 - marg_prob[i]])
            if val:
                sampled_nodes.append(str(int(i/2)))
    return " ".join(sampled_nodes), sampled_nodes


def process_samples(samples, n_dim):
    def one_hot_encode_single(indices, n_dim):
        # Initialize a list of zeros
        encoded = [0] * n_dim

        # Set the specified indices to 1
        for index in indices:
            if 0 <= int(index) < n_dim:
                encoded[int(index)] = 1
            else:
                raise ValueError(f"Index {index} out of range for dimension {n_dim}")

        return encoded

        # Encode each list of indices
    return [one_hot_encode_single(indices, n_dim) for indices in samples]


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
    """in the output file, every line contains the list of nodes with disaster happened."""
    n_samples = 1000

    gibbs_samples = sn.sample_disaster_via_gibbs_sampling(n_samples, "./sampled_output_")
    processed_gibbs_samples = process_samples(gibbs_samples, sn.num_dedges)

    marginal_prob = np.array(processed_gibbs_samples).sum(axis=0)
    marginal_prob = marginal_prob / n_samples
    print(marginal_prob)

    # exact_samples = sn.disaster_exact_LazyPropagation(10, './sampled_output_')
    # sampledoutput = sn.sample_disaster_via_importance_sampling(10, "./sampled_output_")
    # sampledoutput = sn.sample_disaster_via_weighted_sampling(10, "./sampled_output_")
    # sampledoutput = sn.sample_disaster_via_loopy_weighted_sampling(10, "./sampled_output_")
    # sampledoutput = sn.sample_disaster_via_loopy_gibbs_sampling(10, "./sampled_output_")
    # sampledoutput = sn.sample_disaster_via_loopy_belief_propagation(10, "./sampled_output_")
    # sampledoutput = sn.sample_disaster_via_loopy_importance_sampling(10, "./sampled_output_")

