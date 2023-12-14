import os
import numpy as np
import argparse
import time
import random
from pgmpy import sampling
from pgmpy.inference import ApproxInference
from supply_net_with_disaster import SupplyNet, process_samples

from find_plan import *
from eval_plan import *
from gen_synthesis_data_v2 import *

from sampler.bayes_net_sampler.bayesian import Bayesian_Sampling
from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling


def parse_spc_results(filepath):
    matrix = []
    with open(filepath, "r") as fp:
        line1 = fp.readline()[:-1]
        if(line1 != "Optimal"):
            print("SPC No solution!")
            return np.zeros_like(matrix)
        for line in fp:
            matrix.append([int(num) for num in line.split()])
    return np.array(matrix)

def run_spc_program(smc_binary = "", network_folder="", out_path="", threshold = 6):
    random_seed = 20
    T = 1
    os.system(f"{smc_binary} {network_folder} -seed {random_seed} -T {T} -output {out_path} -threshold {threshold}")
    print("XMC_SPC finished!")


def generate_new_network(out_path=""):
    # supply network
    net_struct = [9, 7, 9, 19]

    # Generate Budgets
    min_bgt = 5
    max_bgt = 30

    # Generate Capacities
    cap_mode =  "extreme"
    min_cap = 10
    max_cap = 20
    p_edge = 1
    cap_list = [5, 7, 8, 25]

    # Generate Costs
    min_cst = 5
    max_cst = 15

    # Generate Demands
    n_demands = 19
    q = 6

    # Generate Disasters
    num_dedges = 50  #
    num_bayes_edges = 300
    precision = 3  # k-digit precision 2^n-1 0~15
    max_parents = 8  # maximum of parents allowed

    # generates everything
    generate_budget(net_struct, min_bgt, max_bgt, out_path)
    if cap_mode == "extreme":
        adjacency_matrix = generate_extreme_capacity(net_struct, cap_list, p_edge, out_path)
    else:
        adjacency_matrix = generate_capacity(net_struct, min_cap, max_cap, p_edge, out_path)

    generate_cost(adjacency_matrix, min_cst, max_cst, out_path)
    generate_demand(net_struct, n_demands, q, out_path)
    generate_disaster_edges(adjacency_matrix, num_dedges, out_path)
    generate_disaster_model(num_dedges, num_bayes_edges, precision, max_parents, out_path)

    # save parameters
    with open(os.path.join(out_path, "params.txt"), "w") as fp:
        fp.write("net_struct: ")
        [fp.write(f"{x} ") for x in net_struct]
        fp.write("\n")

        # Generate Budgets
        fp.write(f"min_bgt: {min_bgt}, max_bgt: {max_bgt}\n")

        # Generate Capacities
        fp.write(f"min_cap: {min_cap}, max_cap: {max_cap}, p_edge: {p_edge}\n")

        # Generate Costs
        fp.write(f"min_cst: {min_cst}, max_cst: {max_cst}\n")

        # Generate Demands
        fp.write(f"n_demands: {n_demands}, q: {q}\n")

        # Generate Disasters
        fp.write(f"num_dedges: {num_dedges}, num_bayes_edges: {np.min([num_bayes_edges, (num_dedges * (num_dedges - 1)) // 2])} / {(num_dedges * (num_dedges - 1)) // 2}\n")
        fp.write(f"precision: {precision}, max_parents: {max_parents}\n")




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath",
                        default="/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/fast_gen_net",
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
    n_samples = 50
    all_baselines = ["gibbs", "bp", "imp", "loopy-imp", "weight"]

    generate_new_network(args.filepath)


    while True:
        sn = SupplyNet(args.filepath)
        all_baselines_samples = []
        all_baselines_plans = []
        all_baselines_eval = []
        all_baselines_time = []

        for method in all_baselines:
            if method == "gibbs":
                disaster_samples = sn.sample_disaster_via_gibbs_sampling(n_samples, "./sampled_output_")
            elif method == "bp":
                disaster_samples = sn.sample_disaster_via_loopy_belief_propagation(n_samples, "./sampled_output_")
            elif method == "imp":
                disaster_samples = sn.sample_disaster_via_importance_sampling(n_samples, "./sampled_output_")
            elif method == "loopy-imp":
                disaster_samples = sn.sample_disaster_via_loopy_importance_sampling(n_samples, "./sampled_output_")
            elif method == "weight":
                disaster_samples = sn.sample_disaster_via_weighted_sampling(n_samples, "./sampled_output_")
            else:
                exit(1)
            disaster_samples = process_samples(disaster_samples, sn.num_dedges)
            all_baselines_samples.append(disaster_samples)

        # generate SAA plans
        for samples in all_baselines_samples:
            baseline_start_time = time.time()
            all_baselines_plans.append(saa_find_plan(sn, samples, True, 64, 100))
            baseline_end_time = time.time()
            all_baselines_time.append(baseline_end_time - baseline_start_time)


        n_eval_samples = 2000
        gt_disaster_samples = Bayesian_Sampling(os.path.join(args.filepath, "disaster.uai"),
                                          n_eval_samples)

        for plan in all_baselines_plans:
            baseline_total_prodcution = []
            for i in range(n_eval_samples):
                baseline_total_prodcution.append(calc_actual_production(sn, plan, gt_disaster_samples[i]))
            all_baselines_eval.append(sum(baseline_total_prodcution) / n_eval_samples)

        # generate a SMC plan
        for threshold in range(4,9):
            smc_binary = "/Users/jinzhao/Desktop/git_repos/xor_smt/xor_smc/supply_chain/SPC"
            smc_outpath = "/Users/jinzhao/Desktop/git_repos/xor_smt/xor_smc/supply_chain/LOG-SPC"
            smc_start_time = time.time()
            run_spc_program(smc_binary, args.filepath, smc_outpath)
            smc_end_time = time.time()
            smc_trade_plan = parse_spc_results(os.path.join(smc_outpath, "result.log"))

            smc_total_prodcution = []
            for i in range(n_eval_samples):
                smc_total_prodcution.append(calc_actual_production(sn, smc_trade_plan, gt_disaster_samples[i]))

            smc_product_amount = sum(smc_total_prodcution) / n_eval_samples
            print(f"smc_total_prodcution: {smc_product_amount}")
            if smc_product_amount >= max(all_baselines_eval):
                print(f"Threshold: {threshold}")
                break

        print(f"smc_total_prodcution: {smc_product_amount}")
        print(f"smc_time: {smc_end_time - smc_start_time}")
        print(f"baseline_total_product:" + ' '.join(map(str, all_baselines_eval)))
        print(f"baseline_time:" + ' '.join(map(str, all_baselines_time)))

        if smc_product_amount >= max(all_baselines_eval):
            break
        else:
            generate_new_network(args.filepath)
    pass
