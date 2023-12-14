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

def run_spc_program(smc_binary = "", network_folder="", out_path=""):
    random_seed = 20
    T = 1
    threshold = 5
    os.system(f"{smc_binary} {network_folder} -seed {random_seed} -T {T} -output {out_path} -threshold {threshold}")
    print("XMC_SPC finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath",
                        default="/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/fast_gen_net_save3",
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
        all_baselines_plans.append(saa_find_plan(sn, samples, True, 250, 100))
        baseline_end_time = time.time()
        all_baselines_time.append(baseline_end_time - baseline_start_time)

    # generate a SMC plan
    smc_binary = "/Users/jinzhao/Desktop/git_repos/xor_smt/xor_smc/supply_chain/SPC"
    smc_outpath = "/Users/jinzhao/Desktop/git_repos/xor_smt/xor_smc/supply_chain/LOG-SPC"
    smc_start_time = time.time()
    run_spc_program(smc_binary, args.filepath, smc_outpath)
    smc_end_time = time.time()
    smc_trade_plan = parse_spc_results(os.path.join(smc_outpath, "result.log"))

    n_eval_samples = 10000
    gt_disaster_samples = Bayesian_Sampling(os.path.join(args.filepath, "disaster.uai"),
                                            n_eval_samples)

    smc_total_prodcution = []
    for i in range(n_eval_samples):
        smc_total_prodcution.append(calc_actual_production(sn, smc_trade_plan, gt_disaster_samples[i]))

    smc_product_amount = sum(smc_total_prodcution) / n_eval_samples
    print(f"smc_total_prodcution: {smc_product_amount}")

    for plan in all_baselines_plans:
        baseline_total_prodcution = []
        for i in range(n_eval_samples):
            baseline_total_prodcution.append(calc_actual_production(sn, plan, gt_disaster_samples[i]))
        all_baselines_eval.append(sum(baseline_total_prodcution) / n_eval_samples)

    print(f"smc_total_prodcution: {smc_product_amount}")
    print(f"smc_time: {smc_end_time - smc_start_time}")
    print(f"baseline_total_product:" + ' '.join(map(str, all_baselines_eval)))
    print(f"baseline_time:" + ' '.join(map(str, all_baselines_time)))