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
            exit(0)
        for line in fp:
            matrix.append([int(num) for num in line.split()])
    return np.array(matrix)

def run_spc_program(smc_binary = "", network_folder="", out_path=""):
    random_seed = 20
    T = 1
    os.system(f"{smc_binary} {network_folder} -seed {random_seed} -T {T} -output {out_path}")
    print("XMC_SPC finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath",
                        default="/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/large_sized_network_medium_distribution",
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

    # disaster_samples = process_samples(disaster_samples, sn.num_dedges)


    # generate a MIP plan
    n_samples = 1
    baseline_trade_plans = []
    for i in range(n_samples):
        baseline_trade_plans.append(mip_find_plan(sn, disaster_samples[i], False))


    # generate a SMC plan
    smc_binary = "/home/jinzhao/jinzhao/xor_smt/xor_smc/supply_chain/SPC"
    smc_outpath = "/home/jinzhao/jinzhao/xor_smt/xor_smc/supply_chain/LOG-SPC"
    run_spc_program(smc_binary, args.filepath, smc_outpath)
    smc_trade_plan = parse_spc_results(os.path.join(smc_outpath, "result.log"))


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
