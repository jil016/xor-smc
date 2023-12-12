import sys
import time
import math
import os
import numpy as np
from collections import Counter
from pysat.formula import CNF
from pysat.solvers import Solver

from sampler.nelson.random_sat import Monte_Carlo_sampler
from sampler.xor_sampling.xor_sampler import XOR_Sampling
from sampler.gibbs_sampler.gibbs_mrf import Gibbs_Sampling
from cnf2uai import cnf_to_uai


def draw_from_xor_sampling(num_samples, input_file, cnf_instance, prob):
    cnf_to_uai(cnf_instance, prob, input_file + ".weight.uai")
    returned_samples = XOR_Sampling(input_file + ".weight.uai", num_samples)

    return returned_samples


def draw_from_gibbs_sampling(num_samples, input_file, cnf_instance, prob):
    cnf_to_uai(cnf_instance, prob, input_file + ".weight.uai")
    returned_samples = Gibbs_Sampling(input_file + ".weight.uai", num_samples)

    return returned_samples

def draw_from_unigen(num_samples, input_file):
    tmpfile = '/tmp/unigen.txt'
    cmd = """./scripts/sampler/uniformSATSampler/unigen --input {} --samples {} --sampleout {} > /tmp/tmp.txt""".format(input_file,
                                                                                                                        num_samples,
                                                                                                                        tmpfile)

    os.system(cmd)
    sampled_assignments = []
    with open(tmpfile, 'r') as fr:
        for li in fr.read().split("\n"):
            one_sol = li.strip().split(" ")[:-1]
            # print(one_sol)
            if len(li) <= 1:
                continue
            sampled_assig = [1 if int(x) > 0 else 0 for x in one_sol]
            # print(sampled_assig)
            sampled_assignments.append(sampled_assig)
    return sampled_assignments


def draw_from_weightgen(num_samples, input_file, cnf_instance, prob, device='cuda'):
    cnf_instance.to_file(input_file + ".weight")
    with open(input_file + ".weight", "a") as fw:
        for xi in range(cnf_instance.nv):
            fw.write("w {} {} 0\n".format(xi + 1, 1.0 - prob[xi]))
            fw.write("w -{} {} 0\n".format(xi + 1, prob[xi]))
    kappa = 0.4
    timeout = 72000
    satTimeout = 3000
    epsilon = 0.8
    delta = 0.2
    tilt = 5
    pivotAC = 2 * math.ceil(4.4817 * (1 + 1 / epsilon) * (1 + 1 / epsilon))

    numIterations = int(math.ceil(35 * math.log((3 * 1.0 / delta), 2)))

    pivotUniGen = math.ceil(4.03 * (1 + 1 / kappa) * (1 + 1 / kappa))
    st = time.time()
    tmp_file = "/tmp/randksat.weightgen.log"
    cmd = """./src/lll/sampler/weightedSATSampler/weightgen --samples={} --kappa={} --pivotUniGen={} --maxTotalTime={} \
                --startIteration=0 --maxLoopTime={} --tApproxMC=17 --pivotAC=46 --gaussuntil=400 \
                --verbosity=0 --ratio={} {} {}""".format(num_samples, kappa, pivotUniGen, timeout,
                                                         satTimeout, tilt, input_file + ".weight", tmp_file)

    os.system(cmd)
    sampled_assignments = []
    with open(tmp_file, 'r') as fr:
        for li in fr:
            if len(li) <= 1:
                continue
            one_sol = li.strip().split(" ")[1:-1]
            sampled_assig = [1 if int(x) > 0 else 0 for x in one_sol]
            sampled_assignments.append(sampled_assig)
    return sampled_assignments


def draw_from_quicksampler(num_samples, input_file):
    sampled_assignments = []
    formula = CNF(input_file)
    solver = Solver(bootstrap_with=formula.clauses)
    iter = 0
    while len(sampled_assignments) < num_samples:
        iter += 1
        cmd = """./src/lll/sampler/uniformSATSampler/quicksampler -n {} -t 180.0 {} >/tmp/tmp.log""".format(num_samples, input_file)

        os.system(cmd)
        print(iter, len(sampled_assignments), end="\r")
        sys.stdout.flush()
        lines = []
        with open(input_file + '.samples', 'r') as fr:
            for li in fr.read().split("\n"):
                one_sol = li.strip().split(" ")[-1].strip()
                if len(one_sol) <= 1:
                    continue
                lines.append(one_sol)
        os.remove(input_file + '.samples')
        cnt = 0
        for k in lines:
            assignment = []
            for idx, x in enumerate(k):
                if x == '0':
                    assignment.append(-idx - 1)
                else:
                    assignment.append(idx + 1)
            if solver.solve(assignment):
                sampled_assignments.append([1 if int(x) > 0 else 0 for x in k])
            else:
                cnt += 1

    return sampled_assignments


def draw_from_kus(num_samples, input_file):
    tmpfile = "/tmp/randksat.kus.txt"
    cmd = "python3  ./src/lll/sampler/uniformSATSampler/KUS.py --samples {} --outputfile {} {} >/tmp/tmp.log".format(num_samples,
                                                                                                                     tmpfile,
                                                                                                                     input_file)

    os.system(cmd)
    sampled_assignments = []
    with open(tmpfile, 'r') as fr:
        for li in fr.read().split("\n"):
            one_sol = li.strip().split(" ")

            if len(li) <= 1:
                continue
            sampled_assig = [0, ] * len(one_sol)
            for x in one_sol:
                if int(x) > 0:
                    sampled_assig[abs(int(x)) - 1] = 1
            sampled_assignments.append(sampled_assig)
    return sampled_assignments


def draw_from_cmsgen(num_samples, input_file):
    tmpfile = "/tmp/randksat.cmsgen.txt"
    cmd = "./src/lll/sampler/uniformSATSampler/cmsgen --samples {} --samplefile {} {} >/tmp/tmp.log".format(num_samples,
                                                                                                            tmpfile,
                                                                                                            input_file)
    os.system(cmd)
    sampled_assignments = []
    with open(tmpfile, 'r') as fr:
        for li in fr.read().split("\n"):
            one_sol = li.strip().split(" ")[:-1]
            if len(li) <= 1:
                continue
            # print(li)
            sampled_assig = [1 if int(x) > 0 else 0 for x in one_sol]
            sampled_assignments.append(sampled_assig)
    # print(len(sampled_assignments), sampled_assignments[0].shape)
    return sampled_assignments

