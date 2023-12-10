import numpy as np
from pgmpy.models import BayesianNetwork
from pgmpy.sampling import BayesianModelSampling

import uuid
import re
import os


class UaiFile(object):
    def __init__(self, filename):
        self.filename = filename

        self.inst_type = ""
        self.n_var = 0
        self.dims = []
        self.n_cliques = 0
        self.cliques = []
        self.factors = []
        self.raw_factors = []
        self.readUai(filename)
        return

    def readFileByTokens(self, path, specials=[]):
        spliton = '([\s' + ''.join(specials) + '])'
        with open(path, 'r') as fp:
            for line in fp:
                tok = [t.strip() for t in re.split(spliton, line) if t and not t.isspace()]
                for t in tok: yield t

    def readUai(self, filename):
        dims = []  # store dimension (# of states) of the variables
        cliques = []  # cliques (scopes) of the factors we read in
        factors = []  # the factors themselves

        gen = self.readFileByTokens(filename, '(),')  # get token generator for the UAI file
        inst_type = next(gen)
        n_var = int(next(gen))  # get the number of variables
        dims = [int(next(gen)) for i in range(n_var)]  # and their dimensions (states)
        n_cliques = int(next(gen))  # get the number of cliques / factors
        cliques = [None] * n_cliques
        for c in range(n_cliques):
            c_size = int(next(gen))  # (size of clique)
            cliques[c] = [int(next(gen)) for i in range(c_size)]

        factors = []
        raw_factors = []
        for c in range(n_cliques):  # now read in the factor tables:
            t_size = int(next(gen))  # (# of entries in table = # of states in scope)
            factor_size = tuple(dims[v] for v in cliques[c]) if len(cliques[c]) else (1,)
            f_table = np.empty(t_size)
            for i in range(t_size):
                f_table[i] = float(next(gen))

            f_table_organized = np.empty(factor_size)
            for index in np.ndindex(factor_size):
                # Convert index to binary and reverse it
                b_str = ''.join(str(bit) for bit in index)
                f_table_organized[index] = f_table[int(b_str, 2)]

            factors.append(np.array(f_table_organized, dtype=float))
            raw_factors.append(f_table)

        self.inst_type = inst_type
        self.n_var = n_var
        self.dims = dims
        self.n_cliques = n_cliques
        self.cliques = cliques
        self.factors = factors
        self.raw_factors = raw_factors

        return factors, n_var

    def writePgmpyUai(self, filename):
        with open(filename, 'w') as fp:
            fp.write(self.inst_type + "\n")
            fp.write("{:d}\n".format(self.n_var))  # number of variables in model
            fp.write(" ".join(map(str, self.dims)) + "\n")  # write dimensions of each variable
            fp.write("{:d}\n".format(self.n_cliques));  # number of factors
            for clique in self.cliques:
                fp.write(str(len(clique)) + " " + " ".join(map(str, clique)))
                fp.write("\n")
            fp.write("\n")
            for factor in self.raw_factors:
                fp.write(str(factor.size) + "\n")
                line1 = factor[::2]
                line2 = factor[1::2]
                [fp.write(f"{v} ") for v in line1]
                fp.write("\n")
                [fp.write(f"{v} ") for v in line2]
                fp.write("\n\n")


def Bayesian_Sampling(filename, n_samples):
    uai = UaiFile(filename)
    random_filename = f"tmp_bayesnet_sampling_{uuid.uuid4()}.uai"
    uai.writePgmpyUai(random_filename) # convert to pgmpy type

    bn_model = BayesianNetwork.load(random_filename, filetype='uai')
    sampler = BayesianModelSampling(bn_model)
    samples = sampler.forward_sample(size=n_samples, show_progress=False)
    probs = np.sum(samples.values, axis=0) / n_samples

    print("P(X) = 1: ", probs)

    os.remove(random_filename)
    return samples



if __name__ == "__main__":
    filename = "/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network/disaster.uai"
    Bayesian_Sampling(filename, 20000)
