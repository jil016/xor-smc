import numpy as np
from pgmpy.models import BayesianNetwork
from pgmpy.sampling import BayesianModelSampling

import re


class UaiFile(object):
    def __init__(self, filename):
        self.filename = filename

        self.inst_type = ""
        self.n_var = 0
        self.dims = []
        self.n_cliques = 0
        self.cliques = []
        self.factors = []
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

        factors = [None] * n_cliques
        for c in range(n_cliques):  # now read in the factor tables:
            t_size = int(next(gen))  # (# of entries in table = # of states in scope)
            factor_size = tuple(dims[v] for v in cliques[c]) if len(cliques[c]) else (1,)
            f_table = np.empty(t_size)
            for i in range(t_size):
                f_table[i] = float(next(gen))

            f_table_T = np.empty(factor_size)
            for index in np.ndindex(factor_size):
                # Convert index to binary and reverse it
                b_str = ''.join(str(bit) for bit in index)
                f_table_T[index] = f_table[int(b_str, 2)]

            factors[c] = np.array(f_table_T, dtype=float)

        self.inst_type = inst_type
        self.n_var = n_var
        self.dims = dims
        self.n_cliques = n_cliques
        self.cliques = cliques
        self.factors = factors

        return factors, n_var

    def writeUai(self, filename):
        pass


def Bayesian_Sampling(filename, n_samples):
    bn_model = BayesianNetwork.load(filename, filetype='uai')
    sampler = BayesianModelSampling(bn_model)
    samples = sampler.forward_sample(size=n_samples)

    probs = np.sum(samples.values, axis=0) / n_samples

    print(samples.iloc[0])
    print("P(X) = 1: ", probs)


def my_Bayesian_Sampling(filename, n_samples):
    uai = UaiFile(filename)
    is_fixed = [0] * uai.n_var

    pass


if __name__ == "__main__":
    uai = UaiFile("/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network/disaster.uai")

    my_Bayesian_Sampling("disaster.uai", 20000)
