import re
import numpy as np
import copy
import warnings
warnings.filterwarnings("ignore")


class UaiFile_old(object):
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
            f_table = f_table.reshape(factor_size)

            f_table_T = np.transpose(f_table, tuple(np.argsort(cliques[c])))
            # factors[c] = tuple((cliques[c], np.array(f_table_T, dtype=float)))
            factors[c] = np.array(f_table_T, dtype=float)

        self.inst_type = inst_type
        self.n_var = n_var
        self.dims = dims
        self.n_cliques = n_cliques
        self.cliques = cliques
        self.factors = factors

        return factors, n_var

    def writeUai(self, filename):
        with open(filename, 'w') as fp:
            fp.write(self.inst_type + "\n")
            fp.write("{:d}\n".format(self.n_var))  # number of variables in model
            fp.write(" ".join(map(str, self.dims)) + "\n")  # write dimensions of each variable
            fp.write("{:d}\n".format(self.n_cliques))  # number of factors
            for clique in self.cliques:
                fp.write(str(len(clique)) + " " + " ".join(map(str, clique)))
                fp.write("\n")
            fp.write("\n")
            for factor in self.factors:
                fp.write(str(factor.size) + "\n")
                fp.write(str(factor).replace(' [', '').replace('[', '').replace(']', ''))
                fp.write("\n\n")


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

def Gibbs_Sampling(filename, n_samples, initial=None, burnin_time=50):
    # Takes too much time
    uai = UaiFile(filename)
    if initial == None:
        term = np.random.binomial(1, 0.5, size=(n_samples, uai.n_var))
    else:
        term = copy.deepcopy(initial)

    # Vectorize
    for j in range(burnin_time):
        for x_idx in range(uai.n_var):
            # sample each element x
            p = np.ones([n_samples, 2], dtype=float)

            # compute marginal distribution
            for f_idx in range(uai.n_cliques):
                # each factor
                arr = np.array(uai.cliques[f_idx])
                loc = np.where(arr == x_idx)[0]
                if len(loc) == 0:
                    continue

                label = term[:, tuple(uai.cliques[f_idx])]
                label[:, loc[0]] = 0
                p[:, 0] *= uai.factors[f_idx][tuple(label.T)]
                label[:, loc[0]] = 1
                p[:, 1] *= uai.factors[f_idx][tuple(label.T)]

            p /= np.sum(p, axis=1)[:, None]

            rand_ = np.random.uniform(0.0, 1.0, size=(n_samples, 1))

            term[(rand_.squeeze() < p[:, 0].squeeze()), x_idx] = 0
            term[(rand_.squeeze() >= p[:, 0].squeeze()), x_idx] = 1


    return term


if __name__ == "__main__":
    n_samples = 2000
    filename = "/Users/jinzhao/Desktop/git_repos/xor_smt/data/supply_chain/network/disaster.uai"
    samples = Gibbs_Sampling(filename,
                             n_samples,
                             None,
                             100)

    probs = samples.sum(axis=0)
    probs = probs / n_samples
    print(probs)