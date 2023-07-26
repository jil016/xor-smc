"""
   pip install z3-solver
"""
import numpy as np
from z3 import *
import scipy.sparse as scisparse
import time
def transtive_closure_BMM(A):
    """
    Let A be the adjacency matrix of some graph G. Then (A + I)^n is the transitive closure of G.
    https://people.csail.mit.edu/virgi/6.890/lecture3.pdf
    """
    temp = A
    S = A
    old_S = None
    # st = time.time()
    for i in range(A.shape[0]):
        S = np.logical_or(S, temp)
        if np.sum(S != old_S) == 0:
            break
        old_S = S
        temp = temp @ A
    # print("time used for Logical MM:", time.time() - st)
    scc = np.logical_and(S, np.transpose(S))
    return S, scc
    # seq_len = A.shape[0]
    # I = np.eye(seq_len, dtype=bool)
    # A_plus_I = np.logical_or(A, I)
    # S = A_plus_I
    # old_S = None
    # st = time.time()
    # iterations = int(np.log2(seq_len))+1
    # for i in range(iterations):
    #     S = S @ S
    #     if np.sum(S != old_S) == 0:
    #         break
    #     old_S = S
    # print("time used for TC:", time.time() - st)
    # S = S.astype(int)
    #
    # scc = np.logical_and(S, np.transpose(S))
    # return S, scc

def expression_construction(n,m):
    # all paths from one location to another location
    path = np.random.randint(low=0, high=2, size=(n, n))
    # C[t, :] all the locations that is needed by factory t.
    C = np.random.randint(low=0, high=2, size=(m, n))
    # threshold
    threshold_exp2_q = np.exp2(np.random.randint(low=0,high=5, size=m))



    c = BitVec('c', m)

    I_st = BitVecVal('I', m)


    print("vector c", c)
    expression = Sum([Implies(c[i], I_st[i]) for i in range(m)])
    print(expression)

if __name__ == '__main__':
    expression_construction(n=20, m=4)