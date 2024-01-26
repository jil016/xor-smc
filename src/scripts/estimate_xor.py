"""
this example code estimate how many number of XOR can be "not satisfiable"
"""
import os.path

import sympy.printing as printing

from sympy.logic.boolalg import And, Or, Not, Xor, to_cnf
from sympy.logic.inference import satisfiable
from sympy import symbols, Symbol

import numpy as np

np.random.seed(10086)


def sympy_formula(n=256):
    K = n // 2
    variables = [str(i+1) for i in range(n)]
    xor_assignments = np.random.randint(2, size=(K, n))
    list_of_xor_terms = []
    for k in range(K):
        temp = []
        for i in range(n):
            if xor_assignments[k, i] == 1:
                temp.append(variables[i])
        list_of_xor_terms.append(temp)
        print(list_of_xor_terms)
        header = """p cnf {} {}\n""".format(n, len(list_of_xor_terms))
        write_to_file(list_of_xor_terms, "./", header)


def write_to_file(xor_formula, basepath, header):
    with open(os.path.join(basepath, str(len(xor_formula))+".cnf"), 'w') as fw:
        fw.write(header)
        for line in xor_formula:
            fw.write("x" + " ".join(line) + " 0\n")
        fw.close()


if __name__ == '__main__':
    sympy_formula()
