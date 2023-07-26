import time

def sympy_formula():
    from sympy import symbols, Symbol

    import numpy as np

    np.random.seed(10086)
    import sympy.printing as printing

    from sympy.logic.boolalg import And, Or, Not, Xor, to_cnf
    from sympy.logic.inference import satisfiable

    m = 256
    K = m//2
    y = [Symbol(f'y_{i}') for i in range(m)]
    xor_assignments = np.random.randint(2, size=(K, m))
    print(xor_assignments)
    list_of_xor_terms = []
    for k in range(K):
        temp = []
        for i in range(m):
            if xor_assignments[0, i] == 1:
                temp.append(y[i])
        list_of_xor_terms.append(Xor(*temp))

    st = time.time()
    out = satisfiable(to_cnf(list_of_xor_terms[0]))
    used = time.time() - st
    print("one XOR terms", out, 'Time usage',  used)

    st = time.time()
    out = satisfiable(And(to_cnf(list_of_xor_terms[0]), to_cnf(list_of_xor_terms[1])))
    used = time.time() - st
    print("two XOR terms", out, 'Time usage',  used)

    st = time.time()
    out = satisfiable(And(to_cnf(list_of_xor_terms[0]), to_cnf(list_of_xor_terms[1]), to_cnf(list_of_xor_terms[2])))
    used = time.time() - st
    print("three XOR terms", out, 'Time usage',  used)

    st = time.time()
    out = satisfiable(
        And(to_cnf(list_of_xor_terms[0]), to_cnf(list_of_xor_terms[1]), to_cnf(list_of_xor_terms[2]), to_cnf(list_of_xor_terms[3])))
    used = time.time() - st
    print("four XOR terms", out, 'Time usage',  used)

    st = time.time()
    out = satisfiable(
        And(to_cnf(list_of_xor_terms[0]), to_cnf(list_of_xor_terms[1]), to_cnf(list_of_xor_terms[2]), to_cnf(list_of_xor_terms[3]),
            to_cnf(list_of_xor_terms[4]), to_cnf(list_of_xor_terms[5]), to_cnf(list_of_xor_terms[6]), to_cnf(list_of_xor_terms[7]),
            to_cnf(list_of_xor_terms[8]), to_cnf(list_of_xor_terms[9])
            ))
    used = time.time() - st
    print("10 XOR terms", out,"time used", used)
    phi = Symbol('phi')
    Pkm1, Pkm2, Pkm3, Pkm4 = symbols("P_{k-1} P_{k-2} P_{k-3} P_{k-4}")
    formula = phi
    component = Or(Not(Pkm1), list_of_xor_terms[0])
    formula = And(formula, component)
    print(printing.latex(to_cnf(formula)))
    print(satisfiable(to_cnf(formula)))

    component = Or(And(Not(Pkm1), Not(Pkm2)), list_of_xor_terms[1])
    formula = And(formula, component)
    print(to_cnf(formula))

    component = Or(And(Not(Pkm1), Not(Pkm2), Not(Pkm3)), list_of_xor_terms[2])
    formula = And(formula, component)
    print(printing.latex(to_cnf(formula)))
    print(to_cnf(formula))
    component = Or(And(Not(Pkm1), Not(Pkm2), Not(Pkm3), Not(Pkm4)), list_of_xor_terms[3])
    formula = And(formula, component)
    print(printing.latex(to_cnf(formula)))
    print(to_cnf(formula))


def cryptominisat_formula():
    """
    git clone https://github.com/msoos/cryptominisat
    python setup.py build
    python setup.py install
    """
    from pycryptosat import Solver
    s = Solver()
    s.add_clause([1])
    s.add_clause([-2])
    s.add_clause([-1, 2, 3])
    s.add_xor_clause([1,2,3])
    sat, solution = s.solve()
    print(sat, solution)

def cdlc_crypto_formula():
    # https://sites.google.com/view/crypto-sat/
    pass

if __name__ == "__main__":
    sympy_formula()
