import numpy as np
from sympy import *
from sympy.logic import And, Not, Xor, Or
import itertools

np.random.seed(10086)


def half_adder(a: bool, b: bool):
    save = Xor(a, b)
    carry = And(a, b)
    return save, carry


def full_adder(a: bool, b: bool, c: bool):
    save = Xor(a, b, c)
    carry = Or(And(c, Or(a, b)), And(a, b))
    return save, carry


def bit_add(a: list, b: list, n: int) -> list:
    a = a[::-1]
    b = b[::-1]
    c = np.zeros(n + 1, dtype=bool)
    r = np.zeros(n + 1, dtype=bool)
    r[0], c[1] = half_adder(a[0], b[0])
    for k in range(1, n):
        r[k], c[k + 1] = full_adder(a[k], b[k], c[k])

    val = 0
    for k in range(n + 1):
        val += r[k] * (2 ** k)
    r = r[::-1]
    c = c[::-1]
    print("r", r.astype(int))
    print("c", c.astype(int))
    print(val)
    return val


def test_addition():
    def addition(a: list, b: list, n: int) -> list:
        A = 0
        B = 0
        a = a[::-1]
        b = b[::-1]
        for k in range(n):
            A += a[k] * (2 ** k)
            B += b[k] * (2 ** k)
        C = A + B
        print(" a={} binary format {}\n b={} binary format:{}\n a+b={} binary format:{}".format(A, bin(A), B, bin(B), C, bin(C)))
        return C

    n = 8
    for i in range(100):
        a = np.random.randint(low=0, high=2, size=n, dtype=bool)
        b = np.random.randint(low=0, high=2, size=n, dtype=bool)
        a[0] = 0
        b[0] = 0
        print("a:{}".format(a.astype(int)))
        print("b:{}".format(b.astype(int)))
        C = addition(a, b, n)
        bit_C = bit_add(a, b, n)
        assert C == bit_C, "they are not equal!"


def save_carry_trick(a: list) -> tuple:
    """
    an array of length three. a0 + a1 + a2 = d + e.
    :param a:
    :return:
    """
    save = Xor(a[0], a[1], a[2])
    carry = Or(And(a[0], a[1]), And(a[0], a[2]), And(a[1], a[2]))
    return save, carry


def bit_mul(a: list, b: list):
    """
    https://ttu-ir.tdl.org/bitstream/handle/2346/21613/31295013250955.pdf?sequence=1
    https://en.wikipedia.org/wiki/Binary_multiplier
    :param a:
    :param b:
    :return:
    """
    a = a[::-1]
    b = b[::-1]
    pairwise_terms = []
    pairwise_terms_symbols = []
    for i, ai in enumerate(a):
        pairwise_terms.append([And(ai, bi) for bi in b])
        pairwise_terms_symbols.append([Symbol(f"p{i}{j}") for j in range(len(b))])
        print(pairwise_terms[-1])
        print(pairwise_terms_symbols[-1])
    print('-' * 30)
    output = []

    for i in range(len(a) + len(b)):
        tmp = []
        # if i < len(a):
        # tmp = [pairwise_terms_symbols[i][0], ]
        for j in range(i + 1):
            if i - j >= len(pairwise_terms_symbols) or j >= len(pairwise_terms_symbols[i - j]):
                continue
            tmp.append(pairwise_terms_symbols[i - j][j])
        output.append(tmp)
        print(i, ':', tmp)
    ### TODO: add them up
    print('-' * 20)
    carries = [[] for i in range(len(a) + len(b))]
    for i in range(1, len(a)):
        print('i=', i, output[i])
        ai = output[i].pop(0)
        bi = output[i].pop(0)
        save, carry = half_adder(ai, bi)
        carries[i + 1].insert(0, carry)
        output[i].insert(0, save)
    print("output P:")
    for i, oi in enumerate(output):
        print(i, len(oi), oi)
    print("carries:")
    for i, ci in enumerate(carries):
        print(i, ci)
    for layer in range(2, len(a)):
        print(f'-----{layer}--------')
        for i in range(len(a) - 1):
            ai = output[i + layer].pop(0)
            bi = output[i + layer].pop(0)
            ci = carries[i + layer].pop(0)
            print("i={},{},\t{},\t{}".format(i, ai, bi, ci))
            # print(i,bi)
            save, carry = full_adder(ai, bi, ci)
            carries[i + layer + 1].append(carry)
            output[i + layer].insert(0, save)
        print('_' * 20)
        print("output P:")
        for i, oi in enumerate(output):
            print(i, len(oi), oi)
        print("carries:")
        for i, ci in enumerate(carries):
            print(i, ci)
        print('_' * 20)
    #

    ai = output[len(a)].pop(0)
    bi = carries[len(a)].pop(0)
    save, carry = half_adder(ai, bi)
    carries[len(a) + 1].append(carry)
    output[len(a)].append(save)
    #
    for i in range(len(a) + 1, len(a) + len(b)-1):
        ai = output[i].pop(0)
        bi = carries[i].pop(0)
        ci = carries[i].pop(0)
        print("i={}, {},\t{},\t{}".format(i, ai, bi, ci))
        # print(i,bi)
        save, carry = full_adder(ai, bi, ci)
        carries[i + 1].append(carry)
        output[i].insert(0, save)
    final_carry=carries[len(a)+len(b)-1].pop()
    output[len(a)+len(b)-1]=final_carry

    print("final output")
    print("output P:")
    for i, oi in enumerate(output):
        print(i, oi)


def binary_majority(boolean_values: list, n: int):
    """
    https://electronics.stackexchange.com/questions/436849/n-bit-majority-digital-logic
    :param boolean_values:
    :param n: the length of boolean_values array
    the expression will be very long.
    """
    print(boolean_values)
    k = n // 2 + 1
    tmp = []

    for selected_variables in itertools.combinations(boolean_values, k):
        print(selected_variables)
        tmp.append(Not(And(*selected_variables)))
    val = Not(And(*tmp))
    print(val)
    return val


def test_bit_majority():
    n = 8
    for i in range(100):
        a = np.random.randint(low=0, high=2, size=n, dtype=bool)
        print("a:{}".format(a.astype(int)))

        C = binary_majority(a.astype(bool), n)

        real_C = np.sum(a) * 2 >= n
        print(C, real_C)
        assert C == real_C, "they are not equal!"


def test_bit_mul(n=4):
    a = [Symbol(f'a{i}') for i in range(n)]
    b = [Symbol(f'b{i}') for i in range(n)]
    bit_mul(a, b)


if __name__ == '__main__':
    test_bit_majority()
    test_bit_mul()

