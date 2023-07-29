import numpy as np

import boolexpr as bx
import itertools


def half_adder(a, b):
    save = bx.xor(a, b)
    carry = bx.and_(a, b)
    return save, carry


def full_adder(a, b, c):
    save = bx.xor(a, b, c)
    carry = bx.or_(bx.and_(c, bx.or_(a, b)), bx.and_(a, b))
    return save, carry


def int2binlist(x: int):
    x_list = [(int)(i) for i in str(bin(x))[2:]]
    x_list = x_list[::-1]
    return x_list


def crop_bits(a: list, a_min: int, out_min: int, n: int):
    # ================ Input Number ===============#
    # a = [a_0, ..., a_k] representing 
    #   (a_0 * 2^{a_min+0} + ... + a_k * 2^{a_min+k})
    # ================ Output Number ===============#
    # [b_0, ..., b_n] representing 
    #   (b_0 * 2^{out_min+0} + ... + b_n * 2^{out_min+n})

    if (out_min == None):  # not sure about the min bit
        out_min = a_min

    if (out_min < a_min):  # need extra zeros on the left
        a = [False] * (a_min - out_min) + a
        a_min = out_min

    if (n < 0):  # return all meaningful bits
        return a

    if (out_min - a_min + n > len(a)):  # need extra zeros on the right
        a = a + [False] * (out_min - a_min + n - len(a))

    return a[out_min - a_min:out_min - a_min + n]


def bin_add_int(a: list, b: list, n=-1) -> list:
    # ================ Input Numbers ===============#
    # a = [a_0, ..., a_k] representing (a_0 * 2^0 + ... + a_k * 2^k)
    # b = [b_0, ..., b_t] representing (b_0 * 2^0 + ... + b_t * 2^t)
    # ================ Output Number ===============#
    # n: number of output bits; if(n < 0), return all meaningful bits.

    if (n < 0):
        n = max([len(a), len(b)]) + 1
    carries = [0]
    saves = []

    if (n > len(a)):
        a = a + [0] * (n - len(a))
    if (n > len(b)):
        b = b + [0] * (n - len(b))

    for k in range(n):
        s, c = full_adder(a[k], b[k], carries[k])
        saves.append(s)
        carries.append(c)

    return saves


def bin_add_float(a: list, b: list,
                  a_min: int, b_min: int,
                  out_min: int, n: int) -> list:
    # ================ Input Numbers ===============#
    # a = [a_0, ..., a_k] representing 
    #   (a_0 * 2^{a_min+0} + ... + a_k * 2^{a_min+k})
    # b = [b_0, ..., b_t] representing 
    #   (b_0 * 2^{b_min+0} + ... + b_t * 2^{b_min+t})
    # a_min: lowest significant digit of a
    # b_min: lowest significant digit of b
    # ================ Output Number ===============#
    # out_min: lowest significant digit of the output 
    # n: number of TOTAL output bits if (n<0), return all meaningful digits

    # align the lowest significant digit
    a = [0] * (a_min - b_min) + a
    b = [0] * (b_min - a_min) + b

    ab_min = min(a_min, b_min)
    res = bin_add_int(a, b, -1)

    # crop saves for output
    res = crop_bits(res, ab_min, out_min, n)
    return res


def bin_mul_int(a: list, b: list, n: int):
    # ================ Input Numbers ===============#
    # a = [a_0, ..., a_k] representing (a_0 * 2^0 + ... + a_k * 2^k)
    # b = [b_0, ..., b_t] representing (b_0 * 2^0 + ... + b_t * 2^t)
    # ================ Output Number ===============#
    # n: number of output bits; if(n < 0), return all meaningful bits.

    pairwise_terms = []
    # pairwise_terms_symbols = []

    for i, ai in enumerate(a):
        pairwise_terms.append([bx.and_(ai, bj) for bj in b])
        # pairwise_terms_symbols.append([Symbol(f"p{i}_{j}") for j in range(len(b))])  

    ### TODO: add them up
    vertical_accum = pairwise_terms[0]
    res = []
    for i in range(1, len(a)):
        # cut the tail
        res.append(vertical_accum[0])
        vertical_accum = bin_add_int(vertical_accum[1:], pairwise_terms[i], -1)

    res = res + vertical_accum
    if n < 0:
        return res  # return all digits
    elif n > len(res):
        return res + [False] * (n - len(res))  # add 0 to the right
    else:
        return res[:n]  # return n digits


def bin_mul_float(a: list, b: list,
                  a_min: int, b_min: int,
                  out_min: int, n: int) -> list:
    # ================ Input Numbers ===============#
    # a = [a_0, ..., a_k] representing 
    #   (a_0 * 2^{a_min+0} + ... + a_k * 2^{a_min+k})
    # b = [b_0, ..., b_t] representing 
    #   (b_0 * 2^{b_min+0} + ... + b_t * 2^{b_min+t})
    # a_min: lowest significant digit of a
    # b_min: lowest significant digit of b
    # ================ Output Number ===============#
    # out_min: lowest significant digit of the output 
    # n: number of TOTAL output bits

    ab_min = a_min + b_min

    res = bin_mul_int(a, b, -1)

    # crop result for output
    res = crop_bits(res, ab_min, out_min, n)
    return res


def bin_mat_mul(a: list, b: list,
                a_min: int, b_min: int,
                mid_min: int, n_mid: int):
    # ================ Input Matrices  ===============#
    # a = [[a0_0, a0_1, ..., a0_k],
    #      [a1_0, a1_1, ..., a1_k],
    #      ...,
    #      [ad_0, ad_1, ..., ad_k]]
    # 
    # b = [[b0_0, b0_1, ..., b0_t],
    #      [b1_0, b1_1, ..., b1_t],
    #      ...,
    #      [bd_0, bd_1, ..., bd_t]]
    #
    # We also need: 
    # - mid_min as the lowest digit in middle steps
    # - n_mid as the number of maximum digits in middle steps

    assert len(a) == len(b), f"Dimension of a doesn't match with b, got {len(a)} != {len(b)}"

    res = bin_mul_float(a[0], b[0], a_min, b_min, mid_min, n_mid)
    for i in range(1, len(a)):
        temp = bin_mul_float(a[i], b[i], a_min, b_min, mid_min, n_mid)
        res = bin_add_float(res, temp,
                            mid_min, mid_min,
                            mid_min, n_mid)

    return res


def bit_g(a, b):
    # whether a > b
    return bx.and_(a,bx.not_(b))


def bit_l(a, b):
    # whether a < b
    return bit_g(b, a)


def bit_geq(a, b):
    return bx.not_(bit_l(a, b))


def bit_leq(a, b):
    return bit_geq(b, a)


def bit_eq(a, b):
    return bx.not_(bx.xor(a, b))


def bin_geq_int(a: list, b: list):
    # assume lowest digits are aligned already
    # align the highest digit
    a = a + [bx.ZERO] * (len(b) - len(a))
    b = b + [bx.ZERO] * (len(a) - len(b))

    res = True
    for i in range(max(len(a), len(b))):
        res = res & bit_eq(a[i], b[i])
        res = res | bit_g(a[i], b[i])
    return res


def bin_leq_int(a: list, b: list):
    # assume lowest digits are aligned already
    # align the highest digit
    return bin_geq_int(b, a)


def bin_eq_int(a: list, b: list):
    # assume lowest digits are aligned already
    # align the highest digit
    a = a + [bx.ZERO] * (len(b) - len(a))
    b = b + [bx.ZERO] * (len(a) - len(b))

    res = True
    for i in range(max(len(a), len(b))):
        res = res & bit_eq(a[i], b[i])
    return res


# =============================================================== #
# =========================== Testers =========================== #
# =============================================================== #

def test_bin_add():
    #
    ctx = bx.Context()
    #

    n = 4
    a = [ctx.get_var(f'a{i}') for i in range(n)]
    b = [ctx.get_var(f'b{i}') for i in range(n)]
    print(bin_add_float(a, b, -2, -1, 0, 0))

    a = [0, 0, 1]
    b = [0, 0, 0, 1]
    print(bin_add_float(a, b, -2, -1, -2, 1))


def test_bin_mul():
    a = [0, 0, 1]
    b = [0, 0, 0, 1]
    print(bin_mul_int(a, b, -1))

    a = [0, 0, 0, 1]
    b = [0, 0, 0, 1]
    print(bin_mul_float(a, b, 1, -2, 4, 4))


def test_bin_mat_mal():
    ctx = bx.Context()
    a = [[0, 1, 0], [0, 0, 1]]
    b = [[0, 0, 1], [1, 0, 0]]
    n = 3
    dim = 2
    a = []
    b = []
    for i in range(dim):
        a.append([ctx.get_var(f"a{i}_{j}") for j in range(n)])
    for i in range(dim):
        b.append([ctx.get_var(f"b{i}_{j}") for j in range(n)])
    res = bin_mat_mul(a, b, 0, 0, 0, -1)
    print(bin_mat_mul(a, b, 0, 0, 0, -1))


def test_bit_compare():
    for a in [True, False]:
        for b in [True, False]:
            print(f"bit_greater({a},{b}):", bit_g(a, b))
            print(f"bit_less({a},{b}):", bit_l(a, b))
            print(f"bit_geq({a},{b}):", bit_geq(a, b))
            print(f"bit_leq({a},{b}):", bit_leq(a, b))
            print(f"bit_eq({a},{b}):", bit_eq(a, b))


def test_bin_compare():
    ctx = bx.Context()
    a = np.random.randint(2, size=4).tolist()
    b = np.random.randint(2, size=3).tolist()
    print("a: ", a)
    print("b: ", b)
    print("bin_geq_int({a}, {b}): ", bin_geq_int(a, b))
    print("bin_geq_int({b}, {a}): ", bin_geq_int(b, a))

    a = [ctx.get_var(f'a{i}') for i in range(3)]
    b = [ctx.get_var(f'b{i}') for i in range(3)]
    print("bin_geq_int({a}, {b}): ", bin_geq_int(a, b))
    print("bin_geq_int({b}, {a}): ", bin_geq_int(b, a))


def test_majority():
    ctx = bx.Context()
    p = np.random.randint(2, size=7,dtype=bool).tolist()
    p += [ctx.get_var(f't')]
    print("p: ", p)
    print("Majority of p is 1: ", bx.majority(*p).sat())


if __name__ == '__main__':
    test_majority()
