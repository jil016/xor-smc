import numpy as np
import boolexpr as bx

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
    #================ Input Number ===============#
    # a = [a_0, ..., a_k] representing 
    #   (a_0 * 2^{a_min+0} + ... + a_k * 2^{a_min+k})
    #================ Output Number ===============#
    # [b_0, ..., b_n] representing 
    #   (b_0 * 2^{out_min+0} + ... + b_n * 2^{out_min+n})
    

    if(out_min == None): # not sure about the min bit
        out_min = a_min
    
    if(out_min < a_min): # need extra zeros on the left
        a = [False]*(a_min - out_min) + a
        a_min = out_min

    if (n < 0): # return all meaningful bits
        return a
    
    if(out_min - a_min + n > len(a)): # need extra zeros on the right
        a = a + [False] *(out_min - a_min + n - len(a))
    
    return a[out_min - a_min:out_min - a_min + n]


def bin_add_int(a: list, b: list, n = -1) -> list:
    #================ Input Numbers ===============#
    # a = [a_0, ..., a_k] representing (a_0 * 2^0 + ... + a_k * 2^k)
    # b = [b_0, ..., b_t] representing (b_0 * 2^0 + ... + b_t * 2^t)
    #================ Output Number ===============#
    # n: number of output bits; if(n < 0), return all meaningful bits.

    if(n < 0):
        n = max([len(a), len(b)]) + 1
    carries = [0]
    saves = []

    if(n > len(a)):
        a = a + [0] * (n - len(a))
    if(n > len(b)):
        b = b + [0] * (n - len(b))

    for k in range(n):
        s, c = full_adder(a[k], b[k], carries[k])
        saves.append(s)
        carries.append(c)

    return saves


def majority():
    pass