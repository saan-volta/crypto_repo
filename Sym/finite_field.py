"""
This concerns operations done over GF(2^8), the polynomial finite field used by Rijndael.
Most of it is necessary only for computing the S-box, which is not typically needed as it can be simply stored.
1. S-box entry calculation
    1.1 poly_degree
    1.2 poly_mult
    1.3 poly_div
    1.4 field_inverse
2. Mix columns
    2.1 field_mult

"""
import numpy as np


def poly_degree(x):
    return x.bit_length() - 1


def poly_div(a, b):
    q = 0
    while poly_degree(a) >= poly_degree(b):
        shift = poly_degree(a) - poly_degree(b)
        q ^= (1 << shift)
        a ^= b << shift
    return q, a


# general polynomial mult, no field
def poly_mult(a, b):
    out = 0
    while b > 0:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


# EXT-EUCLID to find inverse in GF(2^8)
def field_inverse(a, mod=0x11B):
    if a == 0: return 0

    r0, r1 = mod, a
    s0, s1 = 0, 1

    while r1 != 0:
        q, _ = poly_div(r0, r1)
        r0, r1 = r1, r0 ^ poly_mult(q, r1)
        s0, s1 = s1, s0 ^ poly_mult(q, s1)

    assert(r0 == 1)
    return s0



def field_mult(a, b, mod=0x11B):
    # x^8+x^4+x^3+x^1+x^0 = 0x11b
    out = 0
    while b > 0:
        if b & 1:
            out ^= a
        a <<= 1
        if a & 0x100: # reduce if overflow
            a ^= mod
        b >>= 1
    return out & 0xFF


# matmul of a matrix and a vector. Used in MixColumns
def field_mat_vec_mult(M, v, mod=0x11b):
    d1, d2 = M.shape
    assert(d2 == len(v))
    out = np.zeros((d1,), dtype=int)
    for i in range(d1):
        for j in range(d2):
            out[i] ^= field_mult(M[i,j], v[j], mod=mod)
    return out


