"""
AES key scheduling and reversal
"""
import numpy as np
from .finite_field import rijndael_sbox, rijndael_inv_sbox

def keyschedule(bK):
    # 1. parse as matrix
    K = np.array(list(bK)).astype(np.uint8).reshape(4, 4).T
    keys = []
    keys.append(K)  # initial

    rcons = np.array([0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36])
    for i in range(10):
        prev = keys[i]
        t = prev[:, 3]
        t = sub_word(rot_vector(t, 1))
        t[0] ^= rcons[i]
        newk = np.zeros((4, 4), dtype=np.uint8)
        newk[:, 0] = prev[:, 0] ^ t
        newk[:, 1] = prev[:, 1] ^ newk[:, 0]
        newk[:, 2] = prev[:, 2] ^ newk[:, 1]
        newk[:, 3] = prev[:, 3] ^ newk[:, 2]
        keys.append(newk)

    return keys


def rot_word(x):
    return np.concatenate([x[1:], np.array([x[0]])])


def rot_vector(x, i):
    return np.concatenate([x[i:], x[:i]])


def sub_word(x):
    x1 = x >> 4
    x2 = x & 0b00001111
    return rijndael_sbox[x1, x2]

def inv_sub_word(x):
    x1 = x >> 4
    x2 = x & 0b00001111
    return rijndael_inv_sbox[x1, x2]


def reverse_round(rk, rcon):
    prevk = np.zeros((4,4), dtype=np.uint8)

    prevk[:, 3] = rk[:, 2] ^ rk[:, 3]
    t = prevk[:, 3]
    t = sub_word(rot_vector(t, 1))
    t[0] ^= rcon

    prevk[:, 2] = rk[:, 1] ^ rk[:, 2]
    prevk[:, 1] = rk[:, 0] ^ rk[:, 1]
    prevk[:, 0] = t ^ rk[:, 0]
    return prevk

def reverse_schedule(rk_last):
    rcons = np.array([0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36])
    rk = rk_last

    for i in range(9, -1, -1):
        rk = reverse_round(rk, rcons[i])

    return rk.T.flatten()
