"""
AES128 main file.
"""
import numpy as np
from .finite_field import *

"""
Helper functions
"""
# left circular shift of n by d
def lcircshift(n, d, bit_width=8):
    n &= (1 << bit_width)-1 # cap
    diff = bit_width - d
    return (n >> diff) | ((n << d) & (1 << bit_width)-1)


def sbox_forward(b):
    b_inv = field_inverse(b)
    s = b_inv ^ lcircshift(b_inv,1) ^ lcircshift(b_inv,2) ^ lcircshift(b_inv,3) ^ lcircshift(b_inv,4) ^ 0x63
    return s

def sbox_backward(s):
    b_inv = lcircshift(s, 1) ^ lcircshift(s, 3) ^ lcircshift(s, 6) ^ 0x05
    b = field_inverse(b_inv)
    return b

"""
AES full implementation
"""

class AES128:
    def __init__(self, key, mode='ECB', sbox=None):
        self.mode = mode
        self.sbox = np.array([
            [99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254,215, 171, 118],
            [202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156,164, 114, 192],
            [183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113,216, 49, 21],
            [4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235,39, 178, 117],
            [9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41,227, 47, 132],
            [83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74,76, 88, 207],
            [208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80,60, 159, 168],
            [81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16,255, 243, 210],
            [205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100,93, 25, 115],
            [96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222,94, 11, 219],
            [224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145,149, 228, 121],
            [231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101,122, 174, 8],
            [186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75,189, 139, 138],
            [112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134,193, 29, 158],
            [225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206,85, 40, 223],
            [140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176,84, 187, 22]
        ], dtype=np.uint8) if sbox is None else sbox

        self.inv_sbox = np.array([
            [ 82,   9, 106, 213,  48,  54, 165,  56, 191,  64, 163, 158, 129,243, 215, 251],
            [124, 227,  57, 130, 155,  47, 255, 135,  52, 142,  67,  68, 196,222, 233, 203],
            [ 84, 123, 148,  50, 166, 194,  35,  61, 238,  76, 149,  11,  66,250, 195,  78],
            [  8,  46, 161, 102,  40, 217,  36, 178, 118,  91, 162,  73, 109,139, 209,  37],
            [114, 248, 246, 100, 134, 104, 152,  22, 212, 164,  92, 204,  93,101, 182, 146],
            [108, 112,  72,  80, 253, 237, 185, 218,  94,  21,  70,  87, 167,141, 157, 132],
            [144, 216, 171,   0, 140, 188, 211,  10, 247, 228,  88,   5, 184,179,  69,   6],
            [208,  44,  30, 143, 202,  63,  15,   2, 193, 175, 189,   3,   1,19, 138, 107],
            [ 58, 145,  17,  65,  79, 103, 220, 234, 151, 242, 207, 206, 240,180, 230, 115],
            [150, 172, 116,  34, 231, 173,  53, 133, 226, 249,  55, 232,  28,117, 223, 110],
            [ 71, 241,  26, 113,  29,  41, 197, 137, 111, 183,  98,  14, 170,24, 190,  27],
            [252,  86,  62,  75, 198, 210, 121,  32, 154, 219, 192, 254, 120,205,  90, 244],
            [ 31, 221, 168,  51, 136,   7, 199,  49, 177,  18,  16,  89,  39,128, 236,  95],
            [ 96,  81, 127, 169,  25, 181,  74,  13,  45, 229, 122, 159, 147,201, 156, 239],
            [160, 224,  59,  77, 174,  42, 245, 176, 200, 235, 187,  60, 131,83, 153,  97],
            [ 23,  43,   4, 126, 186, 119, 214,  38, 225, 105,  20,  99,  85,33,  12, 125]
        ], dtype=np.uint8) # TODO: if custom sbox provided, invert it here.

        self.mix_mat = np.array([
            [0x02, 0x03, 0x01, 0x01],
            [0x01, 0x02, 0x03, 0x01],
            [0x01, 0x01, 0x02, 0x03],
            [0x03, 0x01, 0x01, 0x02]
        ], dtype=np.uint8)

        self.round_keys = self.keyschedule(key)

    def encrypt(self, pt_bytes, iv=None, n_rounds=10):
        read_pt = np.frombuffer(pt_bytes, dtype=np.uint8)
        pt_full = np.copy(read_pt)
        pt_split = [pt_full[i:i+16] for i in range(0, len(pt_full), 16) ]

        out = []

        # assert(len(pt_full) % 16 == 0)
        last = pt_split[len(pt_split)-1]
        if len(last) < 16:
            pad = np.array([16- len(last)] * (16-len(last)))
            pt_split[len(pt_split) - 1] = np.concatenate([last, pad])


        match self.mode:
            case 'ECB':
                for block in pt_split:
                    out.append(self.forward_block(block, n_rounds))

            case 'CBC':
                assert(iv is not None)
                iv = np.copy(np.frombuffer(iv, dtype=np.uint8))
                prev = iv
                for i in range(len(pt_split)):
                    assert (len(prev) == 16)
                    enc = self.forward_block(prev ^ pt_split[i], n_rounds)
                    out.append(enc)
                    prev = out[i]



            case 'CTR':
                assert(iv is not None)
                iv = np.copy(np.frombuffer(iv, dtype=np.uint8))
                for c, block in enumerate(pt_split):
                    counter = np.zeros((16,), dtype=np.uint8)
                    counter[15] = c
                    enc = self.forward_block(iv ^ counter, n_rounds)
                    print(enc.dtype)
                    out.append(enc ^ block)

        print(np.array(out).shape)
        return np.array(out, dtype=np.uint8).flatten().tobytes()




    def forward_block(self, block, n_rounds):
        assert (len(block) == 16)  # single block
        data = np.copy(block).reshape(4, 4).T.astype(np.uint8)
        data ^= self.round_keys[0]

        for i, k in zip(range(n_rounds), self.round_keys[1:]):
            # 1. bytesub
            for r in range(4):
                for c in range(4):
                    b = data[r, c]
                    nib1 = b >> 4
                    nib2 = b & 0b1111
                    data[r, c] = self.sbox[nib1, nib2]

            #             print(f'BYTESUB\n {np.vectorize(hex)(data)} \n ---')

            # 2. shiftrows
            for r in range(4):
                data[r] = self.rot_vector(data[r], r)

            #             print(f'SHIFTROWS\n {np.vectorize(hex)(data)} \n ---')

            # 3. mix columns
            if i < 10 - 1:
                for c in range(4):
                    data[:, c] = field_mat_vec_mult(self.mix_mat, data[:, c])

            # 4. add roundkey
            data ^= k

        return data.T.flatten()


    # decryption
    def backward_block(self, block, n_rounds):
        pass


    # input: key as 16 bytes
    def keyschedule(self, bK):
        # 1. parse as matrix
        K = np.array(list(bK)).astype(np.uint8).reshape(4, 4).T
        keys = []
        keys.append(K)  # initial

        rcons = np.array([0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36])
        for i in range(10):
            prev = keys[i]
            t = prev[:, 3]
            t = self.sub_word( self.rot_vector(t, 1) )
            t[0] ^= rcons[i]
            newk = np.zeros((4, 4), dtype=np.uint8)
            newk[:, 0] = prev[:, 0] ^ t
            newk[:, 1] = prev[:, 1] ^ newk[:, 0]
            newk[:, 2] = prev[:, 2] ^ newk[:, 1]
            newk[:, 3] = prev[:, 3] ^ newk[:, 2]
            keys.append(newk)

        return keys


    def rot_word(self, x):
        return np.concatenate([x[1:], np.array([x[0]])])

    def rot_vector(self, x, i):
        return np.concatenate([x[i:], x[:i]])

    def sub_word(self, x):
        x1 = x >> 4
        x2 = x & 0b00001111
        return self.sbox[x1, x2]



