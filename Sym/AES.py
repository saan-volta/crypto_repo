"""
AES128 main file.
"""
import numpy as np
from .finite_field import *
from .key_schedule import *

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
        self.sbox = rijndael_sbox if sbox is None else sbox

        self.inv_sbox = rijndael_inv_sbox # TODO: if custom sbox provided, invert it here.

        self.mix_mat = np.array([
            [0x02, 0x03, 0x01, 0x01],
            [0x01, 0x02, 0x03, 0x01],
            [0x01, 0x01, 0x02, 0x03],
            [0x03, 0x01, 0x01, 0x02]
        ], dtype=np.uint8)

        self.inv_mix_mat = np.array([
            [0x0e, 0x0b, 0x0d, 0x09],
            [0x09, 0x0e, 0x0b, 0x0d],
            [0x0d, 0x09, 0x0e, 0x0b],
            [0x0b, 0x0d, 0x09, 0x0e]
        ], dtype=np.uint8)

        self.round_keys = keyschedule(key)

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

        # print(np.array(out, dtype=np.uint8).flatten())
        return np.array(out, dtype=np.uint8).flatten().tobytes()


    def decrypt(self, pt, iv=None, n_rounds=10):
        read_pt = np.frombuffer(pt, dtype=np.uint8)
        pt_full = np.copy(read_pt)
        pt_split = [pt_full[i:i + 16] for i in range(0, len(pt_full), 16)]

        out = []

        assert(len(pt_full) % 16 == 0)

        match self.mode:
            case 'ECB':
                for block in pt_split:
                    out.append(self.backward_block(block, n_rounds))

        return np.array(out, dtype=np.uint8).flatten().tobytes()


    def forward_block(self, block, n_rounds):
        assert (len(block) == 16)  # single block
        data = np.copy(block).reshape(4, 4).T.astype(np.uint8)
        data ^= self.round_keys[0]
        # print(f'E InitRK \n {np.vectorize(hex)(data)} \n ---')

        for i, k in zip(range(n_rounds), self.round_keys[1:n_rounds+1]):
            # 1. bytesub
            for r in range(4):
                for c in range(4):
                    b = data[r, c]
                    nib1 = b >> 4
                    nib2 = b & 0b1111
                    data[r, c] = self.sbox[nib1, nib2]

            # print(f'E/{i} ByteSub \n {np.vectorize(hex)(data)} \n ---')

            # 2. shiftrows
            for r in range(4):
                data[r] = rot_vector(data[r], r)

            # print(f'E/{i} ShiftRows \n {np.vectorize(hex)(data)} \n ---')
            # 3. mix columns
            if i < n_rounds - 1:
                for c in range(4):
                    data[:, c] = field_mat_vec_mult(self.mix_mat, data[:, c])

            # print(f'E/{i} MixCols \n {np.vectorize(hex)(data)} \n ---')

            # 4. add roundkey
            data ^= k
            # print(f'E/{i} AddRK \n {np.vectorize(hex)(data)} \n ---')

        return data.T.flatten()


    # decryption
    def backward_block(self, block, n_rounds):
        assert (len(block) == 16)  # single block
        data = np.copy(block).reshape(4, 4).T.astype(np.uint8)

        for i, k in zip(range(n_rounds), reversed(self.round_keys[1:n_rounds+1])):
            # 1. inv addroundkey
            data ^= k

            # print(f'D/{i} invAddRK \n {np.vectorize(hex)(data)} \n ---')

            # 2. inv mixcolumns
            if i > 0:
                for c in range(4):
                    data[:, c] = field_mat_vec_mult(self.inv_mix_mat, data[:, c])

            # print(f'D/{i} invMixCol \n {np.vectorize(hex)(data)} \n ---')

            # 3. inv shiftrows
            for r in range(4):
                data[r] = rot_vector(data[r], -1*r)

            # print(f'D/{i} invShiftRows \n {np.vectorize(hex)(data)} \n ---')

            # 4. inv bytesub
            for r in range(4):
                for c in range(4):
                    b = data[r, c]
                    nib1 = b >> 4
                    nib2 = b & 0b1111
                    data[r, c] = self.inv_sbox[nib1, nib2]

            # print(f'D/{i} invByteSub \n {np.vectorize(hex)(data)} \n ---')

        data ^= self.round_keys[0]
        # print(f'D/ invInitRK \n {np.vectorize(hex)(data)} \n ---')
        return data.T.flatten()





