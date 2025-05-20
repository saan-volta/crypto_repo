from Sym.AES import *


def main():
    # example
    # pt = 'abcdefghabcdefgh'.encode()
    ct = str_literal_to_bytes('<CIPHERTEXT>')

    for i in range(1, 256):
        for j in range(1, 256):
            bcand = bytes([i, j] + [0] * 14)
            cipher = AES128(bcand, mode='ECB')
            pt_cand = cipher.decrypt(ct, n_rounds=2, rk_ixs=[10,1])
            last = pt_cand[15]
            if last < 2:
                continue
            if pt_cand[16-last:16] == bytes([last]*last):
                print(pt_cand[:16-last].hex())



    # cipher = AES128(key, mode='ECB')
    # ct = cipher.encrypt(pt, n_rounds=10)

    # print('-'*10 + 'CIPHERTEXT' + '-'*10)
    # print(ct.hex())
    # print('-' * 10 + 'CIPHERTEXT' + '-' * 10)
    # rev_pt = cipher.decrypt(ct, n_rounds=10)
    # print(rev_pt)

    # output: 0725ca5fac66d39fac3390d488c00431


def str_literal_to_bytes(x):
    return bytes([ int('0x'+x[i:i+2], base=16) for i in range(0, len(x), 2)])


if __name__ == '__main__':
    main()
