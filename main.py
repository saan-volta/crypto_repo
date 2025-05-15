from Sym.AES import *


def main():
    # example
    pt = b'\xa1c\xf2\xa1\xa5\xaf\xba\xcc\xbc\xcc\xdcL\xae\xcbvy'
    key = b'g7\xca\x08\xb3\xbe\xaeN\xcc\xb9\xf9D.\xca\x0e\xdb'

    cipher = AES128(key, mode='ECB')
    ct = cipher.encrypt(pt)

    print('-'*10 + 'CIPHERTEXT' + '-'*10)
    print(ct.hex())
    print('-' * 10 + 'CIPHERTEXT' + '-' * 10)

    # output: 0725ca5fac66d39fac3390d488c00431


def str_literal_to_bytes(x):
    return bytes([ int('0x'+x[i:i+2], base=16) for i in range(0, len(x), 2)])


if __name__ == '__main__':
    main()
