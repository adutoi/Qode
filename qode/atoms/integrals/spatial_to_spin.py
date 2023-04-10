import numpy

def one_electron_blocked(O_raw):
    norb1, norb2 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2))
    for p in range(norb1):
        for q in range(norb2):
            O[p, q]             = O_raw[p, q]
            O[p+norb1, q+norb2] = O_raw[p, q]
    return O

def one_electron_alternating(O_raw):
    norb1, norb2 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2))
    for p in range(norb1):
        for q in range(norb2):
            O[2*p,   2*q  ] = O_raw[p, q]
            O[2*p+1, 2*q+1] = O_raw[p, q]
    return O

def two_electron_blocked(O_raw):
    norb1, norb2, norb3, norb4 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2, 2*norb3, 2*norb4))
    for p in range(norb1):
        for q in range(norb2):
            for r in range(norb3):
                for s in range(norb4):
                    O[p, q, r, s]                         = O_raw[p, q, r, s]
                    O[p+norb1, q, r+norb3, s]             = O_raw[p, q, r, s]
                    O[p, q+norb2, r, s+norb4]             = O_raw[p, q, r, s]
                    O[p+norb1, q+norb2, r+norb3, s+norb4] = O_raw[p, q, r, s]
    return O

def two_electron_alternating(O_raw):
    norb1, norb2, norb3, norb4 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2, 2*norb3, 2*norb4))
    for p in range(norb1):
        for q in range(norb2):
            for r in range(norb3):
                for s in range(norb4):
                    O[2*p,   2*q,   2*r,   2*s  ] = O_raw[p, q, r, s]
                    O[2*p+1, 2*q,   2*r+1, 2*s  ] = O_raw[p, q, r, s]
                    O[2*p,   2*q+1, 2*r,   2*s+1] = O_raw[p, q, r, s]
                    O[2*p+1, 2*q+1, 2*r+1, 2*s+1] = O_raw[p, q, r, s]
    return O
