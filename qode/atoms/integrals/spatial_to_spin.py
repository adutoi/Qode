import numpy

def one_electron_blocked(O_raw, n_orb_map):
    norb1, norb2 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2))
    for p in range(norb1):
        for q in range(norb2):
            O[p, q]             = O_raw[p, q]
            O[p+norb1, q+norb2] = O_raw[p, q]
    return O

def one_electron_alternating(O_raw, n_orb_map):
    norb1, norb2 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2))
    for p in range(norb1):
        for q in range(norb2):
            O[2*p,   2*q  ] = O_raw[p, q]
            O[2*p+1, 2*q+1] = O_raw[p, q]
    return O

def one_electron_alternating_blocks(O_raw, n_orb_map):
    # this function orders as frozen_a, frozen_b, occ_a, occ_b, virt_a and virt_b, which is the ordering adcc densities require
    if "o3" in n_orb_map.keys():
        #print("frozen occs used for blocking")
        orb_type_order = ("o3", "o1", "v1")  # frozen_occ, occ, virt
    else:
        orb_type_order = ("o1", "v1")  # occ, virt
    a_sl = []
    b_sl = []
    r_sl = []  # restricted slices
    for orb_type in orb_type_order:
        n = n_orb_map[orb_type]
        if len(b_sl) == 0:
            a_sl.append((0, n))
            b_sl.append((n, 2 * n))
            r_sl.append((0, n))
        else:
            offset = b_sl[-1][1]
            a_sl.append((offset, n + offset))
            offset = a_sl[-1][1]
            b_sl.append((offset, n + offset))
            offset = r_sl[-1][1]
            r_sl.append((offset, offset + n))
            
    #print(n_orb_map, r_sl, a_sl, b_sl)
    norb1, norb2 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2))
    for p in range(len(orb_type_order)):
        for q in range(len(orb_type_order)):
            O[a_sl[p][0]:a_sl[p][1], a_sl[q][0]:a_sl[q][1]] =\
                O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1]]  # aa
            O[b_sl[p][0]:b_sl[p][1], b_sl[q][0]:b_sl[q][1]] =\
                O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1]]  # bb
    #print("O_raw", O_raw)
    #print("O_final_00", O[:norb1, :norb2])
    #print("O_final_11", O[norb1:, norb2:])
    return O

def two_electron_blocked(O_raw, n_orb_map):
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

def two_electron_alternating(O_raw, n_orb_map):
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

def two_electron_alternating_blocks(O_raw, n_orb_map):
    # this function orders as frozen_a, frozen_b, occ_a, occ_b, virt_a and virt_b, which is the ordering adcc densities require
    if "o3" in n_orb_map.keys():
        orb_type_order = ("o3", "o1", "v1")  # frozen_occ, occ, virt
    else:
        orb_type_order = ("o1", "v1")  # occ, virt
    a_sl = []
    b_sl = []
    r_sl = []  # restricted slices
    for orb_type in orb_type_order:
        n = n_orb_map[orb_type]
        if len(b_sl) == 0:
            a_sl.append((0, n))
            b_sl.append((n, 2 * n))
            r_sl.append((0, n))
        else:
            offset = b_sl[-1][1]
            a_sl.append((offset, n + offset))
            offset = a_sl[-1][1]
            b_sl.append((offset, n + offset))
            offset = r_sl[-1][1]
            r_sl.append((offset, offset + n))
            
    #print(n_orb_map, r_sl, a_sl, b_sl)
    norb1, norb2, norb3, norb4 = O_raw.shape
    O = numpy.zeros((2*norb1, 2*norb2, 2*norb3, 2*norb4))
    for p in range(len(orb_type_order)):
        for q in range(len(orb_type_order)):
            for r in range(len(orb_type_order)):
                for s in range(len(orb_type_order)):
                    O[a_sl[p][0]:a_sl[p][1], a_sl[q][0]:a_sl[q][1], a_sl[r][0]:a_sl[r][1], a_sl[s][0]:a_sl[s][1]] =\
                        O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1], r_sl[r][0]:r_sl[r][1], r_sl[s][0]:r_sl[s][1]]  # aaaa
                    O[a_sl[p][0]:a_sl[p][1], b_sl[q][0]:b_sl[q][1], a_sl[r][0]:a_sl[r][1], b_sl[s][0]:b_sl[s][1]] =\
                        O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1], r_sl[r][0]:r_sl[r][1], r_sl[s][0]:r_sl[s][1]]  # abab
                    O[b_sl[p][0]:b_sl[p][1], a_sl[q][0]:a_sl[q][1], b_sl[r][0]:b_sl[r][1], a_sl[s][0]:a_sl[s][1]] =\
                        O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1], r_sl[r][0]:r_sl[r][1], r_sl[s][0]:r_sl[s][1]]  # baba
                    O[b_sl[p][0]:b_sl[p][1], b_sl[q][0]:b_sl[q][1], b_sl[r][0]:b_sl[r][1], b_sl[s][0]:b_sl[s][1]] =\
                        O_raw[r_sl[p][0]:r_sl[p][1], r_sl[q][0]:r_sl[q][1], r_sl[r][0]:r_sl[r][1], r_sl[s][0]:r_sl[s][1]]  # bbbb
    return O
