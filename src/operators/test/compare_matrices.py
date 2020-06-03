#!env python3

import numpy as np


def read_matrix(fname):
    """Read a sparse matrix in COO format"""
    with open(fname,'r') as fid:
        begin_read = False
        for line in fid:
            if line.strip().startswith('%'):
                continue

            if not begin_read:
                nrows_ncols_count = [int(ent) for ent in line.split()]
                begin_read = True
                mat = dict()
                for i in range(nrows_ncols_count[0]):
                    mat[i] = dict()
            else:
                row_col_val = line.split()
                mat[int(row_col_val[0])-1][int(row_col_val[1])-1] = float(row_col_val[2])
    return mat


def permute(mat, perm):
    return dict(*sorted([(perm[i_row],dict(*sorted([(perm[i_col], v) for (i_col, v) in cols.items()]))) for (i_row, cols) in mat.items()]))

def close(m1, m2, perm, eps=1.e-5):
    for i1 in range(len(m1)):
        i2 = perm[i1]
        col1 = m1[i1]
        col2 = m2[i2]
        if len(col1) != len(col2):
            print(p, "NO")
            return False
        for j2, v2 in col2.items():
            jj2 = perm[j2]
            if jj2 in col1:
                if abs(col1[jj2] - v2) > eps:
                    print(p, "NO")
                    return False
            else:
                print(p, "NO")
                return False
    return True

def brute_force(m1, m2):
    import itertools

    p = range(len(m1))
    if close(m1, m2, p):
        return p

    perm_seed = range(len(m1))
    for p in itertools.permutations(perm_seed):
        if close(m1, m2, p):
            return p

    return False


if __name__ == "__main__":
    import sys
    fname = sys.argv[-1]
    m1 = read_matrix(fname+"_np1.test")
    m2 = read_matrix(fname+"_np3.test")
    p = brute_force(m1, m2)
    if p:
        print(p)
        print("SUCCESS")
        sys.exit(0)
    else:
        print("FAIL")
        sys.exit(1)
        
