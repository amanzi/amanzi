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
            print("len col {} not the same as col {}".format(i1,i2))
            return False
        for j1, v1 in col1.items():
            j2 = perm[j1]
            if j2 in col2:
                if abs(col2[j2] - v1) > eps:
                    print("val ({},{} not the same as {},{}".format(j1,v1, j2,col2[j2]))
                    return False
            else:
                print("col {} not at {}".format(j1, j2))
                return False
    return True

def diag(m1,m2, eps=1.e-5):
    valid = []
    for i1 in m1:
        valid_i1 = []
        for i2 in m2:
            if abs(m1[i1][i1] - m2[i2][i2]) < eps:
                valid_i1.append(i2)
        valid.append(valid_i1)

    for i,v in enumerate(valid):
        print(i, v)
        
        
    def gen_perm(current, leftover):
        if len(leftover) == 0:
            yield current
        else:
            for i in valid[len(current)]:
                if i in leftover:
                    current.append(i)
                    leftover.remove(i)
                    for perm in gen_perm(current, leftover):
                        yield perm
            if len(current) == 0:
                return

            leftover.add(current.pop())
            return

    for perm in gen_perm([], set(range(len(m1)))):
        print(perm)
        if close(m1, m2, perm):
            return True
    return False
    
            
        

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
    m1 = read_matrix(fname+"_np1.gold")

    print(fname+"_np1.test")
    for r, col in m1.items():
        print(r, ''.join(["({0}, {1})".format(*c) for c in col.items()]))

    m2 = read_matrix(fname+"_np3.gold")

    print(fname+"_np3.test")
    for r, col in m2.items():
        print(r, ''.join(["({0}, {1})".format(*c) for c in col.items()]))

    p = diag(m1, m2)
    if p:
        print(p)
        print("SUCCESS")
        sys.exit(0)
    else:
        print("FAIL")
        sys.exit(1)
        
