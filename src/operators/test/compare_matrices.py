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

def read_mesh_info(fname):
    """Read a sparse matrices mesh info file to get coordinates"""
    class Point(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)
            self.coord = np.array([float(self.x),float(self.y),float(self.z)])


    with open(fname, 'r') as fid:
        points = [Point(**dict([(k,v) for (k,v) in
                                zip(['rank', 'etype', 'lid', 'gid', 'smap_lid', 'smap_gid', 'x', 'y', 'z'],
                                    line.split())]))
                  for line in fid if not line.strip().startswith('#')]
    return dict([(int(p.smap_gid), p.coord) for p in points])


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

def calc_permutation(d1,d2, eps=1.e-5):
    perm = dict()
    for gid, coord in d1.items():
        done = False
        for gid2, coord2 in d2.items():
            if np.allclose(coord, coord2, eps):
                perm[gid] = gid2
                done = True
                print('found match {}:{} at {}'.format(gid,gid2,coord))
                break

        if not done:
            raise RuntimeError("Could not find match for point {}".format(coord))
        else:
            d2.pop(gid2)
    return perm


if __name__ == "__main__":
    import sys

    fname1 = sys.argv[-1]
    m1 = read_matrix(fname1+".dat")
    d1 = read_mesh_info(fname1+"_map.dat")
    


    # print(fname1)
    # for r, col in m1.items():
    #     print(r, ''.join(["({0}, {1})".format(*c) for c in col.items()]))

    fname2 = sys.argv[-2]
    m2 = read_matrix(fname2+".dat")
    d2 = read_mesh_info(fname2+"_map.dat")

    # print(fname2)
    # for r, col in m2.items():
    #     print(r, ''.join(["({0}, {1})".format(*c) for c in col.items()]))

    p = calc_permutation(d1,d2)
    
    if close(m1,m2,p):
        print("SUCCESS")
        sys.exit(0)
    else:
        print("FAIL")
        sys.exit(1)
        
