# this is the script used to generate test.h5
import h5py
import numpy as np

f = h5py.File('test.h5', 'w')

v1 = f.create_dataset('vec1', data=np.array([0,1,2], 'd'))
intv1 = f.create_dataset('int_vec1', data=np.array([0,1,2], 'i'))

v2 = f.create_group('vec2')
v2.create_dataset('0', data=np.array([0,1,2], 'd'))
v2.create_dataset('1', data=np.array([3,4,5], 'd'))

m1 = f.create_dataset('mat1', data=np.array([[0,1], [2,3]], 'd'))

m2 = f.create_group('mat2')
m2.create_dataset('0', data=np.array([[0,1], [2,3]], 'd'))
m2.create_dataset('1', data=np.array([[2,3], [4,5]], 'd'))

grp = f.create_group('group1')
grp.create_dataset('vec3', data=np.array([0,1,2], 'd'))

m3 = f.create_dataset('mat3', data=np.array([[0,1,2], [3,4,5]],'d'))

f.close()
