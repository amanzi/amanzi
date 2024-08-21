# this is the script used to generate utils_reader.nc
import netCDF4
import numpy as np

f = netCDF4.Dataset('utils_reader.nc', 'w')

f.createDimension('x', 2)
f.createDimension('y1', 2)
f.createDimension('y2', 3)
f.createDimension('t', None)

v1 = f.createVariable('vec1', 'd', ['y2',])
v1[:] = np.array([0,1,2], 'd')

vi1 = f.createVariable('int_vec1', 'i', ['y2',])
vi1[:] = np.array([0,1,2], 'i')

v2 = f.createVariable('vec2', 'd', ['t', 'y2'])
v2[:] = np.array([[0,1,2], [3,4,5]], 'd')

m1 = f.createVariable('mat1', 'd', ['x', 'y1'])
m1[:] = np.array([[0,1],[2,3]])

m2 = f.createVariable('mat2', 'd', ['t', 'x', 'y1'])
m2[:] = np.array([[[0,1],[2,3]], [[2,3],[4,5]]], 'd')

v3 = f.createVariable('/group1/vec3', 'd', ['y2',])
v3[:] = np.array([0,1,2], 'd')

m3 = f.createVariable('mat3', 'd', ['x', 'y2'])
m3[:] = np.array([[0,1,2],[3,4,5]], 'd')
print(m3.shape)
print(m3[0,1])

f.close()
