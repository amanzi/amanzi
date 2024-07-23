import h5py
import numpy as np

nx = 50
ny = 25
ncells = nx * ny

Lx = 100
Ly = 50

x = np.zeros(nx)
y = np.zeros(ny)

for i in range(nx):
    x[i] = (i + 0.5) * (Lx / nx)

for i in range(ny):
    y[i] = (i + 0.5) * (Ly / ny)

Q = np.zeros((ncells,1))

for i in range(nx):
   for j in range(ny):
       n = i * ny + j
       Q[n][0] = np.sin(4.0 / Lx * x[i] + 1.57 / Ly * y[j]) * 2.5e-8

with h5py.File('flow_richards_source.h5', 'w') as f_out:
    f_out.create_dataset('strain_rate.cell.0//0', data = Q)
    


