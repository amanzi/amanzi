import h5py
import numpy as np

with h5py.File("flow_darcy_well_adaptive_mesh.h5", 'r') as f_in:
    print("List of arrays in this file: \n", f_in.keys())
    x = f_in.get('x')
    x_data =  np.array(x)
    y = f_in.get('y')
    y_data =  np.array(y)

ncells = len(x_data)

np_perm_x =  np.zeros((ncells,1))
np_perm_y =  np.zeros((ncells,1))

for i in range(ncells):
   np_perm_x[i][0] = 0.1 + np.sin(x_data[i]) * 0.02
   np_perm_y[i][0] = 2.0 + np.cos(y_data[i]) * 0.4

with h5py.File('perm.h5', 'w') as f_out:
    f_out.create_dataset('permeability.cell.0', data = np_perm_x)
    f_out.create_dataset('permeability.cell.1', data = np_perm_y)
    


