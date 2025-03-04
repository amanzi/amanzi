import numpy as np
import h5py

hout = h5py.File('aperture_dynamic.h5', 'w')
gout = hout.create_group('fracture-aperture.cell.0')

# times for aperture are hardcoded
ncells = 18
aperture0 = np.zeros((ncells, 1))
aperture1 = np.zeros((ncells, 1))
for i in range(ncells):
    aperture0[i][0] = 1e-5
    aperture1[i][0] = 1.2e-5

times = [0.0, 3.0e+9]

gout.create_dataset('0', data=aperture0)
gout.create_dataset('1', data=aperture1)
hout.create_dataset('time', data=times)
hout.close()
