import sys,os

sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'],'tools','meshing_ats'))
import meshing_ats

import numpy as np
from matplotlib import pyplot as plt


# 1 km long hillslope, 10% slope
x = np.linspace(0,1000,101)
z = np.linspace(100,0,101)
print x, z
print len(x),len(z)
m2 = meshing_ats.Mesh2D.from_Transect(x,z)


#Changing organic layer thickness 
def dz_layer1(s):
    if s<100:
        thickness=0.5
    elif ((100<=s)&(s<=200)):
        thickness=-0.0045*s+0.95
    elif ((200<s)&(s<800)):
        thickness=0.05
    elif ((800<=s)&(s<=900)):
        thickness=0.0025*s-1.95
    else:
        thickness=0.3
    return thickness


# layer extrusion for 2D
layer_types = []
layer_data = []
layer_ncells = []
layer_mat_ids = []

# 50 x 2cm cells, labeled according to dz function
ncells = 50
dz = 0.02

centroid_depths = np.arange(dz/2.0, ncells*dz, dz)

for i in range(ncells):
    layer_types.append('constant') #organic
    layer_data.append(dz)
    layer_ncells.append(1)

    # labeling in the top meter varies across the domain
    layer_mat_ids1 = np.zeros((len(x)-1,),'d')
    for j in range(len(x)-1):
        if centroid_depths[i] < dz_layer1(x[j]):
            layer_mat_ids1[j] = 1001
        else:
            layer_mat_ids1[j] = 1002
    layer_mat_ids.append(layer_mat_ids1)

# layer 2
dz = 0.04
for i in range(25):
    layer_types.append('constant') #mineral
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1002 * np.ones((len(x)-1,)))

# expanding
dz = .04
for i in range(15):
    dz *= 1.2
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1002*np.ones((len(x)-1,)))
        
for i in range(4):
    dz *= 2
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(101*np.ones((len(x)-1,)))
 

print layer_data
print sum(layer_data)

layer_types.append('node')
layer_data.append(45 - sum(layer_data))
layer_ncells.append(2)
layer_mat_ids.append(101*np.ones((len(x)-1,)))

print len(layer_data)

#print layer_data
#print np.array([layer_data, np.cumsum(np.array(layer_data)), layer_mat_ids]).transpose()
#print layer_mat_ids
#print len(layer_mat_ids)
#print sum(layer_ncells)

m3 = meshing_ats.Mesh3D.extruded_Mesh2D(m2, layer_types,layer_data, layer_ncells, layer_mat_ids)
m3.write_exodus("hillslope_organic_layerbyid.exo")

