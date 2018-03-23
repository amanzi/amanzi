import numpy as np
from matplotlib import pyplot as plt

def plot(f, scale=100, log=False):
    with open(f,'r') as fid:
        l = fid.readline()
        l = fid.readline()
        cells = []
        cell_volumes = []
        while not l.startswith("face"):
            l = l.split()
            cell_volumes.append(float(l.pop(-1)))
            cells.append(map(lambda a: float(a), l))
            l = fid.readline()
        l = fid.readline()
        faces = []
        face_areas = []
        while len(l) > 0:
            l = l.split()
            face_areas.append(float(l.pop(-1)))
            faces.append(map(lambda a: int(a),l))
            l = fid.readline()
    cells = np.array(cells)
    cell_volumes = np.array(cell_volumes)
    print cell_volumes
    if log:
        plt.scatter(cells[:,0], cells[:,2], s=scale*(np.log(cell_volumes)+15))
    else:
        plt.scatter(cells[:,0], cells[:,2], s=scale*cell_volumes)
    for f in faces:
        if len(f) == 2:
            plt.plot([cells[f[0]][0],cells[f[1]][0]], [cells[f[0]][2],cells[f[1]][2]], 'k')
