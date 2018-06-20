import numpy as np
from matplotlib import pyplot as plt

def read(f):
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
        while not l.startswith("cell_sets"):
            l = l.split()
            face_areas.append(float(l.pop(-1)))
            faces.append(map(lambda a: int(a),l))
            l = fid.readline()

        l = fid.readline()
        sets = dict()
        while len(l) > 0:
            l = l.split()
            if len(l) > 1:
                name = l[0]
                ents = [int(ent) for ent in l[1:]]
                sets[name] = ents
            l = fid.readline()

    cells = np.array(cells)
    cell_volumes = np.array(cell_volumes)
    return cells, cell_volumes, faces, sets


def color(setname):
    if setname.startswith("tree"):
        return 'b'
    if setname.startswith("shrub"):
        return 'r'


def plot((cells, volumes, faces, sets), scale=100, log=False):
    if log:
        def size(cv):
            return scale * (np.log(cv) + 12) / 10
    else:
        def size(cv):
            return scale * cv
        
    plt.scatter(cells[:,0], cells[:,2], s=size(volumes), c="brown")

    for f in faces:
        if len(f) == 2:
            plt.plot([cells[f[0]][0],cells[f[1]][0]], [cells[f[0]][2],cells[f[1]][2]], 'k')

    for n, cset in sets.iteritems():
        c = color(n)
        plt.scatter(cells[cset,0], cells[cset,2], s=size(volumes[cset]), c=c)
        
