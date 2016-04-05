import sys,os
try:
    import elementtree.ElementTree as etree
except ImportError:
    import xml.etree.ElementTree as etree
import numpy as np

class ETreeParser(object):
    def __init__(self, directory='.', base='visdump_data.h5', suffix='xmf'):
        self._dir = directory
        self._base = base
        self._suffix = suffix

    def _filename(self, key):
        return os.path.join(self._dir, ".".join((self._base,key,self._suffix)))

    def _loadTree(self, key):
        return etree.parse(self._filename(key)).getroot()

    def getTime(self, key):
        return float(self._loadTree(key).getchildren()[0].getchildren()[0].find("Time").get("Value"))


def get_times(directory='.', base='visdump_data.h5'):
    import h5py
    dat = h5py.File(os.path.join(directory, base),'r')
    keys = dat[dat.keys()[0]].keys()
    keys.sort(lambda a,b: int.__cmp__(int(a),int(b)))
    parser = ETreeParser(directory, base)
    return [parser.getTime(key) for key in keys]

def readATS(directory='.', base='visdump_data.h5', inds=None):
    import h5py
    dat = h5py.File(os.path.join(directory,base),'r')
    keys = dat[dat.keys()[0]].keys()
    keys.sort(lambda a,b: int.__cmp__(int(a),int(b)))
    parser = ETreeParser(directory,base)

    if inds is None:
        times = [parser.getTime(key) for key in keys]
        return keys, times, dat
    else:
        times = [parser.getTime(keys[ind]) for ind in inds]
        keys = [keys[ind] for ind in inds]
        return keys, times, dat

def getSurfaceData(keys, dat, name):
    return np.array([dat[name][key][0] for key in keys])

