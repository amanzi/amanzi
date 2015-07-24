#!/usr/bin/env python
# -------------------------------------------------------------
# file: combine_xmf.py
#
# Combines xmf files from restarted runs into one VisIt readable file.
#

import sys,os

# load etree
try:
    import elementtree.ElementTree as etree
except ImportError:
    import xml.etree.ElementTree as etree


class Step(object):
    def __init__(self, dirname, filename):
        self.dirname = dirname
        self.filename = filename
        self.full_filename = os.path.join(dirname, filename)
        self.new_filename = os.path.join(dirname, "combined_"+filename)

        self._tree = None
        self._grid = None
        try:
            self._tree = self._loadTree()
        except:
            pass
        else:
            self._grid = self._tree.getroot().find("Domain").find("Grid")

    def _loadTree(self):
        return etree.parse(self.full_filename)

    def time(self):
        if self._grid is not None:
            return self._grid.find("Time").get("Value")
        else:
            return None

    def fixItem(self, item):
        di = item.find("DataItem")
        if di is not None:
            meshname = di.text.strip().split(":")[0]
            new_meshname = os.path.join(self.dirname, meshname)
            di.text = di.text.replace(meshname, new_meshname)

    def fixStep(self):
        if self._grid is not None:
            for child in self._grid:
                self.fixItem(child)

    def write(self):
        if self._tree is not None:
            self._tree.write(self.new_filename)

class CombinedDirectory(object):
    def __init__(self, dirname, filename):
        self.dirname = dirname
        self.filename = filename
        self.full_filename = os.path.join(dirname, filename)
        self._tree = self._loadTree()
        self._grid = self._tree.getroot().find("Domain").find("Grid")

    def _loadTree(self):
        return etree.parse(self.full_filename)

    def loadSteps(self):
        self._steps = [Step(self.dirname, child.get("href")) for child in self._grid]

    def fixDirectory(self):
        for step in self._steps:
            step.fixStep()

        for (child, step) in zip(self._grid, self._steps):
            child.set("href", step.new_filename)

    def writeSteps(self):
        for step in self._steps:
            step.write()

    def clean(self, time, eps=1.e-16):
        epstime = time - eps
        while self._steps[-1].time() is not None and float(self._steps[-1].time()) >= epstime:
            self._grid.remove(self._grid.getchildren()[-1])
            self._steps.pop()

    def write(self, filename):
        self._tree.write(filename)

def combined_file(dirlist, clean_time):
    assert len(dirlist) > 0

    if clean_time > 0.:
        for lcv in range(len(dirlist)-1):
            endtime = dirlist[lcv+1]._steps[0].time()
            dirlist[lcv].clean(float(endtime), clean_time)

    primary = dirlist.pop(0)
    for el in dirlist:
        for child in el._grid:
            primary._grid.append(child)
    return primary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Combine *.VisIt.xmf files into a single VisIt readable file.")
    parser.add_argument("DIRECTORIES", nargs="+", type=str,
                        help="List of directories to combine.")

    parser.add_argument("--prefix", "-p", type=str,
                        default="visdump_data",
                        help="File prefix, PREFIX.VisIt.xmf")

    parser.add_argument("--clean_time", "-c", type=float,
                        default=-1.,
                        help="Clean the overlap, removing times that are within this tolerance of each other. REQUIRES directories to be specified in order!")

    args = parser.parse_args()
    if len(args.DIRECTORIES) < 1:
        parser.print_help()
        raise RuntimeError("Specify nonzero length list of directories.")

    # create the directory objects
    dirs = [CombinedDirectory(dirname, args.prefix+".VisIt.xmf") for dirname in args.DIRECTORIES]

    # fix the directories and steps and write the steps
    for el in dirs:
        el.loadSteps()
        el.fixDirectory()
        el.writeSteps()

    # combine the index files into one file, cleaning if requested
    combined = combined_file(dirs, args.clean_time)

    # write the combined index
    combined.write(args.prefix+".VisIt.xmf")


