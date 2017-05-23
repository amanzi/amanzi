#!/usr/bin/env python

"""ATS input converter from 0.86 to 0.90"""

import sys, os
try:
    amanzi_xml = os.path.join(os.environ["AMANZI_SRC_DIR"], "tools","amanzi_xml")
except KeyError:
    pass
else:
    if amanzi_xml not in sys.path:
        sys.path.append(amanzi_xml)


from amanzi_xml.utils import search as asearch
from amanzi_xml.utils import io as aio
from amanzi_xml.utils import errors as aerrors

import flatten_pks

def fixEvaluator(xml, name, newname):
    try:
        pd = asearch.childByNamePath(xml, "state/field evaluators/%s"%name)
    except aerrors.MissingXMLError:
        pd.set("name", newname)
    pass
    

def compressibility(xml):
    for c in asearch.generateElementByNamePath(xml, "pore compressibility"):
        c.set("name", "pore compressibility [Pa^-1]")

def diffusion(xml):
    for diff in asearch.generateElementByNamePath(xml, "Diffusion"):
        diff.set("name", "diffusion")
    for diff in asearch.generateElementByNamePath(xml, "Diffusion PC"):
        diff.set("name", "diffusion preconditioner")


def update(xml):
    flatten_pks.flatten_pks(xml)
#    fixEvaluator(xml, "ponded_depth", "surface-ponded_depth")
#    fixEvaluator(xml, "ponded_depth_bar", "surface-ponded_depth_bar")
    compressibility(xml)

        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.86 to 0.9x")
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    xml = aio.fromFile(args.infile)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
