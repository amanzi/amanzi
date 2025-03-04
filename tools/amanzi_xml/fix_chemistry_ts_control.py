#!/usr/bin/env python3
"""Amanzi and ATS input converter to fix timestep controller in chemistry."""

import sys, os

from amanzi_xml.utils import search as asearch
from amanzi_xml.utils import io as aio
from amanzi_xml.utils import errors as aerrors
from amanzi_xml.common import parameter, parameter_list

def addPrimaryVariable(pk_list):
    if not pk_list.isElement("primary variable key") and \
       not pk_list.isElement("primary variable key suffix"):
        pk_list.setParameter("primary variable key suffix", "string", "total_component_concentration")


def fixTSControl(pk_list):
    ts_control = pk_list.sublist("timestep controller")
    if ts_control.isElement("timestep controller type"):
        tsc_type = ts_control.getElement("timestep controller type").getValue()
        return

    elif pk_list.isElement("timestep control method"):
        method = pk_list.pop("timestep control method").getValue()
        if method == "simple":
            method = "standard"
    else:
        method = "fixed" # default in chem PK    

    ts_control.setParameter("timestep controller type", "string", method)
    tsc_params = ts_control.sublist(f"timestep controller {method} parameters")

    # parameters in the fixed timestep controller
    if method == "fixed":
        if pk_list.isElement("initial timestep (s)"):
            tsc_params.setParameter("timestep [s]", "double",
                                  pk_list.pop("initial timestep (s)").getValue())
        else:
            tsc_params.setParameter("timestep [s]", "double", 1.e16);

        # remove dead stuff
        for dead in ["min timestep (s)",
                     "max timestep (s)",
                     "timestep cut threshold",
                     "timestep increase threshold",
                     "timestep cut factor",
                     "timestep increase factor"]:
            if pk_list.isElement(dead):
                pk_list.pop(dead)

    # parameters in the simple timestep controller
    elif method == "standard":
        if pk_list.isElement("initial timestep (s)"):
            tsc_params.setParameter("initial timestep [s]", "double",
                                  pk_list.pop("initial timestep (s)").getValue())
        if tsc_params.isElement("min timestep (s)"):
            tsc_params.setParameter("min timestep [s]", "double",
                                  pk_list.pop("min timestep (s)").getValue())
        if pk_list.isElement("max timestep (s)"):
            tsc_params.setParameter("max timestep [s]", "double",
                                  pk_list.pop("max timestep (s)").getValue())
        if pk_list.isElement("timestep cut threshold"):
            tsc_params.setParameter("max iterations", "int",
                                  pk_list.pop("timestep cut threshold").getValue())
        if pk_list.isElement("timestep increase threshold"):
            tsc_params.setParameter("min iterations", "int",
                                  pk_list.pop("timestep increase threshold").getValue())
        if pk_list.isElement("timestep cut factor"):
            tsc_params.setParameter("timestep reduction factor", "double",
                                  pk_list.pop("timestep cut factor").getValue())
        if pk_list.isElement("timestep increase factor"):
            tsc_params.setParameter("timestep increase factor", "double",
                                  pk_list.pop("timestep increase factor").getValue())
            


def preOrder(xml):
    yield xml
    for child in xml:
        for c in preOrder(child):
            yield c

            
def findAllPKType(xml, pk_type):
    """Amanzi doesn't always include "PK type" in the PKs sublist,
    only in the tree.  This works, but makes it hard to find all
    PKs of a given type.
    """
    flat_pks = xml.sublist("PKs")
    
    pk_tree = asearch.find_path(xml, ["PK tree",])
    for pk in preOrder(pk_tree):
        if isinstance(pk, parameter_list.ParameterList):
            if pk.isElement("PK type"):
                if pk.getElement("PK type").getValue() == pk_type:
                    yield flat_pks.getElement(pk.getName())


def fixAll(xml):
    for pk_list in findAllPKType(xml, "chemistry alquimia"):
        fixTSControl(pk_list)
        addPrimaryVariable(pk_list)
    for pk_list in findAllPKType(xml, "chemistry amanzi"):
        fixTSControl(pk_list)
        addPrimaryVariable(pk_list)
        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    parser.add_argument("-s", "--skip", action="store_true", help="read and write only, no changes")
    args = parser.parse_args()

    # check for orig file
    print("Converting file: %s"%args.infile)
    xml = aio.fromFile(args.infile, True)

    if not args.skip:
        fixAll(xml)

    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)

