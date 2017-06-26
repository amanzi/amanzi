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

def water_energy(xml):
    try:
        water = asearch.childByNamePath(xml, "state/field evaluators/water_content")
    except aerrors.MissingXMLError:
        pass
    else:
        ftype = asearch.childByNamePath(water, "field evaluator type")
        if ftype.value == "permafrost water content":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("three phase water content")
            else:
                ftype.setValue("liquid+ice water content")
                    
        elif ftype.value == "richards water content":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("liquid+gas water content")
            else:
                ftype.setValue("richards water content")

    try:
        energy = asearch.childByNamePath(xml, "state/field evaluators/energy")
    except aerrors.MissingXMLError:
        pass
    else:
        ftype = asearch.childByNamePath(energy, "field evaluator type")
        if ftype.value == "three phase energy":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("three phase energy")
            else:
                ftype.setValue("liquid+ice energy")
                    
        elif ftype.value == "two phase energy":
            if water.isElement("include water vapor") and water.getElement("include water vapor").value:
                ftype.setValue("liquid+gas energy")
            else:
                ftype.setValue("richards energy")

            
def adds_source_units(xml):
    fevals = asearch.childByNamePath(xml, "state/field evaluators")
    try:
        te = asearch.childByName(fevals, "total_energy_source")
    except aerrors.MissingXMLError:
        pass
    else:
        if not te.isElement("mass source units"):
            te.setParameter("mass source units", "string", "mol m^-2 s^-1")

    try:
        te = asearch.childByName(fevals, "surface-total_energy_source")
    except aerrors.MissingXMLError:
        pass
    else:
        if not te.isElement("mass source units"):
            te.setParameter("mass source units", "string", "m s^-1")


def seepage_face_bcs(xml):
    for bclist in asearch.generateElementByNamePath(xml, "boundary conditions"):
        if bclist.isElement("seepage face"):
            agroup = bclist.sublist("seepage face")[0]
            if agroup.isElement("boundary pressure"):
                bclist.sublist("seepage face").set("name", "seepage face pressure")
            elif agroup.isElement("boundary head"):
                bclist.sublist("seepage face").set("name", "seepage face head")
            else:
                raise RuntimeError("unknown seepage condition")
    

            


               

def update(xml):
    flatten_pks.flatten_pks(xml)
#    fixEvaluator(xml, "ponded_depth", "surface-ponded_depth")
#    fixEvaluator(xml, "ponded_depth_bar", "surface-ponded_depth_bar")
    compressibility(xml)
    diffusion(xml)
    water_energy(xml)
    adds_source_units(xml)
    seepage_face_bcs(xml)
        

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
    

    
