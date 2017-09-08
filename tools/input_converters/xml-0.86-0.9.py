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
        pass
    else:
        pd.set("name", newname)
    

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
        etype = asearch.childByNamePath(energy, "field evaluator type")
        if ftype.value == "three phase water content":
            etype.setValue("three phase energy")
        elif ftype.value == "liquid+ice water content":
            etype.setValue("liquid+ice energy")
        elif ftype.value == "liquid+gas water content":
            etype.setValue("liquid+gas energy")
        elif ftype.value == "richards water content":
            etype.setValue("richards energy")

            
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

def primary_variable(xml):
    for pv in asearch.generateElementByNamePath(xml, "primary variable"):
        pv.name = "primary variable key"
    for pv in asearch.generateElementByNamePath(xml, "primary variable key"):
        if pv.value == "ponded_depth":
            pv.setValue("surface-ponded_depth")
            
def mesh_list(xml):
    try:
        xml.pop("Native Unstructured Input")
    except aerrors.MissingXMLError:
        pass

    try:
        xml.pop("grid_option")
    except aerrors.MissingXMLError:
        pass

    # move domain mesh parameters to a sublist
    mesh = asearch.childByName(xml, "mesh")
    domain = mesh.sublist("domain")
    to_pop = []
    for el in mesh:
        if el.get("name") == "surface mesh":
            el.set("name", "surface")
        elif el.get("name") == "column meshes":
            el.set("name", "column")
        elif el.get("name") == "column surface meshes":
            el.set("name", "column surface")
        elif el.get("name") in ["surface", "column", "column surface", "domain", "subgrid"]:
            pass
        elif el.get("name") in ["framework"]:
            to_pop.append(el)
        else:
            domain.append(el)

    if domain.isElement("framework"):
        domain.pop("framework")
            
    for el in domain:
        if mesh.isElement(el.get("name")):
            mesh.pop(el.get("name"))
    for el in to_pop:
        mesh.pop(el.get("name"))        

    # move surface mesh parameters to a sublist
    if mesh.isElement("surface"):
        surf_list = mesh.sublist("surface")
        surf_p_list = surf_list.sublist("surface")
        if surf_list.isElement("surface sideset name"):
            surf_p_list.append(surf_list.pop("surface sideset name"))
        if surf_list.isElement("surface sideset names"):
            surf_p_list.append(surf_list.pop("surface sideset names"))
        
    # make sure all left are mesh sublists, add a mesh type parameter
    valid_types = ["read mesh file", "generate mesh", "logical mesh",
                   "aliased", "surface", "column", "column surface", "subgrid"]    
    for el in mesh:
        if not el.isElement("mesh type"):
            assert(el.get("type") == "ParameterList")
            found = False
            for valid_type in valid_types:
                if el.isElement(valid_type):
                    print "setting type: ", valid_type
                    el.setParameter("mesh type", "string", valid_type)
                    asearch.childByName(el, valid_type).set("name", valid_type+" parameters")
                    found = True
                    continue
            assert(found)
        
    return


def vis(xml):
    if xml.isElement("visualization") and not asearch.childByName(xml, "visualization").isElement("domain"):
        vis_domain = xml.pop("visualization")
        vis_domain.set("name", "domain")
    
        vis_list = xml.sublist("visualization")
        vis_list.append(vis_domain)
        if xml.isElement("visualization surface"):
            vis_surf = xml.pop("visualization surface")
            vis_surf.set("name", "surface")
            vis_list.append(vis_surf)

        if xml.isElement("visualization columns"):
            vis_col = xml.pop("visualization columns")
            vis_col.set("name", "column_*")
            vis_list.append(vis_col)

def update(xml):
    flatten_pks.flatten_pks(xml)
    fixEvaluator(xml, "ponded_depth", "surface-ponded_depth")
    fixEvaluator(xml, "ponded_depth_bar", "surface-ponded_depth_bar")
    fixEvaluator(xml, "manning_coefficient", "surface-manning_coefficient")
    fixEvaluator(xml, "unfrozen_fraction", "surface-unfrozen_fraction")
    fixEvaluator(xml, "unfrozen_effective_depth", "surface-unfrozen_effective_depth")
    compressibility(xml)
    diffusion(xml)
    water_energy(xml)
    adds_source_units(xml)
    seepage_face_bcs(xml)
    primary_variable(xml)
    mesh_list(xml)
    vis(xml)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix a number of changes from ATS input spec 0.86 to 0.9x")
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    print "Converting file: %s"%args.infile
    xml = aio.fromFile(args.infile, True)
    update(xml)
    if args.inplace:
        aio.toFile(xml, args.infile)
    else:
        aio.toFile(xml, args.outfile)
    sys.exit(0)
    

    
