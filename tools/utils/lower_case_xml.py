"""
Converts native-spec xml files from old capitalized form to new lower-case form.

Usage:

export AMANZI_SRC_DIR=/path/to/my/amanzi/repo
python lower_case_xml.py my_native_file.xml ... my_other_native_file.xml

NOTE: all work is done in-place!
"""


import sys,os
sys.path.append(os.path.join(os.environ["AMANZI_SRC_DIR"],'tools', 'amanzi_xml'))

import amanzi_xml.utils.errors as xml_errors
import amanzi_xml.utils.search as xml_search
import amanzi_xml.utils.io as xml_io




changes_regions = ["Coordinate",
                   "Low Coordinate",
                   "High Coordinate",
                   "Location",
                   "Direction",
                   "Tolerance",
                   "Entity",
                   "File",
                   "Format",
                   "Label",
                   "Value",
                   "Number of points",
                   "Points",
                   "Union",
                   "Entity GIDs",]

changes_mesh = ["Number of Cells",
                "Domain Low Corner",
                "Domain High Corner",
                "Domain Low Coordinate",
                "Domain High Coordinate",
                "Framework",
                "Verify Mesh",
                "File",
                "Format"]


def lower_case(xml):
    """Converts an xml object, in-place"""

    regions = xml_search.getElementByNamePath(xml, "Regions")
    for label in changes_regions:
        for param in xml_search.generateElementByNamePath(regions, label):
            param.set('name', label.lower())

    for param in xml_search.generateElementByNamePath(xml, "Region: Plane/location"):
        param.set('name', 'point')
    for param in xml_search.generateElementByNamePath(xml, "Region: Plane/direction"):
        param.set('name', 'normal')
        
        
    mesh = xml_search.getElementByNamePath(xml, "Mesh")
    for label in changes_mesh:
        for param in xml_search.generateElementByNamePath(mesh, label):
            param.set('name', label.lower())

    


if __name__ == "__main__":
    fnames = [f for f in sys.argv if f.endswith('.xml')]
    for f in fnames:
        try:
            xml = xml_io.fromFile(f)
        except xml_errors.NotNativeSpecError:
            pass
        except RuntimeError,msg:
            print "Attempted to convert invalid xml file:", f
            print "  Error:", msg

        else:        
            lower_case(xml)
            xml_io.toFile(xml, f)





