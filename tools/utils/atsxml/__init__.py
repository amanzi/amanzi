#from atsxml import ATSXML, replace_regions, get_root, get_regions
from atsxml import *
try:
    from xml_functions import *
except:
    print "Error: Unable to locate amanzi xml_functions module"
    print "Add $AMANZI_SRC_DIR/tools/amanzi_xml to your PYTHONPATH environment variable"
    sys.exit()
#__xall__ = ['ATSXML','Region','replace_regions']
