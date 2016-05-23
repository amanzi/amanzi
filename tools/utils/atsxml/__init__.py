import sys,os
if os.environ.has_key("AMANZI_SRC_DIR"):
    sys.path.append(os.path.join(os.environ["AMANZI_SRC_DIR"], "tools", "amanzi_xml"))

from atsxml import *
from amanzi_xml.utils.search import *
