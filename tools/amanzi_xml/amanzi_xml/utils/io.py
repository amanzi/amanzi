import parser, errors
try:
    import elementtree.ElementTree as etree
except ImportError:
    import xml.etree.ElementTree as etree


def fromFile(file_or_filename):
    """Reads a amanzi-xml hierarchy from a file or file handle"""
    elem = etree.parse(file_or_filename)

    try:
        return parser.fromElement(elem.getroot())
    except:
        return elem.getroot()
    # except RuntimeError, msg:
    #     if "amanzi_input" in msg.__str__():
    #         raise errors.NotNativeSpecError()
    #     else:
    #         raise RuntimeError(msg)

def fromString(string):
    """Reads a amanzi-xml hierarchy from a string"""
    elem = etree.fromstring(string)
    return parser.fromElement(elem)

def toString(plist):
    """Writes a amanzi-xml hierarchy to a string"""
    plist.indent(0)
    return etree.tostring(plist)

def toFile(plist, file_or_filename):
    """Writes a amanzi-xml hierarchy to a file"""
    plist.indent(0)
    tree = etree.ElementTree(plist)
    tree.write(file_or_filename)

def extractDoxygenXML(filename, example_header="Native Spec Example"):
    """Extremely simple parser to pull out an xml example from source."""
    with open(filename,'r') as fid:
        line = fid.readline()

        # find the start of the native spec section
        try:
            while not "Native Spec Example" in line:
                line = fid.readline()
        except StopIteration:
            return None

        # strip any blank lines
        line = fid.readline()
        while line.strip() == "":
            line = fid.readline()

        if not line.startswith("    "):
            return None

        # start reading
        lines = []
        while line.startswith("    "):
            lines.append(line)
            line = fid.readline()

        return fromString("".join(lines))
            
    
