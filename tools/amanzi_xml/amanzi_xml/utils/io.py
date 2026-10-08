from . import parser, errors
import xml.etree.ElementTree as etree


def _serialize_amanzi_xml(write, elem, qnames, namespaces,
                          short_empty_elements, **kwargs):
    """This is a serializer, based on etree._serialize_xml(), that
    canonicalizes according to how Teuchos::XMLObject writes XML.

    This allows reading/writing with Teuchos and reading/writing with
    amanzi_xml to result in the same output.
    """
    tag = elem.tag
    text = elem.text
    if tag is etree.Comment:
        write("<!--%s-->" % text)
    elif tag is etree.ProcessingInstruction:
        write("<?%s?>" % text)
    else:
        tag = qnames[tag]
        if tag is None:
            if text:
                write(etree._escape_cdata(text))
            for e in elem:
                _serialize_amanzi_xml(write, e, qnames, None,
                                      short_empty_elements=short_empty_elements)
        else:
            write("<" + tag)
            items = list(elem.items())
            if items or namespaces:
                if namespaces:
                    for v, k in sorted(namespaces.items(),
                                       key=lambda x: x[1]):  # sort on prefix
                        if k:
                            k = ":" + k
                        write(" xmlns%s=\"%s\"" % (
                            k,
                            etree._escape_attrib(v)
                            ))
                for k, v in items:
                    if isinstance(k, etree.QName):
                        k = k.text
                    if isinstance(v, etree.QName):
                        v = qnames[v.text]
                    else:
                        v = etree._escape_attrib(v)
                    write(" %s=\"%s\"" % (qnames[k], v))
            if text or len(elem) or not short_empty_elements:
                write(">")
                if text:
                    write(etree._escape_cdata(text))
                for e in elem:
                    _serialize_amanzi_xml(write, e, qnames, None,
                                          short_empty_elements=short_empty_elements)
                write("</" + tag + ">")
            else:
                # NOTE: this is the only current change relative to _serialize_xml()
                # write(" />")  # serialize_xml versions
                write("/>") # Teuchos version
    if elem.tail:
        write(etree._escape_cdata(elem.tail))


etree._serialize['amanzi_xml'] = _serialize_amanzi_xml


def fromFile(file_or_filename, ensure_is_plistable=False):
    """Reads a amanzi-xml hierarchy from a file or file handle"""
    etree_parser = etree.XMLParser(target=etree.TreeBuilder(insert_comments=True))
    elem = etree.parse(file_or_filename, etree_parser)

    if ensure_is_plistable:
        return parser.fromElement(elem.getroot())
    else:
        try:
            return parser.fromElement(elem.getroot())
        except Exception:
            if ensure_is_plistable:
                raise errors.NotNativeSpecError()
            else:
                return elem.getroot()


def fromString(string):
    """Reads a amanzi-xml hierarchy from a string"""
    etree_parser = etree.XMLParser(target=etree.TreeBuilder(insert_comments=True))
    elem = etree.fromstring(string, etree_parser)
    return parser.fromElement(elem)


def toString(plist):
    """Writes a amanzi-xml hierarchy to a string"""
    plist.indent(0)
    return etree.tostring(plist, encoding='unicode', method='amanzi_xml')


def toFile(plist, file_or_filename):
    """Writes a amanzi-xml hierarchy to a file"""
    plist.indent(0)
    tree = etree.ElementTree(plist)
    tree.write(file_or_filename, method='amanzi_xml')

    
def extractDoxygenXML(filename, example_header="Native Spec Example"):
    """Extremely simple parser to pull out an xml example from source."""
    with open(filename,'r') as fid:
        line = fid.readline()

        # find the start of the native spec section
        try:
            while "Native Spec Example" not in line:
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
            
    
