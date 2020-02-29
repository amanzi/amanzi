"""Recovers the input .xml file used in an ATS run from the logfile
generated when the run was made.

This can be useful to recover from accidental deletion of the input
file!
"""

import xml.etree.ElementTree as ElementTree

startstring = "======================> dumping parameter list <======================"
endstring = "======================> done dumping parameter list. <================"


def read_lines(fid):
    """Reads the raw plist from the file into a string"""
    xml_lines = []
    line = fid.readline()
    while line is not None and not line.strip().startswith(startstring):
        line = fid.readline()

    if line is not None:
        line = fid.readline() # pops the startstring
        while line is not None and not line.strip().startswith(endstring):
            xml_lines.append(line)
            line = fid.readline()

    if line is None:
        raise RuntimeError("Malformed logfile: cannot find the starting string and ending string.")

    return "".join(xml_lines)

def get_xml(fname):
    """Reads the file and extracts an XML"""
    with open(fname, 'r') as fid:
        xml_string = read_lines(fid)

    print(xml_string)
    elem = ElementTree.fromstring(xml_string)
    return elem

        
def cleanup(xml):
    """Strips a bunch of extra stuff trilinos adds."""
    to_clean = ['isUsed', 'isDefault', 'id', 'docString']
    def gen(xml):
        for e in xml:
            for x in gen(e):
                yield x
            yield e

    for e in gen(xml):
        for k in to_clean:
            if k in e.keys():
                e.attrib.pop(k)

    xml.attrib['name'] = 'Main'

def write(xml, filename):
    """Writes the xml to file."""
    tree = ElementTree.ElementTree(xml)
    tree.write(filename)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("LOG_FILE", type=str,
                        help="Logfile to recover from.")
    parser.add_argument("OUTPUT_FILE", type=str,
                        help="XML filename to write.")
    args = parser.parse_args()

    xml = get_xml(args.LOG_FILE)
    cleanup(xml)
    write(xml, args.OUTPUT_FILE)
    exit(0)
    
