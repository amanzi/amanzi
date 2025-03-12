import sys,os
import xml.etree.ElementTree as etree
from amanzi_xml.utils import io as aio

def fromFile(file_or_filename):
    etree_parser = etree.XMLParser(target=etree.TreeBuilder(insert_comments=True))
    return etree.parse(file_or_filename, etree_parser)

def toFile(xml, file_or_filename):
    xml.write(file_or_filename, method='amanzi_xml')



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="input filename")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true", help="fix file in place")
    group.add_argument("-o", "--outfile", help="output filename")

    args = parser.parse_args()

    # check for orig file
    print("Converting file: %s"%args.infile)
    xml = fromFile(args.infile)

    if args.inplace:
        toFile(xml, args.infile)
    else:
        toFile(xml, args.outfile)
    sys.exit(0)

