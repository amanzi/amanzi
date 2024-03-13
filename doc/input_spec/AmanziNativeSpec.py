#!/usr/bin/env python

"""Generates the native spec documentation by sucking in .hh file documentation into a template"""

import sys,os
import string

sys.path.append(os.path.join(os.environ['AMANZI_SRC_DIR'],
                             "tools", "amanzi_xml"))
import amanzi_xml.utils.io

def readFileDoc(filename):
    """Parses provided file for doc strings, returning a header that starts with "//!" and any blocks that start with "/*!" """
    header = None
    doc = []
    
    with open(filename) as fid:
        # find the start indicator
        line = fid.readline()

        done = False
        while not done:
            if not line:
                done = True
                continue

            if header is None and line.strip().startswith("//!"):
                header = line.strip()[3:]

            elif line.strip().startswith("/*!"):
                example = [line.strip()[3:].split("*/")[0]+"\n",]
                if "*/" not in line:
                    line = fid.readline()
                
                    while not "*/" in line:
                        example.append(line)
                        line = fid.readline()
                        if line is None:
                            done = True
                            continue
                    example.append(line.split("*/")[0]+"\n")
                doc.extend(example)

            line = fid.readline()

    # if header is not None:
    #     doc.insert(0, header)
        
    return "".join(doc)


def findAndRead(filename):
    """Searches in AMANZI_SRC_DIR for filename, then reads that."""
    filename = filename + ".hh"
    
    for root, dirnames, filenames in os.walk(os.environ["AMANZI_SRC_DIR"]):
        if filename in filenames:
            return readFileDoc(os.path.join(os.environ["AMANZI_SRC_DIR"], root, filename))

    raise RuntimeError("file not found")


def parseTemplate(filename):
    """Parses a template, replacing with filename docs as provided."""
    fname_in = filename+".pre"
    fname_out = filename
    print("Reading template: '{}'".format(fname_in))

    with open(fname_in,'r') as fid:
        template = fid.read()

    missed_keys = []
    parsed = string.Formatter().parse(template)
    keys = [par[1] for par in parsed if par[1] is not None]

    template_dict = dict()
    for key in keys:
        try:
            template_dict[key] = findAndRead(key.strip())
        except RuntimeError:
            missed_keys.append(key)
            template_dict[key] = "** DOC GENERATION ERROR: file not found '{}' **".format(key)
                
    print("  Unable to find:", missed_keys)
    with open(fname_out,'w') as fid:
        fid.write(template.format(**template_dict))

    print("Wrote: '{}'".format(fname_out))
        

if __name__ == "__main__":
    parseTemplate(sys.argv[-1].strip(".pre"))

