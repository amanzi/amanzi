import xml.etree.ElementTree as ET
import sys

#
#  Update input files from 2.1.1 to 2.2.0
#
#  Transformations:
#
#  - Update the version number from 2.1.1 to 2.2.0
#  - Structured AMR Controls: 
#            Update list of subelements name int to a list of integers within the original element
#                  Note: updates refinement_ratio, regrid_interval, blocking_factor, number_error_buffer_cells, max_grid_size
#


def usage():

    usage  = '\n'
    usage += '  UpdateSpec_211-220 Usage:\n'
    usage += '    python UpdateSpec_211-220.py oldfile.xml (newfile.xml)\n'
    usage += '\n'
    usage += '  UpdateSpec_211-220:\n'
    usage += '    - reads a 2.1.1 amanzi input file and updates it to comform to the 2.2.0 input schema\n'
    usage += '    - if the optional argument newfile.xml is not specified, the updated xml will be written out to the oldfile.xml specified\n'
    usage += '    - it is recommended not to overwrite the oldfile.xml so that updates can be compared\n'
    usage += '    - to compare the old and new files use: diff -Eb or diff -w \n'
    usage += '\n'

    return usage

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def read_input_tree(inputfile):

    try:
        tree = ET.parse(inputfile)
    except Exception:
        raise SyntaxError("Error reading inputfile - " + inputfile + "\n\n" + usage())

    return tree

def report(tree, outputfile):

    root = tree.getroot()
    indent(root, 1)
    tree.write(outputfile)

def v211_update(tree):

    # try reading in inputfile
    try:
        root = tree.getroot()
    except Exception:
        raise RuntimeError("Error getting xml tree root\n\n" + usage())

    # check version for 2.1.1
    version = root.get('version')
    if version == '2.1.0':
        # call UpdateSpec_210-211 to bring input format up to 2.1.1 spec
        import UpdateSpec_210to211
        tree = UpdateSpec_210to211.v210_update(tree)
    elif version == '2.0.0':
        #call v200_to_210 before continuing
        raise RuntimeError("Error reading input file, can not currently update from 2.0.0 input format")
    elif (version == '1.2.3' or version == '1.2.2'):
        raise RuntimeError("Error reading input file, can not currently update from 1.2.x input format")

    # update version number
    ### EIB >> update this
    root.set('version','2.2.0')

    print("\n\n  Beginning update from old 2.1.1 to 2.2.0\n\n")

    # continue to updating format

    # get simulation type
    type =  root.get('type')

    ### Output Updates ###
    
    vis = root.find('./output/vis')
    if (vis is not None):
        wr = vis.find('write_regions')
        if (wr is not None):
            print("    Update form of write_regions")
            regs = wr.text
            wr.clear()
            field = ET.SubElement(wr,'field')
            field.set("name","regions")
            field.set("regions",regs)

    ### Structured Updates ###

    if (type == "structured"):
        # if structured: update int subelements to list of integers
        amr = root.find('./numerical_controls/structured_controls/str_amr_controls')
        if (amr is not None):
            options = ['refinement_ratio', 'regrid_interval', 'blocking_factor', 'number_error_buffer_cells', 'max_grid_size']
            for child in amr:
                if child.tag in options:
                    print(("    Update form of arm control: ",child.tag))
                    ints = ""
                    # loop over subelements to get list of values, append to new element
                    for kid in child:
                        ints = ints + kid.text + " "
                    # remove the subelements
                    kids = child.findall('int')
                    for e in kids:
                        child.remove(e)
                    # update current element text
                    child.text = ints

            

    return tree


    
if __name__=="__main__":
  
    try:
        inputfile = sys.argv[1]
        try:
            outputfile = sys.argv[2]
        except Exception:
            outputfile = inputfile
            print(">> no output filename specified, will overwrite inputfile")
    except Exception:
        raise RuntimeError("Error reading arguments\n\n" + usage())

    tree = read_input_tree(inputfile)
    new_tree = v211_update(tree)
    report(new_tree, outputfile)

    print("\n\n  Use 'diff -w' or 'diff -Eb' to compare xml files\n\n")
   
