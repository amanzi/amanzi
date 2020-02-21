#!/usr/bin/env python

#
#  Update input files from 1.2.0 to 1.2.1
#
#  Transformations:
#
#  - Update the version number from 1.2.0 to 1.2.1
#  - Adjust Time Integration Mode to "Transient with Static Flow" 
#    * This is triggered when Flow Model is "Off", Transport is "On",
#      and the Time Integraion Mode is "Transient".
#  - Remove Concentration Units
#  - Checkpoint Data 
#       File Name Digits string -> int
#  - Visualization Data:
#       File Name Digist string -> int
#       Cycle Macro -> Cycle Macros; string -> Array(string); "string" -> "{string}"
#       Time Macro -> Time Macros
#
#  * Still need to add:
#       regions -> Write Regions
#       write partitions -> Write Partitions
#

import sys
import re, optparse


def CheckVersion(ncp_pre_lines):

    version_string=''

    for line in ncp_pre_lines:

        if ( "Amanzi input format version".lower() in line.lower() ):
            version_string=line.split('"')[-2]
            version_major, version_minor, version_patch = version_string.split('.')
            if ( int(version_major) == 1 and int(version_minor) == 2 or int(version_patch) == 0 ):
                break
            else:
                print('Error: This scripts translates version 1.2.0 to 1.2.1')
                print('       But this input file is version', version_string)
                sys.exit()

    if ( version_string == '' ):
        print('Error: The version string is empty!')
        sys.exit()

    return version_string

def WriteHeader(xml_output,Header_lines,version_string):

    for line in Header_lines:
       
        AmanziVersionName="Amanzi Input Format Version"
        if ( AmanziVersionName.lower() in line.lower() ):
            # case insensitive replacement
            line=re.sub("(?i)"+"amanzi input format version",AmanziVersionName,line,re.IGNORECASE)
            xml_output.write("%s\n" % ( line.replace(version_string,'1.2.1' ) ) )
        else:
            xml_output.write("%s\n" % ( line ) )


    return

def WriteModel(xml_output,Model_lines):

    # Do we need to switch time integration modes:
    FlowOff=False
    TransportOn=False
    TIModeTransient=False
    TIMode_tag=False

    for line in Model_lines:
        if ( 'name="Flow Model"' in line and 'value="Off"' in line ):
            FlowOff=True
        elif ( 'name="Transport Model"' in line and 'value="On"' in line ):
            TransportOn=True
        elif ( 'name="Time Integration Mode"' in line ) :
            TIMode_tag=True
        elif ( TIMode_tag ):
            if ( 'name="Transient"' in line ):
                TIModeTransient=True
            elif ('</ParameterList>' in line):
                TIMode_tag=False
    
    # Transient with Static
    TransientWithStatic_fix=FlowOff and TransportOn and TIModeTransient
    print('TransientWithStatic_fix ', TransientWithStatic_fix, FlowOff, TransportOn, TIModeTransient)
    
    TIMode_tag=False
    for line in Model_lines:
        if ( 'name="Concentration Units"' in line ):
            continue
        elif ( TransientWithStatic_fix ):
            if ( 'name="Flow Model"' in line ):
                xml_output.write("%s\n" % ( re.sub(r'value="Off"','value="Single Phase"', line)) )
            elif ( 'name="Time Integration Mode"' in line ):
                TIMode_tag=True
                xml_output.write("%s\n" % ( line ) )
            elif ( TIMode_tag ):
                if ( 'name="Transient"' in line ):
                    xml_output.write("%s\n" % ( re.sub(r'name="Transient"','name="Transient with Static Flow"', line)) )
                else:
                    xml_output.write("%s\n" % ( line ) )
            elif ('</ParameterList>' in line):
                TIMode_tag=False
                xml_output.write("%s\n" % ( line ) )
            else:
                xml_output.write("%s\n" % ( line ) )
        else:
            xml_output.write("%s\n" % ( line ) )
                        
    return


def WriteOutput(xml_output,Output_lines):

    Visualization_tag=False
    Checkpoint_tag=False

    for line in Output_lines:

        if ( 'name="Visualization Data"' in line ):
            Visualization_tag=True
            xml_output.write("%s\n" % ( line ) )
        elif ( 'name="Checkpoint Data"' in line ):
            Checkpoint_tag=True
            xml_output.write("%s\n" % ( line ) )
        elif ( Visualization_tag ):
            if ( "</ParameterList>" in line ):
                Visualization_tag=False
                xml_output.write("%s\n" % ( line ) )
            elif ( "File Name Digits" in line ):
                xml_output.write("%s\n" % ( re.sub(r'type=.*" ','type="int" ', line)) )
            elif ( "Cycle Macro" in line ):
                new_line = re.sub(r'name="Cycle Macro"','name="Cycle Macros"', line) 
                new_line = re.sub(r'type="string"', 'type="Array(string)"', new_line)
                xml_output.write("%s\n" % new_line )
            elif ( "Time Macro" in line ):
                xml_output.write("%s\n" % ( re.sub(r'name="Time Macro"','name="Time Macros"', line)) )
            else:
                xml_output.write("%s\n" % ( line ) )
        elif ( Checkpoint_tag ):
            if ( "</ParameterList>" in line ):
                Checkpoint_tag=False
                xml_output.write("%s\n" % ( line ) )
            elif ( "File Name Digits" in line ):
                xml_output.write("%s\n" % ( re.sub(r'type=.*" ','type="int" ', line)) )
            elif ( "Time Macro" in line ):
                xml_output.write("%s\n" % ( re.sub(r'name="Time Macro"','name="Time Macros"', line)) )
            else:
                xml_output.write("%s\n" % ( line ) )
        else:
            xml_output.write("%s\n" % ( line ) )

    return


# Create Parser, add file options 

p = optparse.OptionParser()
p.add_option("--xml_file",action="store", type='string', dest="xml_input")
p.add_option("--stdout",action="store_true", dest="xml_stdout", default=False)
p.add_option("--overwrite",action="store_true", dest="xml_overwrite", default=False)

(opt, args) = p.parse_args()
p.set_defaults(xml_input="amanzi.xml")

#
# Read the input file
#
xml_input = open(opt.xml_input)
xml_lines = xml_input.readlines()
xml_input.close()

#
# Determine level spacing:
#
lev_spc=[]
cl=0
for line in xml_lines:
    if ( "<ParameterList" in line ):
        cl=cl+1
        spaces=len(line.expandtabs()) - len(line.expandtabs().lstrip())
        if ( len(lev_spc) < cl ):
            lev_spc.append(spaces)
    elif ( "</ParameterList" in line ):
        cl=cl-1

#
# We maybe adding one new level of indentation
#
lev_spc.append(lev_spc[len(lev_spc)-1]+2)

print('Debugging: Level indentation = ', lev_spc)

#
# Break up file into sections
#
Header_tag=True
Model_tag=False
Output_tag=False

Header_lines=[]
Model_lines=[]
Output_lines=[]

for line in xml_lines:
    if ( Header_tag ):
        Header_lines.append(line.rstrip("\n"))
        if ( "name" in line and "Execution Control" in line ):
            Header_tag=False
            Model_tag=True
        #endif
    elif ( Model_tag ):
        Model_lines.append(line.rstrip("\n"))
        if ( "ParameterList" in line and "Output" in line ):
            Model_tag=False
            Output_tag=True
        #endif
    elif ( Output_tag ):
        Output_lines.append(line.rstrip("\n"))


#
# Check the version
#  
version_string=CheckVersion(Header_lines)

# Overwrite or strip xml suffix and add version 

if ( opt.xml_stdout ):
    xml_output = sys.stdout
elif ( opt.xml_overwrite ):
    xml_output = open(opt.xml_input,'w')
else:
    xml_output_base=opt.xml_input.rstrip(".xml")
    xml_output = open(xml_output_base+"_121.xml",'w')
    

WriteHeader(xml_output,Header_lines,version_string)
WriteModel(xml_output,Model_lines)
WriteOutput(xml_output,Output_lines)

