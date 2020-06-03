#!/usr/bin/env python

#
#  Update input files from 1.2.2 to 1.2.3
#
#  Transformations:
#
#  - Update the version number from 1.2.2 to 1.2.3
#  - Observations: "Time Macro" -> Time Macros
#                  Updates value=" " to avlu
#

import sys
import re, optparse


def CheckVersion(ncp_pre_lines):

    version_string=''

    for line in ncp_pre_lines:

        if ( "Amanzi input format version".lower() in line.lower() ):
            version_string=line.split('"')[-2]
            version_major, version_minor, version_patch = version_string.split('.')
            if ( int(version_major) == 1 and int(version_minor) == 2 or int(version_patch) == 2 or int(version_patch) == 3 ):
                break
            else:
                print('Error: This scripts translates version 1.2.2 to 1.2.3')
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
            xml_output.write("%s\n" % ( line.replace(version_string,'1.2.3' ) ) )
        else:
            xml_output.write("%s\n" % ( line ) )


    return

def WriteOutput(xml_output,Output_lines):

    Observation_Data_tag=False
    Observation_Entry_tag=False

    for line in Output_lines:

        if ( 'name="Observation Data"' in line ):
            Observation_Data_tag=True
            xml_output.write("%s\n" % ( line ) )
        elif ( Observation_Data_tag ):
            if ( 'name="Observation Output Filename"' in line ):
                xml_output.write("%s\n" % ( line ) )
            elif ( '<ParameterList name=' in line ):
                Observation_Entry_tag=True
                xml_output.write("%s\n" % ( line ) )
            elif ( Observation_Entry_tag ):
                if ( 'name="Time Macro"' in line ):
                    # extract the value
                    new_value= 'value="{'+re.findall(r'value="([^"]*)"', line)[0]+'}"'
                    # update the value
                    line=re.sub(r'value=".*"', new_value, line)
                    # update the type
                    line=re.sub(r'type="string"', 'type="Array(string)"', line)
                    xml_output.write("%s\n" % ( re.sub(r'name="Time Macro"','name="Time Macros"', line)) )
                elif ( '</ParameterList>' in line ):
                    Observation_Entry_tag=False
                    xml_output.write("%s\n" % ( line ) )
                else:
                    xml_output.write("%s\n" % ( line ) )
            elif ( '</ParameterList>' in line ):
                Observation_Data_tag=False
                xml_output.write("%s\n" % ( line ) )
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
Output_tag=False

Header_lines=[]
Output_lines=[]

for line in xml_lines:
    if ( Header_tag ):
        Header_lines.append(line.rstrip("\n"))
        if ( "name" in line and "Execution Control" in line ):
            Header_tag=False
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
    xml_output = open(xml_output_base+"_123.xml",'w')
    

WriteHeader(xml_output,Header_lines,version_string)
WriteOutput(xml_output,Output_lines)

