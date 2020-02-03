#!/usr/bin/env python

#
#  Translates old-style "Unstructured Algorithm" monolithic list
#  into sublists:
#
#  Transport Process Kernel
#  Chemistry Process Kernel
#  Steady-State Implicit Time Integration
#  Transient Implicit Time Integration
#  Steady-State Psuedo-Time Implicit Solver
#  Linear Solver
#  Preconditioners
#    Trilinos ML
#    Hypre AMG
#    Block ILU
#
#  Commit 4231:4c89f903ce2d on December 20, 2012 marked the 
#  transition to the new sublist format
#

import sys
import re, optparse


def PackUnstructuredLists(ncp_lines,KeyList, SubLists):


    # List of deprecated control parameters 

    deprecated=["steady error rel tol",
                "steady error abs tol",
                "transient error rel tol",
                "transient error abs tol",
                ]

    # Default preconditioner

    preconditioner = 'Trilinos ML'

    for key in KeyList:

        for line in ncp_lines:
            
            drop=False
            for d in deprecated:
                if (d in line.lower() ):
                    drop=True
                    
            if ( key.lower() in line.lower() and not drop ):
                if ( key == "linear" and "nonlinear" in line ):
                    break
                elif ( key == "linear" and "pseudo" in line ):
                    break
                elif ( key == "transport" and "chemistry" in line ):
                    break
                elif ( "use Hypre AMG" in line and 'value="true"' in line ):
                    preconditioner='Hypre AMG'
                elif ( "use Block ILU" in line and 'value="true"' in line ):
                    preconditioner='Block ILU'
                else:
                    SubLists[key].append(line.rstrip().lstrip())

    if ( len(SubLists["steady"]) > 1 and preconditioner != 'Trilinos ML' ):
        SubLists["steady"].append('<Parameter name="steady preconditioner" type="string" value="'+preconditioner+'"/>')

    if ( len(SubLists["transient"]) > 1 and preconditioner != 'Trilinos ML' ):
        SubLists["transient"].append('<Parameter name="transient preconditioner" type="string" value="'+preconditioner+'"/>')

    if ( len(SubLists["pseudo"]) > 1 and preconditioner != 'Trilinos ML' ):
        SubLists["pseudo"].append('<Parameter name="pseudo time integrator preconditioner" type="string" value="'+preconditioner+'"/>')

    return KeyList, SubLists


def WriteUnstructuredLists(xml_output,KeyList,SubLists,lev_spc):

    PreconHeader=False
    PreconFooter=False
    extra_level=0

    # need to use >>> text = "abcdef"
    #             >>> print "<%*s>" % (len(text)+2,text)
    # <  abcdef>

    text='<ParameterList name="Unstructured Algorithm">'
    xml_output.write ("%*s\n" % ( len(text)+lev_spc[3],text ) )

    for key in KeyList:

        if ( len(SubLists[key]) > 1 ):
            
            if ( key in ['ML','Hypre','Block'] and not PreconHeader ):
                xml_output.write ("          %s\n" % ( '<ParameterList name="Preconditioners">' ))
                PreconHeader=True
                extra_level=1

            line='<ParameterList name="'+SubLists[key][0]+'">'
            xml_output.write ("%*s\n" % ( len(line)+lev_spc[4+extra_level], line ) )

            for line in SubLists[key][1:]:
                xml_output.write("%*s\n" % (len(line)+lev_spc[5+extra_level], line) )
                
            line='</ParameterList>'
            xml_output.write ("%*s\n" % ( len(line)+lev_spc[4+extra_level], line ) )
                    
        if ( key is 'Block' and PreconHeader ):
            extra_level=0
            line='</ParameterList>'
            xml_output.write ("%*s\n" % ( len(line)+lev_spc[4], line ) )

    line='</ParameterList>'
    xml_output.write ("%*s\n" % ( len(line)+lev_spc[3], line ) )

    return 0


def CheckVersion(ncp_pre_lines):

    version_string=''

    for line in ncp_pre_lines:

        if ( "Amanzi input format version".lower() in line.lower() ):
            version_string=line.split('"')[-2]
            version_major, version_minor, version_patch = version_string.split('.')
            if ( int(version_major) < 1 ):
                break
            elif ( int(version_major) == 1 and int(version_minor) == 0 ):
                break
            else:
                print('Error: This scripts translates versions <= 1.0.0!')
                print('       But this input file is version', version_string)
                sys.exit()

    if ( version_string == '' ):
        print('Error: The version string is empty!')
        sys.exit()

    return version_string

def WritePreNCP(xml_output,ncp_pre_lines,version_string):

    for line in ncp_pre_lines:
        
        AmanziVersionName="Amanzi Input Format Version"
        if ( AmanziVersionName.lower() in line.lower() ):
            # case insensitive replacement
            line=re.sub("(?i)"+"amanzi input format version",AmanziVersionName,line,re.IGNORECASE)
            xml_output.write("%s\n" % ( line.replace(version_string,'1.1.0' ) ) )
        elif ( "Flow Model" in line and "Steady State Saturated" in line ):
            xml_output.write("%s\n" %( line.replace("Steady State Saturated","Single Phase") ) )
        else:
            xml_output.write("%s\n" % ( line ) )



def WritePostNCP(xml_output,ncp_pre_lines):

    for line in ncp_pre_lines:
        xml_output.write("%s\n" % ( line ) )


def CheckValidValues(xml_input):

    i=1

    MultiLine=False
    TestValue=False

    for line in xml_input:
       
        i=i+1
        # print i, line

        if ( "type=" in line and "value=" not in line ):
            MultiLine=True
            data_type=line.split('type="')[1].rstrip('" \n')
            # print "In MultiLine type = ", data_type
        elif ( "value" in line ):
            TestValue=True
            if ( MultiLine ):
                data_value=line.split('value="')[1].strip('\n"/>{}').split(',')
                # print "In Multiline value = ", data_value
            else:
                data_type=line.split('"')[3]
                data_value=line.split('"')[5].strip('{}').split(",")
                # print "In SingleLine type = ", data_type
                # print "In SingleLine value = ", data_value

        if ( TestValue ):

            # print  "data_type= ", data_type
            # print "data_value=", data_value

            TestValue=False
            MultiLine=False

            if ( "double" in data_type ):
                for v in data_value:
                    try:
                        d=float(v)
                    except:
                        print("Invalid double in line", i)
                        print(line)

            elif ("int" in data_type ):

                for v in data_value:
                    try:
                        d=int(v)
                    except:
                        print("Invalid integer in line", i)
                        print(line)



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
# Define booleans and lists for each of the three sections:
#  
#   ncp      - Numerical Control Parameters
#   ncp_pre  - Lines preceeding Numerical Control Parameters
#   ncp_post - Lines following Numerical Control Parameters
#
ncp_pre=True
ncp=False
ncp_post=False

ncp_pre_lines=[]
ncp_lines=[]
ncp_post_lines=[]

for line in xml_lines:
    if ( ncp_pre ):
        if ( "Unstructured Algorithm" in line ):
            ncp_pre=False
            ncp=True
            ncp_lines.append(line.rstrip())
        else:
            ncp_pre_lines.append(line.rstrip("\n"))
        #endif
    elif ( ncp ):
        ncp_lines.append(line.rstrip())
        if ( "ParameterList" in line ):
            ncp=False
            ncp_post=True
        #endif
    elif ( ncp_post ):
        ncp_post_lines.append(line.rstrip("\n"))
    #endif

# List of possible sublists for the "Unstructured Algorithm" sublist

KeyList  = [ "transport", "chemistry", "steady", "transient",
             "pseudo", "linear", "ML", "Hypre", "Block", 
             ]

SubLists = {"transport": ["Transport Process Kernel",] ,
            "chemistry": ["Chemistry Process Kernel",] ,
            "steady":    [ "Steady-State Implicit Time Integration", ] ,
            "transient": ["Transient Implicit Time Integration", ] , 
            "pseudo":    ["Steady-State Pseudo-Time Implicit Solver", ] ,
            "linear":    ["Linear Solver", ] , 
            "ML":        ["Trilinos ML", ] ,
            "Hypre":     ["Hypre AMG", ] , 
            "Block":     ["Block ILU", ] , 
            }


version_string=CheckVersion(ncp_pre_lines)
valid_numbers=CheckValidValues(xml_lines)

# Overwrite or strip xml suffix and add a sublist indicator

if ( opt.xml_stdout ):
    xml_output = sys.stdout
elif ( opt.xml_overwrite ):
    xml_output = open(opt.xml_input,'w')
else:
    xml_output_base=opt.xml_input.rstrip(".xml")
    xml_output = open(xml_output_base+"_UAsublists.xml",'w')
    

WritePreNCP(xml_output,ncp_pre_lines,version_string)

PackUnstructuredLists(ncp_lines,KeyList, SubLists)
WriteUnstructuredLists(xml_output,KeyList,SubLists,lev_spc)

WritePostNCP(xml_output,ncp_post_lines)
