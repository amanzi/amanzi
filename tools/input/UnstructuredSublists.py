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

def PackUnstructuredLists(KeyList, SubLists):

    for key in KeyList:
        for line in ncp_lines:
            if ( key.lower() in line.lower() ):
                if ( key == "linear" and "nonlinear" in line ):
                    break
                else:
                    SubLists[key].append(line.rstrip().replace("\t","        "))

    return KeyList, SubLists


def WriteUnstructuredLists(xml_output,KeyList,SubLists):

    PreconHeader=False
    PreconFooter=False
    TwoSpaces=''

    xml_output.write ("        %s\n" % ( '<ParameterList name="Unstructured Algorithm">' ))

    for key in KeyList:

        if ( len(SubLists[key]) > 1 ):
            
            if ( key in ['ML','Hypre','Block'] and not PreconHeader ):
                xml_output.write ("          %s\n" % ( '<ParameterList name="Preconditioners">' ))
                PreconHeader=True
                TwoSpaces='  '
            xml_output.write ("%s          %s%s%s\n" % ( TwoSpaces, '<ParameterList name="', SubLists[key][0], '">' ))

            for line in SubLists[key][1:]:
                xml_output.write("%s  %s \n" % (TwoSpaces, line) )
                
            xml_output.write ("%s          %s\n" % (TwoSpaces, '</ParameterList>' ) )
                    
        if ( key is 'Block' and PreconHeader ):
            TwoSpaces=''
            xml_output.write ("%s          %s\n" % (TwoSpaces, '</ParameterList>' ) )
            
    xml_output.write ("        %s\n" % ( '</ParameterList>' ) )

    return 0


def WritePreNCP(xml_output,ncp_pre_lines):

    for line in ncp_pre_lines:
        xml_output.write("%s\n" % ( line ) )

def WritePostNCP(xml_output,ncp_pre_lines):

    for line in ncp_pre_lines:
        xml_output.write("%s\n" % ( line ) )


# Create Parser, add file options 

p = optparse.OptionParser()
p.add_option("--xml_file",action="store", type='string', dest="xml_input")
p.add_option("--stdout",action="store_false", dest="xml_stdout", default=False)

(opt, args) = p.parse_args()
p.set_defaults(xml_input="amanzi.xml")

# strip xml suffix and add a sublist indicator

if ( opt.xml_stdout ):
    xml_output = sys.stdout
else:
    xml_output_base=opt.xml_input.rstrip(".xml")
    xml_output = open(xml_output_base+"_UAsublists.xml",'w')
    
#
# Read the input file
#
xml_input = open(opt.xml_input)
xml_lines = xml_input.readlines()
xml_input.close()

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
            ncp_pre_lines.append(line.rstrip())
        #endif
    elif ( ncp ):
        ncp_lines.append(line.rstrip())
        if ( "ParameterList" in line ):
            ncp=False
            ncp_post=True
        #endif
    elif ( ncp_post ):
        ncp_post_lines.append(line.rstrip())
    #endif

# List of possible sublists for the "Unstructured Algorithm" sublist

KeyList  = [ "transport", "chemistry", "steady", "transient",
             "psuedo", "linear", "ML", "Hypre", "Block", 
             ]

SubLists = {"transport": ["Transport Process Kernel",] ,
            "chemistry": ["Chemistry Process Kernel",] ,
            "steady":    [ "Steady-State Implicit Time Integration", ] ,
            "transient": ["Transient Implicit Time Integration", ] , 
            "psuedo":    ["Steady-State Psuedo-Time Implicit Solver", ] ,
            "linear":    ["Linear Solver", ] , 
            "ML":        ["Trilinos ML", ] ,
            "Hypre":     ["Hypre AMG", ] , 
            "Block":     ["Block ILU", ] , 
            }


WritePreNCP(xml_output,ncp_pre_lines)

PackUnstructuredLists(KeyList, SubLists)
WriteUnstructuredLists(xml_output,KeyList,SubLists)

WritePostNCP(xml_output,ncp_post_lines)
