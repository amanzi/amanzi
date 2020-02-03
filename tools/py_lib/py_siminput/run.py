#!/usr/bin/env python
#
# Copyright 2010 Los Alamos National Security, LLC
# Written by Robert B. Lowrie (CCS-2)
#
'''
Functions for running a simulation.
'''

_usage = '''
Usage: %s [options] input_file

where,

  input_file :  The input file.  If the extension on the file name is
  ".xml", then the executable reads this file directly.  Otherwise,
  the file is assumed to be a Python front-end input file, which is
  then converted to an xml file, after which the executable is run.

Options [default]:

  -d, --doc
          Dump input variable documentation to stdout and exit.
          Unavailable if input_file is an xml file. [false]

  -h, --help
          Print this message and exit. [false]

  -n, --procs
          Run on this many processors. [1]

  -p, --python
          Dump input variables, in python format, to stdout and exit.
          Some editing may be required in order to use this dump as an
          input file. Unavailable if input_file is an xml file. [false]
          
  -x, --xml
          Convert input_file to xml and exit. [false]

  Options -d and -x may be combined; if either of these 2 options are
  specified, the executable is not run.

  The -d and -p options do not document all possible input
  combinations.  For example, a linear solver type could allow CG or
  GMRES.  If the solver type is set to CG in input_file, only CG
  related input variables will be documented with the -d option.
  There is no way from the -d or -p options to know that GMRES is
  available as a solver option.

Example:

  %s --procs=4 input.py
'''

import subprocess, getopt, sys, os
from py_siminput.input_interface import Dump_XML, Dump_Doc, Dump_Python

# global variables

prefix_default = 'time %s'
script_name    = os.path.basename(sys.argv[0])

build_type = 'host'
        
def Usage(mesg):
    '''
    Prints documentation and exits.
    '''
    print(mesg)
    print(_usage % (script_name, script_name))
    sys.exit(0)

def system(command, debug):
    '''
    Runs command on the underlying operating system.
    '''
    if debug:
        print("      " + command) # for debugging
    else:
        r = os.system(command)
        if r != 0:  # alternatively, "if not os.WIFEXITED(r)"
            print("Error running: " + command)
            sys.exit(1)

def append_to_file(file):
    '''
    Can be used to append items to a file, in a system() call.
    '''
    return " 1>> " + file + " 2>&1"

def mpi_only_run_string(num_pes,exe,inputfile):
    '''
    This function may be used to define a default run string for
    mpi runs, dependent on the platform.
    '''
    machine = subprocess.getoutput('uname')

    if machine == "AIX":
        # IBM AIX (with poe)
        s = "poe %s -procs %d -rmpool 2 -retry 1" % (exe,num_pes,inputfile)
    elif machine == "OSF1":
        # COMPAQ
        s= "prun -n %d %s %s" % (num_pes,exe,inputfile)
    else:
        # Set the default
        s = "mpirun -np %d %s %s" % (num_pes,exe,inputfile)

    return s

def mcmpi_run_string(num_pes,exe1,exe2,inputfile):
    '''
    This function defines a run string for mcmpi runs. only has default.
    platforms should be added a la mpi_run_string as needed.
    '''
    s = "mpirun -np %d %s %s : -np %d %s" % (num_pes,exe2,inputfile,
                                             num_pes,exe1)
    return s

def mpidacs_run_string(num_pes,exe1,exe2,inputfile):
    '''
    This function defines a run string for mpidacs runs. only has default.
    platforms should be added a la mpi_run_string as needed.
    '''
    s = "mpiexec --mca mpi_paffinity_alone 1 -np %d %s %s %s" % (num_pes,exe2,
                                                                 exe1,
                                                                 inputfile)

    return s

def run_string(num_pes,exe1,exe2,inputfile):
    if build_type == 'host' and num_pes > 1:
        return mpi_only_run_string(num_pes,exe1,inputfile)
    elif build_type == 'host' and num_pes == 1:
        return exe1 + " " + inputfile
    elif build_type == 'mcmpi':
        return mcmpi_run_string(num_pes,exe1,exe2,inputfile)
    elif build_type == 'ppe':
        return mpidacs_run_string(num_pes,exe1,exe2,inputfile)

def unit_test_args(exe1,exe2,input,args):
    '''
    This is a utility function that parses the command-line arguments
    from dracos launchtest and returns the pair (number_of_procs,
    run_string).

    Typically, args is sys.argv.
    '''
    # serial default
    rs = "%s"
    num_pes = 1

    # check for --procs, i.e., a parallel run.
    if len(args) > 1:
        try:
            optlist, args = getopt.getopt(args[1:], 'n:',
                                          ['procs='])
        except getopt.error as val:
            print('ERROR ' + val.msg)
            sys.exit(1)

        for o, a in optlist:
            if o in ('-n', '--procs'):
                num_pes = int(a)

    rs = run_string(num_pes,exe1,exe2,input)

    # all done
    return (num_pes, rs)

def get_root(input):
    input_dict = {}
    exec(compile(open(input).read(), input, 'exec'), input_dict)
    
    if 'root' not in input_dict:
        Usage("No 'root' variable defined in input file %s" % (input))

    return input_dict['root']

def run_root(root, input_xml, command, debug = 0):
    '''
    This is a low-level function, used by the run() function.  Also,
    it may be useful in more sophisticated scripts.
    
    root      = Root of input tree.
    input_xml = File name where xml is dumped.
    command   = Command string.
    debug     = Print the action to be taken, but do not actually do it.
    '''
    Dump_XML(root, input_xml)
    system(command, debug)

def run(exe1,exe2):
    '''
    Runs exes.  See _usage string for more information on the actions of
    this function.
    '''
    prefix      = prefix_default
    debug       = 0
    dump_doc    = 0
    dump_python = 0
    dump_xml    = 0
    input       = 'input.xml'
    num_pes     = 1

    # Parse command line options

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'dhn:px',
                                      ['doc',
                                       'help',
                                       'procs=',
                                       'python',
                                       'xml'])
    except getopt.error as val:
        Usage('ERROR: ' + val.msg)

    for o, a in optlist:
        if o in ('-d', '--doc'):
            dump_doc = 1
        elif o in ('-h', '--help'):
            Usage("");
        elif o in ('-p', '--python'):
            dump_python = 1
        elif o in ('-n', '--procs'):
            num_pes = int(a)
        elif o in ('-x', '--xml'):
            dump_xml = 1

    if len(args) > 1:
        Usage("Too many arguments.")

    if len(args) < 1:
        Usage("Must supply input filename.")

    input = args[0]

    # Make sure input exists

    if not os.path.isfile(input):
        Usage("Unable to find input file: %s" % (input))

    # Create the xml filename

    input_base = os.path.basename(input)

    input_is_xml = False

    if input_base[-4:] == ".xml":
        input_is_xml = True
        input_xml = input_base
    elif len(input_base) < 3 or input_base[-3:] != ".py":
        input_xml = input_base + ".xml"
    else:
        input_xml = input_base[0:-3] + ".xml"

    # Filter input file through python.  The input file must generate
    # an input interface named 'root'.

    if input_is_xml:
        if dump_xml:
            Usage('--xml specified when input is already xml.')
        if dump_python:
            Usage('--python unavailable when input is xml.')
        if dump_doc:
            Usage('--doc unavailable when input is xml.')
        command = run_string(num_pes,exe1,exe2,input_xml)
        print('Running: %s ...' % (command))
        system(command,debug)
    else:
        print('Generating %s from %s ...' % (input_xml, input))
        root = get_root(input)

        run_sim = 1 # if true, run simulation.

        if dump_xml:
            run_sim = 0

        if dump_doc:
            # Dump documentation
            Dump_Doc(root)
            run_sim = 0

        if dump_python:
            # Dump python
            d = Dump_Python(1)
            d.dump(root, 'root')
            run_sim = 0
        
        if run_sim:
            # Run the solver
            command = run_string(num_pes,exe1,exe2,input_xml)
            print('Running: %s ...' % (command))
            run_root(root, input_xml, command, debug)
        elif dump_xml:
            # just dump the xml
            Dump_XML(root, input_xml)
