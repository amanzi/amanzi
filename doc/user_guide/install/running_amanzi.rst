============================================================
Running Amanzi 
============================================================

Introduction
++++++++++++

After building Amanzi, the executable ``amanzi`` is installed in ``<prefix>/bin/``.  Amanzi can be run in serial or parallel.  The input file is specified using the ``--xml_file=`` commandline option.  If the input file uses a current XML schema, the schema must be specified using the option ``--xml_schema=``.  The current XMl schema is located in the Amanzi source at ``amanzi/doc/input_spec/schema``.  Input files using previous Teuchos ParameterLists may skip this option. Other commandline options can be listed using the option ``--help`` and are listed below.  ::

   Usage: amanzi [options]
     options:
       --help                          Prints this help message
       --pause-for-debugging           Pauses for user input to allow attaching a debugger
       --echo-command-line             Echo the command-line but continue as normal
       --xml_file              string  XML options file
                                      (default: --xml_file="")
       --xml_schema            string  XML Schema File
                                       (default: --xml_schema="")
       --print_version         bool    Print version number and exit.
       --no_print_version              (default: --no_print_version)
       --print_tplversions     bool    Print version numbers of third party libraries and exit.
       --no_print_tplversions          (default: --no_print_tplversions)
       --print_all             bool    Print all pre-run information.
       --no_print_all                  (default: --no_print_all)
       --print_paths           bool    Print paths of the xml input file and the xml schema file.
       --no_print_paths                (default: --no_print_paths)

