#! /bin/sh
# -------------------------------------------------------------
# file: splitexo.sh
#
# This script is used to split ExodusII files for use with Amanzi when
# the Mesh_maps_stk class is used.  I tries to encapuslate the use of
# SEACAS utilities nem_slice and nem_spread and avoid their
# idosyncrancies.  You need to have SEACAS installed in order 
#
# Usage Example: 
#
# shprompt> setenv ACCESS /the/path/to/seacas
# shprompt> ls fbasin_unstr_400_V02.*
# fbasin_unstr_400_V02.exo
# shprompt> sh splitexo.sh fbasin_unstr_400_V02.exo 4
# shprompt> ls fbasin_unstr_400_V02.*
# fbasin_unstr_400_V02.exo
# fbasin_unstr_400_V02.par.4.0
# fbasin_unstr_400_V02.par.4.1
# fbasin_unstr_400_V02.par.4.2
# fbasin_unstr_400_V02.par.4.3
#
# A temporary directory is required. The default is /tmp. If another
# location is desired, set the TMPDIR environment variable.
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created November 15, 2010 by William A. Perkins
# Last Change: Wed Dec 15 08:46:58 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
# -------------------------------------------------------------

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------

# This script requires the SEACAS utilities nem_slice and
# nem_spread. SEACAS source can be obtained from
# http://sourceforge.net/projects/seacas or built through
# the SuperBuild TPL build. Once installed, set the environment
# variable ACCESS to the path where SEACAS was installed.

seacas=${ACCESS-/files0/ascem/src/seacas}
nem_slice="$seacas/bin/nem_slice"
nem_spread="$seacas/bin/nem_spread"

if [ ! -x "$nem_slice" ]; then
    echo "$0: error: Unable to find required utility $nem_slice" 1>&2 
    echo "set $ACCESS environment variable to SEACAS installation location" 1>&2
    exit 2
fi
if [ ! -x "$nem_spread" ]; then
    echo "$0: error: Unable to find required utility $nem_spread" 1>&2 
    echo "set $ACCESS environment variable to SEACAS installation location" 1>&2
    exit 2
fi

# -------------------------------------------------------------
# cleanup
# -------------------------------------------------------------
cleanup() {
    rm -f nem_slice.inp nem_spread.inp
    if [ -n ${base} ]; then
        rm -f "${base}.nemI"
    fi
    if [ -n ${tmpdir} ]; then
        rm -rf "$tmpdir"
    fi
}

# try to clean up after exit
trap cleanup 0 1 2 15

# -------------------------------------------------------------
# handle command line
# -------------------------------------------------------------
progname=`basename $0`
usage="Usage: $progname file.exo numpart"
if [ $# -lt 2 ]; then
    echo "$usage" 1>&2
    exit 2
fi

infile=$1
nproc=$2

if [ ! -f $infile ]; then
    echo "$progname: error: cannot find input file \"$infile\"" 1>&2 
    echo "$usage" 1>&2
    exit 2
fi

if ! (expr "$nproc" + 0 2>&1) > /dev/null; then
    echo "$progname: error: number of partitions is not numeric (got $nproc)" 1>&2 
    echo "$usage" 1>&2
    exit 2
fi

if [ "$nproc" -le 1 ]; then
    echo "$progname: error: number of partitions must be more than 1 (got $nproc)" 1>&2 
    echo "$usage" 1>&2
    exit 2
fi

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

base=`expr "$infile" : '\(.*\)\.[^\.]*$' `
if [ -z "$base" ]; then
    base="$infile"
fi

# make a temporary place to put the output; nem_spread needs some
# arcane/stupid directory structure in which to put the files; this
# needs to be made;

thetmpdir="${TMPDIR-/tmp}"

if [ ! -w "$thetmpdir" ]; then
    echo "$0: error: \"$thetmpdir\": unable to write to temporary directory" 1>&2
    exit 3
fi

tmpdir="$thetmpdir/start1"

if ! mkdir "$tmpdir"; then
    echo "$0: error: \"$tmpdir\": unable to create temporary directory" 1>&2
    exit 3
fi

# for nem_slice
# The "INPUT EXODUSII FILE" line does not seem to be doing much
# We will leave it in there but specify the input file on the command line

outpath="$thetmpdir/start"

cat > nem_slice.inp <<EOF
INPUT EXODUSII FILE             = $infile
OUTPUT NEMESISI FILE            = ${base}.nemI
GRAPH TYPE                      = ELEMENTAL
MACHINE DESCRIPTION             = MESH=$nproc
DECOMPOSITION METHOD            = INERTIAL,KL, NUM_SECTS=1, CNCTD_DOM
OUTPUT VISUALIZATION FILE       = false
EOF

"$nem_slice" -a nem_slice.inp $infile


cat > nem_spread.inp <<EOF
Input FEM file          = $infile
LB file                 = ${base}.nemI
Debug                   = 4
Parallel Disk Info      = number=1
Parallel file location  = root=${outpath}, subdir=.
EOF

"$nem_spread" 

cp "$tmpdir/${base}"* .

