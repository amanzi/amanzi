#! /bin/sh
# -------------------------------------------------------------
# file: splitem.sh
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created November 15, 2010 by William A. Perkins
# Last Change: Mon Nov 15 09:41:29 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
# -------------------------------------------------------------

set -xue

seacas=${ACCESS-/files0/ascem/src/seacas}
nem_slice="$seacas/bin/nem_slice"
nem_spread="$seacas/bin/nem_spread"

# -------------------------------------------------------------
# splitexo
# -------------------------------------------------------------
splitexo() {
    infile=$1
    nproc=$2

    if [ ! -f $infile ]; then
        echo "cannot find input file \"$infile\"" 1>&2 
        return 2
    fi

    base=`expr "$infile" : '\(.*\)\.[^\.]*$' `
    if [ -z "$base" ]; then
        base="$infile"
    fi

    cat > nem_slice.inp <<EOF
INPUT EXODUSII FILE             = $infile
OUTPUT NEMESISI FILE            = ${base}.nemI
GRAPH TYPE                      = ELEMENTAL
MACHINE DESCRIPTION             = MESH=$nproc
DECOMPOSITION METHOD            = INERTIAL,KL, NUM_SECTS=1, CNCTD_DOM
OUTPUT VISUALIZATION FILE       = false
EOF

    "$nem_slice" -a nem_slice.inp 


    # nem_spread needs some arcane/stupid directory structure in which
    # to put the files; this needs to be made; 

    mkdir "./split1" || true

    cat > nem_spread.inp <<EOF
Input FEM file          = $infile
LB file                 = ${base}.nemI
Debug                   = 4
Parallel Disk Info      = number=1
Parallel file location  = root=./split, subdir=.
EOF
    
    "$nem_spread" 

    rm -f ${base}.nemI nem_slice.inp nem_spread.inp

    return 0
}




# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

files=" \
    hex_11x11x11_ss.exo     4 \
    hex_11x11x11_ss.exo     3 \
    hex_11x11x11_ss.exo     2 \
    hex_4x4x4_ss.exo        3 \
    hex_4x4x4_ss.exo        2 \
    twoblktet_ss.exo        4 \
    twoblktet_ss.exo        3 \
    twoblktet_ss.exo        2 \
    htc_rad_test-random.exo 4 \
    htc_rad_test-random.exo 3 \
    htc_rad_test-random.exo 2 \
"

set $files

while [ $# -ge 2 ]; do
    splitexo $1 $2
    shift 2
done
