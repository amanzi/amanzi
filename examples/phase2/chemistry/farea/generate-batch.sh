#!/bin/bash

#
# Script to copy pflotran processed constraints into batch_chem cfg files
#

PFLOTRAN=pflotran
PFLOTRAN_FILE=ascem-2012-1d-farea-full.in

PFLOTRAN_BACKGROUND=pflotran-background.txt
BATCHCHEM_BACKGROUND=ascem-2012-farea-background.cfg

PFLOTRAN_SEEPAGE=pflotran-seepage.txt
BATCHCHEM_SEEPAGE=ascem-2012-farea-seepage.cfg

PFLOTRAN_BGD=pflotran.bgd
AMANZI_BGD=ascem-2012-farea.bgd

#
# run pflotran to process the constrains
#

${PFLOTRAN} -pflotranin ${PFLOTRAN_FILE}

#
# convert the background constraint to cfg
#

cat > ${BATCHCHEM_BACKGROUND} <<EOF
[simulation parameters]
description = "Savanna River F-Basin, background condition, 2012 ASCEM demo; based on description and toughreact input file from Sergio Bea. PFloTran preprocessed basis and constraints"

verbosity = verbose
comparison_model = pflotran

database_type = simple
database_file = ascem-2012-farea.bgd
activity_model = debye-huckel
porosity = 0.5
saturation = 1.0
volume = 1.0
# dt = 10 days
delta_time = 864000.0
num_time_steps = 100
output_interval = 1

# all component values must be in the same order as the database file

EOF

cat ${PFLOTRAN_BACKGROUND} >> ${BATCHCHEM_BACKGROUND}

#
# copy the seepage constraint to cfg
#

cat > ${BATCHCHEM_SEEPAGE} <<EOF
[simulation parameters]
description = "Savanna River F-Basin, seepage condition, 2012 ASCEM demo; based on description and toughreact input file from Sergio Bea. PFloTran preprocessed basis and constraints"

verbosity = verbose
comparison_model = pflotran

database_type = simple
database_file = ascem-2012-farea.bgd
activity_model = debye-huckel
porosity = 0.5
saturation = 1.0
volume = 1.0
# dt = 10 days
delta_time = 864000.0
num_time_steps = 100
output_interval = 1

# all component values must be in the same order as the database file

EOF

cat ${PFLOTRAN_SEEPAGE} >> ${BATCHCHEM_SEEPAGE}

#
# rename the processed database file
#

mv ${PFLOTRAN_BGD} ${AMANZI_BGD}

#
# we have to add a mineral name to the ion exchange site, even though it isn't used...
#
perl -w -i -p -e "s/X- ; -1.0 ; /X- ; -1.0 ; Bulk/" ${AMANZI_BGD}

#
# Can't get pflotran to write a reaction rate constant when a prefactor is used...?
#
perl -w -i -p -e "s/-Infinity/-12.967/" ${AMANZI_BGD}

#
# change the mineral surface areas manually...
#
