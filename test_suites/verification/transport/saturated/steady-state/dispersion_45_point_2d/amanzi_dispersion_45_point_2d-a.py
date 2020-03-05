import matplotlib.pyplot as plt
import numpy
import prettytable 
import os, sys
import utils_dispersion_45_point_2d

#
# Amanzi Input File:
#
input_filename = "amanzi_dispersion_45_point_2d-u.xml"
    
#
# Overwrite:  True|False
#
# If Amanzi or AT123D-AT results exist, force them to be recomputed, and hence, overwritten
#
overwrite=False

[obs_slices, subtests, analytic, table_layout] = utils_dispersion_45_point_2d.SetupTests()
[obs_xml, obs_data, obs_scatter] = utils_dispersion_45_point_2d.AmanziResults(input_filename,subtests,obs_slices,overwrite)
analytic_soln = utils_dispersion_45_point_2d.AnalyticSolutions(analytic,obs_slices,overwrite)

#
# One of two plots rotated 'centerline' (y=x) and rotated 'x=0.0' (y=-x)
# 
slice='centerline'

# Plot the data:
fig1 = plt.figure(figsize=(10,6))
axes1 = fig1.add_axes([0.15,0.15,0.80,0.80])
axes1.set_yscale('log')
axes1.set_xlim(obs_slices[slice]['domain'][0],obs_slices[slice]['domain'][1])
axes1.set_ylim(obs_slices[slice]['range'][0],obs_slices[slice]['range'][1])

# Plot rotated centerline (x=y) 
utils_dispersion_45_point_2d.PlotObservations(obs_scatter,slice,subtests,axes1)
utils_dispersion_45_point_2d.PlotAnalyticSoln(analytic_soln,analytic,slice,obs_slices,axes1)

axes1.legend(loc='lower right',fontsize=14)
axes1.set_xlabel('Position along the Plume Centerline, x[m]',fontsize=14)
axes1.set_ylabel('Concentration [$kg/m^3$]',fontsize=14)
axes1.set_title('Concentration along $y=0$ at $t=1440$ days',fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)

# Create/Write the table:
utils_dispersion_45_point_2d.MakeTableCols(table_layout,slice,obs_scatter,subtests,analytic_soln,analytic)


