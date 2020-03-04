import matplotlib.pyplot as plt
import numpy
import prettytable 
import os, re
import utils_dispersion_aligned_point_2d

#
# Amanzi Input File:
#
input_filename = "amanzi_dispersion_aligned_point_2d-u.xml"
    
#
# Overwrite:  True|False
#
# If Amanzi or AT123D-AT results exist, force them to be recomputed, and hence, overwritten
#
overwrite=False

[obs_slices, subtests, analytic, table_layout] = utils_dispersion_aligned_point_2d.SetupTests()
[obs_xml, obs_data, obs_scatter] = utils_dispersion_aligned_point_2d.AmanziResults(input_filename,subtests,obs_slices,overwrite)
analytic_soln = utils_dispersion_aligned_point_2d.AnalyticSolutions(analytic,obs_slices,overwrite)

#
# One of three plots 'centerline', 'x=0.0', and 'x=424.0'
# 
slice='x=424.0'

# Plot the data:
fig1 = plt.figure(figsize=(10,6))
axes1 = fig1.add_axes([0.15,0.15,0.80,0.80])
axes1.set_yscale('log')
axes1.set_xlim(obs_slices[slice]['domain'][0],obs_slices[slice]['domain'][1])

# Plot centerline (y=0) 
utils_dispersion_aligned_point_2d.PlotObservations(obs_scatter,slice,subtests,axes1)
utils_dispersion_aligned_point_2d.PlotAnalyticSoln(analytic_soln,analytic,slice,axes1)

#
axes1.legend(loc='lower left',fontsize=14)
axes1.set_xlabel('Transverse Distance from Plume Centerline, y[m]',fontsize=14)
axes1.set_ylabel('Concentration [$kg/m^3$]',fontsize=14)
axes1.set_title('Concentration at $x=420$ and $t=1440$ days',fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)

#plt.show()
utils_dispersion_aligned_point_2d.MakeTableCols(table_layout,slice,obs_scatter,subtests,analytic_soln,analytic)


