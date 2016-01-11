import h5py
import math
import numpy
import struct
from matplotlib import pyplot as plt

fname_area_weights = 'firetech_weights_for_prism_v3_mesh_384.dat'

nx_sq = 1500
ny_sq = 1500

indata_names = ['aspen_drymass', 'conifer_drymass', 'grass_drymass', 'litter_drymass', 'pipo_drymass', 'spruce_drymass']
outdata_fname = 'drymass_mesh_v3_384.h5'
val_matches_fname = 'drymass_value_matches.txt'

in_data = []
for imass in range(len(indata_names)):
  fname = indata_names[imass] + '.dat'
  fdata = open(fname, 'rb').read()
  struct.unpack("@i", fdata[:4])
  cur_data = struct.unpack("f" * ((len(fdata) - 8) // 4), fdata[4:-4])
  cur_data = numpy.reshape(cur_data, (nx_sq, ny_sq), order='F')
  in_data.append(cur_data)
  struct.unpack("@i", fdata[-4:])
  #  plt.figure()
#  plt.imshow(cur_data, origin='lower left')

#plt.show()

#for imass in range(len(indata_names)):
#  print in_data[imass][0][139]

matches_outf = open(val_matches_fname, 'w')
nmatches = 0
for icol in range(nx_sq):
  for irow in range(ny_sq):
    iunchecked = range(0, len(indata_names), 1)
    cell_same_data = [[] for x in range(len(indata_names) - 1)]
    any_matches = False
    for imass in range(len(indata_names) - 1):
      if len(iunchecked) == 0:
        break
      
      if imass != iunchecked[0]:
        continue

      if in_data[imass][icol][irow] == 0.0:
        del iunchecked[0]
        continue
          
      cmp_val = numpy.float64(in_data[imass][icol][irow])
      del iunchecked[0]
      nunchecked = len(iunchecked)
      nmatched = 0
      for icmpmass in range(nunchecked):
        iumass = icmpmass - nmatched
        icurmass = iunchecked[iumass]
        if abs(numpy.float64(in_data[icurmass][icol][irow]) - cmp_val) < 1.0e-15:
          cell_same_data[imass].append(icurmass)
          del iunchecked[iumass]
          nmatched += 1
      if nmatched > 0:
        any_matches = True
    
    if any_matches:
      matches_outf.write('Square [' + repr(icol) + ', ' + repr(irow) + ']:\t')
      for imass in range(len(indata_names) - 1):
        if len(cell_same_data[imass][:]) > 0:
          matches_outf.write('{' + indata_names[imass])
          for imatchedmass in range(len(cell_same_data[imass][:])):
            matches_outf.write(', ' + indata_names[cell_same_data[imass][imatchedmass]])
          matches_outf.write('}\t')
      matches_outf.write('\n')
      nmatches += 1

matches_outf.close()
print('Number of cells with the same nonzero values for different classes: ' + repr(nmatches))


fid = open(fname_area_weights, 'r')
cur_str = fid.readline()
ntris = int(cur_str)
cur_str = fid.readline()
str_values = cur_str.split()
sqmesh_min_coord = [float(str_values[0]), float(str_values[1])]
sqmesh_step = float(str_values[2])

out_data = [[0.0 for x in range(ntris)] for x in range(len(indata_names))]
for itri in range(ntris):
  cur_str = fid.readline()
  str_values = cur_str.split()
  if itri != int(str_values[0]):
    raise RuntimeError("Expected index of a triangle is not on the beginning of the line!")
  if (len(str_values) - 1) % 3 == 0:
    nweights = (len(str_values) - 1) // 3
  else:
    raise RuntimeError("Unknown weights file format!")

  for iw in range(nweights):
    irow = int(str_values[iw*3 + 1])
    icol = int(str_values[iw*3 + 2])
    wsq = float(str_values[iw*3 + 3])
    for imass in range(len(indata_names)):
      out_data[imass][itri] += wsq*numpy.float64(in_data[imass][icol][irow])

fid.close

h5f = h5py.File(outdata_fname, 'w')
for imass in range(len(indata_names)):
  h5f.create_dataset(indata_names[imass], data=out_data[imass][:])

h5f.close()
