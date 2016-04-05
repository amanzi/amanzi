import sys,os
try:
    from atsxml import ATSXML
except:
    try:
        sys.path.append(os.sep.join(['..','..','..']))
        from atsxml import ATSXML
    except:
        print "Unable to locate atsxml module"
        sys.exit()

x3d = ATSXML('3d_run_from_Ethan.xml')
x2d = ATSXML('geoinv7.xml')
# Replace mesh regions
mapping = {'computational domain':'computational domain',
           'surface domain':'surface domain',
           'top part 1':'computational domain peat',
           'top part 2':'computational domain upper mineral',
           'surface':'surface',
           'bottom face':'bottom face',
           'south edge':'front,back'
          }
x3d.replace_regions(x2d,mapping=mapping)

# Replace BC files and names
x3d.find_value('Snow').set('value','/Ps')
x3d.find_value('Rain').set('value','/Pr')
x3d.find_value('/Windspeed').set('value','/Us')
x3d.find_value('/Relative_humidity').set('value','/RH')
x3d.find_value('/Temperature').set('value','/Ta')
x3d.find_value('/Shortwave').set('value','/QswIn')
x3d.replace_file('/scratch/tundra/dharp/arctic/geophysics/run_newmesh02/10yr-WinterStartRAW-heavystart.h5')

#x3d.replace_by_value('/scratch/tundra/dharp/arctic/geophysics/mesh/from_lucia/arctic_siteb_2d.exo','/scratch/tundra/dharp/arctic/geophysics/mesh/from_lucia/arctic_siteb_2d.par')

# Remove longwave radiation BC
x3d.remove(x3d.find_name('incoming_longwave_radiation'))

# Replace initialization with restart file
e_ic2d = x2d.parent_map[x2d.find_name('restart file')]
for e in x3d.findall_name('initialize from 1D column'):
    x3d.replace_elem(x3d.parent_map[e],e_ic2d)    

#x3d.replace_restart_file('/scratch/tundra/dharp/arctic/geophysics/run_newmesh01/checkpoint10644.h5')

# Change bottom temp bc
x3d.replace_by_value('267.15','263.55')

# Remove wallclock duration
x3d.remove(x3d.find_name('wallclock duration [hrs]'))

x3d.write('geoinv7_serial.xml')
