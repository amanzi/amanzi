import sys,os
try:
    from atsxml import *
except:
    try:
        sys.path.append(os.sep.join(['..','..','..']))
        from atsxml import *
    except:
        print "Unable to locate atsxml module"
        sys.exit()

x3d = get_root('3d_run_from_Ethan.xml')
x2d = get_root('geoinv7.xml')
# Replace mesh regions
mapping = {'computational domain':'computational domain',
           'surface domain':'surface domain',
           'top part 1':'computational domain peat',
           'top part 2':'computational domain upper mineral',
           'surface':'surface',
           'bottom face':'bottom face',
           'south edge':'front,back'
          }
# New functional replace regions
replace_regions(x3d,x2d,mapping=mapping)

# Replace BC files and names
find_value(x3d,'Snow').set('value','/Ps')
find_value(x3d,'Rain').set('value','/Pr')
find_value(x3d,'/Windspeed').set('value','/Us')
find_value(x3d,'/Relative_humidity').set('value','/RH')
find_value(x3d,'/Temperature').set('value','/Ta')
find_value(x3d,'/Shortwave').set('value','/QswIn')
replace_by_name(x3d,'../../Data/Barrow-Anna-2010-2013.h5','/scratch/tundra/dharp/arctic/geophysics/run_newmesh02/10yr-WinterStartRAW-heavystart.h5')

# Remove longwave radiation BC
remove(x3d,find_name(x3d,'incoming_longwave_radiation'))

# Replace initialization with restart file
e_ic2d = create_parent_map(x2d)[findall_name(x2d,'restart file')[0]]
for e in findall_name(x3d,'initialize from 1D column'):
    replace_elem(x3d,create_parent_map(x3d)[e],e_ic2d)    

# Change bottom temp bc
replace_by_value(x3d,'267.15','263.55')

# Remove wallclock duration
remove(x3d,find_name(x3d,'wallclock duration [hrs]'))

write_xml(x3d,'geoinv7_serial.xml')
