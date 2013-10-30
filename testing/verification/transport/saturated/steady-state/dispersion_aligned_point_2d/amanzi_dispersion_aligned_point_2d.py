import matplotlib.pyplot as plt
import numpy
#import model_dispersion_aligned_point_2d
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import prettytable 
import os, re

# load input xml file
#  -- create an ObservationXML object
def loadInputXML(filename):
    Obs_xml = ObsXML(filename)
    return Obs_xml
            
# load the data file
#  -- use above xml object to get observation filename
#  -- create an ObservationData object
#  -- load the data using this object
def loadDataFile(Obs_xml,directory):
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA(os.path.join(directory,output_file))
    Obs_data.getObservationData()
    coords = Obs_xml.getAllCoordinates()

    for obs in Obs_data.observations.itervalues():
        region = obs.region
        obs.coordinate = coords[region]
    return Obs_data


def CollectObservations(Obs_xml, Obs_data, Obs_lines):

    # Create dictionary for scatter plots
    Obs_scatter={}
    for key in Obs_lines:
        Obs_scatter[key]={}
        Obs_scatter[key]['distance']=[]
        Obs_scatter[key]['Tc99']=[]


    # Collect observations in scatter_data
    for key in Obs_lines:

        if (Obs_lines[key]['slice'][0] is 'x'):
            slice_coord=0
        elif (Obs_lines[key]['slice'][0] is 'y'):
            slice_coord=1
            
        if ( Obs_lines[key]['vary'] is 'x' ):
            vary_coord=0
        elif ( Obs_lines[key]['vary'] is 'y'):
            vary_coord=1

        for obs in Obs_data.observations.itervalues(): 

            print obs.coordinate, obs.data
            print obs.coordinate[slice_coord], Obs_lines[key]['slice'][1]
            
            if ( obs.coordinate[slice_coord] == Obs_lines[key]['slice'][1] ):
                Obs_scatter[key]['distance'].append(obs.coordinate[vary_coord])
                Obs_scatter[key]['Tc99'].append(obs.data)

    return Obs_scatter


def PlotObservations(Obs_scatter,slice_name,subtests,axes1):

    # Plot the observations
    for st in Obs_scatter:
        plot_props=subtests[st]['plot_props']
        axes1.scatter(Obs_scatter[st][slice_name]['distance'],Obs_scatter[st][slice_name]['Tc99'],c=plot_props['color'],marker=plot_props['marker'],s=25,label=plot_props['label'])

    return


def CollectAnalyticSolutions(input_file,directory):

    at123d_setup=os.path.join(directory,os.path.splitext(input_file)[0]+"_setup.out")
    at123d_soln=os.path.join(directory,os.path.splitext(input_file)[0]+"_soln.out")

    solution = {}
    solution['source']={}

    try:
        with open(at123d_setup,'r') as f_setup:
            for line in f_setup:
                if ( "NO. OF POINTS IN X-DIRECTION" in line ):
                    solution['nx']=int(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ( "NO. OF POINTS IN Y-DIRECTION" in line ):
                    solution['ny']=int(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ( "NO. OF POINTS IN Z-DIRECTION" in line ):
                    solution['nz']=int(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ("BEGIN POINT OF X-SOURCE LOCATION" in line):
                    solution['source']['x']=float(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ("END POINT OF X-SOURCE LOCATION" in line):
                    solution['source']['x']=(solution['source']['x']+float(re.sub(r'.*. ','', line).rstrip().lstrip()))/2.0
                elif ("BEGIN POINT OF Y-SOURCE LOCATION" in line):
                    solution['source']['y']=float(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ("END POINT OF Y-SOURCE LOCATION" in line):
                    solution['source']['y']=(solution['source']['y']+float(re.sub(r'.*. ','', line).rstrip().lstrip()))/2.0
                elif ("BEGIN POINT OF Z-SOURCE LOCATION" in line):
                    solution['source']['z']=float(re.sub(r'.*. ','', line).rstrip().lstrip())
                elif ("END POINT OF Z-SOURCE LOCATION" in line):
                    solution['source']['z']=(solution['source']['z']+float(re.sub(r'.*. ','', line).rstrip().lstrip()))/2.0
                else:
                    pass

    except IOError:
        raise RuntimeError("Unable to open file"+at123d_setup+", it does not exist OR permissions are incorrect")


    read_x=True
    read_y=True
    read_z=True
    read_t=True
    solution['x']=[]
    solution['y']=[]
    solution['z']=[] 
    solution['c']=[]

    try:
        with open(at123d_soln,'r') as f_soln:
            for line in f_soln:
                if ( read_x ):
                    solution['x']=solution['x']+(' '.join(line.split()).split())
                    if ( len(solution['x']) == solution['nx'] ):
                        read_x=False
                elif ( read_y ):
                    solution['y']=solution['y']+(' '.join(line.split()).split())
                    if ( len(solution['y']) == solution['ny'] ):
                        read_y=False
                elif ( read_z ):
                    solution['z']=solution['z']+(' '.join(line.split()).split())
                    if ( len(solution['z']) == solution['nz'] ):
                        read_z=False
                elif ( read_t ):
                    NL=len(line.split(" "))-1
                    solution['t']=line.split(" ")[NL]
                    read_t=False
                else:
                    solution['c']=solution['c']+(' '.join(line.split()).split())
                    
    except IOError:
        raise RuntimeError("Unable to open file"+at123d_soln+", it does not exist OR permissions are incorrect")
                    
    return solution


def PlotAnalyticSoln(solution, analytic, slice, obs_slices, axes1):

    soln = solution[slice]
    obs_slice = obs_slices[slice]

    # Set the key for the horizontal axis
    coord = obs_slice['vary']
    
    # Convert horizontal axis to float and shift to align the source
    hv = [float(i) - float(soln['source'][coord]) for i in soln[coord]]
    vv = [float(i) for i in soln['c']]
    axes1.plot(hv,vv)


if __name__ == "__main__":

    import os
    import run_amanzi, run_at123d_at

    input_filename = "amanzi_dispersion_aligned_point_2d.xml"

    # Collect slices of concentration from the observations
    #
    # Slice should be all fixed quantities, i.e., Time is fixed as well 
    # Should include observation type (or return slices sorted with observation type as a key).
    obs_slices = { 'centerline'  : {'slice': [ 'y', 0.0 ], 'vary': 'x', 'domain': [-270.0,960.0] },
                   'x=0.0'       : {'slice': [ 'x', 0.0 ], 'vary': 'y', 'domain': [-5.0,270.0] },
                   'x=424.0'     : {'slice': [ 'x', 424.0 ], 'vary': 'y', 'domain': [-5.0,270.0] },
    }

    #
    # Create emtpy dictionaries
    #
    obs_scatter={}
    obs_data={}
    obs_xml={}

    subtests = { 'amanzi_first' : 
                 { 'directory'  : 'amanzi-output-first-order',
                   'parameters' : { 'Transport Integration Algorithm': 'Explicit First-Order' },
                   'plot_props' : { 'marker':'s', 'color':'r', 'label': 'Amanzi 1st' } 
                 },
                 'amanzi_second' : 
                 { 'directory'  : 'amanzi-output-second-order',
                   'parameters' : { 'Transport Integration Algorithm': "Explicit Second-Order" },
                   'plot_props' : { 'marker':'o', 'color':'b', 'label': 'Amanzi 2nd' }
                 },
             }
   

    analytic = { 'centerline'   : 
                 { 'directory'  : 'at123d-at',
                   'input_file' : 'at123d-at_centerline.list',
                   'plot_props' : { 'color': 'blue', 'linestyle' : '-' }
                 },
                 'x=0.0'      : 
                 { 'directory' : 'at123d-at', 
                   'input_file' : 'at123d-at_slice_x=0.list', 
                   'plot_props' : { 'color': 'blue', 'linestyle' : '-' }
                 },
                 'x=424.0'      :   
                 { 'directory'  : 'at123d-at', 
                   'input_file' : 'at123d-at_slice_x=420.list',
                   'plot_props' : { 'color': 'blue', 'linestyle' : '-' }
                 },
               }


    analytic_soln = {} 

    try: 

        for st in subtests:
            print st, subtests[st]['directory']
            run_amanzi.run_amanzi(input_filename, subtests[st]['directory'], subtests[st]['parameters'])
            obs_xml[st]=loadInputXML(input_filename)
            obs_data[st]=loadDataFile(obs_xml[st],subtests[st]['directory'])

            # Collect observations to plot
            obs_scatter[st]=CollectObservations(obs_xml[st], obs_data[st], obs_slices)
            
    finally: 
        pass


    try:

        overwrite=True
        for a in analytic:
            run_at123d_at.run_at123d(analytic[a]['input_file'], analytic[a]['directory'],overwrite)
            analytic_soln[a]=CollectAnalyticSolutions(analytic[a]['input_file'],analytic[a]['directory'])

    finally:
        pass



    try:
        
        slice='centerline'
        #slice='x=0.0'
        #slice='x=424.0'

        # Plot the data:
        fig1 = plt.figure()
        axes1 = fig1.add_axes([0.15,0.15,0.80,0.80])
        axes1.set_yscale('log')
        axes1.set_xlim(obs_slices[slice]['domain'][0],obs_slices[slice]['domain'][1])
        
        # Plot centerline (y=0) 
        PlotObservations(obs_scatter,slice,subtests,axes1)
        PlotAnalyticSoln(analytic_soln,analytic,slice,obs_slices,axes1)
        
        plt.show()
        # Plot the data:
        #cmap = plotTestObservations(obs_xml,obs_data,axes1)

    finally:
        pass

