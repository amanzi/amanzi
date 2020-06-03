import matplotlib.pyplot as plt
import numpy
import model_linear_flux_head_1d
from amanzi_xml.observations.ObservationXMLv2 import ObservationXMLv2 as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import prettytable 
import os 

# load input xml file
#  -- create an ObservationXML object
def loadInputXML(filename):
    Obs_xml = ObsXML(filename)
    return Obs_xml
            
# load the data file
#  -- use above xml object to get observation filename
#  -- create an ObservationData object
#  -- load the data using this object

def loadDataFile(Obs_xml):
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA("amanzi-output/"+output_file)
    Obs_data.getObservationData()
    coords = Obs_xml.getAllCoordinates()

    for obs in Obs_data.observations.values():
        region = obs.region
        obs.coordinate = coords[region]
    return Obs_data

def plotTestObservations(Obs_xml, Obs_data, axes1):

    # === SPECIAL CODE ==== for linear flow problems
    # Collect the z-values from observations
    z_vals = [coord[2] for coord in Obs_xml.coordinates.values()]
    z_vals.sort()
    z_vals = set(z_vals)
    colors = ['b','g','r']
    cmap = dict((zval,color) for (zval,color) in zip(z_vals, colors))

    # Create dictionary for scatter plots
    scatter_data={}
    for key in list(cmap.keys()):
        scatter_data[key]={}
        scatter_data[key]['x']=[]
        scatter_data[key]['pressure']=[]

    # Collect observations in scatter_data
    for obs in Obs_data.observations.values(): 
        scatter_data[obs.coordinate[2]]['x'].append(obs.coordinate[0])
        scatter_data[obs.coordinate[2]]['pressure'].append(obs.data)

    # Plot the observations
    for key in list(cmap.keys()):
        axes1.scatter(scatter_data[key]['x'],scatter_data[key]['pressure'],c=cmap[key],marker='s',s=25,label='Amanzi')

    # Set labels and title
    axes1.set_xlabel('x-coordinate [m]')
    axes1.set_ylabel('Pressure [Pa]')
    axes1.set_title('Aqueous Pressure vs Distance')
    
    return cmap

def plotTestModel(filename, cmap, axes1, Obs_xml, Obs_data):

    # Instantiate the analytic solution
    mymodel = model_linear_flux_head_1d.createFromXML(filename)
    table_values = []

    # Create a set of points to plot the solution
    x = numpy.linspace(mymodel.x_0,mymodel.x_1,11)
    coords = numpy.zeros((11,2))
    coords[:,0] = x

    # Plot a line for each z-coordinate in the observations
    for (z_val, color) in cmap.items():
        coords[:,1] = z_val
        pressure = mymodel.pressure(coords)
        axes1.plot(x,pressure,color,label='$z = %0.2f $'%z_val)
        axes1.legend(loc="upper right" , fancybox = True , shadow = True)
          
def MakeTable(Obs_data,Obs_xml,filename):

    pressure_amanzi = []
    coordinates = []
    mymodel = model_linear_flux_head_1d.createFromXML(filename)

    for obs in Obs_data.observations.values():
        coordinates.append([obs.coordinate[0], obs.coordinate[2]])
        pressure_amanzi.append(str(obs.data).rstrip(']').lstrip('['))

    pressure_analytic = list(mymodel.pressure(numpy.array(coordinates)))

    x = prettytable.PrettyTable(["x [m]", "z [m]", "Analytic [Pa]","Amanzi [Pa]"])
    x.padding_width = 1
    x.hrules = 1
    for coords, p_analytic, p_amanzi in zip(coordinates,pressure_analytic,pressure_amanzi):
        x.add_row([coords[0],coords[1],"%.4f" % float(p_analytic),"%.4f" % float(p_amanzi)])
        
    if os.path.exists("table_values.txt"):
        os.remove("table_values.txt")

    table_file = open("table_values.txt", "w+")
    table_file.write('.. tabularcolumns:: ' + '|R|C|C|C|' + '\n\n')
    table_file.write(x.get_string(sortby="x [m]"))
    table_file.write('\n')
    table_file.close()
        
if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file = "amanzi_linear_flux_head_1d-u.xml"
    run_dir = "amanzi-output"
    try: 
        run_amanzi_standard.run_amanzi(input_file, 2, [input_file], run_dir)
        obs_xml=loadInputXML(input_file)
        obs_data=loadDataFile(obs_xml)

        fig1 = plt.figure()
        axes1=fig1.add_axes([.15,.15,.80,.80])
       
        MakeTable(obs_data,obs_xml,input_file)

        cmap = plotTestObservations(obs_xml,obs_data,axes1)
        plotTestModel(input_file,cmap,axes1,obs_xml,obs_data)
        # plt.show()

    finally:
        pass 


