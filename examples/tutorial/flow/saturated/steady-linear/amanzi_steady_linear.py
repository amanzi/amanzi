import matplotlib.pyplot as plt
import numpy
import model_steady_linear as model
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search

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

    for obs in Obs_data.observations.itervalues():
        region = obs.region
        obs.coordinate = coords[region]
    return Obs_data

def plotExampleObservations(Obs_xml, Obs_data, axes1):
    #=== SPECIAL CODE ==== ONLY EXAMPLE 1 
    # Collect the used z-values
    z_vals = [coord[2] for coord in Obs_xml.coordinates.itervalues()]
    z_vals.sort()
    z_vals = set(z_vals)
    colors = ['r','b','g']
    cmap = dict((zval,color) for (zval,color) in zip(z_vals, colors))
        
    # -- grab the right coordinate and value
    # -- call ObservationPlotter to plot with format="bo"
    coords_list=[]
    pressure_value = []
    for obs in Obs_data.observations.itervalues(): 
        pressure_value.append(obs.data)
        coords_list.append(obs.coordinate)
        
    axes1.scatter([coords_list[0][0],coords_list[2][0],coords_list[3][0]],[pressure_value[0],pressure_value[2],pressure_value[3]],c=cmap[coords_list[0][2]],marker='s',s=50, label = '$Amanzi$')
    axes1.scatter([coords_list[1][0]],[pressure_value[1]],c = cmap[coords_list[1][2]],marker = 's',s = 50,label = '$Amanzi$')
    axes1.scatter([coords_list[4][0]],[pressure_value[4]],c = cmap[coords_list[4][2]],marker = 's',s = 50,label = '$Amanzi$')

    axes1.set_xlabel('Distance in x-direction [meter]')
    axes1.set_ylabel('Pressure [Pa]')
    axes1.set_title('Aqueous Pressure vs Distance')
    
    return cmap

def plotExampleModel(filename, cmap, axes1,Obs_xml, Obs_data):
    mymodel = model.createFromXML(filename)
    table_values = []

    x = numpy.linspace(mymodel.x0,mymodel.x1,11)
    coords = numpy.zeros((11,2))
    coords[:,0] = x

    for (z_val, color) in cmap.iteritems():
        coords[:,1] = z_val
        
        pres = mymodel.run(coords)
        axes1.plot(x,pres,color,label='$z = %0.2f $'%z_val)
        axes1.legend(loc="upper right" , fancybox = True , shadow = True)
          
    all_coordinates = Obs_xml.coordinates
    coordinates = []
    pressure_value = [] 
   
    for obs in Obs_data.observations.itervalues():
        coordinates.append([obs.coordinate[0], obs.coordinate[2]])
        pressure_value.append(obs.data)
    
    coordinates = numpy.array(coordinates)
    pres_analytic = mymodel.run(coordinates)
    pressure_analytic = list(pres_analytic)
    
    for coord,press_amanzi, press_analytic in zip(coordinates,pressure_value,pressure_analytic):
        table_values.append([coord[0],coord[1],press_amanzi,press_analytic])
        
    col_labels = ['x [m]', 'z [m]' ,'Pressure(Amanzi) [Pa]','Pressure (analytic)[Pa]']
    ax = plt.subplot(111, frame_on = False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    the_table = plt.table(cellText = table_values ,  colWidths = [.1 , .1 , .4 , .5], colLabels = col_labels,  cellLoc = 'center', colLoc = 'center', loc = 'best') 
    properties = the_table.properties()
    cells = properties['child_artists']
    for cell in cells: 
        cell.set_height(.06)
    
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(14)

if __name__ == "__main__":

    import os
    import run_amanzi

    input_filename = "amanzi_steady_linear.xml"
    try: 
        run_amanzi.run_amanzi("../"+input_filename)
        obs_xml=loadInputXML(input_filename)
        obs_data=loadDataFile(obs_xml)

        fig1= plt.figure()
        fig2 = plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        cmap = plotExampleObservations(obs_xml,obs_data, axes1)
        plotExampleModel(input_filename, cmap, axes1,obs_xml, obs_data)
        
        plt.show()

    finally:
        pass 


