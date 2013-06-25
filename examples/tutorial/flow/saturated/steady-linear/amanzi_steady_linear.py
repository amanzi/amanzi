from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search

from pylab import *


# load input xml file
#  -- create an ObservationXML object
def loadInputXML(filename):
    Obs_xml = ObsXML(filename)
    Obs_xml.getAllCoordinates()
    return Obs_xml
  

# load the data file
#  -- use above xml object to get observation filename
#  -- create an ObservationData object
#  -- load the data using this object

def loadDataFile(Obs_xml):
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA(output_file)
    Obs_data.getObservationData()
    return Obs_data

def plotExample1(Obs_xml, Obs_data):
#=== SPECIAL CODE ==== ONLY EXAMPLE 1 
#-- create a matplotlib figure
    #fig = plt.figure()

#-- create an axis
   # ax = fig.add_subplot(111)

#loop over all lines
# -- grab the right coordinate and value
# -- call ObservationPlotter to plot with format="bo"
    i = 1
    x = []
    z = []
    while i < len(Obs_xml.coordinates):
        x.append(Obs_xml.coordinates[i][0])
        z.append(Obs_xml.coordinates[i][2])
        i += 1
   
       

    x1=x[0:5:4]
    x1.append(x[2])
    p1=Obs_data.pressure_value[0:5:4]
    p1.append(Obs_data.pressure_value[2])

    x2=x[1]
    p2 = float(Obs_data.pressure_value[1])

    x3=x[3]
    p3=float(Obs_data.pressure_value[3])


    fig1= plt.figure()
    axes1=fig1.add_axes([.1,.1,.8,.8])

    axes1.plot(x1,p1,color="red",marker='+',markersize=8,label='Constant z when z = 49.5 meters')
    axes1.plot(x2,p2,color="purple",marker='s',markersize=8)
    axes1.plot(x3,p3,color="green",marker='o',markersize=8)
    axes1.set_xlabel('Distance in x-direction [meter]')
    axes1.set_ylabel('Pressure [Pa]')
    axes1.legend(loc=2)
    axes1.set_title('Aqueous Pressure vs Distance')
    
    show()

if __name__ == "__main__":

    import os
    import run_amanzi


    CWD = os.getcwd()

    # set up the run directory and cd into it
    run_directory = os.path.join(CWD,"output")
    if os.path.isdir(run_directory):
        [os.remove(f) for f in os.listdir(run_directory)]
    else:
        os.mkdir(run_directory)
    os.chdir(run_directory)

    try:
        run_amanzi.run_amanzi('../amanzi_steady_linear.xml')
        
        obs_xml=loadInputXML('../amanzi_steady_linear.xml')
        obs_data=loadDataFile(obs_xml)
        plotExample1(obs_xml,obs_data)    
    finally:
        os.chdir(CWD)

