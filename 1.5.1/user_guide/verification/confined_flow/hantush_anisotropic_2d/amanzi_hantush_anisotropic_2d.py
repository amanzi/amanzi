import matplotlib.pyplot as plt
import numpy
from numpy import linalg as linalg
import math
from amanzi_xml.observations.ObservationXMLv2 import ObservationXMLv2 as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import model_hantush_anisotropic_2d
import prettytable

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
    output_file = Obs_xml.getObservationFilename()
    Obs_data = ObsDATA("amanzi-output/"+output_file)
    Obs_data.getObservationData()
    coords = Obs_xml.getAllCoordinates()

    for obs in Obs_data.observations.values():
        region = obs.region
        obs.coordinate = coords[region]
    
    return Obs_data

def plotHantushObservations(Obs_xml, Obs_data, axes1):
    colors = ['g','r','b']

    i = 0
    for obs in Obs_data.observations.values():
        drawdown = numpy.array([obs.data])
        axes1.scatter(obs.times, drawdown, marker='s', s=25, c=colors[i])
        i = i + 1
        
    return colors

def plotHantushAnalytic(filename, axes1, Obs_xml, Obs_data):
    colors = ['g','r','b']

    mymodel = model_hantush_anisotropic_2d.createFromXML(filename)
    tindex = numpy.arange(118)
    times = []
    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [s]')
    axes1.set_title('Drawdown vs Time after Pumping')
    
    for i in tindex:
        times.append(1 + math.exp(float(i) * (i+1) / (10.29*len(tindex))))
    
    i = 0
    for obs in Obs_data.observations.values():
        x = obs.coordinate[0]
        y = obs.coordinate[1]
        drawdown=mymodel.runForFixedPoint(times, x, y)
        axes1.plot(times, drawdown, label='$x=%0.0f m$ y=%0.0f m'%(x, y), c=colors[i])
        i = i + 1
    
    axes1.legend(title='Hantush Solution',loc='lower right', fancybox=True, shadow=True)


if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file = os.path.join("amanzi_hantush_anisotropic_2d-u.xml")
    run_dir = "amanzi-output"

    cwd = os.getcwd()
    try: 
        run_amanzi_standard.run_amanzi(input_file, 1, ["porflow4_6.exo",input_file], run_dir)
        obs_xml=loadInputXML(input_file)
        obs_data=loadDataFile(obs_xml)

        fig1= plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        plotHantushObservations(obs_xml,obs_data,axes1)
        plotHantushAnalytic(input_file,axes1,obs_xml,obs_data)
        plt.savefig("hantush_anisotropic_2d.png",format="png")
        # plt.show()

    finally:
        os.chdir(cwd)


