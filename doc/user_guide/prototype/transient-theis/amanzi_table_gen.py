import model_theis 
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import numpy
from prettytable import PrettyTable

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

def MakeTable(Obs_data,Obs_xml,filename):
    #### ==== Initial height of water table is set to 20.0 meters ==== ####
    drawdown_amanzi = []
    coordinates = []
    mymodel = model_theis.createFromXML(filename)
    
  
    for obs in Obs_data.observations.itervalues():
        if obs.coordinate[0] == -55.5:
            coordinates.append([obs.coordinate[0], obs.coordinate[2]])
            pres0  = 101325 -9806.65 * (obs.coordinate[2] - 20)
            pres_drop = (pres0 - numpy.array([obs.data]))
            drawdown = pres_drop / 9806.65
            list(drawdown)
            print "first drawdown",drawdown
            drawdown_amanzi.append(drawdown)
            
    print drawdown_amanzi
   


    for coord in coordinates:
        if coord[0] == -55.5:
            pres_analytic = mymodel.runForFixedRadius(mymodel.times,coord[0])
            drawdown_analytic = list(pres_analytic)
        else:
            pass
   
    x = PrettyTable(["r [m]", "z [m]", "Analytic [m]","Amanzi [m]"])
    x.padding_width = 1

    print "this is coordinates", coordinates
    print "this is p_analytic", drawdown_analytic
    print "amanzi" , drawdown_amanzi

    for p_analytic, p_amanzi in zip(drawdown_analytic, drawdown_amanzi):
        x.add_row([coordinates[2][0],coordinates[2][1],"%.3f" % p_analytic,"%.3f" % p_amanzi[0]])
    
    print x

if __name__=="__main__":
    
    input_filename= "amanzi_transient_theis.xml"
    obs_xml=loadInputXML(input_filename)
    obs_data=loadDataFile(obs_xml)
    MakeTable(obs_data,obs_xml,input_filename)

    
        
