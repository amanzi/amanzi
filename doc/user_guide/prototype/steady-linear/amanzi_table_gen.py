import model_steady_linear as model
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
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
    pressure_amanzi = []
    coordinates = []
    mymodel = model.createFromXML(filename)

    for obs in Obs_data.observations.itervalues():
        coordinates.append([obs.coordinate[0], obs.coordinate[2]])
        pressure_amanzi.append(str(obs.data).rstrip(']').lstrip('['))

    pres_analytic = mymodel.run(coordinates)
    pressure_analytic = list(pres_analytic)
    x = PrettyTable(["x [m]", "z [m]", "Analytic [Pa]","Amanzi [Pa]"])
    x.padding_width = 1

    for coords, p_analytic, p_amanzi in zip(coordinates,pressure_analytic,pressure_amanzi):
        x.add_row([coords[0],coords[1],"%.3f" % float(p_analytic),"%.3f" % float(p_amanzi)])
    
    print x

if __name__=="__main__":
    
    input_filename= "amanzi_steady_linear.xml"
    obs_xml=loadInputXML(input_filename)
    obs_data=loadDataFile(obs_xml)
    MakeTable(obs_data,obs_xml,input_filename)

    
        
