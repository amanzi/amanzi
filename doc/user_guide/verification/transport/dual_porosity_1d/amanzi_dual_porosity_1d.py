import matplotlib.pyplot as plt
import numpy
import math
from amanzi_xml.observations.ObservationXMLv2 import ObservationXMLv2 as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
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
    output_file =  Obs_xml.getObservationFilename()
    Obs_data = ObsDATA("amanzi-output/"+output_file)
    Obs_data.getObservationData()
    coords = Obs_xml.getAllCoordinates()

    for obs in Obs_data.observations.values():
        region = obs.region
        obs.coordinate = coords[region]
    
    return Obs_data

def plotDualPorosityObservations(Obs_xml, Obs_data, axes1):
    for obs in Obs_data.observations.values():
        pass

    axes1.scatter(obs.times, obs.data, marker='+', s=25, c='b', label="Amanzi-U")

def plotSudickyAnalytic(filename, axes1, Obs_xml, Obs_data):
    times=[]
    exact=[]

    try:
        f = open(filename,'r')
        for _ in range(34):
            next(f)

        for line in f:
            data = line.split()
            times.append(float(data[0]))
            exact.append(float(data[3]))

    except:
        pass

    axes1.plot(times, exact, c='r', label="Analytic: Sudicky et al")

    # for press1, press2 in zip(press_analytic,press_amanzi):
    #     table_values.append([press1[0],press1[1],press2])
        

if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file =os.path.join("amanzi_dual_porosity_1d-u.xml")
    run_dir = "amanzi-output"

    cwd = os.getcwd()
    try: 
        max_np = 1
        run_amanzi_standard.run_amanzi(input_file, max_np, [input_file], run_dir)
        obs_xml=loadInputXML(input_file)
        obs_data=loadDataFile(obs_xml)

        fig1 = plt.figure()
        axes1 = fig1.add_axes([0.1,0.1,0.86,0.84])
       
        plotDualPorosityObservations(obs_xml,obs_data,axes1)
        plotSudickyAnalytic("sudicky/tracer_conc.txt",axes1,obs_xml,obs_data)

        axes1.set_ylabel('Normalized concentration [-]')
        axes1.set_xlabel('Simulation time [y]')
        axes1.set_title('Normalized concentation at the fracture outlet')
        axes1.legend(loc="lower right", fancybox=True, shadow=True)

        # plt.show()
        # use Theis as an example of making table
        # MakeTable(obs_data,obs_xml,input_file)

    finally:
        os.chdir(cwd)


