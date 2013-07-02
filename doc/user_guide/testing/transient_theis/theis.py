import matplotlib.pyplot as plt
import numpy

from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search




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



















# if __name__ == "__main__":

#     import os
#     import run_amanzi

#     input_filename = os.path.join("..","theis.xml")

#     CWD = os.getcwd()

#     # set up the run directory and cd into it
#     run_directory = os.path.join(CWD,"output")
#     if os.path.isdir(run_directory):
#          [os.remove(os.path.join(run_directory,f)) for f in os.listdir(run_directory)]
#     else:
#          os.mkdir(run_directory)
         
#     os.chdir(run_directory)

#     try:
#         run_amanzi.run_amanzi('../theis.xml')
#         obs_xml=loadInputXML(input_filename)
#         obs_data=loadDataFile(obs_xml)

#         fig1= plt.figure()
#         fig2 = plt.figure()
#         fig2.set_size_inches(3, 2)
#         axes1=fig1.add_axes([.1,.1,.8,.8])
       
#         cmap = plotExampleObservations(obs_xml,obs_data, axes1)
#         plotExampleModel(input_filename, cmap, axes1, obs_xml, obs_data)
        
#         plt.show()

#     finally:
#         os.chdir(CWD)


