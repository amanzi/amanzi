import sys
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

# load the data from output file, which is in amanzi-output
def load_amanzi_obs():
    output_file = obs_xml.getObservationFilename()
    obs_data = ObsDATA("amanzi-output/"+output_file)
    obs_data.getObservationData()   # get all data
    return obs_data

# load data from analytic solution
def load_ana_solution():
    ana_data =numpy.genfromtxt('analytic/drdn.res')
    return ana_data

def plottest(axes1, obstimes, obsdata, ana_data):

    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [days]')
    axes1.set_title('Drawdown at Two Observation Wells with Distance r1  = 24m and r2 = 100m')
    
    ntime1 = int(len(ana_data[:,0])/2)
    ntime2 = int(len(ana_data[:,0]))

    ntime3 = int(len(obstimes[:,0]))
    ntime4 = int(len(obstimes[:,1]))

    axes1.plot(numpy.log10(ana_data[0:ntime1,0]), ana_data[0:ntime1,1], '-r', label='Butler Strip Solution: r=24m')
    axes1.plot(numpy.log10(obstimes[1:ntime4,1]), obsdata[1:ntime4,1], 'ro', label='Amanzi: r=24m')
    axes1.plot(numpy.log10(ana_data[ntime1+1:ntime2,0]), ana_data[ntime1+1:ntime2,1], '-b', label='Butler Strip Solution: r=100m')
    axes1.plot(numpy.log10(obstimes[1:ntime3,0]), obsdata[1:ntime3,0], 'bo', label='Amanzi: r=100m')

    axes1.legend(loc='upper left')

if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file =os.path.join("amanzi_butler_strip_2d-u.xml")
    run_dir = "amanzi-output"

    CWD = os.getcwd()
    try: 
        run_amanzi_standard.run_amanzi(input_file, 16, [input_file], run_dir)
        obs_xml = loadInputXML(input_file)
        obs_data = load_amanzi_obs()

        obsdata = []
        obstimes = []
        for obs in obs_data.observations.values():
            obsdata.append(obs.data)
            obstimes.append(obs.times)

        obsdata  = numpy.transpose(numpy.array(obsdata))
        obstimes  = numpy.transpose(numpy.array(obstimes))
        obstimes[:] =  obstimes[:]/24./3600.

        ana_data=load_ana_solution()
        ana_data[:,0] =  ana_data[:,0]/24./3600.

        fig1= plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        plottest(axes1,obstimes, obsdata, ana_data)
        # plt.show()

    finally:
        os.chdir(CWD)


