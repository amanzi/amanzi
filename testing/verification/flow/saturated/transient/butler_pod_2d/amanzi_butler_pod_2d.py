import sys
#sys.path.append('.')
sys.path.append('../../../../../../tools/amanzi_xml')
sys.path.append('../../../../../../tools/prettytable')
import matplotlib.pyplot as plt
import numpy
import math
from amanzi_xml.observations.ObservationXML import ObservationXML as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import prettytable

# load the data from output file, which is in amanzi-output
def load_amanzi_obs():
    obs_data_file = open('amanzi-output/observations.out', 'r')
    obs_data_file.readline()
    obs_data_file.readline()
    obs_data = [[float(line.split(',')[4]),float(line.split(',')[5])] for line in obs_data_file]
    return obs_data

# load data from analytic solution
def load_ana_solution():
    ana_data =numpy.genfromtxt('analytic/drdn.res')
    return ana_data

def plottest(axes1, obs_data, ana_data):

    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [days]')
    axes1.set_title('Drawdown at Two Observation Wells')
    
    ntime1 = len(ana_data[:,0])/2
    ntime2= len(ana_data[:,0])

    ntime3 = len(obs_data[:,0])/2
    ntime4 = len(obs_data[:,0])

    axes1.plot(numpy.log10(ana_data[0:ntime1,0]), ana_data[0:ntime1,1], '-r', label='Butler Pod Solution: Obs1')
    axes1.plot(numpy.log10(obs_data[1:ntime3,0]), obs_data[1:ntime3,1], 'ro', label='Amanzi: Obs1')
    axes1.plot(numpy.log10(ana_data[ntime1+1:ntime2,0]), ana_data[ntime1+1:ntime2,1], '-b', label='Butler Pod Solution: Obs2')
    axes1.plot(numpy.log10(obs_data[ntime3+1:ntime4,0]), obs_data[ntime3+1:ntime4,1], 'bo', label='Amanzi: Obs2')

    axes1.legend(loc='upper left')

if __name__ == "__main__":

    import os
    import run_amanzi

    input_filename =os.path.join("amanzi_butler_pod_2d.xml")

    CWD = os.getcwd()
    try: 
        run_amanzi.run_amanzi('../'+input_filename)
        obs_data = load_amanzi_obs()
        obs_data = numpy.array(obs_data) 
        obs_data[:,0] =  obs_data[:,0]/24./3600.
        ana_data=load_ana_solution()
        ana_data[:,0] =  ana_data[:,0]/24./3600.

        fig1= plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        plottest(axes1,obs_data, ana_data)
#        plt.show()

    finally:
        os.chdir(CWD)


