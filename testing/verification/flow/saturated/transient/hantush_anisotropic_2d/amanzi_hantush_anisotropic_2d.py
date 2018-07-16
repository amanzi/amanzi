import sys
import matplotlib.pyplot as plt
import numpy
import math
from amanzi_xml.observations.ObservationXMLv2 import ObservationXMLv2 as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import prettytable

import os
import run_amanzi_standard

def load_amanzi_obs():
    #output_file = Obs_xml.getObservationFilename()
    output_file = "observation.out"
    obs_data = ObsDATA("amanzi-output/" + output_file)
    obs_data.getObservationData()
    return obs_data

def load_ana_solution():
    # load data from analytic solution
    ana_data =numpy.genfromtxt('analytic/test_h_tr.dat')
    return ana_data

def plottest(axes1, obstimes, obsdata, ana_data):
    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [days]')
    axes1.set_title('Drawdown at Observation Well with Distance r1  = 24m and r2 = 100m ')
    
    ntime1 = len(ana_data[:,0])/2
    ntime2= len(ana_data[:,0])

    ntime3 = len(obstimes[:,0])
    ntime4 = len(obstimes[:,1])

    axes1.plot(numpy.log10(ana_data[0:ntime1,0]), ana_data[0:ntime1,1], '-r', label='Analytical Solution: r=24m')
    axes1.plot(numpy.log10(obstimes[1:ntime4,1]), obsdata[1:ntime4,1], 'ro', label='Amanzi: r=24m')
    axes1.plot(numpy.log10(ana_data[ntime1+1:ntime2,0]), ana_data[ntime1+1:ntime2,1], '-b', label='Analytical Solution: r=100m')
    axes1.plot(numpy.log10(obstimes[1:ntime3,0]), obsdata[1:ntime3,0], 'bo', label='Amanzi: r=100m')

    axes1.legend(loc='lower right')

if __name__ == "__main__":

    input_file =os.path.join("amanzi_hantush_anisotropic_2d-u.xml")
    run_dir = "amanzi-output"

    try: 
        run_amanzi_standard.run_amanzi(input_file, 1, ["porflow4_6.exo",input_file], run_dir)
        obs_data = load_amanzi_obs()

        obsdata = []
        obstimes = []
        for obs in obs_data.observations.itervalues():
            obsdata.append(obs.data)
            obstimes.append(obs.times)

        obsdata  = numpy.transpose(numpy.array(obsdata))
        obstimes  = numpy.transpose(numpy.array(obstimes))
        obstimes[:] =  obstimes[:] /24. / 3600.

        # analytic solution is not available yet
        ana_data=load_ana_solution()

        fig1 = plt.figure()
        axes1 = fig1.add_axes([.1,.1,.8,.8])
       
        plottest(axes1, obstimes, obsdata, ana_data)
        plt.savefig("hantush_anisotropic_2d.png",format="png")
        # plt.show()
 
    finally:
        pass

