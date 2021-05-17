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
    output_file = Obs_xml.getObservationFilename()
    Obs_data = ObsDATA("amanzi-output/"+output_file)
    Obs_data.getObservationData()
    return Obs_data

def plotSingleFractureObservations(Obs_xml, Obs_data, axes1):
    regions = []
    for obs in Obs_data.observations.values():
        regions.append(obs.region)

    colors = ['b','g','r']
    cmap = dict((region,color) for (region,color) in zip(regions, colors))

    legends = ["bottom region","fracture","Amanzi"]
    lmap = dict((region,legend) for (region,legend) in zip(regions, legends))

    for obs in Obs_data.observations.values():
        color = cmap[obs.region]
        legend = lmap[obs.region]
        # a few pictures could be made later
        if legend == "Amanzi" :
            axes1.scatter(obs.times, numpy.log10(obs.data), marker='o', s=25, c=color, label=legend)


if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file = os.path.join("amanzi_single_fracture.xml")
    path_to_amanzi = "amanzi-output"

    cwd = os.getcwd()
    try: 
        max_np = 1
        run_amanzi_standard.run_amanzi(input_file, max_np,
                                       ["single_fracture.exo", input_file], path_to_amanzi)
        obs_xml=loadInputXML(input_file)
        obs_data=loadDataFile(obs_xml)

        fig1 = plt.figure()
        axes1 = fig1.add_axes([0.12,0.1,0.86,0.84])
       
        plotSingleFractureObservations(obs_xml,obs_data,axes1)

        # first benchmark (2 figures are planned for future)
        data = []
        data = open("benchmark/UiB_RT0.txt").read().split();
        for i in range(0, len(data)): 
            data[i] = data[i].split(',')
        data = numpy.array(data, dtype=numpy.double)

        colors = ['k','b','g','g']
        for i in range(1, 4):
            x = data[:,0]
            y = numpy.log10(data[:,i])
            if i == 3 :
                axes1.plot(x, y, linestyle='-', c=colors[i], label="UiB RT0")

        # second benchmark (2 figures are planned for future)
        data = []
        data = open("benchmark/USTUTT_MPFA.txt").read().split();
        for i in range(0, len(data)): 
            data[i] = data[i].split(',')
        data = numpy.array(data, dtype=numpy.double)

        colors = ['k','b','g','k']
        for i in range(1, 4):
            x = data[:,0]
            y = numpy.log10(data[:,i])
            if i == 3 :
                axes1.plot(x, y, linestyle='-', c=colors[i], label="USTUTT MPFA")

        axes1.set_xlabel('Simulation time [s]')
        axes1.set_ylabel('Integrated solute flux at the outlet')
        axes1.legend(loc="lower right", fancybox=True, shadow=True)

        # plt.show()

    finally:
        os.chdir(cwd)


