import matplotlib.pyplot as plt
import numpy
import math
from amanzi_xml.observations.ObservationXMLv2 import ObservationXMLv2 as ObsXML
from amanzi_xml.observations.ObservationData import ObservationData as ObsDATA
import amanzi_xml.utils.search as search
import model_theis_isotropic_1d
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

def plotTheisObservations(Obs_xml, Obs_data, axes1):
    #=== SPECIAL CODE === Theis EXAMPLE === Water Table is 12.0 m 
    r_vals =[]
    for obs in Obs_data.observations.values():
        r_vals.append(obs.coordinate[0])
    r_vals = set(r_vals)
    r_vals=list(r_vals)
    r_vals.sort()
    colors = ['b','g','r']
    cmap = dict((rval,color) for (rval,color) in zip(r_vals, colors))
    #print cmap

    for obs in Obs_data.observations.values():
        color = cmap[obs.coordinate[0]]
        pres0 = 101325 -9806.65 * (obs.coordinate[2] - 20)
        pres_drop = (pres0 - numpy.array([obs.data]))
        drawdown = pres_drop / 9806.65
        axes1.scatter(obs.times, drawdown, marker='s', s=25, c=color)
        
    return cmap

def plotTheisAnalytic(filename, cmap, axes1, Obs_xml ,Obs_data):
    mymodel = model_theis_isotropic_1d.createFromXML(filename)
    tindex = numpy.arange(125)
    times = []
    table_values = []
    press_amanzi = []
    press_analytic=[]
    axes1.set_ylabel('Drawdown [m]')
    axes1.set_xlabel('Time after pumping [s]')
    axes1.set_title('Drawdown vs Time after Pumping')
    
    for i in tindex:
        times.append(1+math.exp(float(i)*(i+1)/(10.3*len(tindex))))
    
    for rad in mymodel.r:
        r= abs(rad)
        drawdown=mymodel.runForFixedRadius(times,r)
        axes1.plot(times,drawdown, label='$r=%0.1f m$'%r)
        press_analytic.append([r,drawdown])
    
    for obs in Obs_data.observations.values():
        press_amanzi.append(obs.data)
    
    for press1, press2 in zip(press_analytic,press_amanzi):
        table_values.append([press1[0],press1[1],press2])
        
    axes1.legend(title='Theis Solution',loc='lower right', fancybox=True, shadow=True)

def MakeTable(Obs_data,Obs_xml,filename):
    #### ==== Initial height of water table is set to 20.0 meters             ==== ####
    #### ==== The Table is set to generate at a coordinate of (55,y,z) *only* ==== ####
    #### ==== Modify if statement of another coordinate is desired            ==== ####
    drawdown_amanzi = []
    coordinates = []
    mymodel = model_theis_isotropic_1d.createFromXML(filename)
    
    for obs in Obs_data.observations.values():
        if obs.coordinate[0] == -55.0:
           coordinates.append([abs(obs.coordinate[0]), obs.coordinate[2]])
           for press in obs.data:
               pres0  = 101325 -9806.65 * (obs.coordinate[2] - 20)
               pres_drop = (pres0 - press)
               drawdown = pres_drop / 9806.65
               drawdown_amanzi.append(drawdown)

    drawdown_analytic = mymodel.runForFixedRadius(mymodel.times,coordinates[0][0])
    x = prettytable.PrettyTable(["time [s]","r [m]", "z [m]", "Analytic [m]","Amanzi [m]"])
    x.padding_width = 2
    x.hrules = 1

    for time,d_analytic, d_amanzi in zip(mymodel.times,drawdown_analytic, drawdown_amanzi):
        x.add_row([time,coordinates[0][0],coordinates[0][1],"%.10f" % d_analytic,"%.10f" % d_amanzi])
    
    if os.path.exists("table_values_theis.txt"):
        os.remove("table_values_theis.txt")

    table_file = open("table_values_theis.txt", "w+")
    table_file.write(x.get_string())
    table_file.close()

if __name__ == "__main__":

    import os
    import run_amanzi_standard

    input_file =os.path.join("amanzi_theis_isotropic_1d-u.xml")
    run_dir = "amanzi-output"

    cwd = os.getcwd()
    try: 
        max_np = 10
        run_amanzi_standard.run_amanzi(input_file, max_np, [input_file], run_dir)
        obs_xml=loadInputXML(input_file)
        obs_data=loadDataFile(obs_xml)

        fig1= plt.figure()
        axes1=fig1.add_axes([.1,.1,.8,.8])
       
        cmap = plotTheisObservations(obs_xml,obs_data,axes1)
        plotTheisAnalytic(input_file,cmap,axes1,obs_xml,obs_data)
        # plt.show()
        MakeTable(obs_data,obs_xml,input_file)

    finally:
        os.chdir(cwd)


