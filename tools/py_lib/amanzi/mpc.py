################################################################################

import os, sys
import types

from .trilinos import Parameter, ParameterList

################################################################################

class VizBase(ParameterList):

    def __init__(self,label=None):

        ParameterList.__init__(self,'Viz Parameters')
        self.label = label
        self.dt = -1.0
        self.dnc = -1
        self.file = ''
        self.add_parameter('Output Type',self.label)
        self.add_parameter('Dump time frequency', self.dt)
        self.add_parameter('Dump cycle frequecy', self.dnc)
        self.add_parameter('File name', self.file)

    def set_dt(self,value):
        self.dt = value
        node = self.set_parameter('Dump time frequency',value)
        return node

    def set_dnc(self,value):
        self.dnc = value
        node = self.set_parameter('Dump cycle frequency',value)

        return node

    def set_file(self,value):
        self.file = value
        node = self.set_parameter('File name',value)

        return node

        
class CGNS(VizBase):

    def __init__(self,file=None):

        VizBase.__init__(self,'CGNS')
        if file != None:
            self.set_file(file)

class MPC(ParameterList):

    def __init__(self):

        ParameterList.__init__(self,'MPC')
        self.start_time = 0.0
        self.set_parameter('Start Time',0.0)
        self.end_time = 0.0
        self.set_parameter('End Time', 0.0)
        self.end_cycle = -1
        self.set_parameter('End Cycle', -1)
        self.enable_flow = bool(True)
        self.set_parameter('enable Flow', bool(True))
        self.enable_transport = bool(True)
        self.set_parameter('enable Transport', bool(True))
        self.enable_chemistry = bool(True)
        self.set_parameter('enable Chemistry', bool(True))

        self.viz = CGNS('dummy.cgns')
        viz_root = self.viz.getroot()
        self.attach(viz_root)

    def set_start_time(self,value):
        self.start_time = value
        node = self.set_parameter('Start Time',value)
        return node
 
    def set_end_time(self,value):
        self.end_time = value
        node = self.set_parameter('End Time',value)
        return node

    def set_end_cycle(self,value):
	self.end_cycle = value
	node = self.set_parameter('End Cycle',value)
	return node

###############################################################################

if __name__ == '__main__':

    mpc = MPC()
    mpc.viz.set_file('fbasin.cgns')
    mpc.viz.set_dt(0.5)
    mpc.set_end_time(3600000.0)

    mpc.dumpXML()


        

       

