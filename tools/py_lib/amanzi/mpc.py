################################################################################

import os, sys
import types

from trilinos import Parameter, ParameterList

################################################################################

class VizBase(ParameterList):

    def __init__(self):

        ParameterList.__init__(self,'Viz Parameters')
        self.dt = -1.0
        self.dnc = -1
        self.file = ''
        self.add_parameter('Dump time frequency', self.dt)
        self.add_parameter('Dump cycle frequecy', self.dnc)
        self.add_parameter('File name', self.file)

    def set_dt(self,value):
        self.dt = value
        node = self.find_parameter('Dump time frequency')
        if node != None:
            node.set('value', value)

        return node

    def set_dnc(self,value):
        self.dnc = value
        node = self.find_parameter('Dump cycle frequency')
        if node != None:
            node.set('value', value)

        return node

    def set_file(self,value):
        self.file = value
        node = self.find_parameter('Dump cycle frequency')
        if node != None:
            node.set('value', value)

        return node

        
class CGNS(VizBase):

    def __init__(self,file=None):

        VizBase.__init__(self)
        if file != None:
            self.set_file(file)

class MPC(ParameterList):

    def __init__(self):

        ParameterList.__init__(self,'MPC')
        self.start_time = 0.0
        self.set_parameter('Start Time',0.0)
        self.end_time = 0.0
        self.set_parameter('End Time', 0.0)
        self.enable_flow = bool(True)
        self.set_parameter('enable Flow', bool(True))
        self.enable_transport = bool(True)
        self.set_parameter('enable Transport', bool(True))
        self.enable_chemistry = bool(True)
        self.set_parameter('enable Chemistry', bool(True))

        viz = self.add_sublist('Viz Parameters')
        viz.add_sublist(CGNS('dummy.cgns'))

        self.viz = viz

    def set_start_time(self,value):
        self.start_time = value
        node = self.set_parameter('Start Time',value)
        return node
 
    def set_end_time(self,value):
        self.end_time = value
        node = self.set_parameter('End Time',value)
        return node

###############################################################################

if __name__ == '__main__':

    mpc = MPC()
    mpc.viz.set_file('fbasin.cgns')
    mpc.set_end_time(3600000.0)

    mpc.dumpXML()


        

       

