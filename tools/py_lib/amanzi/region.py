################################################################################

import os, sys
import types

from .trilinos import ParameterList

################################################################################
class BoxInterface(ParameterList):

    def __init__(self,lo=[],hi=[],file=None):

        ParameterList.__init__(self,"box",file)

        if len(lo) != 0:
            self.add_parameter("lo",lo)

        if len(hi) != 0:
            self.add_parameter("hi",hi)


class RegionsInterface(ParameterList):

    def __init__(self,file=None):

        ParameterList.__init__(self,"Regions",file)

    def add_region(self,label=None,lo=[],hi=[]):

        if label == None:
            raise ValueError('Must define a label to add a region')

        if len(lo) == 0:
            raise ValueError('Must define a lower bound region')

        if len(hi) == 0:
            raise ValueError('Must define an upper bound region')

        new_region = self.add_sublist(label)
        box = BoxInterface(lo,hi)
        new_region.add_sublist(box)

        return new_region


    
################################################################################
def RegionInputList(file=None):

    return RegionsInterface(file)



################################################################################
if __name__ == '__main__':

    # Example based on the deep vadose example

    regions = RegionInputList()

    # Rwia region
    regions.add_region('Rwia region',[0.0,0.0,0.0],[103.2,0.0,6.0])

    # Rlm region
    regions.add_region('Rlm region',[0.0,0.0,6.0],[103.2,0.0,11.4])

    # CCug region
    regions.add_region('CCug region',[0.0,0.0,11.4],[103.2,0.0,18.0])

    # CCuz region
    regions.add_region('CCuz region',[0.0,0.0,18.0],[103.2,0.0,22.2])

    regions.dumpXML()


