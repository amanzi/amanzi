import sys,os
import visit as v
import visit_rcParams as vrc

class VisItMWS:
    def __init__(self, figsize=vrc.rcParams['figsize'], 
                 directory="movie", prefix="movie-",
                 update_func=None):
        self.update_func = update_func
        sw = v.SaveWindowAttributes()
        sw.advancedMultiWindowSave = 1
        sw.resConstraint = sw.NoConstraint

        sw.width = figsize[0]
        sw.height = figsize[1]

        if not os.path.isabs(directory):
            directory = os.path.join(os.getcwd(), directory)
            
        if not os.path.isdir(directory):
            os.mkdir(directory)
        sw.fileName = os.path.join(directory, prefix)
        sw.format = sw.PNG

        self.sw = sw

    def subplots(self, ny, nx):
        """Tiles windows like matplotlib"""
        dx = self.sw.width / nx
        dy = self.sw.height / ny
        x = [i*dx for i in range(nx)]
        y = [j*dy for j in range(ny)]
        y.reverse() #rows go from top to bottom, y coordinate starts at bottom

        for j in range(ny):
            for i in range(nx):
                index = i+j*nx + 1
                winstr = "win%d"%index

                getattr(self.sw.subWindowAtts, winstr).position = (x[i],y[j])
                getattr(self.sw.subWindowAtts, winstr).size = (dx,dy)

    def getWindow(self, i):
        """Gets access to the window to set position and size manually."""
        return getattr(self.sw.subWindowAtts, "win%d"%i)
                
    def save(self, index=None):
        v.SetSaveWindowAttributes(self.sw)
        nstates = v.TimeSliderGetNStates()

        if index is None:
            for index in range(nstates):
                v.SetTimeSliderState(index)
                if self.update_func is not None:
                    self.update_func()
                v.SaveWindow()
        elif type(index) is int:
            assert index < nstates
            v.SetTimeSliderState(index)
            if self.update_func is not None:
                self.update_func()
            v.SaveWindow()
        else:
            for ind in index:
                assert ind < nstates
                v.SetTimeSliderState(ind)
                if self.update_func is not None:
                    self.update_func()
                v.SaveWindow()
                
          
