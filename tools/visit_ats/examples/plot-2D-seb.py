import sys
sys.path.append('../')
import visit_ats
#ats = visit_ats.ATSVis("/lclscratch/ecoon/ats/ats-testsuite-2d-3/test6-2/test6-1",
#                       hostname="pingo.lanl.gov", n_windows=3)
#ats = visit_ats.ATSVis("test6-7", n_windows=3)
ats = visit_ats.ATSVis("/Users/ecoon/mnt/pingo/ats/ats-testsuite-2d-3/test6-2/test6-1", n_windows=3)
ats.loadSources(True)

# window 1: temperatures
ats.setActiveWindow(1)
ats.plotTemperature()
#ats.createSubsurfaceContour("temperature.cell.0", 273.15)
ats.plotSurfaceTemperature()
ats.plotSnowTemperature()
ats.writeTime() # write the time on window 1

# window 2: ice sat + depths
ats.setActiveWindow(2)
ats.plotIceSaturation()
ats.plotPondedDepth()
ats.plotSnowDepth()

# window 3: liquid sat
ats.setActiveWindow(3)
ats.plotLiquidSaturation()

# set the view
view = GetView2D()
view.viewportCoords = (0.16,0.97,0.13,1.4)
view.windowCoords = (14,21.5,8,10)

# draw
for win in ats.windows:
    win.slice()
    win.exaggerateVertical(2.)
    win.view = view
    ats.draw()

# save
import visit_save
mws = visit_save.VisItMWS(directory="movie-2D-seb-zoom", prefix="2D-seb-zoom-",
                                                    update_func=ats.update)

mws.subplots(3,1)
mws.save()

# redo zoomed out
view = GetView2D()
view.viewportCoords = (0.16,0.97,0.13,1.4)
view.windowCoords = (0,23.5,5,11)

for win in ats.windows:
    win.view = view
    ats.draw()

# save
mws = visit_save.VisItMWS(directory="movie-2D-seb", prefix="2D-seb-",
                                                    update_func=ats.update)

mws.subplots(3,1)
mws.save()

sys.exit()
