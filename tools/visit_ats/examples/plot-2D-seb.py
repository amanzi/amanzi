import sys
sys.path.append('../')
import visit_ats

# monthly plots -- update rcParams to format time as only month
import visit_rcParams as vrc
vrc.rcParams['time.format'] = "%b %d %Y"
vrc.rcParams['time.fontheight'] = 0.06
#vrc.rcParams['time.round'] = "month"

ats = visit_ats.ATSVis("/Users/ecoon/tmp/pingo/ats-transects/c38/spinup-np8/", n_windows=3)
#ats = visit_ats.ATSVis("/Users/ecoon/tmp/pingo/ats-transects/c38/CESM8_5-last20-2/", n_windows=3)
ats.loadSources(True)

# window 1: temperatures
# ats.setActiveWindow(3)
# ats.plotTemperature()
# #ats.createSubsurfaceContour("temperature.cell.0", 273.15)
# ats.plotSurfaceTemperature()
# ats.plotSnowTemperature()
# ats.writeTime() # write the time on window 1

# window 2: ice sat + depths
ats.setActiveWindow(1)
ats.plotIceSaturation()
ats.plotPondedDepth()
ats.plotSnowDepth()
ats.writeTime() # write the time on window 1

# window 3: liquid sat
ats.setActiveWindow(2)
ats.plotLiquidSaturation()

# redo zoomed out
view = GetView2D()
view.fullFrameActivationMode = view.On
view.viewportCoords = (0.16,0.97,0.13,0.97)
view.windowCoords = (-0.01,18.5,4,6)

# draw
for win in ats.windows:
    win.slice()
    #win.exaggerateVertical(2.)
    win.view = view
ats.draw()

# save
import visit_save
mws = visit_save.VisItMWS(figsize=(2048,1024),
                          directory="movie-c38-last20", prefix="c38-last20-",
                          update_func=ats.update)
#mws.subplots(3,1)
mws.subplots(2,1)
mws.save([0,])

# sys.exit()
