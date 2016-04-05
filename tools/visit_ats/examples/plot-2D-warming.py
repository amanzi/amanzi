import sys
sys.path.append('../')
import visit_ats

# monthly plots -- update rcParams to format time as only month
import visit_rcParams as vrc
vrc.rcParams['time.format'] = "%b %Y"
vrc.rcParams['time.round'] = "month"

# works through sshfs
ats = visit_ats.ATSVis("/scratch/tundra/arctic/simulations-archive/ats-testsuite-2D/test3/test3-warming/", hostname="pingo.lanl.gov", n_windows=3)

ats.loadSources()

# window 1: temperatures
ats.setActiveWindow(1)
ats.plotTemperature()
ats.writeTime() # write the time on window 1

# window 2: ice sat + depths
ats.setActiveWindow(2)
ats.plotIceSaturation()

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

# # save
import visit_save
mws = visit_save.VisItMWS(directory="movie-2D-warming-zoom", prefix="2D-warming-zoom-",
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
mws = visit_save.VisItMWS(directory="movie-2D-warming", prefix="2D-warming-",
                                                    update_func=ats.update)

mws.subplots(3,1)
mws.save()

sys.exit()
