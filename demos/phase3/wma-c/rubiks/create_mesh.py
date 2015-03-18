import pylagrit

lg = pylagrit.PyLaGriT()

mo = lg.create()

mo.createpts_brick_xyz([4,4,4], [0,0,0],[1,1,1])

mo.dump_exo('rubiks.exo')

