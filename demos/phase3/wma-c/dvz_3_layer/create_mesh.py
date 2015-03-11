import pylagrit

lg = pylagrit.PyLaGriT()

mo = lg.create(mesh='hex')

mo.createpts_brick_xyz([433,2,257], [0.,0.,0.],[216.,1.,107.52])

mo.dump_exo('mesh.exo')

