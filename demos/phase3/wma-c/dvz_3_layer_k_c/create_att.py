import numpy

# Creat perms
ks = numpy.array(list(numpy.ones(41040)*1.9976E-12)+list(numpy.ones(41472)*6.9365E-11)+list(numpy.ones(28080)*2.0706E-09))
# Create total_component_concentrations
rs = numpy.array(list(numpy.ones(41040)*0.)+list(numpy.ones(41472)*1.)+list(numpy.ones(28080)*2.))

fh = open('dvz_att.txt','w')
fh.write('2\n')
fh.write('perm cell -1\n')
for i in range(0, len(ks), 10):
    fh.write(" ".join([str(v) for v in ks[i:i+10]])+'\n')
fh.write('total_component_concentration cell -1\n')
for i in range(0, len(rs), 10):
    fh.write(" ".join([str(v) for v in rs[i:i+10]])+'\n')
fh.close()
