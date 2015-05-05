import numpy

vs = numpy.array(list(numpy.ones(41040)*1.9976E-12)+list(numpy.ones(41472)*6.9365E-11)+list(numpy.ones(28080)*2.0706E-09))

fh = open('dvz_att.txt','w')
fh.write('1\n')
fh.write('perm cell -1\n')
for i in range(0, len(vs), 10):
    fh.write(" ".join([str(v) for v in vs[i:i+10]])+'\n')
fh.close()
