import numpy

vs = numpy.array(list(numpy.ones(41040)*1.9976E-12)+list(numpy.ones(41472)*6.9365E-11)+list(numpy.ones(28080)*2.0706E-09))

fh = open('dvz_elsets.txt','w')
fh.write('3\n')
fh.write('1100 cell 41040\n')
vs = numpy.arange(1,41040+1)
for i in range(0, len(vs), 10):
    fh.write(" ".join([str(v) for v in vs[i:i+10]])+'\n')
fh.write('1200 cell 41472\n')
vs = numpy.arange(41041,41472+41040+1)
for i in range(0, len(vs), 10):
    fh.write(" ".join([str(v) for v in vs[i:i+10]])+'\n')
fh.write('1300 cell 28080\n')
vs = numpy.arange(41472+41040+1,28080+41472+41040+1)
for i in range(0, len(vs), 10):
    fh.write(" ".join([str(v) for v in vs[i:i+10]])+'\n')
fh.close()
