#!/usr/bin/env python
from pylab import *
import matplotlib.numerix.ma as ma
import loadfile
import re, optparse

p = optparse.OptionParser()
p.add_option("--file",action="store", type='string', dest="PlotFile")

(opt, args) = p.parse_args()

params = {'backend':'ps',
 	'axes.labelsize':24,
	'text.fontsize':24,
	'xticks.labelsize':16,
	'yticks.labelsize':16,	
	'text.usetex':True,
	'font.family':'serif',
	'font.size':24,
	'font.serif':['Times','Palatino','serif'],
	'ps.usedistiller':'xpdf',
	'lines.linewidth':3}
	
rcParams.update(params)

f = open(opt.PlotFile,"r")
s = f.readlines()
#f.close()
I1 = len(s)
I2 = 0

first_data_row = -1

for i in range(0,I1):
        if s[i].startswith('#'):
                labels = s[i].split()
        else: 
                if first_data_row == -1:
                        first_data_row = i
                I2= I2 + 1
                J1 = len(s[1].split())

#print labels


data = zeros([J1,I2], float)

i1 = 0;
i2 = 0;
while (i1<I1):
        while(s[i1].startswith('#')):
                        i1 = i1 + 1
        if s[i1]==" ":
                break
        l = s[i1].split()
        i = 0
        for j in range(0,J1):
                data[i][i2] = float(l[j])
                i = i + 1
        i1 = i1 + 1
        i2 = i2 + 1
        
f.close()        

#print data

data[0] = data[0]/(3600*24*365.25)

#print data

figure(1, [18,9])

plot(data[0], data[1], "r-", label=opt.PlotFile)

legend(loc=0,handlelen = 0.07)
xlabel(labels[1] + " (years)")
ylabel(labels[2])

savefig(opt.PlotFile+".png", dpi= 300)

#axis([1, 330, 1e-5, 2])

show()