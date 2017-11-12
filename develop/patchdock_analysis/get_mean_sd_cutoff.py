import sys
from numpy import mean,var,sqrt
from pylab import * 

datafile = sys.argv[1]

fl_list = open(datafile,'r')

f = [float(f.strip('\n')) for f in fl_list.readlines() ]

sd = sqrt(var(f))
mn = mean(f)

cutoff = mn+2*sd

print "The cut off value is ", cutoff



hist(f,normed=1)

axvline(cutoff,linewidth=4, color='r')

show()
