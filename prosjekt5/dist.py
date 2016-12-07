#This script plots the distribution of the agents' money
from matplotlib.pyplot import *
from numpy import *
from math import factorial
def PLOT(filename,alpha):
	file = open(filename,"r")
	#lamb = float(filename[10:13])
	lamb = float(filename[10:11])
	alpha = str(alpha)
	money=[]
	for line in file:
		money.append(float(line))
	money =asarray(money)
	N = len(money)
	bins = arange(0,1000,0.01)#float(400*float(alpha)**2)
	y,x = histogram(money,bins = bins)
	Y = y/float(sum(y))
	loglog(x[0:-1],Y,label=alpha)



for i in [0.5,1,1.5,2]:
	PLOT("exp7lambda0alpha%s"%str(i),i)
legend()
show()
