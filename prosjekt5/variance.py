from matplotlib.pyplot import *
from numpy import *
exp=8

file = open("variance","r")
var=[]
for line in file:
	var.append(float(line))
x = linspace(0,exp,len(var))
plot(x,var)
xlabel(r"Cycles/Transactions [$10^x$]",size=16)
ylabel(r"Variance $\sigma_m^2$ [cash$^2$]",size=16)
show()
