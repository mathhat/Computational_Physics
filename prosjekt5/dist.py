#This script plots the distribution of the agents' money
from matplotlib.pyplot import *
from numpy import *
from math import factorial
filename = "exp6lambda0.500000alpha0.000000"
file = open(filename,"r")
lamb = float(filename[10:14])

money=[]
for line in file:
	money.append(float(line))
money =asarray(money)
N = len(money)
mac = max(money)

bins = 600
hist(money,bins = bins)
show()





y,x = histogram(money,bins = bins)
Y = y/float(sum(y))
plot(x[0:-1],Y)

#### Gibbs dist through gamma functions
x = x#linspace(0,140,bins)
def n(lamb):
	return 1.+ 3.*lamb/(1.-lamb)
def P(n,x):
	return n**n/factorial(n-1)*(0.1*x)**(n-1)*exp(-n*0.1*x)
p = P(n(lamb),x)
plot(x,p/float(sum(p)))
show()
####

plot(x,log(0.1*exp(-x*0.1)),label="omg")

#legend()
show()

