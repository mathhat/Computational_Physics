#This script plots the distribution of the agents' money
from matplotlib.pyplot import *
from numpy import *
from math import factorial
filename = "exp7lambda0.500000"
file = open(filename,"r")
lamb = float(filename[10:14])
print lamb
money=[]
for line in file:
	money.append(float(line))
money =asarray(money)
N = len(money)
mac = max(money)

bins = 600
hist(money,bins = bins)
show()



X = linspace(0,14,N)
def n(lamb):
	return 1.+ 3.*lamb/(1.-lamb)
def P(n,X):
	return n**n/factorial(n)*X**(n-1)*exp(-n*X)
show()


y,x = histogram(money,bins = bins)
Y = y/float(sum(y))
plot(x[0:-1],Y)

#### Gibbs dist through gamma functions
X = linspace(0,140,N)
def n(lamb):
	return 1.+ 3.*lamb/(1.-lamb)
def P(n,X):
	return n**n/factorial(n-1)*(0.1*X)**(n-1)*exp(-n*0.1*X)
p = P(n(lamb),X)
plot(X,p*max(Y)/max(p))
show()
####

plot(x,log(0.1*exp(-x*0.1)),label="omg")

#legend()
show()

