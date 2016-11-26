#This script plots the distribution of the agents' money
from matplotlib.pyplot import *
from numpy import *
from math import factorial
file = open("exp7lambda0.000000","r")
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
X = linspace(0,140,N)
def n(lamb):
	return 1.+ 3.*lamb/(1.-lamb)
def P(n,X):
	return n**n/factorial(n)*0.1*X**(n-1)*exp(-n*0.1*X)

plot(X,P(n(0),X)*10*max(Y))
show()


plot(x,log(0.1*exp(-x*0.1)),label="omg")

#legend()
show()

