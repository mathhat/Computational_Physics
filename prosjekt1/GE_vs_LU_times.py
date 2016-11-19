#Here we plot some of the timetables in benchmarks.txt to investigate how much time the optimized GE and armadillo LU algorithm use.

from matplotlib.pyplot import *
from numpy import loadtxt, linspace
0.066885 #500
time_LU=[0.00016,0.001043,0.066885,0.17482,0.474327,0.747078,621.976]

time_GE=[1e-06,7e-06,6.8e-05,0.000779]

n_GE = [10,100,1000,10000]


n_LU = [10,100,500,700,1000,1200,10000] #more data in lu to see how the graph works


plot(n_LU,time_LU,"r-o",label="LU")
ylabel("time (seconds")
xlabel("gridpoints (n)")
title("Times for LU as a function of n gridpoints")
show()

plot(n_GE,time_GE,"r-o",label="GE")
ylabel("time (seconds")
xlabel("gridpoints (n)")
title("Times for GE as a function of n gridpoints")
show()
