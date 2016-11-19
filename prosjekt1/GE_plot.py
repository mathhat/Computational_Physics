#This script plots either the solution from prosjekt1.cpp or project1specific.cpp, based on which was run previously
from numpy import linspace, loadtxt
from matplotlib.pyplot import *

v_num,v_exc = loadtxt("v.dat")[:,0],loadtxt("v.dat")[:,1]
x = linspace(0,1,len(v_num))

plot(x,v_num,label="numerical")
plot(x,v_exc,label="exact")
legend()
title("solutions for u(x) with a %d x %d matrix" % (len(v_num),len(v_num)))
xlabel("x")
ylabel("u(x)")
show()

