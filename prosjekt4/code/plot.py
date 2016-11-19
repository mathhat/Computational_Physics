import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
n = 60

m = np.loadtxt("mag1.dat")
e = np.loadtxt("energy1.dat")

energy = e[:,0]
#mag = m[:,0]
X = e[:,0]



plt.plot(np.log10(X),energy)
plt.plot(np.log10(X),mag)
plt.show()



#Following comments are the algorithm that plots probability density
'''
ex,y= np.histogram(energy,bins ="auto")
y = y[:-1] + (y[1] - y[0])/2   # convert bin edges to centers
pdf = ex/float(sum(ex))
plt.plot(y,pdf,"ro",label="Data")
f = UnivariateSpline(y, ex, s=2*n)
plt.plot(y, f(y)/max(f(y))*max(pdf),label="Spline interpolation")
'''

'''
m2 = np.loadtxt("mag2.dat")
e2 = np.loadtxt("energy2.dat")
energy2 = e2[:,0]
mag2 = m2[:,0]
ex2,y2= np.histogram(energy2,bins =50)
y2 = y2[:-1] + (y2[1] - y2[0])/2   # convert bin edges to centers
pdf2 = ex2/float(sum(ex2))
plt.plot(y2,pdf2,"ro")#,label="Data")
f2 = UnivariateSpline(y2, ex2, s=2*n)
plt.plot(y2, f2(y2)/max(f2(y2))*max(pdf2))#,label="States of T=2.4")
'''
'''
plt.title('Probability of energy state, T = 1 [kT/J]',size = 16)
plt.xlabel(r'Scaled Energy state [$\hat{E}$ = E/($L^2$J)]',size = 16)
plt.ylabel(r'Probability  $P(\hat{E})$',size = 16)
plt.legend(loc="best")
plt.show()
'''
