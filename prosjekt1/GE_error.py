#plotting the relative error as a function og step length
from numpy import *
from matplotlib.pyplot import *

error = loadtxt("error.dat")
print "For n = %d, the maximum relative error between the exact solution and the relative error equals %f" % (len(error)+1,max(error))
h = []

for i in range(1,7):
    h.append(1./(10**i))
#Plotting error values, collected from each c++ run in the terminal.
plot(h,[-1.100580,-3.079400,-5.079180,-7.079180,-9.080470,-11.259100],'r-o',label="1/10")
plot(h[1],[-3.079400],'b-o',label="h = 1/100")
plot(h[2],[-5.079180],'k-o',label="h = 1/1\'000")
plot(h[3],[-7.079180],'y-o',label="h = 1 /10\'000")
plot(h[4],[-9.080470],'c-o',label="h = 1/100\'000")
plot(h[5],[-11.259100],'r-o',label="h = 1/1000\'000")

title("Logarithmic error based on steplength")
xlabel("h")
ylabel("log10 of Relative error ($\epsilon$)")
grid("on")
legend(loc=4)
show()
