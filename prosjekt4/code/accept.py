from numpy import *
from matplotlib.pyplot import *

y=[]
x=[]
k=0
file=open("1111","r")
for line in file:
	k+=1
	line = line.split()
	x.append(float(line[0]))
	y.append(float(line[1])*4./100)
plot(x,y)
xlabel("Temperature [relative]",size = 16)
ylabel("Accepted States",size = 16)

show()
