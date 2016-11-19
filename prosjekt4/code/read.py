#e m cv chi mcc L
from numpy import *
from matplotlib.pyplot import *
lab = "140"
file1 = open("22","r")
file2 = open("2020","r")
file3 = open("4040","r")
file4 = open("6060","r")
file5 = open("101101","r")
#file.readline()

def t (E): 
	return linspace(2.2,2.3,len(E))


def Plot(File,index,lab):#0 temp, 1 energy, 2 CV, 3 mag, 4 susc
	x = []		
	k =0
	for line in File:
		k+=1
		line = line.split()
		x.append(float(line[index])/4)
	plot(t(x),x,label=lab)




#Plot(file6,4,lab)
Plot(file1,4,"L=2")
Plot(file2,4,"L=20")
Plot(file3,4,"L=40")
Plot(file4,4,"L=60")
Plot(file5,4,"L=100")
xlabel("Temperature",size=16)
ylabel("Mean Magnetic Susceptibility",size=16)
legend(loc="best"	)

#
#show()
#file2 = open("2121","r")
#Plot(file2,3,"Mean Magnetization Per Spin")
#xlabel("MCC")
#ylabel("Mean Spin [normalized]")

show()

#Plot(file3,1,"40")
#Plot(file4,1,"60")
#Plot(file5,1,"100")

#	M.append(line[3])
#	CV.append(line[2])
#	CHI.append(line[4])

