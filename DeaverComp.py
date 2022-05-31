# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 11:47:31 2022

@author: 2175469R
"""

from DeaverFittingWork import *


start_time = time.time()
path = r"C:\Users\2175469R\OneDrive - University of Glasgow\Documents\PhD stuff\Measurements\Zephyrnator\20_03_18_Nb Nanobridges JAC1034\BJ02_3_6.csv"


I,Vdata,R = loaddata(path,0)

#limI = I[-1]
#limV = Vdata[-1]

#for i in range(1,301):
    
#    I.append(limI + i*1e-6)
#    Vdata.append(limV + R*i*1e-6)

Ic = 0.000155
In = 0.000016
R = 117
l = 0.75
#plt.scatter(I,Vdata,label = 'data')

initials = [Ic,In,R,l]


VFit = DeaverNoise(I,Ic,In,R,l)
print("--- %s seconds" % (time.time() - start_time))

#write2csv(I,VFit,path,initials)

plt.figure()  

plt.scatter(I,VFit,label = "fit",s = 2)
plt.scatter(I,Vdata,label = "data",s = 2)
plt.xlabel('Current')
plt.ylabel('Voltage')
plt.title("Ic = " + str(Ic*1e6)+", In = "+str(In*1e6)+", l = " +str(l)+ ", R = " +str(R))
#plt.title('Ic = '+str(Ic*1e6) + 'uA ' + 'Noise = '+str(In*1e6) + 'uA ' + 'R = '+str(round(R)) + 'Ohm ' + 'l = '+str(l))
plt.legend(fontsize = 20)

residuals = []
VCrop = []
ICrop = [] 
SStot = []

for i in range(0,len(I)):
    residuals.append((Vdata[i] - VFit[i])**2)
    SStot.append((Vdata[i] - np.mean(Vdata))**2)    
SSres = sum(residuals)
SStot = sum(SStot)
    
Rsquare = 1-((len(Vdata)-1)/(len(Vdata)-4))*(1-(1-SSres/SStot))
    
print(sum(residuals))
print(Rsquare)

"""for i in range(len(I)):
    I[i] = 1000*I[i]
    
for i in range(len(Vdata)):
    Vdata[i] = 1000*Vdata[i]

plt.figure()
plt.scatter(ICrop,VCrop,label = 'fitted')
plt.scatter(I,Vdata,label = 'data')
plt.legend()

    
"""