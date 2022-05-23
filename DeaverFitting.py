# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 12:03:28 2021

@author: calst
"""


from distortedsin import *
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as opt
import time
import matplotlib.pyplot as plt
import pandas as pd

import os.path
import csv
from datetime import datetime

from distortedsin import *
from Euler import singlevalue

def loaddata(path,datatype):
    if datatype == 0:
        df = pd.read_csv(path)
        current = list(df.values[:,0])
        voltage = list(-df.values[:,1])
    if datatype == 1:
            with open(path + "/voltage.txt","rb") as fp:
                voltage = pickle.load(fp)
            with open(path + "/current.txt","rb") as fp:
                current = pickle.load(fp)
                
    offset = voltage[0]   
    
    pc = []
    pv = []
    r = []
    voltage[:] = [(x-offset) for x in voltage]
    #plt.figure()
    #plt.scatter(list(range(len(current))),current)
    #plt.figure()
    #plt.scatter(list(range(len(voltage))),voltage)
    #plt.figure()
    #plt.scatter(current,voltage)
    
    for i in range(0,len(current)-1):

        if current[i] >= 0 and current[i] < 0.00035 and i < 255:
            pc.append(current[i])
            pv.append(voltage[i])
        if current[i] >= 0.00025:
            r.append((voltage[i]-voltage[i-1])/(current[i]-current[i-1]))

    #plt.figure()
    #plt.scatter(pc,pv)
    #plt.figure()
    #plt.scatter(list(range(len(pc))),pc)
    #plt.figure()
    #plt.scatter(list(range(len(pv))),pv)
            
    pc, pv = (list(t) for t in zip(*sorted(zip(pc, pv))))
    #plt.figure()
    #plt.scatter(list(range(len(pc))),pc)
    r = np.mean(r)
    #print(str(r))
    return pc, pv, r

def write2csv(I,V,path,params):
    
    folder = r"C:\Users\2175469R\OneDrive - University of Glasgow\Documents\PhD stuff\Measurements\Zephyrnator\FittingData\\"
    sample =   path.split('\\')[-1].split('.')[0]
    if not(os.path.exists(folder+sample)):
        os.makedirs(folder + sample)
    filename = folder + sample + '\\' + str(params[0]) + '_' + str(params[1]) + '_' + str(params[2]) + '_' + str(params[3])
    filename = filename.replace('.','p')
    filename = filename + '.csv'
    
    if not(os.path.exists(filename)):
        with open(filename,'w',newline = '') as f:
            writer = csv.writer(f)
            writer.writerow(['I', 'V'])
            for i in range(0,len(I)):
                    writer.writerow([str(I[i]), str(V[i])])
    return
            
        
        

def odesolve(a,l,step,flag,theta):
    phi0 = np.pi
    prev = singlevalue(phi0,l,0)
    dydt = [a-prev,]
    I = [prev]
    phi = [(phi0 + step*dydt[0])]
    tps = 0
    tp = []
    count = 0
    while tps <4 and count<10000:
        I.append(singlevalue(phi[count],l,I[count]))
        dydt.append(a-I[count+1])
        newphi = phi[count] + step*dydt[count+1]
        phi.append(newphi)
        if dydt[count+1] < dydt[count] and dydt[count-1] < dydt[count]:
            tps +=1
            tp.append(count)
        count+=1
    t = np.linspace(0,(count)*step,count+1)
    if count>9000:
        tp = [800,990,990,990]
    return phi, t, dydt,I, tp


def generateIV(Irange,Ic,l,R):
    V = []
    times = []
    phase = []
    actualV = []
    turnpoints = []
    for i in range(0,len(Irange)):
        if Irange[i]<1:
            actualV.append(0)
            turnpoints.append([0])
            V.append([0])
            times.append([0])
            continue
            
        phi,t,vals,I,tp= odesolve(Irange[i],l,0.05,1,2*np.pi)
        V.append(vals)
        times.append(t)
        phase.append(phi)
        turnpoints.append(tp)
        actualV.append(np.mean(vals[tp[0]:tp[3]]))
    #plt.figure()
    #plt.scatter(Irange,actualV,s=2)
    return V, times, phase, actualV,turnpoints

def Gaussian(Irange,IBatt,In):
    
    Noise = []
    
    for i in Irange:
        Noise.append((1/(np.sqrt(2*np.pi)*In))*np.exp(-((i-IBatt)**2/(2*(In)**2))))
    return Noise

def endgrad(V,I,cutoff):
    g = []
    for i in range(cutoff,len(V)):
        g.append((V[i]-V[i-1])/(I[i]-I[i-1]))
    
    g = np.mean(g)
    print(str(g))
    return g

def DeaverNoise(I,eIc,eIn,eR,el):
    
    #eIc = params[0]
    #eIn = params[1]
    #eR = params[2]
    #el = params[3]
    
    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print("Ic = ", eIc, "In = ", eIn, "R = ", eR, "l =", el )
    
    INorm = []
    for i in range(len(I)):
        INorm.append(I[i]/eIc)
    
    ext = np.arange(0.02,6.02,0.02) + max(INorm)

    INorm = INorm + list(ext)
    
    V,times,phase,actualV,turnpoints = generateIV(INorm, eIc,el,eR)
    #endgrad(actualV,INorm,200)
    time1 = time.time()


    IExtended = []
    for i in INorm:
        IExtended.append(-i)
    IExtended = IExtended+INorm
    IExtended.sort()

    VFull = actualV[::-1]
    VFull[:] = [-i for i in VFull]
    VFull = VFull+actualV

    VFinal = []
    for IBatt in IExtended:
        VNoise = Gaussian(IExtended,IBatt,eIn/eIc)
        product = np.multiply(VNoise,VFull)
        VFinal.append(np.trapz(product, x = IExtended))
    
    VFinal = VFinal[len(I)+300:len(IExtended)-300]
    VFinal = [x*(eIc*eR) for x in VFinal]
    #endgrad(VFinal,I,200)
        
    return VFinal


def main():
    path = r"C:\Users\2175469R\OneDrive - University of Glasgow\Documents\PhD stuff\Measurements\Zephyrnator\20_03_18_Nb Nanobridges JAC1034\BJ02_3_6.csv"

    I,Vdata,R = loaddata(path,0)
    
    #g = endgrad(Vdata, I, 300)

    Ic = 0.00016
    In = 0.000015
    l = 0.1

    
    initials = [Ic,In,R,l]

    params,covar = opt.curve_fit(DeaverNoise,I,Vdata,p0 = initials,bounds = ([1e-4,0,R-2,0],[2e-4,5e-5,R+2,1]))

    
    VFit = DeaverNoise(I,params[0],params[1],params[2],params[3])
    
    write2csv(I,VFit,path,params)
    

    plt.figure()
    plt.scatter(I,VFit)
    plt.scatter(I,Vdata)
    
    residuals = []
    VCrop = []
    ICrop = [] 
    
    for i in range(0,len(I)):
        residuals.append((Vdata[i] - VFit[i])**2)
    
    
    for i in range(len(I)):
        I[i] = 1000*I[i]
        
    for i in range(len(Vdata)):
        Vdata[i] = 1000*Vdata[i]
    
    plt.figure()
    plt.scatter(ICrop,VCrop,label = 'fitted')
    plt.scatter(I,Vdata,label = 'data')
    plt.legend()
    return params,VFit
    
if __name__ == "__main__":
    params, VFit = main()



    




    
    
    
    