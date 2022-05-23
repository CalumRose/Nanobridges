# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:44:44 2022

@author: 2175469R
"""
import csv
import pandas as pd
import matplotlib.pyplot as plt
import os
from DeaverFitting import *

def loadandplot(path):
    
    name = path.split("\\")[-1].split(".")[0]
    
    df = pd.read_csv(path)
    current = list(df.values[:,0])
    voltage = list(-df.values[:,1])
    
    plt.figure()
    plt.scatter(current,voltage,s=2)
    #plt.scatter(list(range(len(current))),current,s=2)
    plt.title(name)
    
    return

mydir = r"C:\Users\2175469R\OneDrive - University of Glasgow\Documents\PhD stuff\Measurements\Zephyrnator\20_03_18_Nb Nanobridges JAC1034"
files = os.listdir(mydir)

for file in files:
    if file.endswith(".csv"):
        path = os.path.join(mydir,file)
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
        plt.scatter(I,VFit,label = "fitted")
        plt.scatter(I,Vdata, label = "data")
        

