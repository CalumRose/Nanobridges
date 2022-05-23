# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:25:30 2020

@author: calst
"""


from distortedsin import*
import numpy as np
from scipy.integrate import odeint
from funcJackel import *
import time
import matplotlib.pyplot as plt


"""using functions from the distortedsin program, single value returns the 
value of the CPR closest to the previous value used in the differential
 equation solver"""
def singlevalue(phi,l,prev):
    searchpoints = [-1,1]
    s = np.linspace(-1,1,1000)
    funvals = func(s,phi,l)
    tps = turningpoints(funvals)
    tps = s[tps[:]]
    searchpoints[1:1] = tps
    roots, p = bisection(searchpoints,phi,l)
    
    desiredvalue = roots[min(range(len(roots)), key=lambda i: abs(roots[i]-prev))]
    
    return desiredvalue


"""Here the phase gradient and currnet are determenined at each timestep"""
def model(a,l,phi,prev,flag,theta):
    if flag == True:
        I = singlevalue(phi,l,prev)
    else:
        vals = funcJackel(phi,theta)
        I = vals[min(range(len(vals)), key=lambda i: abs(vals[i] - prev))]
    dydt = a - I
    return dydt, I
    

""""odesolve is the iterative process of solving (15) it can be used with or
 without a coupling current"""
def odesolve(a,l,step,flag,theta):
    phi0 = np.pi
    prev = phi0
    dydt, prev = model(a,l,phi0,0,flag,theta)
    I = [prev]
    dydt = [dydt,]
    phi = [(phi0 + step*dydt[0])]
    tps = 0
    count = 0
    while tps <5 and count<100000:
        gradient, prev = (model(a,l,phi[count],prev,flag,theta))
        #gradient = gradient + 0.1*np.sin(step*i*0.11)
        dydt.append(gradient)
        newphi = phi[count] + step*dydt[count+1]
        phi.append(newphi)
        I.append(prev)
        if dydt[count+1] < dydt[count] and dydt[count-1] < dydt[count]:
            tps +=1 
        count+=1
    t = np.linspace(0,(count+1)*step,count+1)
    return phi, t, dydt,I

def main():
    bias = 1.1
    inductance = 0.1
    timestep = 0.01
    flag = 1
    thetaC = 2*np.pi
    
    phase, time, Voltage, supercurrent = odesolve(bias,inductance,timestep,flag,thetaC)
    plt.figure()
    plt.scatter(time,Voltage)
    
if __name__ == "__main__":
    main()


