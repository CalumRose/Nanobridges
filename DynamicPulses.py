# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:25:30 2020

@author: calst
"""

from distortedsin import *
import numpy as np
from scipy.integrate import odeint
from funcJackel import *
import time
import matplotlib.pyplot as plt
import random

"""using functions from the distortedsin program, single value returns the 
value of the CPR closest to the previous value used in the differential
 equation solver"""


def singlevalue(phi, l, prev):
    searchpoints = [-1, 1]
    s = np.linspace(-1, 1, 1000)
    funvals = func(s, phi, l)
    tps = turningpoints(funvals)
    tps = s[tps[:]]
    searchpoints[1:1] = tps
    roots, p = bisection(searchpoints, phi, l)

    desiredvalue = roots[min(range(len(roots)), key=lambda i: abs(roots[i] - prev))]

    return desiredvalue


"""Here the phase gradient and currnet are determenined at each timestep"""


def model(a, l, phi, prev, flag, theta):
    if flag == True:
        I = singlevalue(phi, l, prev)
    else:
        vals = funcJackel(phi, theta)
        I = vals[min(range(len(vals)), key=lambda i: abs(vals[i] - prev))]
    dydt = a - I
    return dydt, I


""""odesolve is the iterative process of solving (15) it can be used with or
 without a coupling current"""


def odesolve(a, l, step, flag, theta, limit):
    #phi0 = random.randint(0,100)/5*np.pi
    phi0 = 0
    tps = 0
    count = 0

    prev = phi0
    dydt, prev = model(a[count], l, phi0, 0, flag, theta)
    I = [prev]
    dydt = [dydt, ]
    phi = [(phi0 + step * dydt[0])]

    #I removed the tps limit because it was causing issue with the plotting.  I want the output to be the same size
    #as the bias current
    while count < limit:
        gradient, prev = (model(a[count], l, phi[count], prev, flag, theta))
        # gradient = gradient + 0.1*np.sin(step*i*0.11)
        dydt.append(gradient)
        newphi = phi[count] + step * dydt[count + 1]
        phi.append(newphi)
        I.append(prev)
        if dydt[count + 1] < dydt[count] and dydt[count - 1] < dydt[count]:
            tps += 1
        count += 1
    t = np.linspace(0, (count + 1) * step, count + 1)
    return phi, t, dydt, I


def main():
    #units normalised to usual SFQ parameters
    bias0 = 0.7
    bias = []
    inductance = 0.2
    timestep = 1e-2
    flag = 1
    thetaC = 0

    #I've added a resistance here that converts the output of J1 into a current to drive J2, but since it's normalised it's just 1
    #it might be required in the future for something (eg. if Ic1 != Ic2)
    R = 1

    limit = 10000

    #this is where you set the dynamic device bias.  I wanted to see a single pulse so I've given it 0.7 Ic + an RF signal
    #with 0.35 Ic amplitude.  This is kind of a rushed way to do it, we can revisit if we need to.
    
    
    for l in np.arange(0, limit, 1):
        bias.append(bias0 + 0.35*np.sin(l/limit*np.pi))

    phase, time, Voltage, supercurrent = odesolve(bias, inductance, timestep, flag, thetaC, limit)

    plt.figure()
    plt.scatter(time,phase,s=2)
    #plt.scatter(time, Voltage, s=2)
    plt.scatter(time, supercurrent, s=2)
    #plt.show()
    
    #I basically just do that again for J2
    for i in range(5):
        bias0 = 0.7
        pBias = []
    
        #Here I set the bias for the second junction: it's the 0.7 Ic DC bias + the voltage pulse from the first junction
        for l in np.arange(0, limit, 1):
            pBias.append(bias0 + Voltage[l] / R)
        
    
        #doing the calculation
        phase, time, Voltage, supercurrent = odesolve(pBias, inductance, timestep, flag, thetaC, limit)
    
        #plotting up the graph, uncomment to see I(t) and V(t) for J2
        plt.figure()
        plt.scatter(time,phase,s=2)
        plt.scatter(time, supercurrent, s=2)
        plt.title("Five self generated pulses for nanobridge with l = 0.2")
        plt.ylabel("Voltage")
        plt.xlabel("Time")
        #plt.scatter(time, supercurrent, s=2)
        #plt.show()
    
if __name__ == "__main__":
    main()