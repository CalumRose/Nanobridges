# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:44:51 2020

@author: calst
"""
import numpy as np

""""funcJackel calculates the CPR for a given value of phase and critical angle
it can return one or more values"""


def funcJackel(phi,theta):
    maxdist = theta/(2*np.pi)
    start = int(np.ceil((phi/(2*np.pi))-maxdist))
    stop = int(np.floor((phi/(2*np.pi))+maxdist))
    vals = []
    for i in range(start,stop+1,1):
        vals.append(phi/theta - 2*np.pi*i/theta)
    return vals

