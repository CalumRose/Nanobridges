# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:51:58 2020

@author: calst

distorted_sin is a method of finding the current values form equation 4 in the
Deaver and Pierce paper

the method finds the turning points of the function at each phi and then uses
these in a bisection root finding algorithm
"""

import numpy as np

def distorted_sin(phi,l):
    current = []
    phases = []
    tprange = np.linspace(-1,1,1000)
    for i in range(0,len(phi)):
        searchpoints = [-1,1]
        funvals = func(tprange,phi[i],l)
        tps = turningpoints(funvals)
        tps = tprange[tps[:]]
        searchpoints[1:1] = tps
        roots, angles = bisection(searchpoints,phi[i],l)    
        current = np.append(current,roots)
        phases = np.append(phases, angles)

    return current, phases

def func(s,phase,l):
    """This is the function as defined in D&P"""
    return s - np.sin(phase-l*2*np.pi*s)

def turningpoints(lst):
    """for a function f, takes f(x) and returns indices of turning points 
    in given list"""
    tps = []
    for i in range(1,len(lst)-1):
        if ((lst[i-1] < lst[i] and lst[i+1] < lst[i])
        or (lst[i-1] > lst[i] and lst[i+1] > lst[i])):
           tps.append(i)
    return tps

def bisection(args,phi,l):
    """This function checks if each of given search values is a root and 
    makes sure this is a root between the two values before using the bisection method"""
    roots = []
    angles = []
    for i in range(0,len(args)-1):
        a = args[i]
        b = args[i+1]
        a_n = func(a,phi,l)
        b_n = func(b,phi,l)
        if a_n == 0:
            roots = roots + [args[i]]
            angles.append(phi)
        elif a_n*b_n < 0: 
            root = [method(a,b,20,phi,l),]
            angles.append(phi)
            roots = roots + root
    return roots, angles

def method(a,b,N,phi,l):
    """this is the bisection method itself"""
    if func(a,phi,l)*func(b,phi,l) >= 0:
       return False
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = func(m_n,phi,l)
        if func(a_n,phi,l)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif func(b_n,phi,l)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            return m_n
    return (a_n+b_n)/2






