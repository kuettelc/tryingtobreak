#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 09:56:37 2018

@author: christian
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import scipy.optimize as spo

prop_torque_coeff = 2.33*10**(-10)
orderfit = 2
vehicleMass = 50 #in grams
hover = vehicleMass/4
meanstart = 10 
meanend = 15 
Rcable = 0.144
files = sorted(glob.glob("*.xlsx"))

def force(V):
    a = meanresistance**2/meankv**2*meanthrustcoeff**2*prop_torque_coeff**2
    b = -(2*V*meanresistance*meanthrustcoeff*prop_torque_coeff/meankv + meanthrustcoeff*meankv**2)
    c = V**2
    F = (-b - np.sqrt(b**2-4*a*c))/2/a
    return F

for idx,xlsx in enumerate(files):
    print(xlsx)
    
    file = pd.read_excel(xlsx)
    
    voltage = file.Voltage
    current = file.Current
    voltage = voltage - current*Rcable
    thrust = file.Thrust
    
    start = np.where(thrust > meanstart)[0][0]
    end = np.where(thrust< meanend)[0][-1]
    
    hover_index = min(range(len(thrust)), key=lambda i: abs(thrust[i]-hover))
    print("Hover voltage: " + str(voltage[hover_index]))
    
    thrust = thrust/1000*9.81
    rpm = file.RPM
    omega = rpm*1000/60*2*np.pi
    omega2 = omega*omega
    
    if(idx == 0):
        meanomega = np.zeros(len(voltage))
        meanthrust = np.zeros(len(voltage))
        allomega = np.zeros((len(files),len(voltage)))
        allomega2 = np.zeros((len(files),len(voltage)))
        allvoltage = np.zeros((len(files),len(voltage)))
        allcurrent = np.zeros((len(files),len(voltage)))
        allthrust = np.zeros((len(files),len(voltage)))
        allthrustcoeff = np.zeros((len(files),len(voltage)))
        allhovervoltage = np.zeros(len(files))
        allhoveromega = np.zeros(len(files))
        
    
    allomega[idx] = omega
    allvoltage[idx] = voltage
    allthrust[idx] = thrust
    allthrustcoeff[idx] = omega2/thrust
    allcurrent[idx] = current
    allhovervoltage[idx] = voltage[hover_index]
    allhoveromega[idx] = omega[hover_index]
    
    
allomega = np.hstack(allomega.tolist())
allvoltage = np.hstack(allvoltage.tolist())
allthrustcoeff = np.hstack(allthrustcoeff.tolist())
allthrust = np.hstack(allthrust.tolist())    
allcurrent = np.hstack(allcurrent.tolist()) 

alltorque = prop_torque_coeff*allomega*allomega
allkv = alltorque/allcurrent

omegarange = allomega[np.where(np.logical_and(allthrust*1000/9.81<meanend,allthrust*1000/9.81>meanstart))]
voltagerange = allvoltage[np.where(np.logical_and(allthrust*1000/9.81<meanend,allthrust*1000/9.81>meanstart))]
thrustcoeffrange = allthrustcoeff[np.where(np.logical_and(allthrust*1000/9.81<meanend,allthrust*1000/9.81>meanstart))]
kvrange = allkv[np.where(np.logical_and(allthrust*1000/9.81<meanend,allthrust*1000/9.81>meanstart))]

meanthrustcoeff = np.mean(thrustcoeffrange)
meankv = np.mean(kvrange)

allresistance = 1/allcurrent*(allvoltage-meankv*allomega)
resistancerange = allresistance[np.where(np.logical_and(allthrust*1000/9.81<meanend,allthrust*1000/9.81>meanstart))]
meanresistance = np.mean(resistancerange)

rows = 4
fig1 = plt.figure(1)
fig1.clf()

fig1.add_subplot(rows,1,1)
plt.plot(allvoltage,allkv,'bo',label="Cloud of kv",markersize = 3)
plt.hlines(meankv,min(allvoltage),max(allvoltage),'r',label='Mean kv')
plt.vlines(min(voltagerange),min(allkv),max(allkv),'k',label = 'Range for mean')
plt.vlines(max(voltagerange),min(allkv),max(allkv),'k')
plt.ylabel('kv [Nm/A]')
plt.xlabel('Voltage [V]')

fig1.add_subplot(rows,1,2)
plt.plot(allvoltage,allresistance,'bo',label="Cloud of Resistance",markersize = 3)
plt.hlines(meanresistance,min(allvoltage),max(allvoltage),'r',label='Mean Resistance')
plt.vlines(min(voltagerange),min(allresistance),max(allresistance),'k',label = 'Range for mean')
plt.vlines(max(voltagerange),min(allresistance),max(allresistance),'k')
plt.ylabel('Resistance [Ohm]')
plt.xlabel('Voltage [V]')

fig1.add_subplot(rows,1,3)
plt.plot(allomega,allthrustcoeff,'bo',label = "Cloud of coeffs", markersize = 3)
plt.hlines(meanthrustcoeff,min(allomega),max(allomega),'r',label = "Average Thrustcoeff (Force to w²)")
plt.vlines(min(omegarange),min(allthrustcoeff),max(allthrustcoeff),label = "Average bounds")
plt.vlines(max(omegarange),min(allthrustcoeff),max(allthrustcoeff))
plt.ylabel('Thrust coefficient [(rad²/s²)/N]')
plt.xlabel('Omega [rad/s]')

voltagespace = np.arange(min(allvoltage),max(allvoltage),0.01)

fig1.add_subplot(rows,1,4)
plt.plot(allvoltage,allthrust*1000/9.81,'bo',label = "Cloud of thrusts", markersize = 3)
plt.plot(voltagespace,force(voltagespace)*1000/9.81,'r',label = "Model", markersize = 3)
plt.ylabel('Thrust [g]')
plt.xlabel('Voltage [V]')

print('\nAverage hover voltage: ' + str(np.mean(allhovervoltage)))
print('Average hover speed: ' + str(np.mean(allhoveromega)))
print('Mean thrust coefficient: ' + str(meanthrustcoeff))
print('Mean resistance: ' + str(meanresistance))
print('Mean kv (in V/(rad/s)): ' + str(meankv))
print('Mean kv (in RPM/V): ' + str(1/meankv*60/2/np.pi))


for a in fig1.axes:
    a.grid(True)
    a.legend()
