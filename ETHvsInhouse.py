#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:28:14 2018

@author: christian
"""
import load_cell_code
import LoadCellDataMuncher as lcdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob

Rcable = 0.144

files = sorted(glob.glob("*.xlsx"))


for idx,xlsx in enumerate(files):
    print(xlsx)
    
    file = pd.read_excel(xlsx)
    
    voltage = file.Voltage
    current = file.Current
    voltage = voltage - current*Rcable
    thrust = file.Thrust
    
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
    
allomega = np.hstack(allomega.tolist())
allvoltage = np.hstack(allvoltage.tolist())
allthrustcoeff = np.hstack(allthrustcoeff.tolist())
allthrust = np.hstack(allthrust.tolist())
allcurrent = np.hstack(allcurrent.tolist())

base_path = '/home/'

directory_list = ('/home/christian/Desktop/load_cell/lucie/Load Cell Data/20160901_ManualMiniMotNewSetup/20160901t125034_0720-14_Parrot-BSamp_manual',)
labels = ('0720-14_Parrot',)

muncher = lcdm.LoadCellDataMuncher(directory_list, base_path, labels)

thrustETH = muncher.data[0].thrust
torqueETH = muncher.data[0].torque
currETH = muncher.data[0].curr
voltETH = muncher.data[0].volt
voltETH = voltETH

rpmETH = np.array([ 5238,
        9496,
        12623,
        14989,
        18234,
        20695,
        22425,
        23142
        ])
omegaETH = rpmETH/60*2*np.pi

densityETH = np.mean(muncher.data[0].density)
densityStandard = 1.204
ETHtoStandard = densityStandard/densityETH
InhousetoStandard = densityStandard/1.1894260886195733

rows = 4
fig1 = plt.figure(1)
fig1.clf()
fig1.add_subplot(rows,1,1)
plt.plot(allvoltage,allthrust*1000/9.81*InhousetoStandard,'bo',label="Cloud of Thrust (In-house), subtracted 0.144 Ohm",markersize = 3)
plt.plot(voltETH,thrustETH*1000/9.81*ETHtoStandard,'ro',label="Thrust (ETH), subtracted nothing",markersize = 3)
plt.xlabel('Voltage [V]')
plt.ylabel('Thrust [g]')
plt.legend()



allvoltage = allcurrent*Rcable + allvoltage

fig1.add_subplot(rows,1,2)
plt.plot(allvoltage,allthrust*1000/9.81*InhousetoStandard,'bo',label="Cloud of Thrust (In-house), subtracted nothing",markersize = 3)
plt.plot(voltETH,thrustETH*1000/9.81*ETHtoStandard,'ro',label="Thrust (ETH), subtracted nothing",markersize = 3)
plt.xlabel('Voltage [V]')
plt.ylabel('Thrust [g]')
plt.legend()

voltETH = sorted(voltETH)

allvoltage = -allcurrent*Rcable + allvoltage

fig1.add_subplot(rows,1,3)
plt.plot(allvoltage,allomega,'bo',label="Cloud of Speed (In-house), subtracted 0.144 Ohm",markersize = 3)
plt.plot(voltETH,omegaETH,'ro',label="Speed (ETH), subtracted nothing",markersize = 3)
plt.xlabel('Voltage [V]')
plt.ylabel('Speed [rad/s]')
plt.legend()

allvoltage = allcurrent*Rcable + allvoltage

fig1.add_subplot(rows,1,4)
plt.plot(allvoltage,allomega,'bo',label="Cloud of Speed (In-house), subtracted nothing",markersize = 3)
plt.plot(voltETH,omegaETH,'ro',label="Speed (ETH), subtracted nothing",markersize = 3)
plt.ylabel('Speed [rad/s]')
plt.xlabel('Voltage [V]')
plt.legend()

fit = np.polyfit(omegaETH**2,sorted(torqueETH),1)

#fig1.add_subplot(rows,1,5)
#plt.plot(omegaETH**2,sorted(torqueETH)/omegaETH/omegaETH,'bo')
