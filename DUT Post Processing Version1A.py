"""Created on Thu Sep 7 13:26:20 2017"""
"""***************************************************************************************************************************************************
Data: 02-11-2017                                                                                                                                     *
author: N.Mkhabela                                                                                                                                   *                                                                
Version: 1A                                                                                                                                          *                                                                                                                 *
This script serves as assistance in obtaining the calibrated power radiated from device tested in the SKA SA reverevertion chamber in dBm and dBuV/m *  
***************************************************************************************************************************************************"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#Filenames
DUTname = 'TwispCigerette.csv'
BASEname = 'TwispCigerette_Baseline.csv'
Cablename = 'CableLoss_Pinelands_6a_6b.csv'
CCFname = 'IL_03_06_15.csv'
 
"""Frequency band """
NumberPoints = 10401
StartFreq = 80000000
StopFreq  = 6000000000
StepSize = (StopFreq - StartFreq)/(NumberPoints -1)
frequency = range(StartFreq,StopFreq+1,int(StepSize))
Frequency=(np.array(frequency))/1e6

"""All losses and Gains during the reverberation chamber measurements"""
CableLoss = np.genfromtxt(Cablename, unpack = True, delimiter = ',', skip_header = 4)      # Cable losses
NewCableLoss = np.interp(frequency,CableLoss[0],CableLoss[1])                              # interpalating the values to correspond with above frequency
InsertionLoss = np.genfromtxt(CCFname, unpack = True, delimiter = ',', skip_header = 8)    # chamber calibration factor(CCF)/Insertion losses
NewInsertionLoss = np.interp(frequency,InsertionLoss[0],InsertionLoss[1])                  # interpalating the values to correspond with above frequency                                                          
RChamber = 10*np.log10(NewInsertionLoss)
Gantenna = 10*np.log10(0.75)                                                               # 0.75 is the Antenna's efficiency  

""" Calibrated Power transmitted from device""" 
Data = np.genfromtxt(DUTname,unpack = True,skip_header = 73 )                              # Data from  Deive
data = np.genfromtxt(BASEname,unpack = True,skip_header = 73)                              # Data from reverberation chamber's environment                 
PtxDevice = Data + NewCableLoss - Gantenna - RChamber                                      # calibrated power  transmitted by device     
PtxBaseline = data + NewCableLoss - Gantenna - RChamber                                    # calibrated power from reverberation chamber's environment           
"""Changing the units from dBm to dBuV/m"""                                                # conversion formular Efield(dBmuV/m) = Transmitted Power(dBm) + 10log(Zo/4pir^2) + 90
r =10                                                                                      # Separation distance (taken as 10 meters away from receiving antenna)
Zo =377                                                                                    # Free Space Impendance
ConvFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90                                       # conversion factor (from dBm to dBuV/m)
EField_Device = PtxDevice + ConvFactor                                                     # Electric field radiated by device
EField_Baseline = PtxBaseline  + ConvFactor                                                # Electric field radiated by reverberation chamber's environment

"""Figure will display the Power(dBm)&Efield(dBuV/m) of the device vs. the reverberation chamber's at different frequencies""" 
plt.figure() 
ax1=plt.subplot(2,1,1)                              
plt.plot(Frequency,PtxDevice,'k',Frequency,PtxBaseline,'b')
plt.title('Twisp Aero X E-Cigarette',fontsize=12),plt.xlabel('Frequency [MHz]',fontsize=8),plt.ylabel('Power [dBm]',fontsize=8)	
plt.legend(['Twisp E-Cigarette radiated emissions','Reverberation Chamber/Baseline'],loc =1,fontsize=8),plt.grid(True,which="both",ls="--"),plt.ylim(-130,-60)
ax2=plt.subplot(2,1,2) 
plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b')
plt.xlabel('Frequency [MHz]',fontsize=8),plt.ylabel('Electric Field Strength\n [dBuV/m]',fontsize=8)
plt.grid(True,which="both",ls="--"),plt.ylim(-40,20),plt.subplots_adjust(hspace=0.4)
plt.savefig('Twisp Aero X E-Cigarette.png',bbox_inches='tight')
plt.show()
     

