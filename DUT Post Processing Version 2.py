"""**************************************************************************************************************************************************************************************************************
 Date : 30/10/2017                                                                                                                                                 
 Author: nmkhabela                                                                                                                  
 Version: 2                                                                                                                         
 This script is for post processing the DUT's data from the reverberation chamber's spectrum analyszer and for required path losses 
 and the shielding effectiveness using CISPR and SARAS levels.                                                                      
**************************************************************************************************************************************************************************************************************"""
import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import messagebox
import csv

#Filenames
DUTname = 'UAV transmitter.csv'
BASEname = 'BASELI.csv'
Cablename = 'CableLoss_Pinelands_6a_6b.csv'
CCFname = 'IL_03_06_15.csv'
 
print('Script to process RFI data starting...')

"""Frequency band """
NumberPoints = 10401
StartFreq = 80000000
StopFreq  = 6000000000
StepSize = (StopFreq - StartFreq)/(NumberPoints -1)
frequency = range(StartFreq,StopFreq+1,int(StepSize))
Frequency=(np.array(frequency))/1e6

"""All losses and Gains during the reverberation chamber measurements"""
CableLoss = np.genfromtxt(Cablename, unpack = True, delimiter = ',', skip_header = 4)     # Cable losses csv file
NewCableLoss = np.interp(frequency,CableLoss[0],CableLoss[1])                             # interpalating values to correspond with above frequency
InsertionLoss = np.genfromtxt(CCFname, unpack = True, delimiter = ',', skip_header = 8)   # chamber calibration factor ((CCF)/Insertion losses csv file) 
NewInsertionLoss = np.interp(frequency,InsertionLoss[0],InsertionLoss[1])                                                                        
RChamber = 10*np.log10(NewInsertionLoss)                                                  # changing interpolated values from linear to log                                                                                            
Gantenna = 10*np.log10(0.75)                                                              # 0.75 is the Antenna's efficiency

""" Calibrated Data (Power) transmitted from device"""
Data = np.genfromtxt(DUTname,unpack = True,skip_header = 77 )                             # Raw data from Deive
data = np.genfromtxt(BASEname,unpack = True,skip_header = 77)                             # Raw data from reverberation chamber's environment                      
PtxDevice = Data + NewCableLoss - Gantenna - RChamber                                     # calibrated data  transmitted by device     
PtxBaseline = data + NewCableLoss - Gantenna - RChamber                                   # calibrated data from reverberation chamber's environment           

"""Changing the units from dBm to dBuV/m"""                                                # conversion formular Efield = Transmitted Power + 10log(Zo/4pir^2) + 90
r =10                                                                                      # Separation distance (taken as 10 meters away from receiving antenna)
Zo =377                                                                                    # Free Space Impendance
ConvFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90                                       # conversion factor (from dBm to dBuV/m)
EField_Device = PtxDevice + ConvFactor                                                     # Electric field radiated by device
EField_Baseline = PtxBaseline  + ConvFactor                                                # Electric field radiated by reverberation chamber's environment


'''Figure will display the chamber calibration factor (CCF) and cable losses'''
plt.figure(1)
plt.subplot(221)
plt.plot(Frequency,NewCableLoss,'k',Frequency,RChamber,'r') 
plt.title('CCF and Cable Losses',fontsize=12)
plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8)
plt.ylabel('Losses [dB]',fontweight='bold',fontsize=8)
plt.legend(['Cable Losses','Insertion Losses/Chamber Calibration Factor'],loc =1,fontsize=8)
plt.grid(True,which="both",ls="--")
#plt.ylim(-40,10)
plt.subplots_adjust(hspace=0.3)


'''Figure will display Raw Radiated Emissions of the device and the reverberation chambers environment at different frequencies''' 
plt.subplot(222)
plt.plot(Frequency,Data,'k',Frequency,data,'b')
plt.title('Raw Data of Device Tested',fontsize=12)
plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize = 8)
plt.ylabel('Power [dBm]',fontweight='bold',fontsize=8)
plt.legend(['Device Under Test','Reverberation Chamber/Baseline'],loc=1,fontsize=8,)
plt.grid(True,which="both",ls="--")
plt.ylim(-115,-70),plt.subplots_adjust(hspace=0.3)

'''Figure will display the calibrated Power(dBm) of the device and the reverberation chambers environment at different frequencies''' 
plt.subplot(223)                                           
plt.plot(Frequency,PtxDevice,'k',Frequency,PtxBaseline,'b')
plt.title('Calibrated Data of Device Tested',fontsize=12)
plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8)
plt.ylabel('Power [dBm]',fontweight='bold',fontsize=8)
plt.legend(['Device Under Test','Reverberation Chamber/Baseline'],loc=1,fontsize=8)
plt.grid(True,which="both",ls="--"),plt.ylim(-110,-50)
plt.subplots_adjust(hspace=0.3)


'''Figure will display the electric field strength of the device and the reverberation chambers environment at different frequencies'''
plt.subplot(224) 
plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b')
plt.title('Calibrated Data of Device Tested',fontsize=12)
plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8)
plt.ylabel('Electric Field Strength [dBuV/m]',fontweight='bold',fontsize=8)
plt.legend(['Device Under Test','Reverberation Chamber/Baseline'],loc =1,fontsize=8)
plt.grid(True,which="both",ls="--"),plt.ylim(-25,35),plt.subplots_adjust(hspace=0.3)
plt.show()

F = csv.writer(open('Radiated Emissions.csv','w', newline='\n'),delimiter=';')
Num = zip(*[Frequency,PtxDevice,EField_Device])
F.writerow(['Radiated Emissions from DUT'])
F.writerow(['Freq[MHz]','P[dBm]','EField[dBuV/m]'])
F.writerows(Num)

"""*****************************************************************************************************************************************************************************************************************""" 
"""*****************************************************************************************************************************************************************************************************************"""
#sarasy = input("Plot the SARAS limits (y/n)?")
sarasy = messagebox.askyesno("Title","Do you want to plot the SARAS limits?")
if(sarasy):
	print("SARAS is being calculated and plotted")
else:
	print("SARAS is NOT going to be calculated or plotted")
	
"""SARAS (South African Radio Astronomy Services) protection limit"""
SARAS = []              
for values in Frequency:
    if values < 2000:
        SARASdBmHz = -17.2708*np.log10(values) - 192.0714                 #for f<2GHz,  frequency units in MHz
    else:
        SARASdBmHz = -0.065676*np.log10(values) - 248.8661                #for f>2GHz,  frequency units in MHz                                   
    SARAS.append(SARASdBmHz)
    SARAS_= np.array([SARAS])                                           	                                          

if(sarasy):
	"""Converting SARAS protection limit from dBm/Hz to dBm"""
	RBW = 10*np.log10(0.01*Frequency*(1000000))                                #SKA SA Continuum Detrimental Threshold Levels RBW equal to 1% of observing frequency (for Spectral Detrimental Threshold Levels RBW equal to 0.001% of observing frequency)
	SARASdBm = SARAS_ + RBW                                                    #SARAS(dBm)=SARAS(dBm/Hz) +10log(RBW) for RBW in Hz
	SARASlevel = np.reshape(np.transpose(SARASdBm),NumberPoints)              
	"""Required Path Loss in dB""" 
	PathLossdB = PtxDevice-SARASdBm                                          #PathLoss(dB) = Power(dBm) - SARAS(dBm)
	PathLoss = np.reshape(np.transpose(PathLossdB),NumberPoints)
	""" Free Space Path Loss in dB"""
	d = 0.05                                                                   #Separation distance from transmitting device to nearest antenna in Km
	FSPL = 20*np.log10(Frequency) + 20*np.log10(d) + 32.45                     # freeSpacePathLoss = 20log(f) + 20log(d) + 32.45 for  f in MHz, d in Km
	"""Shielding Required for the device under test"""
	Shielding = np.reshape(np.transpose(PathLossdB - FSPL),NumberPoints)                                              # from the separation distance chosen from above, the Device under test will need shielding of this value in dB 

	'''                      
	"""Figure will display SARAS protection limit"""
	plt.subplot(2,2,1)
	plt.plot(Frequency,SARAS,'k',linewidth =1,label='SARAS Levels')
	plt.title('South African Radio Astronomy Services Protection Levels',fontsize=12),plt.legend(loc=1,fontsize=8),plt.grid(True,which="both",ls="--"),plt.semilogx()
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8),plt.ylabel('Power Spectral Density [dBm/Hz]',fontweight='bold',fontsize=8),plt.ylim(-260,-200),plt.subplots_adjust(hspace=0.3)

	plt.subplot(2,2,2)
	plt.plot(Frequency,SARASlevel,'k',linewidth =1,label='Continuum(SARAS)')
	plt.title('South African Radio Astronomy Services Protection Levels',fontsize=12),plt.legend(loc=1,fontsize=8),plt.grid(True,which="both",ls="--")
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8),plt.ylabel('Power [dBm]',fontweight='bold',fontsize=8),plt.ylim(-180,-160),plt.subplots_adjust(hspace=0.3)
	'''
	
	plt.figure(2)
	plt.subplot(2,1,1)
	plt.plot(Frequency,PathLoss,'k',Frequency,FSPL,'r')
	plt.legend(['Required Path Loss','Free Space Path Loss'],loc=1,fontsize=8)
	plt.title('Path Losses at 50m ',fontsize=12)
	plt.grid(True,which="both",ls="--")
	plt.subplots_adjust(hspace=0.5)
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8)
	plt.ylabel('Losses [dB]',fontweight='bold',fontsize=8)
	#plt.ylim(40,200)

	"""Figure will display shielding required for the Device Under Test"""
	plt.subplot(2,1,2)
	plt.plot(Frequency,Shielding,'k')
	plt.title('Shielding required for the DUT at 50m',fontsize=12)
	#plt.legend(loc=1,fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.subplots_adjust(hspace=0.5)
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8)
	plt.ylabel('Sheilding [dB]',fontweight='bold',fontsize=8),
	#plt.ylim(-25,25)
	plt.show()

	F = csv.writer(open('SARAS levels(DUT).csv','w', newline=''),delimiter=';')
	Values =zip(*[Frequency,SARAS,SARASlevel,PathLoss,FSPL,Shielding])
	F.writerow(['SARAS Protection levels at 50m Separation Distance'])
	F.writerow(['Freq[MHz]','SARAS[dBm/Hz]','SARAS[dBm]','Required PathLoss[dB]','Free Space PathLoss[dB]','Shielding[dB]'])
	F.writerows(Values)


"""*****************************************************************************************************************************************************************************************************************"""
"""*****************************************************************************************************************************************************************************************************************"""             
#cispray = input("Plot the CISPR-22 Class A Radiated Emission Limits(y/n)?")
cispray = messagebox.askyesno("Post Processing RFI data","Do you want to plot the CISPR-22 Class A Radiated Emission Limits")
if(cispray):
	print("CISPR-22 Class A is being calculated and plotted")
else:
	print("CISPR-22 Class A is NOT going to be calculated or plotted")
	
if(cispray):
	"""CISPR-22 Class A  Radiated Emission Limits at r = 10m"""                                #frequencies 30MHz to 230MHz and 230MHz to 1GHz for 10m separation distance     
	A10m = []
	def H(Frequency):
		return (40 if 30<Frequency<230 else 47 if Frequency <=1000 else np.nan)        
	for i in range(len(Frequency)):
		A10m.append(H(Frequency[i]))                   
	"""CISPR-22 Class B  Radiated Emission Limits at r = 10m"""                               
	B10m = []
	def H(Frequency):
		return (30 if 30<Frequency<230 else 37 if Frequency <=1000 else np.nan)        
	for i in range(len(Frequency)):
		B10m.append(H(Frequency[i]))
	"""CISPR-22 Class A  Radiated Emission Limits at r = 3m"""
	A3m = []
	def H(Frequency):
		return (76 if 1000<Frequency<3000 else 80 if Frequency >=3000 else np.nan)        
	for i in range(len(Frequency)):
		A3m.append(H(Frequency[i]))                   
	"""CISPR-22 Class B  Radiated Emission Limits at r = 3m"""                               #frequencies 1GHz to 3GHz and 3GHz to 6GHz for 3m separation distance
	B3m = []
	def H(Frequency):
		return (70 if 1000<Frequency<3000 else 74 if Frequency >=3000 else np.nan)        
	for i in range(len(Frequency)):
		B3m.append(H(Frequency[i]))   
		
	"""figure will display the CISPR-22 standards Class A and Class B at different frequencies and different distances""" 
	#print("Show CISPR-22 plots")
	plt.figure(3)
	plt.subplot(2,1,1)    
	plt.plot(Frequency,A10m,'k',Frequency,B10m,'b') 
	plt.title('CISPR-22 Radiated Emission Limits at 10m',fontsize = 20)
	plt.xlabel('Frequency[MHz]',fontweight='bold',fontsize=12)
	plt.ylabel('Electric Field Strength [dBuV/m]',fontsize=12)
	plt.ylim(20,60)
	plt.subplots_adjust(hspace=0.4) 
	plt.legend(['Class A  Radiated Emission Limits','Class B  Radiated Emission Limits'],loc =1,fontsize=8)
	plt.grid(True,which="both",ls="--")

	plt.subplot(2,1,2)
	plt.plot(Frequency,A3m,'k',Frequency,B3m,'b') 
	plt.title('CISPR-22 Radiated Emission Limits at 3m',fontsize = 20)
	plt.xlabel('Frequency[MHz]',fontweight='bold',fontsize=12)
	plt.ylabel('Electric Field Strength [dBuV/m]',fontsize=12)
	plt.ylim(60,90)
	plt.subplots_adjust(hspace=0.4)
	plt.legend(['Class A  Radiated Emission Limits','Class B  Radiated Emission Limits'],loc =1,fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.show()
	


