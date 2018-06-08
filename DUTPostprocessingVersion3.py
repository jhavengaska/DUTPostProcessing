#!/usr/bin/env python3
"""**************************************************************************************************************************************************************************************************************
 Date : 09/04/2018                                                                                                                                                
 Author: jhavenga                                                                                                                  
 Version: 3                                                                                                                         
 This script is for post processing the DUT's data from the reverberation chamber's spectrum analyszer and for required path losses 
 and the shielding effectiveness using CISPR and SARAS levels.           
 This script builds on the scripts created by nmkhabela 
**************************************************************************************************************************************************************************************************************"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import tkinter
from tkinter import filedialog
from tkinter import messagebox
import matplotlib.pyplot as plt
import sys
import csv



'''def calcEfield(Prsa, Lcable, Linsertion, Glna, AFact):
	""" This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB), the reverberation chamber calibration factor(dB),
	the gain from the LNA(dB) and the antenna factor(dB/m) to caluclate the E-field at
	the receiving antenna."""
	Vrsa = Prsa + 106.9897 #Gives the answer in dBuv	
	Linsertion = Linsertion * -1.00 #Given as a loss in negative dB of the ACF of the Chamber
	Vd = Vrsa + Lcable + Linsertion #Gives the answer in dBuV of the voltage just before the LNA
	Vlna = Vd - Glna #Gives the answer in dBuV just before the antenna and after the LNA
	E = Vlna + AFact #Factor in the Antenna Factor
	return(E)'''

# def GetHigh(FreqMHZ, PWR, EFS, FStart, FStop):
	# """This function takes the frequency dataset, in MHz, the EIRP of the DUT in dBm
	# and the Electric Field Strength in dbUv/m and calculates the highest 5
	# measured radiations starting from FStart, in MHz, and ending in FStop, in MHZ."""
	# N = 5
	# Max = [-20,-20,-20]
	# Max = [-20] * N
	# Maxp = [-200,-200,-200]
	# Maxp = [-200] * N
	# Mfreq = [0,0,0]
	# Mfreq = [0] * N
	# for i, val in enumerate(FreqMHz):
		# if Fstart <= values <= FStop:
			# #This is for HERA
			# if EFS[i] > Max[0]:
				# for s in range(N-1,0,-1):
					# Max[s] = Max[(s-1)]
					# Maxp[s] = Maxp[(s-1)]
					# Hfreq[s] = Hfreq[(s-1)]
				# Max[0] = EFS[i]
				# Maxp[0] = PWR[i]
				# Hfreq[0] = val
			# elif EField_Device[ind] > Hmax[1]:
				# Hmax[2] = Hmax[1]
				# Hmaxp[2] = Hmaxp[1]
				# Hfreq[2] = Hfreq[1]
				# Hmax[1] = EField_Device[ind]
				# Hmaxp[1] = PtxDevice[ind]
				# Hfreq[1] = values
			# elif EField_Device[ind] > Hmax[2]:
				# Hmax[2] = EField_Device[ind]
				# Hmaxp[2] = PtxDevice[ind]
				# Hfreq[2] = values
				
# def max_indices(arr, k):
    # """
    # Returns the indices of the k first largest elements of arr
    # (in descending order in values)
    # """
    # assert k <= arr.size, 'k should be smaller or equal to the array size'
    # arr_ = arr.astype(float)  # make a copy of arr
    # max_idxs = []
    # for _ in range(k):
        # max_element = np.max(arr_)
        # if np.isinf(max_element):
            # break
        # else:
            # idx = np.where(arr_ == max_element)
        # max_idxs.append(idx)
        # arr_[idx] = -np.inf
    # return max_idxs

def calcEfieldPow(Prsa, Lcable, Linsertion, Glna, Aeff, r):
	""" This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB), the reverberation chamber calibration factor(dB),
	the gain from the LNA(dB) and the frequencies of operation(Hz) to caluclate the E-field at
	the receiving antenna."""
	eps0 = 8.854E-12
	mu0 = 4*np.pi*1E-7
	#c0 = 1./np.sqrt(eps0 * mu0)
	Zo= np.sqrt(mu0/eps0)
	#Zo =377                                                                                    # Free Space Impendance
	Linsertion = (Linsertion - 10*np.log(Aeff)) * -1.00 #Given as a loss in negative dB of the ACF of the Chamber
	Plessloss = Prsa + Lcable + Linsertion - Glna 
	CFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90   #Conversion factor (from dBm to dBuV/m)
	#print(CFactor)
	Efield = Plessloss + CFactor
	return(Efield, Plessloss)

# def EftoEIRP(Ef, d):
	# eirp = Ef + 20*np.log10(d) - 104.68
	# return(eirp)

# >>>>>>>>>>>>>> THIS IS THE PART YOU MUST/CAN EDIT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Remember to edit these according to what you are testing
DeviceName = "Getac Tablet"
outputfilename = "Getac Tablet"
NumberPoints = 64001 #Maybe enter this as a GUI?
StartFreq = 80000000 #Maybe enter this as a GUI?
StopFreq  = 6000000000 #Maybe enter this as a GUI?
BASEname = "Baseline 100kHz 30.csv"	#The filename of the CSV that contains the baseline measurements
DUTname = "Getac Tablet 100kHz 30.csv"	#The filename of the CSV that contains the DUT measurements
AntEff = 0.8 #Antenna Efficiency, used for E-Field calculations
Dist = 10 #Distance from DUT for E-Field calculations (m)

#File names for cable, antenna and chamber factors
Cablename = 'CableLoss_Pinelands_6a_6b.csv'
CCFname = 'ACF_03_06_15.csv'
AntName = 'Antenna_MESA_GLPDA.csv'
LNAName = 'No_LNA.csv'
#LNAName = 'LNA_3.csv'
#sarasy = 1 #Plot SARAS limits (1 = Yes, 0 = No)
#cispray = 1 #Plot CISPR limits (1 = Yes, 0 = No)

# >>>>>>>>>>>>>> THIS IS THE PART YOU MUST/CAN EDIT ENDS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#"""ADD the DUT and Baseline files with an file dialog"""
#root = tkinter.Tk()
#root.withdraw()
#BASEname = filedialog.askopenfilename(title = "Select BASELINE file")
#if(BASEname == ""):
#	print("BASELINE file not selected")
#	sys.exit()
#DUTname = filedialog.askopenfilename(title = "Select DUT file")
#if(DUTname == ""):
#	print("DUT file not selected")
#	sys.exit()
	
savefilename = outputfilename+"_dbm"+".png"
savefilename2 = outputfilename+"_efield"+".png"
csvfilename = outputfilename+".csv"
	
"""Frequency band """
StepSize = (StopFreq - StartFreq)/(NumberPoints -1)
frequency = range(StartFreq,StopFreq+1,int(StepSize))
Frequency=(np.array(frequency))/1e6

"""All losses and Gains during the reverberation chamber measurements"""
CableLoss = np.genfromtxt(Cablename, unpack = True, delimiter = ',', skip_header = 4)     # Cable losses csv file
NewCableLoss = np.interp(frequency,CableLoss[0],CableLoss[1])                             # interpalating values to correspond with above frequency
InsertionLoss = np.genfromtxt(CCFname, unpack = True, delimiter = ',', skip_header = 8)   # chamber calibration factor ((CCF)/Insertion losses csv file) 
NewInsertionLoss = np.interp(frequency,InsertionLoss[0],InsertionLoss[1])                                                                        
RChamber = 10*np.log10(NewInsertionLoss)                                                  # changing interpolated values from linear to log                                                                                            
Antenna = np.genfromtxt(AntName, unpack = True, delimiter = ',', skip_header = 2)	  #Read antenna csv file
AF = np.interp(frequency,Antenna[0], Antenna[3])					  #Interpolate to get Antenna Factor
AGain = np.interp(frequency,Antenna[0], Antenna[1])					  #Interpolate to get Antenna Facto
LNAGainf = np.genfromtxt(LNAName, unpack = True, delimiter = ',', skip_header = 2)
LNAGain = np.interp(frequency,LNAGainf[0], LNAGainf[1])
Gantenna = 10*np.log10(0.75)                                                              # 0.75 is the Antenna's efficiency
#print(RChamber.shape) 
'''
plt.figure()        
plt.plot(frequency,RChamber,'k')
plt.title("Insertion loss vs frequency",fontsize=12)
plt.xlabel('Frequency [Hz]',fontsize=8)
plt.ylabel('Insertion Loss [dB]',fontsize=8)	


plt.figure()        
plt.plot(frequency,NewCableLoss,'k')
plt.title("Cable loss vs frequency",fontsize=12)
plt.xlabel('Frequency [Hz]',fontsize=8)
plt.ylabel('Cable Loss [dB]',fontsize=8)	

#print(str(Gantenna))

plt.figure()        
plt.plot(frequency,AF,'k')
plt.title("Antenna Factor vs frequency",fontsize=12)
plt.xlabel('Frequency [Hz]',fontsize=8)
plt.ylabel('Antenna Factor [dB/m]',fontsize=8)	

plt.figure()        
plt.plot(frequency,AGain,'k')
plt.title("Antenna Gain vs frequency",fontsize=12)
plt.xlabel('Frequency [Hz]',fontsize=8)
plt.ylabel('Antenna Gain [dB]',fontsize=8)	

plt.figure()        
plt.plot(frequency,LNAGain,'k')
plt.title("LNA Gain vs frequency",fontsize=12)
plt.xlabel('Frequency [Hz]',fontsize=8)
plt.ylabel('LNA Gain [dB]',fontsize=8)	
'''

""" Calibrated Data (Power) transmitted from device"""
Data = np.genfromtxt(DUTname,unpack = True, skip_header = 77)            # Raw data from Deive
data = np.genfromtxt(BASEname,unpack = True, skip_header = 77)            # Raw data from reverberation chamber's environment                
#PtxDevice = Data + NewCableLoss - Gantenna - RChamber                                     # calibrated data  transmitted by device     
#PtxBaseline = data + NewCableLoss - Gantenna - RChamber                                   # calibrated data from reverberation chamber's environment           

# """ Calibrated data (Voltage) """
# EVolDevice = calcEfield(Data, NewCableLoss, RChamber, LNAGain, AF)

# plt.figure()        
# plt.plot(frequency,EVolDevice,'k')
# plt.title("E-Field according to voltage method vs frequency",fontsize=12)
# plt.xlabel('Frequency [Hz]',fontsize=8)
# plt.ylabel('E-Field [dBuV/m]',fontsize=8)

#"""Calibrated data (Power) """
#EPowDevice, PowPowDevice = calcEfieldPow(Data, NewCableLoss, RChamber, LNAGain)

# plt.figure()        
# plt.plot(frequency,EPowDevice,'k')
# plt.title("E-Field according to power method vs frequency",fontsize=12)
# plt.xlabel('Frequency [Hz]',fontsize=8)
# plt.ylabel('E-Field [dBuV/m]',fontsize=8)


# """Changing the units from dBm to dBuV/m"""                                                # conversion formular Efield = Transmitted Power + 10log(Zo/4pir^2) + 90
# r =10                                                                                      # Separation distance (taken as 10 meters away from receiving antenna)
# Zo =377                                                                                    # Free Space Impendance
# ConvFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90                                       # conversion factor (from dBm to dBuV/m)
#EField_Device = PtxDevice + ConvFactor                                                     # Electric field radiated by device
#EField_Baseline = PtxBaseline  + ConvFactor                                                # Electric field radiated by reverberation chamber's environment

"""Calibrated data (Power) """
EField_Device, PtxDevice = calcEfieldPow(Data, NewCableLoss, RChamber, LNAGain, AntEff, Dist)
EField_Baseline, PtxBaseline = calcEfieldPow(data, NewCableLoss, RChamber, LNAGain, AntEff, Dist)

#Get min values for plots
minE = min(EField_Baseline)
if(min(EField_Device) < minE):
	minE = min(EField_Device)
minE = minE - 7	

minP = min(data)
if(min(Data) < minP):
	minP = min(Data)
minP = minP - 7	

minTxP = min(PtxBaseline)
if(min(PtxDevice) < minTxP):
	minTxP = min(PtxDevice)
minTxP = minTxP - 7

#Get max values for plots
maxE = max(EField_Baseline)
if(max(EField_Device) > maxE):
	maxE = max(EField_Device)
maxE = maxE + 7

maxP = max(data)
if(max(Data) > maxP):
	maxP = max(Data)
maxP = maxP + 20	

maxTxP = max(PtxBaseline)
if(max(PtxDevice) > maxTxP):
	maxTxP = max(PtxDevice)
maxTxP = maxTxP + 20
#print(str(minE)+" and " + str(maxE))
#print(str(minP)+" and " + str(maxP))

"""Figure will display the Power(dBm)&Efield(dBuV/m) of the device vs. the reverberation chamber's at different frequencies""" 
plt.figure() 
ax1=plt.subplot(2,1,1)                              
titel = DeviceName + " Spectrum Analyzer Raw Data"
plt.plot(Frequency,Data,'k',Frequency,data,'b')
plt.title(titel,fontsize=12)
plt.xlabel('Frequency [MHz]',fontsize=8)
plt.ylabel('Power [dBm]',fontsize=8)	
plt.legend([DeviceName,'Reverberation Chamber/Baseline'],loc =1,fontsize=8)
plt.grid(True,which="both",ls="--")
plt.ylim(minP,maxP)
ax2=plt.subplot(2,1,2) 
plt.plot(Frequency,PtxDevice,'k',Frequency,PtxBaseline,'b')
titel = DeviceName + " EIRP"
plt.title(titel,fontsize=12)
plt.xlabel('Frequency [MHz]',fontsize=8)
plt.ylabel('Power [dBm]',fontsize=8)	
plt.legend([DeviceName,'Reverberation Chamber/Baseline'],loc =1,fontsize=8)
plt.grid(True,which="both",ls="--")
plt.ylim(minTxP,maxTxP)
plt.subplots_adjust(hspace=0.4)
plt.savefig(savefilename,bbox_inches='tight')
plt.show()

titel = DeviceName + " Electric Field Strength @ "+str(Dist)+"10m"
plt.figure()
plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b')
plt.title(titel,fontsize=12)
plt.xlabel('Frequency [MHz]',fontsize=8)
plt.ylabel('Electric Field Strength\n [dBuV/m]',fontsize=8)
plt.grid(True,which="both",ls="--")
plt.ylim(minE,maxE),                        
plt.savefig(savefilename2,bbox_inches='tight')
plt.show()

#Write CSV file with radiated files
#F = csv.writer(open('Radiated Emissions.csv','w', newline='\n'),delimiter=';')
#Num = zip(*[Frequency,PtxDevice,EField_Device])
#F.writerow(['Radiated Emissions from DUT'])
#F.writerow(['Freq[MHz]','P[dBm]','EField[dBuV/m]'])
#F.writerows(Num)


'''
"""*****************************************************************************************************************************************************************************************************************""" 
"""*****************************************************************************************************************************************************************************************************************"""
#sarasy = input("Plot the SARAS limits (y/n)?")
#print("Ask about SARAS")
#root = tkinter.Tk()
#root.withdraw()
#sarasy = messagebox.askyesno("Title","Do you want to plot the SARAS limits?")
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
	#"""Converting SARAS protection limit from dBm/Hz to dBm"""
	RBW = 10*np.log10(0.01*Frequency*(1000000))                                #SKA SA Continuum Detrimental Threshold Levels RBW equal to 1% of observing frequency (for Spectral Detrimental Threshold Levels RBW equal to 0.001% of observing frequency)
	#print("The RBW = "+str(RBW))
	#print("The SARAS = "+str(SARAS_))
	SARASdBm = SARAS_ + RBW                                                    #SARAS(dBm)=SARAS(dBm/Hz) +10log(RBW) for RBW in Hz
	SARASlevel = np.reshape(np.transpose(SARASdBm),NumberPoints)              
	#"""Required Path Loss in dB""" 
	PathLossdB = PtxDevice-SARASdBm                                          #PathLoss(dB) = Power(dBm) - SARAS(dBm)
	PathLoss = np.reshape(np.transpose(PathLossdB),NumberPoints)
	#""" Free Space Path Loss in dB"""
	d = 0.05                                                                   #Separation distance from transmitting device to nearest antenna in Km
	FSPL = 20*np.log10(Frequency) + 20*np.log10(d) + 32.45                     # freeSpacePathLoss = 20log(f) + 20log(d) + 32.45 for  f in MHz, d in Km
	#"""Shielding Required for the device under test"""
	Shielding = np.reshape(np.transpose(PathLossdB - FSPL),NumberPoints)       # from the separation distance chosen from above, the Device under test will need shielding of this value in dB 

	                      
	"""Figure will display SARAS protection limit"""
	plt.subplot(2,2,1)
	plt.plot(Frequency,SARAS,'k',linewidth =1,label='SARAS Levels')
	plt.title('South African Radio Astronomy Services Protection Levels',fontsize=12),plt.legend(loc=1,fontsize=8),plt.grid(True,which="both",ls="--"),plt.semilogx()
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8),plt.ylabel('Power Spectral Density [dBm/Hz]',fontweight='bold',fontsize=8),plt.ylim(-260,-200),plt.subplots_adjust(hspace=0.3)

	plt.subplot(2,2,2)
	plt.plot(Frequency,SARASlevel,'k',linewidth =1,label='Continuum(SARAS)')
	plt.title('South African Radio Astronomy Services Protection Levels',fontsize=12),plt.legend(loc=1,fontsize=8),plt.grid(True,which="both",ls="--")
	plt.xlabel('Frequency [MHz]',fontweight='bold',fontsize=8),plt.ylabel('Power [dBm]',fontweight='bold',fontsize=8),plt.ylim(-180,-160),plt.subplots_adjust(hspace=0.3)
	
	
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

	#"""Figure will display shielding required for the Device Under Test"""
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

	#F = csv.writer(open('SARAS levels(DUT).csv','w', newline=''),delimiter=';')
	#Values =zip(*[Frequency,SARAS,SARASlevel,PathLoss,FSPL,Shielding])
	#F.writerow(['SARAS Protection levels at 50m Separation Distance'])
	#F.writerow(['Freq[MHz]','SARAS[dBm/Hz]','SARAS[dBm]','Required PathLoss[dB]','Free Space PathLoss[dB]','Shielding[dB]'])
	#F.writerows(Values)


#"""*****************************************************************************************************************************************************************************************************************"""
#"""*****************************************************************************************************************************************************************************************************************"""             
#cispray = input("Plot the CISPR-22 Class A Radiated Emission Limits(y/n)?")
#root = tkinter.Tk()
#root.withdraw()
#cispray = messagebox.askyesno("Post Processing RFI data","Do you want to plot the CISPR-22 Class A Radiated Emission Limits")
if(cispray):
	print("CISPR-22 Class A is being calculated and plotted")
else:
	print("CISPR-22 Class A is NOT going to be calculated or plotted")
	
if(cispray):
	#"""CISPR-22 Class A  Radiated Emission Limits at r = 10m"""                                #frequencies 30MHz to 230MHz and 230MHz to 1GHz for 10m separation distance     
	A10m = []
	def H(Frequency):
		return (40 if 30<Frequency<230 else 47 if Frequency <=1000 else np.nan)        
	for i in range(len(Frequency)):
		A10m.append(H(Frequency[i]))                   
	#"""CISPR-22 Class B  Radiated Emission Limits at r = 10m"""                               
	B10m = []
	def H(Frequency):
		return (30 if 30<Frequency<230 else 37 if Frequency <=1000 else np.nan)        
	for i in range(len(Frequency)):
		B10m.append(H(Frequency[i]))
	#"""CISPR-22 Class A  Radiated Emission Limits at r = 3m"""
	A3m = []
	def H(Frequency):
		return (76 if 1000<Frequency<3000 else 80 if Frequency >=3000 else np.nan)        
	for i in range(len(Frequency)):
		A3m.append(H(Frequency[i]))                   
	#"""CISPR-22 Class B  Radiated Emission Limits at r = 3m"""                               #frequencies 1GHz to 3GHz and 3GHz to 6GHz for 3m separation distance
	B3m = []
	def H(Frequency):
		return (70 if 1000<Frequency<3000 else 74 if Frequency >=3000 else np.nan)        
	for i in range(len(Frequency)):
		B3m.append(H(Frequency[i]))   
		
	#"""figure will display the CISPR-22 standards Class A and Class B at different frequencies and different distances""" 
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
'''
#Add part that automatically gets the 5 top frequencies in each band(HERA, UHF, L, S, X and C-Bass)
Hmax = [-15,-15,-15]
Hmaxp = [-200,-200,-200]
Hfreq = [0,0,0]
N1max = [-15,-15,-15]
N1maxp = [-200,-200,-200]
N1freq = [0,0,0]
Umax = [-15,-15,-15]
Umaxp = [-200,-200,-200]
Ufreq = [0,0,0]
Lmax = [-15,-15,-15]
Lmaxp = [-200,-200,-200]
Lfreq = [0,0,0]
Smax = [-15,-15,-15]
Smaxp = [-200,-200,-200]
Sfreq = [0,0,0]
Cmax = [0,0,0]
Cmaxp = [-200,-200,-200]
Cfreq = [0,0,0]
N2max = [0,0,0]
N2maxp = [-200,-200,-200]
N2freq = [0,0,0]
Xmax = [0,0,0]
Xmaxp = [-200,-200,-200]
Xfreq = [0,0,0]

# indices = max_indices(EField_Device, 3)
# print(indices)
# print(str(EField_Device[indices[0][0]]))

for ind, values in enumerate(Frequency):
	if 100 <= values <= 200:
		#This is for HERA
		if EField_Device[ind] > Hmax[0]:
			Hmax[2] = Hmax[1]
			Hmaxp[2] = Hmaxp[1]
			Hfreq[2] = Hfreq[1]
			Hmax[1] = Hmax[0]
			Hmaxp[1] = Hmaxp[0]
			Hfreq[1] = Hfreq[0]
			Hmax[0] = EField_Device[ind]
			Hmaxp[0] = PtxDevice[ind]
			Hfreq[0] = values
		elif EField_Device[ind] > Hmax[1]:
			Hmax[2] = Hmax[1]
			Hmaxp[2] = Hmaxp[1]
			Hfreq[2] = Hfreq[1]
			Hmax[1] = EField_Device[ind]
			Hmaxp[1] = PtxDevice[ind]
			Hfreq[1] = values
		elif EField_Device[ind] > Hmax[2]:
			Hmax[2] = EField_Device[ind]
			Hmaxp[2] = PtxDevice[ind]
			Hfreq[2] = values
	if 200 < values < 580:
		#This is for None 1 (200 - 580)
		if EField_Device[ind] > N1max[0]:
			N1max[2] = N1max[1]
			N1maxp[2] = N1maxp[1]
			N1freq[2] = N1freq[1]
			N1max[1] = N1max[0]
			N1maxp[1] = N1maxp[0]
			N1freq[1] = N1freq[0]
			N1max[0] = EField_Device[ind]
			N1maxp[0] = PtxDevice[ind]
			N1freq[0] = values
		elif EField_Device[ind] > N1max[1]:
			N1max[2] = N1max[1]
			N1maxp[2] = N1maxp[1]
			N1freq[2] = N1freq[1]
			N1max[1] = EField_Device[ind]
			N1maxp[1] = PtxDevice[ind]
			N1freq[1] = values
		elif EField_Device[ind] > N1max[2]:
			N1max[2] = EField_Device[ind]
			N1maxp[2] = PtxDevice[ind]
			N1freq[2] = values
	if 580 <= values <= 1015:
		#This is for UHF
		if EField_Device[ind] > Umax[0]:
			Umax[2] = Umax[1]
			Umaxp[2] = Umaxp[1]
			Ufreq[2] = Ufreq[1]
			Umax[1] = Umax[0]
			Umaxp[1] = Umaxp[0]
			Ufreq[1] = Ufreq[0]
			Umax[0] = EField_Device[ind]
			Umaxp[0] = PtxDevice[ind]
			Ufreq[0] = values
		elif EField_Device[ind] > Umax[1]:
			Umax[2] = Umax[1]
			Umaxp[2] = Umaxp[1]
			Ufreq[2] = Ufreq[1]
			Umax[1] = EField_Device[ind]
			Umaxp[1] = PtxDevice[ind]
			Ufreq[1] = values
		elif EField_Device[ind] > Umax[2]:
			Umax[2] = EField_Device[ind]
			Umaxp[2] = PtxDevice[ind]
			Ufreq[2] = values
	if 900 <= values <= 1670:
		#This is for L-Band
		if EField_Device[ind] > Lmax[0]:
			Lmax[2] = Lmax[1]
			Lmaxp[2] = Lmaxp[1]
			Lfreq[2] = Lfreq[1]
			Lmax[1] = Lmax[0]
			Lmaxp[1] = Lmaxp[0]
			Lfreq[1] = Lfreq[0]
			Lmax[0] = EField_Device[ind]
			Lmaxp[0] = PtxDevice[ind]
			Lfreq[0] = values
		elif EField_Device[ind] > Lmax[1]:
			Lmax[2] = Lmax[1]
			Lmaxp[2] = Lmaxp[1]
			Lfreq[2] = Lfreq[1]
			Lmax[1] = EField_Device[ind]
			Lmaxp[1] = PtxDevice[ind]
			Lfreq[1] = values
		elif EField_Device[ind] > Lmax[2]:
			Lmax[2] = EField_Device[ind]
			Lmaxp[2] = PtxDevice[ind]
			Lfreq[2] = values
	if 1750 <= values <= 3500:
		#This is for S-Band
		if EField_Device[ind] > Smax[0]:
			Smax[2] = Smax[1]
			Smaxp[2] = Smaxp[1]
			Sfreq[2] = Sfreq[1]
			Smax[1] = Smax[0]
			Smaxp[1] = Smaxp[0]
			Sfreq[1] = Sfreq[0]
			Smax[0] = EField_Device[ind]
			Smaxp[0] = PtxDevice[ind]
			Sfreq[0] = values
		elif EField_Device[ind] > Smax[1]:
			Smax[2] = Smax[1]
			Smaxp[2] = Smaxp[1]
			Sfreq[2] = Sfreq[1]
			Smax[1] = EField_Device[ind]
			Smaxp[1] = PtxDevice[ind]
			Sfreq[1] = values
		elif EField_Device[ind] > Smax[2]:
			Smax[2] = EField_Device[ind]
			Smaxp[2] = PtxDevice[ind]
			Sfreq[2] = values
	if (4500 <= values <= 5500):
		#This is for C-BASS
		if(EField_Device[ind] > Cmax[0]):
			Cmax[2] = Cmax[1]
			Cmaxp[2] = Cmaxp[1]
			Cfreq[2] = Cfreq[1]
			Cmax[1] = Cmax[0]
			Cmaxp[1] = Cmaxp[0]
			Cfreq[1] = Cfreq[0]
			Cmax[0] = EField_Device[ind]
			Cmaxp[0] = PtxDevice[ind]
			Cfreq[0] = values
		elif(EField_Device[ind] > Cmax[1]):
			Cmax[2] = Cmax[1]
			Cmaxp[2] = Cmaxp[1]
			Cfreq[2] = Cfreq[1]
			Cmax[1] = EField_Device[ind]
			Cmaxp[1] = PtxDevice[ind]
			Cfreq[1] = values
		elif(EField_Device[ind] > Cmax[2]):
			Cmax[2] = EField_Device[ind]
			Cmaxp[2] = PtxDevice[ind]
			Cfreq[2] = values
	if 5500 < values < 8000:
	#This is for None 2 (5500 - 8000)
		if EField_Device[ind] > N2max[0]:
			N2max[2] = N2max[1]
			N2maxp[2] = N2maxp[1]
			N2freq[2] = N2freq[1]
			N2max[1] = N2max[0]
			N2maxp[1] = N2maxp[0]
			N2freq[1] = N2freq[0]
			N2max[0] = EField_Device[ind]
			N2maxp[0] = PtxDevice[ind]
			N2freq[0] = values
		elif EField_Device[ind] > N2max[1]:
			N2max[2] = N2max[1]
			N2maxp[2] = N2maxp[1]
			N2freq[2] = N2freq[1]
			N2max[1] = EField_Device[ind]
			N2maxp[1] = PtxDevice[ind]
			N2freq[1] = values
		elif EField_Device[ind] > N2max[2]:
			N2max[2] = EField_Device[ind]
			N2maxp[2] = PtxDevice[ind]
			N2freq[2] = values	
	if (8000 <= values <= 14500):
		#This is for X-Band
		if(EField_Device[ind] > Xmax[0]):
			Xmax[2] = Xmax[1]
			Xmaxp[2] = Xmaxp[1]
			Xfreq[2] = Xfreq[1]
			Xmax[1] = Xmax[0]
			Xmaxp[1] = Xmaxp[0]
			Xfreq[1] = Xfreq[0]
			Xmax[0] = EField_Device[ind]
			Xmaxp[0] = PtxDevice[ind]
			Xfreq[0] = values
		elif(EField_Device[ind] > Xmax[1]):
			Xmax[2] = Xmax[1]
			Xmaxp[2] = Xmaxp[1]
			Xfreq[2] = Xfreq[1]
			Xmax[1] = EField_Device[ind]
			Xmaxp[1] = PtxDevice[ind]
			Xfreq[1] = values
		elif(EField_Device[ind] > Xmax[2]):
			Xmax[2] = EField_Device[ind]
			Xmaxp[2] = PtxDevice[ind]
			Xfreq[2] = values	

myFile = open(csvfilename, 'w')
myFile.write("Band"+","+"Frequency"+"," + "E-Field"+"," + "EIRP"+'\n')
print("Results for "+DeviceName);
print("Maximum radiated emssions for HERA");
for i in range(3):
	print(str(i+1)+" >> "+str(Hfreq[i]) + " >> " + str(Hmax[i]) + " >> " + str(Hmaxp[i]))
	myFile.write("HERA"+","+str(Hfreq[i]) + "," + str(Hmax[i]) + "," + str(Hmaxp[i])+'\n')
print("Maximum radiated emssions for None 1(200 - 580)");
for i in range(3):
	print(str(i+1)+" >> "+str(N1freq[i]) + " >> " + str(N1max[i]) + " >> " + str(N1maxp[i]))
	myFile.write("None1 (200 - 580)"+","+str(N1freq[i]) + "," + str(N1max[i]) + "," + str(N1maxp[i])+'\n')
print("Maximum radiated emssions for UHF");
for i in range(3):
	print(str(i+1)+" >> "+str(Ufreq[i]) + " >> " + str(Umax[i]) + " >> " + str(Umaxp[i]))
	myFile.write("UHF Band"+","+str(Ufreq[i]) + "," + str(Umax[i]) + "," + str(Umaxp[i])+'\n')
print("Maximum radiated emssions for L-Band");
for i in range(3):
	print(str(i+1)+" >> "+str(Lfreq[i]) + " >> " + str(Lmax[i]) + " >> " + str(Lmaxp[i]))
	myFile.write("L-Band"+","+str(Lfreq[i]) + "," + str(Lmax[i]) + "," + str(Lmaxp[i])+'\n')
print("Maximum radiated emssions for S-Band");
for i in range(3):
	print(str(i+1)+" >> "+str(Sfreq[i]) + " >> " + str(Smax[i]) + " >> " + str(Smaxp[i]))
	myFile.write("S-Band"+","+str(Sfreq[i]) + "," + str(Smax[i]) + "," + str(Smaxp[i])+'\n')
print("Maximum radiated emssions for None2(5500 - 8000)");
for i in range(3):
	print(str(i+1)+" >> "+str(N2freq[i]) + " >> " + str(N2max[i]) + " >> " + str(N2maxp[i]))
	myFile.write("None2 (5500 - 8000)"+","+str(N2freq[i]) + "," + str(N2max[i]) + "," + str(N2maxp[i])+'\n')
print("Maximum radiated emssions for C-BASS");
for i in range(3):
	print(str(i+1)+" >> "+str(Cfreq[i]) + " >> " + str(Cmax[i]) + " >> " + str(Cmaxp[i]))
	myFile.write("C-BASS"+","+str(Cfreq[i]) + "," + str(Cmax[i]) + "," + str(Cmaxp[i])+'\n')
print("Maximum radiated emssions for X-Band");
for i in range(3):
	print(str(i+1)+" >> "+str(Xfreq[i]) + " >> " + str(Xmax[i]) + " >> " + str(Xmaxp[i]))
	myFile.write("X-Band"+","+str(Xfreq[i]) + "," + str(Xmax[i]) + "," + str(Xmaxp[i])+'\n')
myFile.close() 
