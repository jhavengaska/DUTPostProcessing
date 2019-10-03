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
#matplotlib.use('WxAgg')
#matplotlib.rcParams['font.size'] = 16
#import tkinter
#from tkinter import filedialog
#from tkinter import messagebox
import matplotlib.pyplot as plt
import RFIEq as rfi
import sys
import csv

# >>>>>>>>>>>>>> THIS IS THE PART YOU MUST/CAN EDIT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Remember to edit these according to what you are testing
DeviceName = "Fusion_Splicer" 
outputfilename = "Fusion_Splicer"
NumberPoints = 64001 #Points
StartFreq = 80e6 #In Hz
StopFreq  = 6e9 #In Hz
#BASEname = "HERA_BG_NO_ATT.csv"
BASEname = "FusionSplicer_BG_10dBmRef_100k.csv"	#The filename of the CSV that contains the baseline measurements
DUTname = "FusionSplicer_ChargingandUsing_10dBmRef_100k.csv"	#The filename of the CSV that contains the DUT measurements
AntEff = 0.8 #Antenna Efficiency, used for E-Field calculations
Dist = 10 #Distance from DUT for E-Field calculations (m)
Dist2 = 3
PLOTSARAS = 0
CISPR = 0
SHIELD = 0
CISPRCSV = "CSPRFailures.txt"
SARASCSV = "SARASFailures.txt"
CSV_EIRP = 0
NORMAL = 1
MAXBANDEMISSIONS = 1

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
	
	
"""Frequency band """
StepSize = (StopFreq - StartFreq)/(NumberPoints -1)
frequency = np.arange(StartFreq,(StopFreq+StepSize),(StepSize))
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
NoLNA = np.zeros(frequency.size)
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

"""
titel = DeviceName + " Spectrum Analyzer Raw Data"
plt.plot(Frequency,Data,'k',Frequency,data,'b')
plt.title(titel,fontsize=12)
plt.xlabel('Frequency [MHz]',fontsize=8)
plt.ylabel('Power [dBm]',fontsize=8)	
plt.legend([DeviceName,'Reverberation Chamber/Baseline'],loc=1,fontsize=8)
plt.grid(True,which="both",ls="--")
plt.ylim(minP,maxP)
"""


"""Figure will display the Power(dBm)&Efield(dBuV/m) of the device vs. the reverberation chamber's at different frequencies""" 
#plt.figure(figsize=(10,5)) 
#ax1=plt.subplot(2,1,1)                              
#titel = DeviceName + " Spectrum Analyzer Raw Data"
#plt.plot(Frequency,Data,'k',Frequency,data,'b')
#plt.title(titel,fontsize=12)
#plt.xlabel('Frequency [MHz]',fontsize=8)
#plt.ylabel('Power [dBm]',fontsize=8)	
#plt.legend([DeviceName,'Reverberation Chamber/Baseline'],loc=1,fontsize=8)
#plt.grid(True,which="both",ls="--")
#plt.ylim(minP,maxP)
#ax2=plt.subplot(2,1,2) 
#plt.plot(Frequency,PtxDevice,'k',Frequency,PtxBaseline,'b')
#plt.plot(Frequency,PtxDevice,'k')
#titel = DeviceName + " EIRP"
#plt.title(titel,fontsize=12)
#plt.xlabel('Frequency [MHz]',fontsize=8)
#plt.ylabel('EIRP [dBm]',fontsize=8)	
#plt.legend([DeviceName,'RC Background'],loc =1,fontsize=8)
#plt.grid(True,which="both",ls="--")
#plt.ylim(minTxP,maxTxP)
#plt.subplots_adjust(hspace=0.4)
#plt.savefig(savefilename,bbox_inches='tight')
#plt.show()

if(CISPR):
	B = []
	def H(Frequency):
		return (30 if 30<Frequency<230 else 37 if 230<=Frequency <=1000 else 50 if 1000<Frequency <=3000 else 54 if 3000<Frequency <=6000 else np.nan)        
	for i in range(len(Frequency)):
		B.append(H(Frequency[i]))
		
	B3m = []
	def H(Frequency):
		return (70 if 1000<Frequency<3000 else 74 if Frequency >=3000 else np.nan)        
	for i in range(len(Frequency)):
		B3m.append(H(Frequency[i]))   
		
	"""Calibrated data (Power) """
	EField_CISPR_Device, PtxCISPRDevice = rfi.calcEfPowRCCISPR(Data, NewCableLoss, RChamber, LNAGain, AntEff,Frequency)
	EField_CISPR_Baseline, PtxCISPRBaseline = rfi.calcEfPowRCCISPR(data, NewCableLoss, RChamber, LNAGain, AntEff,Frequency)

	maxFail = 0;
	maxFailF = 0;
	
	myFile = open(CISPRCSV, 'a')
	myFile.write("/******************************************************************************************\\\n")
	myFile.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"+DeviceName+"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
	for i,f in enumerate(Frequency):
		if EField_CISPR_Device[i] > B[i]:
			myFile.write("The device emits above the CISPR threshold level\n")
			print("The device emits above the CISPR threshold level")
			Edif = EField_CISPR_Device[i] - B[i]
			if(Edif > maxFail):
				maxFail = Edif
				maxFailF = f
			myFile.write(str(f) + " MHz >> " + str(EField_CISPR_Device[i]) + "(Efield) >> " + str(B[i]) + "(CISPR Level), difference = " + str(Edif)+"\n")
			print(str(f) + " MHz >> " + str(EField_CISPR_Device[i]) + "(Efield) >> " + str(B[i]) + "(CISPR Level), difference = " + str(Edif))
	myFile.write("Max Failure: "+str(maxFail)+" dB @ " +str(maxFailF)+" MHz\n")
	myFile.write("\******************************************************************************************/\n")
	myFile.close() 
	
	
	#Get min values for plots
	minE = min(EField_CISPR_Baseline)
	#minE = min(EField_Device)
	if(min(EField_CISPR_Device) < minE):
		minE = min(EField_CISPR_Device)
	minE = minE - 7	

	minP = min(data)
	#minP = min(Data)
	if(min(Data) < minP):
		minP = min(Data)
	minP = minP - 7	

	minTxP = min(PtxCISPRBaseline)
	#minTxP = min(PtxDevice)
	if(min(PtxCISPRDevice) < minTxP):
		minTxP = min(PtxCISPRDevice)
	minTxP = minTxP - 7

	#Get max values for plots
	maxE = 54
	#maxE = max(EField_Device)
	#maxE = max(EField_Baseline)
	if(max(EField_CISPR_Device) > maxE):
		maxE = max(EField_CISPR_Device)
	maxE = maxE + 7

	maxP = max(data)
	#maxP = max(Data)
	if(max(Data) > maxP):
		maxP = max(Data)
	maxP = maxP + 20	

	maxTxP = max(PtxCISPRBaseline)
	#maxTxP = max(PtxDevice)
	if(max(PtxCISPRDevice) > maxTxP):
		maxTxP = max(PtxCISPRDevice)
	maxTxP = maxTxP + 20
	#print(str(minE)+" and " + str(maxE))
	#print(str(minP)+" and " + str(maxP))

	
	titel = DeviceName + "\n Electric Field Strength @ 10 m for Freq < 1 GHz and 3 m for Freq > 1 GHz"
	plt.figure(figsize=(10,5))
	#plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b')
	#plt.plot(Frequency,EField_Device,'-k',Frequency,B,'-r')
	plt.plot(Frequency,EField_CISPR_Device,'-k',Frequency,EField_CISPR_Baseline,'b',Frequency,B,'-r')
	plt.legend(["DUT", "Baseline","CISPR-22 (Class B)"],loc=4,fontsize=8)
	plt.title(titel,fontsize=12)
	plt.xlabel('Frequency [MHz]',fontsize=8)
	plt.ylabel('Electric Field Strength\n [dBuV/m]',fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.ylim(minE,maxE),            
	savefilename = outputfilename+"_CISPR_B_EField"+".png"
	plt.savefig(savefilename,bbox_inches='tight', dpi=600)
	plt.show()

if(PLOTSARAS):
	Atten = 80.0	
	D = 10.0
	
	"""Calibrated data (Power) """
	EField_Device, PtxDevice = rfi.calcEfPowRC(Data, NewCableLoss, RChamber, LNAGain, AntEff, D)
	EField_Baseline, PtxBaseline = rfi.calcEfPowRC(data, NewCableLoss, RChamber, LNAGain, AntEff, D)

	#EField_Device = EField_Device - Atten
	#EField_Baseline = EField_Baseline - Atten

	freq, E_cont_threshold, E_spec_threshold, P_cont_threshold, P_spec_threshol =\
		rfi.get_SARAS(Dl=500.0, Dh=D)
		
	#plt.figure(figsize=(15, 10))
	#plt.plot(freq,E_spec_threshold,'r',freq,E_cont_threshold,'c')
	#plt.show()

	#for i,n in enumerate(E_spec_threshold):
	#	print(str(freq[i]) + "-->" + str(n))

	E_spec_threshold = np.interp(Frequency,freq, E_spec_threshold)
	E_cont_threshold = np.interp(Frequency,freq, E_cont_threshold)
	
	E_spec_threshold = E_spec_threshold + Atten
	E_cont_threshold = E_cont_threshold + Atten
	

	if(SHIELD):
		Shield_spec = EField_Device - E_spec_threshold
		Shield_cont = EField_Device - E_cont_threshold

		#print(str(len(E_spec_threshold)))

		#plt.figure(figsize=(15, 10))
		plt.figure()
		plt.plot(Frequency,Shield_spec,'r',Frequency,Shield_cont,'b')
		plt.title("Shielding needed to pass thresholds",fontsize=12)
		plt.xlabel('Frequency [MHz]',fontsize=8)
		plt.ylabel('Shielding [dB]',fontsize=8)
		plt.legend(["Spectral line threshold","Continuum line threshold"],loc=4,fontsize=8)
		plt.savefig(outputfilename+"ShieldingNeeded.png",bbox_inches='tight')
		plt.show()

	#Get min values for plots
	minE = min(EField_Baseline)
	#minE = min(EField_Device)
	if(min(EField_Device) < minE):
		minE = min(EField_Device)
	minE = minE - 7	

	minP = min(data)
	#minP = min(Data)
	if(min(Data) < minP):
		minP = min(Data)
	minP = minP - 7	

	minTxP = min(PtxBaseline)
	#minTxP = min(PtxDevice)
	if(min(PtxDevice) < minTxP):
		minTxP = min(PtxDevice)
	minTxP = minTxP - 7

	#Get max values for plots
	maxE = max(E_cont_threshold)
	#maxE = max(EField_Device)
	#maxE = max(EField_Baseline)
	if(max(EField_Device) > maxE):
		maxE = max(EField_Device)
	maxE = maxE + 7

	maxP = max(data)
	#maxP = max(Data)
	if(max(Data) > maxP):
		maxP = max(Data)
	maxP = maxP + 20	

	maxTxP = max(PtxBaseline)
	#maxTxP = max(PtxDevice)
	if(max(PtxDevice) > maxTxP):
		maxTxP = max(PtxDevice)
	maxTxP = maxTxP + 20
	#print(str(minE)+" and " + str(maxE))
	#print(str(minP)+" and " + str(maxP))
	
	
	titel = DeviceName + " Electric Field Strength @ "+str(Dist)+"m"
	plt.figure(figsize=(10,5))
	plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b',Frequency,E_spec_threshold,'-r',Frequency,E_cont_threshold,'--r')
	plt.legend([DeviceName, "RC Background","Spectral line threshold @ "+str(D)+" m","Continuum line threshold @ "+str(D)+" m"],loc=4,fontsize=8)
	plt.title(titel,fontsize=12)
	plt.xlabel('Frequency [MHz]',fontsize=8)
	plt.ylabel('Electric Field Strength\n [dBuV/m]',fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.ylim(minE,maxE)
	savefilename = outputfilename+"_SARAS_EField"+".png"
	plt.savefig(savefilename,bbox_inches='tight', dpi=600)
	plt.show()
	
	maxFail = 0;
	maxFailF = 0;
	myFile = open(SARASCSV, 'a')
	myFile.write("/******************************************************************************************\\\n")
	myFile.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"+DeviceName+"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
	for i,f in enumerate(Frequency):
		if EField_Device[i] > E_spec_threshold[i]:
			myFile.write("The device emits above the SARAS Spectral Line threshold level\n")
			print("The device emits above the SARAS Spectral Line threshold levels")
			Edif = EField_Device[i] - E_spec_threshold[i]
			if(Edif > maxFail):
				maxFail = Edif
				maxFailF = f
			myFile.write(str(f) + " MHz >> " + str(EField_Device[i]) + "(Efield) >> " + str(E_spec_threshold[i]) + "(SARAS Spectral-Line Level), difference = " + str(Edif)+"\n")
			print(str(f) + " MHz >> " + str(EField_Device[i]) + "(Efield) >> " + str(E_spec_threshold[i]) + "(SARAS Spectral-Line Level), difference = " + str(Edif))
	myFile.write("Max Failure: "+str(maxFail)+" dB @ " +str(maxFailF)+" MHz\n")
	myFile.write("\******************************************************************************************/\n")
	myFile.close() 
	
	if(CSV_EIRP):
		myFile = open(outputfilename+"_EIRP_120kHz.csv", 'w+')
		myFile.write("Frequency(MHz), EIRP(dBm)\n")
		for i,f in enumerate(Frequency):
			myFile.write(str(f)+','+str(PtxDevice[i])+'\n')
		myFile.close()
		
if(NORMAL):
	D = 10.0
	
	"""Calibrated data (Power) """
	EField_Device, PtxDevice = rfi.calcEfPowRC(Data, NewCableLoss, RChamber, LNAGain, AntEff, D)
	EField_Baseline, PtxBaseline = rfi.calcEfPowRC(data, NewCableLoss, RChamber, LNAGain, AntEff, D)

	
	#plt.figure(figsize=(15, 10))
	#plt.plot(freq,E_spec_threshold,'r',freq,E_cont_threshold,'c')
	#plt.show()

	#for i,n in enumerate(E_spec_threshold):
	#	print(str(freq[i]) + "-->" + str(n))

	#Get min values for plots
	minE = min(EField_Baseline)
	#minE = min(EField_Device)
	if(min(EField_Device) < minE):
		minE = min(EField_Device)
	minE = minE - 7	

	minP = min(data)
	#minP = min(Data)
	if(min(Data) < minP):
		minP = min(Data)
	minP = minP - 7	

	minTxP = min(PtxBaseline)
	#minTxP = min(PtxDevice)
	if(min(PtxDevice) < minTxP):
		minTxP = min(PtxDevice)
	minTxP = minTxP - 7

	#Get max values for plots
	maxE = -100.00
	#maxE = max(EField_Device)
	#maxE = max(EField_Baseline)
	if(max(EField_Device) > maxE):
		maxE = max(EField_Device)
	maxE = maxE + 7

	maxP = max(data)
	#maxP = max(Data)
	if(max(Data) > maxP):
		maxP = max(Data)
	maxP = maxP + 20	

	maxTxP = max(PtxBaseline)
	#maxTxP = max(PtxDevice)
	if(max(PtxDevice) > maxTxP):
		maxTxP = max(PtxDevice)
	maxTxP = maxTxP + 20
	#print(str(minE)+" and " + str(maxE))
	#print(str(minP)+" and " + str(maxP))
	
	titel = DeviceName + " Electric Field Strength @ "+str(Dist)+"m"
	plt.figure(figsize=(10,5))
	plt.plot(Frequency,EField_Device,'k',Frequency,EField_Baseline,'b')
	plt.legend([DeviceName, "RC Background"],loc=4,fontsize=8)
	plt.title(titel,fontsize=12)
	plt.xlabel('Frequency [MHz]',fontsize=8)
	plt.ylabel('Electric Field Strength\n [dBuV/m]',fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.ylim(minE,maxE)
	savefilename = outputfilename+"_EField"+".png"
	plt.savefig(savefilename,bbox_inches='tight', dpi=600)
	plt.show()	
	
	titel = DeviceName + " EIRP"
	plt.figure(figsize=(10,5))
	plt.plot(Frequency,PtxDevice,'k',Frequency,PtxBaseline,'b')
#plt.plot(Frequency,PtxDevice,'k')
	titel = DeviceName + " EIRP"
	plt.title(titel,fontsize=12)
	plt.xlabel('Frequency [MHz]',fontsize=8)
	plt.ylabel('EIRP [dBm]',fontsize=8)	
	plt.legend([DeviceName,'RC Background'],loc =1,fontsize=8)
	plt.grid(True,which="both",ls="--")
	plt.ylim(minTxP,maxTxP)
	#plt.xlim(65,80)
	#plt.subplots_adjust(hspace=0.4)
	savefilename = outputfilename+"_EIRP"+".png"
	plt.savefig(savefilename,bbox_inches='tight', dpi=600)
	plt.show()

	if(MAXBANDEMISSIONS):
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

		csvfilename = outputfilename+"_maxband.csv"
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
		myFile.write("*****,******,******,******"+'\n')	
		for i,n in enumerate(Frequency):
			myFile.write("----"+","+str(n) + "," + str(EField_Device[i]) + "," + str(PtxDevice[i])+'\n')
		myFile.close() 