"""**************************************************************************************************************************************************************************************************************
 Date : 28/02/2019                                                                                                                                                
 Author: jhavenga                                                                                                                  
 Version: 1                                                                                                                         
All of the helpfull RFI equations for quick scripting to help me has been placed inside this python file.
**************************************************************************************************************************************************************************************************************"""
import numpy as np

def dBm2dBuv(dBm, Zrf):
	"""
	# Function to convert between dBm and dBuV. The input parameters are the dBm (power) and
	# the Zrf (input impedance)[usually 50 ohm] of the measurement equipment.
	"""
	dBuv = dBm + 10*np.log10(Zrf) - 20*np.log10(1e-6) + 10*np.log10(1e-3)
	return(dBuv)
	
def dBuV2dBm(dBuV, Zrf):
	"""
	Function to convert between dBuV and dBm. The input parameters are the dBuV (voltage) and
	the Zrf (input impedance)[usually 50 ohm] of the measurement equipment.
	"""
	dBm = dBuV - 10*np.log10(Zrf) + 20*np.log10(1e-6) - 10*np.log10(1e-3)
	return(dBm)

#Function that returns the SARAS levels
#Created by Dr. A.J. Otto
def get_SARAS(Dl=500, Dh=1, Dm=1):
	"""
	SARAS Continuum and Spectral Line
	Inputs  :   Dl=500m
		Dh=1m (1m Threshold)
	Returns :   [0] Frequency
		[1] E Continuum
		[2] E Spectral Line
		[3] PSD Continuum
		[4] PSD Spectral Line 
	"""
	N = 10
	freq_low_low = np.linspace(50, 300, N)
	freq_low = np.linspace(301, 2000, N)
	freq_high = np.linspace(2001, 25000, N)

	''' CONTINUUM SPECIFICATION '''
	SARAS_cont_low_low = -17.2708 * np.log10(freq_low_low)-192.0714
	SARAS_cont_low = -17.2708 * np.log10(freq_low)-192.0714
	SARAS_cont_high = -0.065676 * np.log10(freq_high)-248.8661

	''' SPECTRAL LINE SPECIFICATION '''
	SARAS_spec_low_low = SARAS_cont_low_low + 15.
	SARAS_spec_low = SARAS_cont_low + 15.
	SARAS_spec_high = SARAS_cont_high + 15.

	''' RBW's '''
	RBW_cont_low_low = 10.*np.log10((1./100.) * freq_low_low * 1.E6)
	RBW_cont_low = 10.*np.log10((1./100.) * freq_low * 1.E6)
	RBW_cont_high = 10.*np.log10((1./100.) * freq_high * 1.E6)

	RBW_spec_low_low = 10.*np.log10((0.001/100.) * freq_low_low * 1.E6)
	RBW_spec_low = 10.*np.log10((0.001/100.) * freq_low * 1.E6)
	RBW_spec_high = 10.*np.log10((0.001/100.) * freq_high * 1.E6)

	''' PATH LOSS '''
	c0 = 3E8
	if Dl != 0:
		D = Dl
		pathloss_low_low = 10 * np.log10(((4*np.pi*D)/(c0/(freq_low_low * 1.E6)))**2)
	else:
		pathloss_low_low = 0

	if Dh != 0:
		D = Dh
		pathloss_low = 10 * np.log10(((4*np.pi*D)/(c0/(freq_low * 1.E6)))**2)
		pathloss_high = 10 * np.log10(((4*np.pi*D)/(c0/(freq_high * 1.E6)))**2)

		''' PSD THRESHOLD LEVELS '''
		PSD_cont_thresh_low_low = SARAS_cont_low_low + pathloss_low_low
		PSD_cont_thresh_low = SARAS_cont_low + pathloss_low
		PSD_cont_thresh_high = SARAS_cont_high + pathloss_high

		PSD_spec_thresh_low_low = SARAS_spec_low_low + pathloss_low_low
		PSD_spec_thresh_low = SARAS_spec_low + pathloss_low
		PSD_spec_thresh_high = SARAS_spec_high + pathloss_high
	else:
		''' PSD THRESHOLD LEVELS '''
		PSD_cont_thresh_low_low = SARAS_cont_low_low + pathloss_low_low
		PSD_cont_thresh_low = SARAS_cont_low 
		PSD_cont_thresh_high = SARAS_cont_high

		PSD_spec_thresh_low_low = SARAS_spec_low_low + pathloss_low_low
		PSD_spec_thresh_low = SARAS_spec_low
		PSD_spec_thresh_high = SARAS_spec_high


	''' E-FIELD THRESHOLD LEVELS '''
	''' E-field at distance Dm '''

	Dm = 10.
	E_cont_low_low = 20.*np.log10(np.sqrt(((10.**((PSD_cont_thresh_low_low)/10.)*0.001) * ((1./100.) * freq_low_low * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)
	E_cont_low = 20.*np.log10(np.sqrt(((10.**((PSD_cont_thresh_low)/10.)*0.001) * ((1./100.) * freq_low * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)
	E_cont_high = 20.*np.log10(np.sqrt(((10.**((PSD_cont_thresh_high)/10.)*0.001) * ((1./100.) * freq_high * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)

	E_spec_low_low = 20.*np.log10(np.sqrt(((10.**((PSD_spec_thresh_low_low)/10.)*0.001) * ((0.001/100.) * freq_low_low * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)
	E_spec_low = 20.*np.log10(np.sqrt(((10.**((PSD_spec_thresh_low)/10.)*0.001) * ((0.001/100.) * freq_low * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)
	E_spec_high = 20.*np.log10(np.sqrt(((10.**((PSD_spec_thresh_high)/10.)*0.001) * ((0.001/100.) * freq_high * 1.E6) *377.) / (4*np.pi*Dm**2.)) / 1E-6)

	freq = []
	freq.extend(freq_low_low)
	freq.extend(freq_low)
	freq.extend(freq_high)

	E_cont_threshold = []
	for a in range(0, len(E_cont_low_low)):
		E_cont_threshold.append(E_cont_low_low[a])
	for a in range(0, len(E_cont_low)):
		E_cont_threshold.append(E_cont_low[a])
	for a in range(0, len(E_cont_high)):
		E_cont_threshold.append(E_cont_high[a])

	E_spec_threshold = []
	for a in range(0, len(E_spec_low_low)):
		E_spec_threshold.append(E_spec_low_low[a])
	for a in range(0, len(E_spec_low)):
		E_spec_threshold.append(E_spec_low[a])
	for a in range(0, len(E_spec_high)):
		E_spec_threshold.append(E_spec_high[a])

	P_cont_threshold = []
	for a in range(0, len(RBW_cont_low_low)):
		P_cont_threshold.append(PSD_cont_thresh_low_low[a] + RBW_cont_low_low[a])
	for a in range(0, len(RBW_cont_low)):
		P_cont_threshold.append(PSD_cont_thresh_low[a] + RBW_cont_low[a])
	for a in range(0, len(RBW_cont_high)):
		P_cont_threshold.append(PSD_cont_thresh_high[a] + RBW_cont_high[a])

	P_spec_threshold = []
	for a in range(0, len(RBW_spec_low_low)):
		P_spec_threshold.append(PSD_spec_thresh_low_low[a] + RBW_spec_low_low[a])
	for a in range(0, len(RBW_spec_low)):
		P_spec_threshold.append(PSD_spec_thresh_low[a] + RBW_spec_low[a])
	for a in range(0, len(RBW_spec_high)):
		P_spec_threshold.append(PSD_spec_thresh_high[a] + RBW_spec_high[a])

	return freq, E_cont_threshold, E_spec_threshold, P_cont_threshold, P_spec_threshold

def calcEfPowOS(Prsa, Lcable, Glna, Aeff, r):
	"""
	Calculate E-Field and EIRP for Open Site measurement
	Calculate the E-Field at distance [r] and the EIRP of the open site test measurements made at distance [r].
	This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB),
	the gain from the LNA(dB) and the frequencies of operation(Hz) to caluclate the E-field at
	the receiving antenna.

	"""
	Vrsa = dBm2dBuv(Prsa,50)
	Vlessloss = Prsa + Lcable - Glna	
	Efield = Vlessloss + Aeff
	EIRP = Efield + 20*np.log10(r) - 104.68
	return(Efield, EIRP)
	
def calcEfPowRC(Prsa, Lcable, Linsertion, Glna, Aeff, r):
	""" 
	Calculate E-Field and Power for Reverberation chamber measurement.
	Calculate the E-Field at distance [r] and the EIRP of the DUT tested in the reverberation chamber.
	This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB), the reverberation chamber calibration factor(dB),
	the gain from the LNA(dB) and the frequencies of operation(Hz) to caluclate the E-field at
	the receiving antenna.
	"""
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
	
def calcEfPowRCCISPR(Prsa, Lcable, Linsertion, Glna, Aeff, Freq):
	""" 
	Calculate E-Field and Power for Reverberation chamber measurement to apply CISPR.
	Calculate the E-Field at distance of 10 m for freq up to 1 GHz and for a distance 
	of 3 m for freqyency between 1 GHz and 6 GHz and the EIRP of the DUT tested in the reverberation chamber.
	This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB), the reverberation chamber calibration factor(dB),
	the gain from the LNA(dB) and the frequencies of operation(Hz) to caluclate the E-field at
	the receiving antenna.
	"""
	eps0 = 8.854E-12
	mu0 = 4*np.pi*1E-7
	#c0 = 1./np.sqrt(eps0 * mu0)
	Zo= np.sqrt(mu0/eps0)
	#Zo =377                                                                                    # Free Space Impendance
	Linsertion = (Linsertion - 10*np.log(Aeff)) * -1.00 #Given as a loss in negative dB of the ACF of the Chamber
	Plessloss = Prsa + Lcable + Linsertion - Glna
	Efield = np.ones(Plessloss.size)
	for i,F in enumerate(Freq):
		if F < 1000:
			r = 10
			CFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90   #Conversion factor (from dBm to dBuV/m)
			Efield[i] = Plessloss[i] + CFactor
		else:
			r = 3
			CFactor = 10*np.log10(Zo/(4*(np.pi)*(r**2))) + 90   #Conversion factor (from dBm to dBuV/m)
			Efield[i] = Plessloss[i] + CFactor
	#print(CFactor)
	return(Efield, Plessloss)
	
def calcEfPowAC(Prsa, Lcable, Glna, Gant, freq, r):
	""" 
	Calculate E-Field and Power for anechoic chamber measurement.
	Calculate the E-Field at distance [r] and the EIRP of the DUT tested in the reverberation chamber.
	This function takes the measured power(dBm) from the spectrum analyzer
	along with the cable losses (dB), the reverberation chamber calibration factor(dB),
	the gain from the LNA(dB) and the frequencies of operation(Hz) to caluclate the E-field at
	the receiving antenna.
	"""
	eps0 = 8.854E-12
	mu0 = 4*np.pi*1E-7
	c0 = 1./np.sqrt(eps0 * mu0)
	Zo= np.sqrt(mu0/eps0)                                                                               # Free Space Impendance
	Efield = Prsa + 20.0*np.log10(freq) - Gant + Lcable - Glna + 10.0*np.log10((4.0*np.pi*Zo)/(c0**2)) + 90
	EIRP = Efield + 20*np.log10(r) - 104.68
	return(Efield, EIRP)