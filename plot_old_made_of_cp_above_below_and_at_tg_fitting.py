# Date: 2019
#
# Description: The purpose of this file is to plot Polystyrene (PS) Thermodynamics Properties
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
#from matplotlib.ticker import AutoMinorLocator
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from findVectors import findVectors
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from wrapperFlexibilityFunctions import calculateSpecificHeat
from isListOrNpyArray import *
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
from optimizeResidualFunctions import pureEOSResidual,pureChemicalPotentialResidual
from loadSpecificHeatExperimentalData import *

def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.0000000000000001,0.9999999999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

def plot(P,T,R,M,**kwargs):     #S/m*k_B  and **kwargs must contain "three" characteristic parameters and "three" flexibility parameters.
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	Pratio=Pstarstar/Pstar
	Tratio=Tstarstar/Tstar
	Rratio=Rstarstar/Rstar
	epsilon_2=kB*Tstarstar
	z=Rratio+2
	# print Pratio
	# print Tratio
	# print epsilon_2
	# print Rratio+2
	A0=(1/T)*((1+(Ptilde/Rtilde**2))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))
	# A0=((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))
	E0=(Pstar/Rstar)*(-Rtilde+((Tratio*Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	H0=(Pstar/Rstar)*(-Rtilde+((Tratio*Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))+(Ptilde/Rtilde))
	S0=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	S0_term1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde))
	S0_term2=(Pstar/(Rstar*Tstar))*(-((ln(Rtilde))/r))
	S0_term3=(Pstar/(Rstar*Tstar))*(((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde)))))
	S0_term4=(Pstar/(Rstar*Tstar))*(((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))

	#((Tratio**2)/(Pratio*Ttilde))
	#3*(1+(Ptilde/(Rtilde**2)))-3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r))+(((Ttilde*T*A0)/(Rtilde))*(((Rtilde**2)/(1-Rtilde)**2)-(1/r)))
	F0=(Rratio*exp(-Tratio**2/(Ttilde*Pratio)))/(1+Rratio*exp(-Tratio**2/(Ttilde*Pratio)))
	#wrong # d2S_dT2_P=(Pstar/(Rstar*Tstar))*(Pratio/(Tratio*(T**2)))*(((Tratio**2)/(Pratio*Ttilde))*F0*(1-F0)*(((Tratio**2)/(Pratio*Ttilde))*(1-2*F0)-3)+((Tratio)/(Pratio*Ttilde))*((T*A0)**2)*Rtilde*(3*(1+(Ptilde/(Rtilde**2)))-3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r))+(((Ttilde*T*A0)/(Rtilde))*(((Rtilde**2)/(1-Rtilde)**2)-(1/r)))))
	d2S_dT2_P=(Pstar/(Rstar*Tstar*(T**2)))*((((epsilon_2/(kB*T))**2)*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T))))*(1-((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))*(((epsilon_2/(kB*T))*(1-(2*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))))-3))+((epsilon_2/(kB*T))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	)))

	Cv0=(Pstar/(Rstar*Tstar))*(((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))
	Deriv_S_const_V=(1/T)*((Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
	Cp0=((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	Deriv_S_const_P=(1/T)*(((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
	#C2_with v equal v0=(Pstar/(Rstar*Tstar))*((((((Tratio**2)*Rratio)/(Ttilde**2))*(exp(-(Tratio/Ttilde))))/((1+(Rratio*(exp(-(Tratio/Ttilde)))))**2)))
	#C2=(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	#C2=(Pstar/(Rstar*Tstar))*(((((   (Tratio**2)*r*Rratio/N   )/(Ttilde**2))*(exp(-(((r*Tratio)/N)/Ttilde))))/((1+(Rratio*(exp(-(((r*Tratio)/N)/Ttilde)))))**2)))
	#Following is Entropy of Condo for Pure Polymer:
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/(Rratio+2))-1)/r)-((r-2)/r)*(ln(1-F0)-(F0*Tratio/Ttilde)))
	dP_dT=5.12
	dPtilde_dT=dP_dT/Pstar
	dS_dT_condo=((1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-(((((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(1+((1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde))


	return S0,F0,Cv0,Deriv_S_const_V,Cp0,Deriv_S_const_P,E0,H0,S_condo,d2S_dT2_P,S0_term1,S0_term2,S0_term3,S0_term4,dS_dT_condo


#Condo Random Parameters for PMMA
Pstar = 503.0
Tstar = 696.0
Rstar = 1.269

# Kier PS Parameters
# Pstar = 421.76
# Tstar = 687.78
# Rstar = 1.118

# Below Tg:

Tstarstar1=895.344
Rstarstar1=3.0*Rstar
Tratio=Tstarstar1/Tstar
Pstarstar1=Pstar*Tratio

# Pstarstar1=0.00249935
# Tstarstar1=11.9757060
# Rstarstar1=1.1013e+40

#Above Tg:
# Pstarstar2=1.2855e-05
# Tstarstar2=1.58456913
# Rstarstar2=3.999e+141

#Excluding Tg:
# Pstarstar3=1.7274e-04
# Tstarstar3=4.36981412
# Rstarstar3=1.3155e+80

# P0=0.0001
# M0=258000

P0=0.101325
M0=2580000000


#Initializing the array of densities.
T0 = npy.linspace(100,600,300)
S1=npy.zeros(len(T0))		#Total Entropy below Tg
S2=npy.zeros(len(T0))		#Total Entropy
S3=npy.zeros(len(T0))		#Total Entropy
R0=npy.zeros(len(T0))

S0_term1=npy.zeros(len(T0))		
S0_term2=npy.zeros(len(T0))		
S0_term3=npy.zeros(len(T0))		
S0_term4=npy.zeros(len(T0))		

S1_condo=npy.zeros(len(T0))		#Total Entropy below Tg
S2_condo=npy.zeros(len(T0))		#Total Entropy
S3_condo=npy.zeros(len(T0))		#Total Entropy

dS_dT_condo=npy.zeros(len(T0))

d2S_dT2_P1=npy.zeros(len(T0))
d2S_dT2_P2=npy.zeros(len(T0))
d2S_dT2_P3=npy.zeros(len(T0))

F1=npy.zeros(len(T0))		#Fraction of Excited Segments
E1=npy.zeros(len(T0))	
H1=npy.zeros(len(T0))	
Cv1=npy.zeros(len(T0))
Cp1=npy.zeros(len(T0))
dS_v1=npy.zeros(len(T0))	#dS/dT at constant V
dS_p1=npy.zeros(len(T0))	#dS/dT at constant P

F2=npy.zeros(len(T0))		#Fraction of Excited Segments
E2=npy.zeros(len(T0))	
H2=npy.zeros(len(T0))	
Cv2=npy.zeros(len(T0))
Cp2=npy.zeros(len(T0))
dS_v2=npy.zeros(len(T0))	#dS/dT at constant V
dS_p2=npy.zeros(len(T0))	#dS/dT at constant P

F3=npy.zeros(len(T0))		#Fraction of Excited Segments
E3=npy.zeros(len(T0))	
H3=npy.zeros(len(T0))	
Cv3=npy.zeros(len(T0))
Cp3=npy.zeros(len(T0))
dS_v3=npy.zeros(len(T0))	#dS/dT at constant V
dS_p3=npy.zeros(len(T0))	#dS/dT at constant P


#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'plot_condo'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

for i in range(0,len(T0)):
	R0[i]=density(P0,T0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T0)):
	S1[i],F1[i],Cv1[i],dS_v1[i],Cp1[i],dS_p1[i],E1[i],H1[i],S1_condo[i],d2S_dT2_P1[i],S0_term1[i],S0_term2[i],S0_term3[i],S0_term4[i],dS_dT_condo[i] = plot(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar1,Tstarstar=Tstarstar1,Rstarstar=Rstarstar1)
	#S2[i],F2[i],Cv2[i],dS_v2[i],Cp2[i],dS_p2[i],E2[i],H2[i],S2_condo[i],d2S_dT2_P2[i] = plot(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar2,Tstarstar=Tstarstar2,Rstarstar=Rstarstar2)
	#S3[i],F3[i],Cv3[i],dS_v3[i],Cp3[i],dS_p3[i],E3[i],H3[i],S3_condo[i],d2S_dT2_P3[i] = plot(P0,T0[i],R0[i],M0,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar3,Tstarstar=Tstarstar3,Rstarstar=Rstarstar3)

# d2S_p1=(npy.diff(dS_p1)/npy.diff(T0))
# d2S_v1=npy.diff(dS_v1)/npy.diff(T0)
# dF_dT=npy.diff(F1)/npy.diff(T0)
#random=R0*S1

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

T1 = npy.linspace(100,600,299)

# plt.plot(T0,R0,'k',color='y',lw=linewidth,ls='-',label='R0')

# plt.plot(T1,d2S_p1,'k',color='b',lw=linewidth,ls='-',label='d2S_p1')
# plt.plot(T1,d2S_v1,'k',color='g',lw=linewidth,ls='-',label='d2S_v1')
# plt.plot(T0,d2S_dT2_P1,'k',color='r',lw=linewidth,ls='-',label='d2S_dT2_P1_Formula')
plt.plot(T0,S1_condo,'k',color='b',lw=linewidth,ls='-',label='S_condo')
# plt.plot(T0,S2_condo,'k',color='b',lw=linewidth,ls='-',label='S_condo_Above_Tg')
# plt.plot(T0,S3_condo,'k',color='g',lw=linewidth,ls='-',label='S_condo_Exluding Tg')

plt.plot(T0,dS_dT_condo,'k',color='g',lw=linewidth,ls='-',label='dS_dT_condo')

# plt.plot(T0,S0_term1,'k',color='b',lw=linewidth,ls='-',label='S_term1')
# plt.plot(T0,S0_term2,'k',color='g',lw=linewidth,ls='-',label='S_term2')
# plt.plot(T0,S0_term3,'k',color='r',lw=linewidth,ls='-',label='S_term3')
# plt.plot(T0,S0_term4,'k',color='m',lw=linewidth,ls='-',label='S_term4')

# plt.plot(T0,S1,'k',color='k',lw=linewidth,ls='-',label='S_Below_Tg')
# plt.plot(T0,S2,'k',color='b',lw=linewidth,ls='-',label='S_Above_Tg')
# plt.plot(T0,S3,'k',color='g',lw=linewidth,ls='-',label='S_Exluding Tg')

# plt.plot(T0,F1,'k',color='g',lw=linewidth,ls='-',label='F1_Below_Tg')
# plt.plot(T0,F2,'k',color='b',lw=linewidth,ls='-',label='F2_Above_Tg')
# plt.plot(T0,F3,'k',color='g',lw=linewidth,ls='-',label='F3_Exluding Tg')

# plt.plot(T0,E1,'k',color='g',lw=linewidth,ls='-',label='E1_Below_Tg')
# plt.plot(T0,E2,'k',color='b',lw=linewidth,ls='-',label='E2_Above_Tg')
# plt.plot(T0,E3,'k',color='g',lw=linewidth,ls='-',label='E3_Exluding Tg')

# plt.plot(T0,H1,'k',color='g',lw=linewidth,ls='-',label='H1_Below_Tg')
#plt.plot(T0,H2,'k',color='b',lw=linewidth,ls='-',label='H2_Above_Tg')
#plt.plot(T0,H3,'k',color='g',lw=linewidth,ls='-',label='H3_Exluding Tg')

# plt.plot(T0,Cv1,'k',color='r',lw=linewidth,ls='-',label='Cv1_Below_Tg')
# plt.plot(T0,Cv2,'k',color='b',lw=linewidth,ls='-',label='Cv2_Above_Tg')
# plt.plot(T0,Cv3,'k',color='g',lw=linewidth,ls='-',label='Cv3_Exluding Tg')

# plt.plot(T0,dS_v1,'k',color='m',lw=linewidth,ls='-',label='dS_v1_Below_Tg')
# plt.plot(T0,dS_v2,'k',color='b',lw=linewidth,ls='-',label='dS_v2_Above_Tg')
# plt.plot(T0,dS_v3,'k',color='g',lw=linewidth,ls='-',label='dS_v3_Exluding Tg')

# plt.plot(T0,Cp1,'k',color='c',lw=linewidth,ls='-',label='Cp1_Below_Tg')
# plt.plot(T0,Cp2,'k',color='b',lw=linewidth,ls='-',label='Cp2_Above_Tg')
# plt.plot(T0,Cp3,'k',color='g',lw=linewidth,ls='-',label='Cp3_Exluding Tg')

# plt.plot(T0,dS_p1,'k',color='y',lw=linewidth,ls='-',label='dS_p1_Below_Tg')
# plt.plot(T0,dS_p2,'k',color='b',lw=linewidth,ls='-',label='dS_p2_Above_Tg')
# plt.plot(T0,dS_p3,'k',color='g',lw=linewidth,ls='-',label='dS_p3_Exluding Tg')

#plt.axvline(x=380,lw=0.5,color='k', linestyle='-.')
#plt.axvline(x=370,lw=0.5,color='k', linestyle='-.')
plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'Entropy',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_infinitekilo_below_Tg_Enthalpy'+img_extension,dpi=img_dpi)

plt.show()

