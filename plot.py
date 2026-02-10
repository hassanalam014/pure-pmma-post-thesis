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

def ThermodynamicPlots(P,T,R,M,**kwargs):
	
	print T

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	
	Pratio=Pstarstar/Pstar
	Tratio=Tstarstar/Tstar
	Rratio=Rstarstar/Rstar

	#Following are My Theory Equations with, in general, v!=v_0:
	F=(Rratio*exp(-Tratio**2/(Ttilde*Pratio)))/(1+Rratio*exp(-Tratio**2/(Ttilde*Pratio)))
	A=(1/T)*((1+(Ptilde/Rtilde**2))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))
	E=(Pstar/Rstar)*(-Rtilde+((Tratio*Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	H=(Pstar/Rstar)*(-Rtilde+((Tratio*Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))+(Ptilde/Rtilde))
	S=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	S_term1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde))
	S_term2=(Pstar/(Rstar*Tstar))*(-((ln(Rtilde))/r))
	S_term3=(Pstar/(Rstar*Tstar))*(((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde)))))
	S_term4=(Pstar/(Rstar*Tstar))*(((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	Cp=((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
	Cv=(Pstar/(Rstar*Tstar))*(((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))
	dS_dT_p=(1/T)*(((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
	dS_dT_v=(1/T)*((Pstar/(Rstar*Tstar))*((((((Tratio**3)*Rratio/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+(Rratio*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
	d2S_dT2_p=(Pstar/(Rstar*Tstar*(T**2)))*((((Tstarstar/(T))**2)*((Rratio*exp(-Tstarstar/(T)))/(1+Rratio*exp(-Tstarstar/(T))))*(1-((Rratio*exp(-Tstarstar/(T)))/(1+Rratio*exp(-Tstarstar/(T)))))*(((Tstarstar/(T))*(1-(2*((Rratio*exp(-Tstarstar/(T)))/(1+Rratio*exp(-Tstarstar/(T)))))))-3))+((Tstarstar/(T))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	)))

	#Following is Equations of Condo Paper for Pure Polymer:
	dP_dT=5.12
	dPtilde_dT=dP_dT/Pstar
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((Rratio)+2))-1)/r)-((r-2)/r)*(ln(1-(((Rratio)*exp(-epsilon_2/(kB*T)))/(1+(Rratio)*exp(-epsilon_2/(kB*T)))))-((((Rratio)*exp(-epsilon_2/(kB*T)))/(1+(Rratio)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))
	# S_condo_again=(Pstar/(Rstar*Tstar))*(-(((1-Rtilde)/Rtilde)*(ln(1-Rtilde)))-((ln(Rtilde))/r)+((ln(r))/r)-1-(((ln(2/(Rratio+2)))-1)/r)-(((r-2)/r)*((ln(1-((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))))-((((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T)))))
	dS_dTg_condo=(1+(ln(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(ln(1-Rtilde)+Rtilde))-(((((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(((1-(((Rratio)*exp(-epsilon_2/(kB*T)))/(1+((Rratio)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
	# dS_dTg_condo_again=(((r-2)/r)*((epsilon_2/(kB*T))**2)*(((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T))))))*(1-((Rratio*(exp(-epsilon_2/(kB*T))))/(1+(Rratio*(exp(-epsilon_2/(kB*T)))))))/T)*((Rtilde)**2)*(((Ttilde/Rtilde)*((1/r)+(Rtilde/(1-Rtilde))))-2))+((((ln(1-Rtilde))/(Rtilde))+1-(1/r))*((dPtilde_dT)+((1/Tstar)*((ln(1-Rtilde))+((1-(1/r))*Rtilde)))))

	return F,A,E,H,S,S_term1,S_term2,S_term3,S_term4,Cp,Cv,dS_dT_p,dS_dT_v,d2S_dT2_p,S_condo,dS_dT_condo


#Condo Random Parameters for PMMA
Pstar = 503.0
Tstar = 696.0
Rstar = 1.269

# Kier PS Parameters
# Pstar = 421.76
# Tstar = 687.78
# Rstar = 1.118

#Condo Flexibiity Parameters
Vratio=1.0    #Vratio = v/v_0
epsilon_2=7444.0
z=5.0
Rratio=z-2.0
Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio
Pstarstar=Pratio*Pstar
Rstarstar=Rratio*Rstar

P=0.101325
M=2580000000.0

#Initializing the array of densities.
T=npy.linspace(100,600,300)
# T_numerical = npy.linspace(100,600,299)

R=npy.zeros(len(T))
F=npy.zeros(len(T))
A=npy.zeros(len(T))
E=npy.zeros(len(T))	
H=npy.zeros(len(T))	
S=npy.zeros(len(T))
S_term1=npy.zeros(len(T))		
S_term2=npy.zeros(len(T))		
S_term3=npy.zeros(len(T))		
S_term4=npy.zeros(len(T))
Cp=npy.zeros(len(T))
Cv=npy.zeros(len(T))
dS_dT_p=npy.zeros(len(T))	#dS/dT at constant P
dS_dT_v=npy.zeros(len(T))	#dS/dT at constant V
d2S_dT2_p=npy.zeros(len(T))
S_condo=npy.zeros(len(T))
dS_dT_condo=npy.zeros(len(T))

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
output_folder = 'Plots'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

for i in range(0,len(T)):
	R[i]=density(P,T[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

for i in range(0,len(T)):
	F[i],A[i],E[i],H[i],S[i],S_term1[i],S_term2[i],S_term3[i],S_term4[i],Cp[i],Cv[i],dS_dT_p[i],dS_dT_v[i],d2S_dT2_p[i],S_condo[i],dS_dT_condo[i] = ThermodynamicPlots(P,T[i],R[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar,Pstarstar=Pstarstar,Tstarstar=Tstarstar,Rstarstar=Rstarstar)

# d2S_dT2_p_numerical=(npy.diff(dS_dT_p)/npy.diff(T))
# d2S_dT2_v_numerical=npy.diff(dS_dT_v)/npy.diff(T)
# dF_dT_numerical=npy.diff(F)/npy.diff(T)

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

plt.plot(T,R,'k',color='y',lw=linewidth,ls='-',label='R')
# plt.plot(T,F,'k',color='y',lw=linewidth,ls='-',label='F')
# plt.plot(T,E,'k',color='g',lw=linewidth,ls='-',label='E')
# plt.plot(T,H,'k',color='m',lw=linewidth,ls='-',label='H')
# plt.plot(T,S,'k',color='k',lw=linewidth,ls='-',label='S')
# plt.plot(T,S_term1,'k',color='b',lw=linewidth,ls='-',label='S_term1')
# plt.plot(T,S_term2,'k',color='g',lw=linewidth,ls='-',label='S_term2')
# plt.plot(T,S_term3,'k',color='r',lw=linewidth,ls='-',label='S_term3')
# plt.plot(T,S_term4,'k',color='m',lw=linewidth,ls='-',label='S_term4')
# plt.plot(T,Cp,'k',color='c',lw=linewidth,ls='-',label='Cp')
# plt.plot(T,Cv,'k',color='r',lw=linewidth,ls='-',label='Cv')
# plt.plot(T,dS_dT_p,'k',color='y',lw=linewidth,ls='-',label='dS_dT_p')
# plt.plot(T,dS_dT_v,'k',color='m',lw=linewidth,ls='-',label='dS_dT_v')
# plt.plot(T,d2S_dT2_p,'k',color='r',lw=linewidth,ls='-',label='d2S_dT2_p')
# plt.plot(T,S_condo,'k',color='b',lw=linewidth,ls='-',label='S_condo')
# plt.plot(T,dS_dT_condo,'k',color='g',lw=linewidth,ls='-',label='dS_dT_condo')

# plt.plot(T_numerical,d2S_dT2_p_numerical,'k',color='b',lw=linewidth,ls='-',label='d2S_dT2_p_numerical')
# plt.plot(T_numerical,d2S_dT2_v_numerical,'k',color='g',lw=linewidth,ls='-',label='d2S_dT2_v_numerical')

plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'TD Property',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()

