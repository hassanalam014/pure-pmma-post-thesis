# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
#from matplotlib.ticker import AutoMinorLocator
# from all_p_params import *
# from loadSpecificHeatExperimentalData import *
from lmfit import minimize, Parameters, report_fit
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
from sympy import Symbol, nsolve
import sympy
import mpmath


def density(P,T,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	phi = bisect(pureEOSResidual,0.000000001,0.9999999999,args=(P,T,M,Pstar,Tstar,Rstar))
	
	R = phi*Rstar
		
	return R

#Programs Specific for Condo Paper 

def pureEntropyResidual(T,P,M,z,epsilon_2,Pstar,Tstar,Rstar):    #S/nN
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Ptilde=P/Pstar
	Rtilde=R/Rstar
	vtilde=1/Rtilde

	# f=((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))
	
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T)))))-((((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))
	# d2S_dT2_P=(Pstar/(Rstar*Tstar))*(1/((T**2)))*(((epsilon_2)/(kB*T))*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T))))*(1-((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))*(((epsilon_2)/(kB*T))*(1-2*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))-3)+((1/(Ttilde))*((T*((1/T)*((1+(Ptilde/Rtilde**2))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*(3*(1+(Ptilde/(Rtilde**2)))-3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r))+(((Ttilde*T*((1/T)*((1+(Ptilde/Rtilde**2))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))/(Rtilde))*(((Rtilde**2)/(1-Rtilde)**2)-(1/r))))))



	return S_condo

def glassTempFromEntropy(P,M,z,epsilon_2,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(pureEntropyResidual,45,1000000000,args=(P,M,z,epsilon_2,Pstar,Tstar,Rstar))
	
	return Tg

def pureDerivativeResidual(T,P,M,dP_dT,z,epsilon_2,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	vtilde=1/Rtilde
	dPtilde_dT=dP_dT/Pstar
	
	#f=((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))
	# res=(1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-((f*epsilon_2)/(kB*T**2))*(1+((1-f)*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
	res=(1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-(((((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(1+((1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)

	return res

def glassTempFromDeriv(P,M,dP_dT,z,epsilon_2,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(pureDerivativeResidual,100,10000,args=(P,M,dP_dT,z,epsilon_2,Pstar,Tstar,Rstar))
	
	return Tg


# def plots(P,T,M,**kwargs):
	
# 	for key,value in kwargs.items():
# 		exec "%s=%s" % (key,value)
	
# 	r = (Pstar*M)/(kB*Tstar*Rstar)
	
# 	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

# 	Ptilde=P/Pstar
# 	Ttilde=T/Tstar
# 	Rtilde=R/Rstar

# 	#Pstarstar=epsilon_2/v
# 	Tstarstar=epsilon_2/kB
# 	Rstarstar=z*Rstar
	
# 	# Pratio=Pstarstar/Pstar
# 	Pratio=Tratio
# 	Tratio=Tstarstar/Tstar
# 	Rratio=Rstarstar/Rstar
	
# 	E0=(Pstar/Rstar)*(-Rtilde+((Tratio*(z)*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+(z)*exp(-(Tratio**2)/(Pratio*Ttilde)))))
# 	H0=(Pstar/Rstar)*(-Rtilde+((Tratio*(z)*exp(-(Tratio**2)/(Pratio*Ttilde)))/(1+(z)*exp(-(Tratio**2)/(Pratio*Ttilde))))+(Ptilde/Rtilde))
# 	S0=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*(z)*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+(z)*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+(z)*exp(-(Tratio**2)/(Pratio*Ttilde)))))
# 	F0=((z)*exp(-Tratio**2/(Ttilde*Pratio)))/(1+(z)*exp(-Tratio**2/(Ttilde*Pratio)))
# 	Cv0=(Pstar/(Rstar*Tstar))*(((((Tratio**3)*(z)/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+((z)*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))
# 	Deriv_S_const_V=(1/T)*((Pstar/(Rstar*Tstar))*((((((Tratio**3)*(z)/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+((z)*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
# 	Cp0=((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*(z)/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+((z)*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2)))
# 	Deriv_S_const_P=(1/T)*(((Pstar/(Rstar*Tstar))*((((1+(Ptilde/(Rtilde**2)))**2)/(((Ttilde/Rtilde)*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))))+(Pstar/(Rstar*Tstar))*((((((Tratio**3)*(z)/Pratio)/(Ttilde**2))*(exp(-(((Tratio**2)/Pratio)/Ttilde))))/((1+((z)*(exp(-(((Tratio**2)/Pratio)/Ttilde)))))**2))))
	
# 	#Following is Entropy of Condo for Pure Polymer:
# 	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-F0)-(F0*Tratio/Ttilde)))

# 	return S_condo,S0,F0,Cv0,Deriv_S_const_V,Cp0,Deriv_S_const_P,E0,H0


### Programs as per my theory

def glassTransitionCriteria(T,P,M,z,epsilon_2,v_ratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar
	v_0=kB*Tstar/Pstar
	v=v_ratio*v_0

	Tstarstar=epsilon_2/kB
	Tratio=Tstarstar/Tstar
	Pratio=Tratio
	Rratio=z
	
	#MY Theory
	# S_myTheory=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde)))))
	# Own_criteria_2=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-v_ratio*(((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))))
	# Own_criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-v_ratio-(((v_ratio*Pratio)/Tratio)*ln(1+Rratio)))

	# d2S_dT2_P=(Pstar/(Rstar*Tstar*(T**2)))*(((((v*epsilon_2)/(v_0*kB*T))**2)*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T))))*(1-((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))*((((v*epsilon_2)/(v_0*kB*T))*(1-(2*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))))-3))+((v/(v_0*Ttilde))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	)))
	# dCp_dT_P=(Pstar/(Rstar*Tstar*(T)))*(((((v*epsilon_2)/(v_0*kB*T))**2)*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T))))*(1-((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))*((((v*epsilon_2)/(v_0*kB*T))*(1-(2*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))))-2))+((v/(v_0*Ttilde))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(2*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))-2)))
	# d2S_dT2_P_1=(Pstar/(Rstar*Tstar*(T**2)))*(((((v*epsilon_2)/(v_0*kB*T))**2)*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T))))*(1-((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))*((((v*epsilon_2)/(v_0*kB*T))*(1-(2*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/(1+z*exp(-(v*epsilon_2)/(v_0*kB*T)))))))-3)))
	# d2S_dT2_P_2=(Pstar/(Rstar*Tstar*(T**2)))*((v/(v_0*Ttilde))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	))
	# dS_dT_P=(Pstar/(Rstar*Tstar))*((1/T)*((((1+(Ptilde/((Rtilde)**2)))**2)/((Ttilde/(Rtilde))*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+((v_0/v)*(((v*epsilon_2)/(v_0*kB*T))**2)*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/((1+(z*exp(-(v*epsilon_2)/(v_0*kB*T))))**2)))))
	# dS_dT_V=(Pstar/(Rstar*Tstar))*((1/T)*(((v_0/v)*(((v*epsilon_2)/(v_0*kB*T))**2)*((z*exp(-(v*epsilon_2)/(v_0*kB*T)))/((1+(z*exp(-(v*epsilon_2)/(v_0*kB*T))))**2)))))

	#Condo Theory
	#dS_dT_condo=((1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-(((((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(1+((1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde))
	S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T)))))-((((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))

	res=S_condo

	return res

def glassTemp(P,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,z,epsilon_2,v_ratio,Pstar,Tstar,Rstar))
	
	return Tg

def ResidualArray(params,P,Tg,M):
	
	Pstar = params['Pstar'].value
	Tstar = params['Tstar'].value
	Rstar = params['Rstar'].value
	epsilon_2 = params['epsilon_2'].value
	z = params['z'].value
	v_ratio = params['v_ratio'].value
	kwargs = {'Pstar':Pstar,'Tstar':Tstar,'Rstar':Rstar,'epsilon_2':epsilon_2,'z':z,'v_ratio':v_ratio}
	print z
	print v_ratio
	print epsilon_2
	residual=npy.zeros(len(P))
	#Tg_calculated=npy.zeros(len(P))

	for j in range(0,len(P)):
		Tg_calculated = glassTemp(P[j],M,**kwargs)
		residual[j] = (Tg[j]-Tg_calculated)
	
	return residual


#For PS

# Pstar = 357
# Tstar = 735
# Rstar = 1.105

# For PMMA
Pstar = 503.0
Tstar = 696.0
Rstar = 1.269

P=0.101325
T=378.0
M=110000000000.0

R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
#epsilon_2=eps_2(P,T,R,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)



r = (Pstar*M)/(kB*Tstar*Rstar)

Ptilde=P/Pstar
Ttilde=T/Tstar
Rtilde=R/Rstar
vtilde=1/Rtilde	
# Pratio=Pstarstar/Pstar
# Tratio=Tstarstar/Tstar
# (z-2)=Rstarstar/Rstar

dP_dT=5.12
dPtilde_dT=dP_dT/Pstar
# dT=50.0

# T_condo=T+dT  #Any random value
# P_condo=dP_dT*(T_condo-T)+P
# print P_condo
# R_condo=density(P_condo,T_condo,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
# R_condotilde=R_condo/Rstar
# P_condotilde=P_condo/Pstar
# T_condotilde=T_condo/Tstar
# R_condotilde=R_condo/Rstar
# v_condotilde=1/R_condotilde	
# f=(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))
# F0=(((z-2)*exp(-Tratio**2/(Ttilde*Pratio)))/(1+(z-2)*exp(-Tratio**2/(Ttilde*Pratio))))
# print R
# print R_condo

#f=((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T))))
# A0=((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))


# d2S_dT2_P=(Pstar/(Rstar*Tstar*(T**2)))*((((epsilon_2/(kB*T))**2)*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T))))*(1-((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))*(((epsilon_2/(kB*T))*(1-(2*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))))-3))+((epsilon_2/(kB*T))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	)))
# d2S_dT2_P_line=(Pstar/(Rstar*Tstar*(T_condo**2)))*((((epsilon_2/(kB*T_condo))**2)*((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo))))*(1-((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo)))))*(((epsilon_2/(kB*T_condo))*(1-(2*((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo)))))))-3))+((epsilon_2/(kB*T_condo))*((T_condo*((1/T_condo)*((1+(P_condotilde/(R_condotilde**2)))/(((T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))-2))))**2)*R_condotilde*((3*(1+(P_condotilde/(R_condotilde**2))))-(3*(T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))+((T_condotilde*T_condo*((1/T_condo)*((1+(P_condotilde/(R_condotilde**2)))/(((T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))-2)))/R_condotilde)*(((R_condotilde/(1-R_condotilde))**2)-(1/r)))	)))

mpmath.mp.dps = 15
z = Symbol('z')
epsilon_2 = Symbol('epsilon_2')
S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T)))))-((((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))
# S_condo_line=(Pstar/(Rstar*Tstar))*(-((1-R_condotilde)*(ln(1-R_condotilde))/R_condotilde)-((ln(R_condotilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-(((z-2)*exp(-epsilon_2/(kB*T_condo)))/(1+(z-2)*exp(-epsilon_2/(kB*T_condo)))))-((((z-2)*exp(-epsilon_2/(kB*T_condo)))/(1+(z-2)*exp(-epsilon_2/(kB*T_condo))))*epsilon_2/(kB*T_condo))))
# d2S_dT2_P=(Pstar/(Rstar*Tstar*(T**2)))*((((epsilon_2/(kB*T))**2)*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T))))*(1-((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))*(((epsilon_2/(kB*T))*(1-(2*((z*exp(-epsilon_2/(kB*T)))/(1+z*exp(-epsilon_2/(kB*T)))))))-3))+((epsilon_2/(kB*T))*((T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2))))**2)*Rtilde*((3*(1+(Ptilde/(Rtilde**2))))-(3*(Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))+((Ttilde*T*((1/T)*((1+(Ptilde/(Rtilde**2)))/(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))/Rtilde)*(((Rtilde/(1-Rtilde))**2)-(1/r)))	)))
# d2S_dT2_P_line=(Pstar/(Rstar*Tstar*(T_condo**2)))*((((epsilon_2/(kB*T_condo))**2)*((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo))))*(1-((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo)))))*(((epsilon_2/(kB*T_condo))*(1-(2*((z*exp(-epsilon_2/(kB*T_condo)))/(1+z*exp(-epsilon_2/(kB*T_condo)))))))-3))+((epsilon_2/(kB*T_condo))*((T_condo*((1/T_condo)*((1+(P_condotilde/(R_condotilde**2)))/(((T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))-2))))**2)*R_condotilde*((3*(1+(P_condotilde/(R_condotilde**2))))-(3*(T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))+((T_condotilde*T_condo*((1/T_condo)*((1+(P_condotilde/(R_condotilde**2)))/(((T_condotilde/R_condotilde)*((R_condotilde/(1-R_condotilde))+(1/r)))-2)))/R_condotilde)*(((R_condotilde/(1-R_condotilde))**2)-(1/r)))	)))
res=(1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-(((((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(1+((1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
print(nsolve((res, S_condo), (z, epsilon_2), (10.00, 7483.0),verify=False))

'''
z=5.79499958093261  #4.81366805048692				#4.86914200453447				#4.81366805059756
epsilon_2=11513.8982105149  #7114.32278870336		#7210.44342406770			#7114.32278892748

# res=(1+(math.log(1-Rtilde)/Rtilde))*(1/Rtilde)*(dPtilde_dT+(1/Tstar)*(math.log(1-Rtilde)+Rtilde))-(((((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T)))))*epsilon_2)/(kB*T**2))*(1+((1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+((z-2)*exp(-epsilon_2/(kB*T))))))*epsilon_2)/(kB*T))*(2*Rtilde-(Ttilde/(1-Rtilde))+Ttilde)
# S_condo=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)-((1/r)*ln(1/r))-1-((ln(2/((z-2)+2))-1)/r)-((r-2)/r)*(ln(1-(((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T)))))-((((z-2)*exp(-epsilon_2/(kB*T)))/(1+(z-2)*exp(-epsilon_2/(kB*T))))*epsilon_2/(kB*T))))
# print S_condo
# print res

P_line = npy.linspace(0.101325,200,15)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))
#d2S_dT2_P_line=npy.zeros(len(P_line))


#Ideal Experimental Straight Line Data
for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)


#								#MyOwnCriteria_1 (Best fit)				#Condo Theory			#My Theory d2S/dT2|p=0 with v!=v_0	##My Theory dCp/dT|p=0 with v!=v_0
# z=5.0							#3.0									#4.915257075699675 		#0.3411518312632549					#0.04394803
# epsilon_2=7443.0			#7444.52718 							#7235.098320856251 		#10818.738735660516					#39795.7441
# v_ratio=1.47042331				#0.37042331								#----------------		#1.0143895679458907					#40.2622941
#Initializing the array of densities.
P0 = npy.linspace(0.101325,200,15)
R0=npy.zeros(len(P0))
Tg_From_S=npy.zeros(len(P0))		
Tg_From_Deriv=npy.zeros(len(P0))
Tg_calculated=npy.zeros(len(P0))

for i in range(0,len(P0)):

	Tg_From_S[i]=glassTempFromEntropy(P0[i],M,z,epsilon_2,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	#Tg_From_Deriv[i]=glassTempFromDeriv(P0[i],M,dP_dT,z,epsilon_2,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	# Tg_calculated[i]=glassTemp(P0[i],M,z=z,epsilon_2=epsilon_2,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	print P0[i]

#Experiment Data of PMMA
P1=[0,30,40,80,120,140,180]
Tg1=[378.0,386.5,386.5,397.5,408.1,408.1,419.7]

# for i in range(0,len(P0)):
# 	R0[i]=density(P0[i],T0[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)	
	
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
output_folder = 'plot_Condo'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

#plt.plot(P0,Tg_From_S,'k',color='g',lw=linewidth,ls='-',label='T Versus P of Pure PMMA Condo from S')
#plt.plot(P0,Tg_From_Deriv,'k',color='r',lw=linewidth,ls='-',label='T Versus P of Pure PMMA Condo from Deriv')

plt.plot(P0,Tg_From_S,'k',color='g',lw=linewidth,ls='-',label='Condo Theory Curve')
plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Pure PMMA Ideal Straight Line')
plt.plot(P1,Tg1,'sk',ms=markersize,label='Exp. Data from Condo Ref 53')

# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()

print z
print epsilon_2
'''
