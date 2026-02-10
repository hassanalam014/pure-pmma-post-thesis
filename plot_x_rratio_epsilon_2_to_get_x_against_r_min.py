# Date: 2019
#
# Description: The purpose of this file is to ..............
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from math import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculateCharacteristicParameters,returnCharacteristicParameters
# from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity
# from wrapperFlexibilityFunctions import calculateSpecificHeat
# from isListOrNpyArray import *
from Parameters_of_Different_Polymers import *
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

### Programs My Theory

def glassTransitionCriteria(T,P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)
	
	R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Pratio=Tratio/Vratio

	Tstarstar=Tratio*Tstar
	Pstarstar=Pratio*Pstar
	Rstarstar=Rratio*Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	res=Own_Criteria_1

	return res

def glassTemp(P,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	Tg = bisect(glassTransitionCriteria,100,10000,args=(P,M,x,Rratio,Tratio,Vratio,Pstar,Tstar,Rstar))
	
	return Tg

def glassTransitionDerivative(epsilon_2,P,T,R,M,Rratio,Vratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	dS_dT_Own_Criteria_1=(Vratio*Ttilde*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((1/Rtilde)*(dPtilde_dTtilde)))/(Rtilde*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+(((Vratio*epsilon_2/(kB*T))**2)*((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))*(1-((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))))
	# print dS_dT_Own_Criteria_1
	res=dS_dT_Own_Criteria_1

	return res

def eps_2(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	epsilon_2 = bisect(glassTransitionDerivative,3000,8000,args=(P,T,R,M,Rratio,Vratio,Pstar,Tstar,Rstar))
	
	return epsilon_2


def OwnCriteria1(x,P,T,R,M,Pratio,Tratio,Rratio,Pstar,Tstar,Rstar):  
	
	r = (Pstar*M)/(kB*Tstar*Rstar)

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	Rtilde=R/Rstar

	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))

	res=Own_Criteria_1

	return res


def x_atm(P,T,R,M,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	epsilon_2 = bisect(OwnCriteria1,0.01,0.99,args=(P,T,R,M,Pratio,Tratio,Rratio,Pstar,Tstar,Rstar))
	
	return epsilon_2

P = P_atm
T=Tg_atm
M=M_infinity
R=density(P,T,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
r = (Pstar*M)/(kB*Tstar*Rstar)
dP_dT_atm=1/dTg_dP_atm

Ptilde=P/Pstar
Ttilde=T/Tstar
Rtilde=R/Rstar
dPtilde_dT=dP_dT_atm/Pstar
dPtilde_dTtilde=dP_dT_atm*Tstar/Pstar

Vratio=1.0
#####################################################################################################

######################################################################################################
x = npy.linspace(0.23,0.7,20)
epsilon_2=npy.zeros(len(x))
Rratio=npy.zeros(len(x))

for i in range(0,len(x)): 
	#Simultaneous Equation Solver
	mpmath.mp.dps = 15
	Rrat = Symbol('Rrat')
	eps_2 = Symbol('eps_2')
	Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rrat*(exp(-((Vratio*eps_2))/(kB*T)))/(1+Rrat*exp(-((Vratio*eps_2))/(kB*T))))+((1/Vratio)*ln(1+Rrat*exp(-(Vratio*eps_2)/(kB*T))))-(x[i])-((x[i]/Vratio)*ln(1+Rrat)))
	dS_dT_Own_Criteria_1=(Vratio*Ttilde*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((1/Rtilde)*(dPtilde_dTtilde)))/(Rtilde*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+(((Vratio*eps_2/(kB*T))**2)*((Rrat*exp(-(Vratio*eps_2)/(kB*T)))/(1+Rrat*exp(-(Vratio*eps_2)/(kB*T))))*(1-((Rrat*exp(-(Vratio*eps_2)/(kB*T)))/(1+Rrat*exp(-(Vratio*eps_2)/(kB*T))))))
	answer=nsolve((Own_Criteria_1, dS_dT_Own_Criteria_1), (Rrat, eps_2), (10.00, 7483.0),verify=False)

	print answer
	Rratio[i]=answer[0]
	epsilon_2[i]=answer[1]
	# print Rratio
	# print epsilon_2



Tstarstar=npy.zeros(len(x))
Tratio=npy.zeros(len(x))
Pratio=npy.zeros(len(x))
Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio

print Rratio
print epsilon_2
print x
# #To Verify:
# Own_Criteria_1=(Pstar/(Rstar*Tstar))*(-((1-Rtilde)*(ln(1-Rtilde))/Rtilde)-((ln(Rtilde))/r)+((1/Ttilde)*Rratio*(exp(-((Tratio)**2)/(Pratio*Ttilde)))/(1+Rratio*exp(-((Tratio)**2)/(Pratio*Ttilde))))+((Pratio/Tratio)*ln(1+Rratio*exp(-(Tratio**2)/(Pratio*Ttilde))))-(x)-((((x)*Pratio)/Tratio)*ln(1+Rratio)))
# dS_dT_Own_Criteria_1=(Vratio*Ttilde*((1-(1/r))+((ln(1-Rtilde))/Rtilde))*(((1-(1/r))+((ln(1-Rtilde))/Rtilde))+((1/Rtilde)*(dPtilde_dTtilde)))/(Rtilde*(((Ttilde/Rtilde)*((Rtilde/(1-Rtilde))+(1/r)))-2)))+(((Vratio*epsilon_2/(kB*T))**2)*((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))*(1-((Rratio*exp(-(Vratio*epsilon_2)/(kB*T)))/(1+Rratio*exp(-(Vratio*epsilon_2)/(kB*T))))))

# print Own_Criteria_1
# print Own_Criteria_11
# print dS_dT_Own_Criteria_1

'''


epsilon_2=eps_2(P,T,R,M,Rratio=Rratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
print epsilon_2
Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio

x=x_atm(P,T,R,M,Pratio=Pratio,Tratio=Tratio,Rratio=Rratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
print x
######################################################################################################

# Vratio=1.0
# x=0.539866552563
# epsilon_2=4952.31845932
# Rratio=2.57080854813003

P_line = npy.linspace(0.101325,200,15)
T_line = npy.zeros(len(P_line))
R_line=npy.zeros(len(P_line))
#d2S_dT2_P_line=npy.zeros(len(P_line))


#Ideal Experimental Straight Line Data
for i in range(0,len(P_line)):
	T_line[i]=((P_line[i]-P)/dP_dT_atm)+T
	#R_line[i]=density(P_line[i],T_line[i],M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)

###############################################################################

###############################################################################

# epsilon_2 = 6319.72664027		#PMMA 140kilo=6251.31324				#PS=6480.99556
# Rratio =	1.19880282		#PMMA 140kilo=1.20630412				#PS=1.85399320
# e=2.718281
# x=0.32105263			#PMMA 140kilo= 0.359891658197082    	#PS=0.381795#1/e

Tstarstar=epsilon_2/kB
Tratio=Tstarstar/Tstar
Pratio=Tratio/Vratio

#Initializing the array of densities.
P = npy.linspace(0.101325,300,20)
R=npy.zeros(len(P))
# Tg_From_S=npy.zeros(len(P))		
# Tg_From_Deriv=npy.zeros(len(P0))
Tg_calculated=npy.zeros(len(P))

for i in range(0,len(P)):

	# Tg_From_S[i]=CondoGlassTempFromEntropy(P0[i],M,z,epsilon_2,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	#Tg_From_Deriv[i]=CondoGlassTempFromDeriv(P0[i],M,dP_dT_atm,z,epsilon_2,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	Tg_calculated[i]=glassTemp(P[i],M=M,x=x,Rratio=Rratio,Tratio=Tratio,Vratio=Vratio,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	print P[i]

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
output_folder = 'plots'

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

plt.plot(P,Tg_calculated,'k',color='g',lw=linewidth,ls='-',label='Fit')
plt.plot(P_line,T_line,'k',color='r',lw=linewidth,ls='-',label='Ideal Straight Line')
plt.plot(Pg_exp,Tg_exp,'sk',ms=markersize,label='Exp. Data')

# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature (K)',fontsize=axis_size)
# plt.axis([0,50,350,400])
plt.legend(loc=4,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

#figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_Tg vs P'+img_extension,dpi=img_dpi)

plt.show()

'''
